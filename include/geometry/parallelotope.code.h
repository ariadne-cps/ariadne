/***************************************************************************
 *            parallelotope.code.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to bouthe Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
#include <iostream>
#include <string>
#include <sstream>

#include <vector>

#include "numeric/exceptions.h"
#include "numeric/rational.h"
#include "numeric/interval.h"

#include "base/stlio.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/qr_matrix.h"

#include "linear_programming/linear_program.h"

#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/box.h"
#include "geometry/list_set.h"
#include "geometry/zonotope.h"
#include "geometry/polyhedron.h"
#include "geometry/polytope.h"

#include "output/logging.h"

#include "parallelotope.h"


namespace {

using namespace Ariadne;




template<class R>
void
compute_orthogonal_over_approximation(Parallelotope<Interval<R>,R>& p, const Zonotope<Interval<R>,R>& z)
{
  //std::cerr << "compute_orthogonal_over_approximation(Parallelotope&,Zonotope)";
  //std::cerr << "z="<<z<<std::endl;
  dimension_type d=z.dimension();
  size_type ng=z.number_of_generators();
  p.resize(z.dimension());
  typedef Interval<R> I;
  const Point<I>& ic=z.centre();
  const Matrix<R>& G=z.generators();
  
  Point<I> c(d);
  Matrix<R> EG(d,d+ng);
  for(uint i=0; i!=d; ++i) {
    c[i]=midpoint(ic[i]);
    for(uint j=0; j!=ng; ++j) {
      EG(i,j)=G(i,j);
    }
    EG(i,ng+i)=radius(ic[i]);
  }
  //std::cerr << "c="<<c<<"\nEG=\n"<<EG<<std::endl;


  // Use QR without column pivoting
  try { 
    Matrix<I> iQ;
    Matrix<I> iR;
    QRMatrix<I> iQR(EG);
    iQ=iQR.Q();
    iR=iQR.R();
    //std::cerr << "Q=\n"<<midpoint(iQ)<<"\nR=\n"<<midpoint(iR)<<std::endl;
    //std::cerr << "QR-G=\n"<<midpoint(iQ*iR-EG)<<std::endl;
    Vector<I> rs(d);
    Matrix<I> iO=iQ;
    for(size_type i=0; i!=d; ++i) {
      for(size_type j=i; j!=d+ng; ++j) {
        rs(i)+=abs(iR(i,j));
      }
      for(size_type k=0; k!=d; ++k) {
        iO(k,i)*=rs(i);
      }
    }
    //std::cerr << "rs="<<midpoint(rs)<<std::endl;
    //std::cerr << "O="<<midpoint(iO)<<std::endl;
    
    //assert(norm(row_norms(Qmid.inverse()*A)+(cmid-c))<=R(1));
    p=Parallelotope<I,R>(c,midpoint(iO));
    //return;
  } 
  catch(DivideByZeroException& e) {
    //std::cerr << "QR(A) with A=" << Matrix<R>(EG) << ": " <<  std::flush;
    p=z.bounding_box();
  }
  



  // Don't use standard QR factorization since this does not rearrange columns
  // Instead use own version which stretches directions
  //std::cerr << "\n\n\n" << "EG="<<EG<<std::endl;
  for(uint k=0; k!=d; ++k) {
    // Find the column with the highest norm
    uint hnc=0;
    R hns=0;
    for(size_type j=k; j!=d+ng; ++j) {
      R ns=0;
      for(size_type i=0; i!=d; ++i) {
        ns=add_approx(ns,mul_approx(EG(i,j),EG(i,j)));
      }
      if(ns>hns) {
        hnc=j;
        hns=ns;
      }
    }
     
    //std::cerr << "k="<<k<<" hnc="<<hnc<<" hns="<<hns<<std::endl;
    // Swap the column with the highest norm with the current column
    if(hnc!=k) {
      R tmp;
      for(size_type i=0; i!=d; ++i) {
        tmp=EG(i,k);
        EG(i,k)=EG(i,hnc);
        EG(i,hnc)=tmp;
      }
    }

    // Normalise the current column, storing the multiplier in a scale factor sf
    R sf=sqrt_up(hns);
    R hnr=div_up(R(1),sqrt_down(hns));
    for(uint i=0; i!=d; ++i) {
      EG(i,k)=mul_approx(EG(i,k),hnr);
    }
    
    // Orthogonalise the remaining columns, adding the absolute value of the inner product to the scale factor sf
    for(uint j=k+1; j!=d+ng; ++j) {
      R ip=0;
      for(uint i=0; i!=d; ++i) {
        ip=add_approx(ip,mul_approx(EG(i,k),EG(i,j)));
      }
      for(uint i=0; i!=d; ++i) {
        EG(i,j)=sub_approx(EG(i,j),mul_approx(EG(i,k),ip));
      }
      sf=add_up(sf,abs(ip));
    }
    
    // Scale the working column with the highest norm
    for(uint i=0; i!=d; ++i) {
      if(EG(i,k)<0) {
        EG(i,k)=mul_down(EG(i,k),sf);
      } else {
        EG(i,k)=mul_up(EG(i,k),sf);
      }
    }

    // End of iteration
  }
        
  //std::cerr << "nEG="<<EG<<std::endl;

  Matrix<R> O(d,d,EG.begin(),EG.row_increment(),EG.column_increment());
  //std::cerr << "O="<<O<<std::endl;
  p=Parallelotope<I,R>(c,O);
  //std::cerr << "\n\n\n" << std::endl;
}

} // namespace





namespace Ariadne {


template<class XC, class XG>
tribool
Parallelotope<XC,XG>::_instantiate() 
{
  //Parallelotope<XC,XG>(*goa)(const Parallelotope< Interval<R> >&) = &over_approximation<R>;
  Box<R> r;
  Parallelotope<XC,XG> p;
  Zonotope<XC,XG> z;
  subset(r,p);
  over_approximation(p);
  orthogonal_over_approximation(z);
  return false;
}





template<class XC, class XG>
void
Parallelotope<XC,XG>::resize(dimension_type d) 
{
  this->Zonotope<XC,XG>::resize(d,d);
}

template<class XC, class XG>
typename Parallelotope<XC,XG>::F
Parallelotope<XC,XG>::volume() const
{
  return determinant(this->generators());
}



template<class XC, class XG>
void 
Parallelotope<XC,XG>::_compute_generators_inverse() const 
{  
  this->_generators_inverse=inverse(this->generators());
}


template<class XC, class XG>
ListSet< Parallelotope<XC,XG> >
Parallelotope<XC,XG>::divide() const 
{
  const Zonotope<XC,XG>& z=*this;
  ListSet< Zonotope<XC,XG> > zls=z.divide();
  return ListSet< Parallelotope<XC,XG> >(zls);
}

template<class XC, class XG>
ListSet< Parallelotope<XC,XG> >
Parallelotope<XC,XG>::subdivide() const 
{
  const Zonotope<XC,XG>& z=*this;
  ListSet< Zonotope<XC,XG> > zls=z.subdivide();
  return ListSet< Parallelotope<XC,XG> >(zls);
}




template<class XC, class XG>
Vector<typename Parallelotope<XC,XG>::F>
Parallelotope<XC,XG>::coordinates(const Point<R>& s) const 
{
  Vector<F> p=s.position_vector();
  Vector<F> c=this->centre().position_vector();
  Vector<F> d=p-c;
  Matrix<F> G=this->generators();
  return solve(G,d);
}




template<class R>
Parallelotope<R,R>
over_approximation(const Parallelotope<R,R>& p)
{
  return p;
}

template<class R>
Parallelotope<R,R>
over_approximation(const Parallelotope< Interval<R>,R >& ep)
{
  typedef Interval<R> I;
  const Zonotope<I,R>& ez=ep;
  Zonotope<R,R> z=over_approximation(ez);
  return static_cast<Parallelotope<R,R>&>(z);
}

template<class R>
Parallelotope<Interval<R>,R>
over_approximation(const Parallelotope< Interval<R> >& ip)
{
  typedef Interval<R> I;
  const Zonotope<I,I>& iz=ip;
  Zonotope<I,R> ez=over_approximation(iz);
  return static_cast<Parallelotope<I,R>&>(ez);
}



template<class R>
Parallelotope<R>
orthogonal_over_approximation(const Zonotope<R>& z)
{
  typedef Interval<R> I;
  Zonotope<I,R> ez(z);
  Parallelotope<I,R> ep;
  compute_orthogonal_over_approximation(ep,ez);
  return over_approximation(ep);
}


/*
template<class R>
Parallelotope<R>
orthogonal_over_approximation(const Zonotope<Interval<R>,R>& ez)
{
  ARIADNE_LOG(5,"Parallelotope<R> orthogonal_over_approximation(Zonotope<I,R>)\n");
  typedef Interval<R> I;
  typedef typename traits<R>::approximate_arithmetic_type A;
  const Point<I>& c=ez.centre();
  const Matrix<R>& G=ez.generators();
  
  size_type d=ez.dimension();
  
  Point<R> cmid=midpoint(c);
  Vector<I> cerr=c-cmid;
  
  Matrix<I> Q;
  try { 
    QRMatrix<I> QR(G);
    Matrix<I> Q=QR.Q();
  } 
  catch(DivideByZeroException& e) {
    std::cerr << "QR(A) with A=" << G << ": " <<  std::flush;
  }
  Matrix<R> Qmid=midpoint(Q);
  Vector<I> Rrwnrm=row_norms(transpose(Qmid)*G);
  for(size_type i=0; i!=d; ++i) {
    R scale=(Rrwnrm(i)+cerr(i)).upper();
    for(size_type j=0; j!=d; ++j) {
      Qmid(i,j)=mul_up(Qmid(i,j),scale);
    }
  }
  
  // Check to make such result is valid.
  //assert(norm(row_norms(Qmid.inverse()*A)+(cmid-c))<=R(1));
  return Parallelotope<R>(cmid,Qmid);
}
*/

template<class R>
Parallelotope<Interval<R>,R>
orthogonal_over_approximation(const Zonotope<Interval<R>,R>& z)
{
  ARIADNE_LOG(5,"Parallelotope<R> orthogonal_over_approximation(Zonotope<I,I>)\n");
  Parallelotope<Interval<R>,R> p;
  ::compute_orthogonal_over_approximation(p,z);
  return p;
}

template<class R>
Parallelotope< Interval<R> >
orthogonal_over_approximation(const Zonotope< Interval<R> >& iz)
{
  ARIADNE_LOG(5,"Parallelotope<R> orthogonal_over_approximation(Zonotope<I,I>)\n");
  Zonotope<Interval<R>,R> z=over_approximation(iz);
  Parallelotope<Interval<R>,R> p;
  ::compute_orthogonal_over_approximation(p,z);
  return p;
}






template<class XC, class XG>
std::ostream&
Parallelotope<XC,XG>::write(std::ostream& os) const
{
  const Parallelotope<XC,XG>& p=*this;
  if(p.dimension() > 0) {
    os << "Parallelotope( centre=" << p.centre()
       << ", directions=" << p.generators()
       << " ) ";
  }
  
  return os;
}


template<class XC, class XG>
std::istream& 
Parallelotope<XC,XG>::read(std::istream& is)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



} // namespace Ariadne
