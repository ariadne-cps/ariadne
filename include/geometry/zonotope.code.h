/***************************************************************************
 *            zonotope.code.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *
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
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
#include <iostream>
#include <vector>
#include <algorithm>

#include "zonotope.h"

#include "base/lvalue.h"
#include "base/array.h"
#include "exceptions.h"
#include "numeric/traits.h"
#include "numeric/error_float.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/diagonal_matrix.h"
#include "linear_algebra/lu_matrix.h"
#include "linear_algebra/qr_matrix.h"

#include "linear_programming/linear_program.h"

#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/rectangle.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"

#include "output/logging.h"



namespace {
  
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;

template<class R> inline
void 
accumulate(R& value, R& error, uint n, const R* aptr, const R* bptr) {
  Interval<R> v=value;
  for(uint i=0; i!=n; ++i) {
    v+=aptr[i]*bptr[i];
  }
  value=v.midpoint();
  error=add_up(error,v.radius());
}

template<class R>
void
write_vector_slice(std::ostream& os, const LinearAlgebra::VectorSlice<R>& vs) {
  os << LinearAlgebra::Vector<R>(vs);
}

template<class R>
void
write_vector_slice(std::ostream& os, const LinearAlgebra::VectorSlice<const R>& vs) {
  os << LinearAlgebra::Vector<R>(vs);
}


template<class R> inline
tribool norm_grtr(const LinearAlgebra::Vector<R>& v1, const LinearAlgebra::Vector<R>& v2) 
{
  return LinearAlgebra::norm(v1)>LinearAlgebra::norm(v2);
}



Rational med_approx(const Rational& ql, const Rational& qu) {
  return (ql+qu)/2;
}

Rational rad_up(const Rational& ql, const Rational& qu) {
  return (ql+qu)/2;
}



template<class R, class Tag>
Zonotope<Rational> rational_zonotope(const Zonotope<R,Tag>& z);

tribool contains(const Zonotope<Rational>& z, const Point<Rational>& r);
tribool disjoint(const Zonotope<Rational>& z, const Rectangle<Rational>& r);
tribool superset(const Zonotope<Rational>& z, const Rectangle<Rational>& r);
tribool subset(const Zonotope<Rational>& z, const Rectangle<Rational>& r);

ListSet< Zonotope<Rational> > subdivide(const Zonotope<Rational>&);

template<class R>
ListSet< Zonotope<R,UniformErrorTag> > subdivide(const Zonotope<R,UniformErrorTag>&);

//template<class R>
//ListSet< Zonotope<R,ExactTag> > subdivide(const Zonotope<R,ExactTag>&);


template<class R> inline
Zonotope<Rational>
rational_zonotope(const Zonotope<R,ExactTag>& z) 
{
  return Zonotope<Rational>(z.centre(),z.generators());
}

template<class R> inline
Zonotope<Numeric::Rational>
rational_zonotope(const Zonotope<R,UniformErrorTag>& z) 
{
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  typedef Interval<Rational> IRational;
  dimension_type d=z.dimension();
  size_type ng=z.number_of_generators();
  const Point<IRational>& zc=z.centre();
  const Matrix<Rational>& zG=z.generators();
  Point<Rational> nc=midpoint(zc);
  Matrix<Rational> nG(d,ng+d);
  MatrixSlice<Rational>(d,ng,nG.data().begin(),nG.row_increment(),nG.column_increment())=zG;
  for(dimension_type i=0; i!=z.dimension(); ++i) {
    nG(i,ng+i)=zc[i].radius();
  }
  return Zonotope<Rational,ExactTag>(nc,nG);
}



template<class R> inline
void
convert(Zonotope<R>& z, const Rectangle<R>& r) {
  z=r;
}


template<class R> void instantiate_zonotope();
template<> void instantiate_zonotope<Numeric::Rational>();

} // namespace 






namespace Ariadne {

extern int Geometry::verbosity;



template<class R, class Tag>
Geometry::Rectangle<R>
Geometry::bounding_box(const Zonotope<R,Tag>& z)
{
  typedef Numeric::Interval<R> I;
  LinearAlgebra::Vector<I> v=z.domain().position_vectors();
  LinearAlgebra::Vector< Interval<R> > b=z.centre().position_vector()+z.generators()*v;
  return Rectangle<R>(b);
}






template<class R>       
Geometry::Zonotope<R,ExactTag>::Zonotope(const Rectangle<R>& r) 
  : _centre(r.dimension()), _generators(r.dimension(),r.dimension())
{
  dimension_type d=r.dimension();
  Point<R>& c=this->_centre;
  LinearAlgebra::Matrix<R>& G=this->_generators;
  for(size_type i=0; i!=d; ++i) {
    c[i]=med_approx(r.lower_bound(i),r.upper_bound(i));
    for(size_type j=0; j!=d; ++j) {
      G(i,j)=0;
    }
    G(i,i)=rad_up(r.lower_bound(i),r.upper_bound(i));
  }
}

template<class R>       
Geometry::Zonotope<R,UniformErrorTag>::Zonotope(const Rectangle<R>& r) 
  : _centre(r.dimension()), _generators(r.dimension(),r.dimension())
{
  *this=r;
}


template<class R>       
Geometry::Zonotope<R,UniformErrorTag>&
Geometry::Zonotope<R,UniformErrorTag>::operator=(const Rectangle<R>& r) 
{
  dimension_type d=r.dimension();
  Point<I>& c=this->_centre;
  LinearAlgebra::Matrix<R>& G=this->_generators;
  c.resize(d);
  G.resize(d,d);
  for(size_type i=0; i!=d; ++i) {
    c[i]=med_approx(r.lower_bound(i),r.upper_bound(i));
    for(size_type j=0; j!=d; ++j) {
      G(i,j)=0;
    }
    G(i,i)=rad_up(r.lower_bound(i),r.upper_bound(i));
  }
  return *this;
}





template<class R>
tribool 
Geometry::contains(const Zonotope<R,ExactTag>& z, const Point<R>& pt)  
{
  return ::contains(::rational_zonotope(z),Point<Rational>(pt));
}

template<class R>
tribool 
Geometry::contains(const Zonotope<R,UniformErrorTag>& z, const Point<R>& pt)  
{
  return ::contains(::rational_zonotope(z),Point<Rational>(pt));
}


template<class R>
tribool
Geometry::disjoint(const Zonotope<R,ExactTag>& z, const Rectangle<R>& r)
{
  return ::disjoint(::rational_zonotope(z),Rectangle<Rational>(r));
}

template<class R>
tribool
Geometry::disjoint(const Zonotope<R,UniformErrorTag>& z, const Rectangle<R>& r)
{
  return ::disjoint(::rational_zonotope(z),Rectangle<Rational>(r));
}


template<class R>
tribool
Geometry::superset(const Zonotope<R,ExactTag>& z, const Rectangle<R>& r)
{
  return ::superset(::rational_zonotope(z),Rectangle<Rational>(r));
}

template<class R>
tribool
Geometry::superset(const Zonotope<R,UniformErrorTag>& z, const Rectangle<R>& r)
{
  return ::superset(::rational_zonotope(z),Rectangle<Rational>(r));
}



template<class R, class Tag>
tribool
Geometry::subset(const Zonotope<R,Tag>& z, const Rectangle<R>& r)
{
  return Geometry::subset(bounding_box(z),r);
}

template<class R, class Tag>
tribool
Geometry::subset(const Zonotope<R,Tag>& z, const Polyhedron<R>& p)
{
  LinearAlgebra::Vector< Interval<R> > im=(p.A()*z.generators())*z.domain().position_vectors()+(p.A()*z.centre().position_vector()+p.b());
  return im>=0;
}

template<class R,class Tag>
Geometry::ListSet< Geometry::Zonotope<R,Tag> >
Geometry::subdivide(const Zonotope<R,Tag>& z)
{
  return ::subdivide(z);
}

template<class R,class Tag>
void 
Geometry::approximate(Zonotope<R>& r, const Zonotope<R,Tag>& z) 
{
  r=Zonotope<R>(approximation(z.centre()),LinearAlgebra::approximation<R>(z.generators()));
}


template<class R,class Tag>
void 
Geometry::over_approximate(Zonotope<R,Tag>& z, const Rectangle<R>& r) 
{
  typedef Numeric::Interval<R> I;
  dimension_type d=r.dimension();
  I* cptr=z.begin();
  R* gptr=cptr+d;
  const R* rptr=r.begin();
  for(size_type i=0; i!=d; ++i) {
    cptr[i]=med_approx(r.lower_bound(i),r.upper_bound(i));
    for(size_type j=0; j!=d; ++j) {
      gptr[i*d+j]=0;
    }
    gptr[(d+1)*i]=div_up(sub_up(r.upper_bound(i),r.lower_bound(i)),2);
  }
}



template<class R> 
void
Geometry::over_approximate(Zonotope<R,ExactTag>& z, 
                           const Zonotope<R,ExactTag>& ez)
{
  z=ez;
}


template<class R> 
void
Geometry::over_approximate(Zonotope<R,ExactTag>& z, const Zonotope<R,UniformErrorTag>& ez)
{
  dimension_type d=ez.dimension();
  size_type ng=ez.number_of_generators();
  Point<R> c=midpoint(ez.centre());
  Matrix<R> G(d,ng+d);
  MatrixSlice<R>(d,ng,G.begin(),G.row_increment(),G.column_increment())=ez.generators();
  for(size_type i=0; i!=d; ++i) {
    G(i,i+ng)=ez.centre()[i].radius();
  }
  z=Zonotope<R,ExactTag>(c,G);
}


template<class R> 
void
Geometry::over_approximate(Geometry::Zonotope<R,ExactTag>& z, const Zonotope<R,IntervalTag>& iz)
{
  dimension_type d=iz.dimension();
  size_type ng=iz.number_of_generators();
  typedef Numeric::Interval<R> I;
  const Point<I>& izc=iz.centre();
  const Matrix<I>& izG=iz.generators();
  Point<R> c=midpoint(izc);
  Matrix<R> G(d,ng+d);
  MatrixSlice<R>(d,ng,G.begin(),G.row_increment(),G.column_increment())=midpoint(iz.generators());
  for(size_type i=0; i!=d; ++i) {
    G(i,i+ng)=izc[i].radius();
    for(size_type j=0; j!=ng; ++j) {
      G(i,i+ng)=add_up(G(i,i+ng),izG(i,j).radius());
    }
  }
  z=Zonotope<R,ExactTag>(c,G);
}


template<class R> 
void
Geometry::over_approximate(Geometry::Zonotope<R,UniformErrorTag>& az, const Zonotope<R,UniformErrorTag>& ez)
{
  az=ez;
}


template<class R> 
void
Geometry::over_approximate(Geometry::Zonotope<R,UniformErrorTag>& ez, const Zonotope<R,IntervalTag>& iz)
{
  dimension_type d=iz.dimension();
  size_type ng=iz.number_of_generators();
  typedef Numeric::Interval<R> I;
  const Point<I>& izc=iz.centre();
  const Matrix<I>& izG=iz.generators();
  Point<I> c=izc;
  Matrix<R> G=midpoint(iz.generators());
  MatrixSlice<R>(d,ng,G.begin(),G.row_increment(),G.column_increment())=midpoint(iz.generators());
  for(size_type i=0; i!=d; ++i) {
    R e=c[i].radius();
    for(size_type j=0; j!=ng; ++j) {
      e=add_up(e,izG(i,j).radius());
    }
    c[i]+=I(-e,e);
  }
  ez=Zonotope<R,UniformErrorTag>(c,G);
}



template<class R> 
void
Geometry::nonsingular_over_approximate(Zonotope<R,ExactTag>& z, const Zonotope<R,IntervalTag>& iz)
{
  assert(z.dimension()==z.number_of_generators());
  typedef Numeric::Interval<R> I;
  typedef Numeric::ErrorFloat<R> F;
  const dimension_type& d=z.dimension();
  const size_type& ng=z.number_of_generators();
  const Point<I>& c=z.centre();
  const Matrix<I>& G=z.generators();
  const Vector<I> e(ng,I(-1,1));
  Point<R> ac=midpoint(c);
  Matrix<R> aG=midpoint(z.generators());
  Matrix<I> aGinv=inverse(aG);
  Vector<I> nd=aGinv*(c-ac)+(aGinv*G)*e;
  for(size_type j=0; j!=ng; ++j) {
    R sf=Numeric::next_up(Numeric::next_up(nd[j].upper()));
    sf=add_approx(sf,R(0.00000000001));
    for(size_type i=0; i!=z.dimension(); ++i) {
      aG(i,j)=mul_approx(aG(i,j),sf);
    }
  }
  
  // Check result
  using namespace std;
  aGinv=inverse(aG);
  nd=aGinv*(c-ac)+(aGinv*G)*e;
  Vector<I> ue(d,I(-1,1));
  cout << "G="<<G <<"\nnG="<<aG<<"\nnGinv="<<aGinv<<"\nnGinv*G="<<aGinv*G<<endl;
  cout <<"nGinv*G*e"<<(aGinv*G)*ue<<endl;
  std::cerr<<"nd="<<nd<<std::endl;
  std::cerr<<"ue="<<ue<<std::endl;
  assert(LinearAlgebra::refines(nd,ue));
  assert(false);
  z=Zonotope<R>(ac,aG);
}  



template<class R> 
void
Geometry::orthogonal_over_approximate(Zonotope<R,UniformErrorTag>& ez, const Zonotope<R,IntervalTag>& iz)
{
  //assert(iz.dimension()==iz.number_of_generators());
  typedef Numeric::Interval<R> I;
  typedef Numeric::ErrorFloat<R> F;
  const dimension_type& d=iz.dimension();
  const Point<I>& c=iz.centre();
  const Matrix<I>& G=iz.generators();
  const Vector<I> e=iz.domain().position_vectors();
  Point<R> ac=midpoint(c);
  Matrix<R> aG=midpoint(G);
  Matrix<R> aQ,aR;
  make_lpair(aQ,aR)=qr_approx(aG);
  //std::cerr << "iG="<<G<<"\nG="<<aG<<"\nQ="<<aQ<<"\nR="<<aR<<std::endl;
  Matrix<I> aQinv=inverse(aQ);
  Vector<I> nd=aQinv*(c-ac)+(aQinv*G)*e;
  DiagonalMatrix<R> aD(radius(nd));
  //std::cerr << "D="<<aD<<std::endl;
  Matrix<R> nG=mul_approx(aQ,aD);
  Vector<I> ne(d,I(-1,1));
  //std::cerr << "Q*D="<<nG<<std::endl;
  Point<I> nc=ac+(aQ*aD-nG)*ne;

  ez=Zonotope<R,UniformErrorTag>(nc,nG);
}  

template<class R> 
void
Geometry::orthogonal_over_approximate(Zonotope<R,UniformErrorTag>& z, const Zonotope<R,UniformErrorTag>& ez)
{
  orthogonal_over_approximate(z,Zonotope<R,IntervalTag>(ez));
}

template<class R> 
void
Geometry::orthogonal_over_approximate(Zonotope<R,IntervalTag>& z, const Zonotope<R,IntervalTag>& iz)
{
  Zonotope<R,UniformErrorTag> ez;
  orthogonal_over_approximate(ez,iz);
  z=ez;
}


template<class R> 
Geometry::Zonotope<R> 
Geometry::over_approximation(const Zonotope<R>& z)
{
  return z;
}




template<class R> 
Geometry::Zonotope<R,R> 
Geometry::approximation(const Zonotope< Numeric::Interval<R>, Numeric::Interval<R> >& z)
{
  return Zonotope<R,R>(midpoint(z.centre()),
                       midpoint(z.generators()));
}


template<class R> 
Geometry::Zonotope<R,R> 
Geometry::approximation(const Zonotope<Numeric::Interval<R>,R>& z)
{
  return Zonotope<R,R>(midpoint(z.centre()),
                       z.generators());
}


template<class R> 
Geometry::Zonotope<R,R> 
Geometry::approximation(const Zonotope<R,R>& z)
{
  return z; 
}


/*
template<class R> 
Geometry::Zonotope<Numeric::Interval<R>,R> 
Geometry::orthogonal_over_approximation(const Zonotope<R,R>& z)
{
  // FIXME: Subdivide in zero order as well!
  static bool warn=true;
  if(warn) {
    std::cerr << std::endl << "WARNING: orthogonal_over_approximation(Zonotope<I,R>) does not over-approximate roundoff errors." << std::endl;
    warn=false;
  }
  Zonotope<R,R> oaz=over_approximation(z);
  
  QRMatrix< Interval<R> > QR(oaz.generators());
  Point< Interval<R> > c(oaz.centre());
  Matrix<R> G(z.dimension(),z.number_of_generators());

  Matrix< Interval<R> > q=QR.Q();
  Matrix< Interval<R> > r=QR.R();
  for(uint i=0; i!=z.dimension();++i) {
    Interval<R> a=0;
    for(uint j=i; j!=z.number_of_generators(); ++j) {
      a+=r(i,j);
    }
    for(uint k=0; k!=z.dimension(); ++k) {
      Interval<R> b=q(k,i)*a;
      G(k,i)=b.midpoint();
      c[k]+=(b-b.midpoint());
    }
  }
  return Zonotope<R,R>(midpoint(c),G);
}

template<class R> 
Geometry::Zonotope<Numeric::Interval<R>,R> 
Geometry::orthogonal_over_approximation(const Zonotope<Numeric::Interval<R>,R>& z)
{
  Zonotope<R,R> oaz=over_approximation(z);
  
  QRMatrix< Interval<R> > QR(oaz.generators());
  Point< Interval<R> > c(oaz.centre());
  Matrix<R> G(z.dimension(),z.number_of_generators());

  Matrix< Interval<R> > q=QR.Q();
  Matrix< Interval<R> > r=QR.R();
  for(uint i=0; i!=z.dimension();++i) {
    Interval<R> a=0;
    for(uint j=i; j!=z.number_of_generators(); ++j) {
      a+=r(i,j);
    }
    for(uint k=0; k!=z.dimension(); ++k) {
      Interval<R> b=q(k,i)*a;
      G(k,i)=b.midpoint();
      c[k]+=(b-b.midpoint());
    }
  }
  return Zonotope<Interval<R>,R>(c,G);
}

template<class R> 
Geometry::Zonotope< Numeric::Interval<R> > 
Geometry::orthogonal_over_approximation(const Zonotope< Numeric::Interval<R> >& z)
{
  Zonotope<R,R> oaz=over_approximation(z);
  
  QRMatrix< Interval<R> > QR(oaz.generators());
  Point< Interval<R> > c(oaz.centre());
  Matrix<R> G(z.dimension(),z.number_of_generators());

  Matrix< Interval<R> > q=QR.Q();
  Matrix< Interval<R> > r=QR.R();
  for(uint i=0; i!=z.dimension();++i) {
    Interval<R> a=0;
    for(uint j=i; j!=z.number_of_generators(); ++j) {
      a+=r(i,j);
    }
    for(uint k=0; k!=z.dimension(); ++k) {
      Interval<R> b=q(k,i)*a;
      G(k,i)=b.midpoint();
      c[k]+=(b-b.midpoint());
    }
  }
  return Zonotope< Interval<R> >(c,G);
}
*/





template<class R, class Tag>
std::ostream&
Geometry::operator<<(std::ostream& os, const Zonotope<R,Tag>& z) 
{
  typedef Numeric::Interval<R> I;
  os << "["<<z.centre();
  for(size_type j=0; j!=z.number_of_generators(); ++j) {
    os << ";";
    ::write_vector_slice(os,z.generators().column(j));
  }
  os << "]";
  return os;
}



template<class R, class Tag>
std::istream& 
Geometry::operator>>(std::istream& is, Zonotope<R,Tag>& z)
{
  Point<R> centre;
  LinearAlgebra::Matrix<R> generators;
  char c0,c1,c2;
  is >> c0 >> centre >> c1 >> generators >> c2;
  z = Zonotope<R,Tag>(centre,generators);
  return is;
}







template<class R>
void
Geometry::Zonotope<R,ExactTag>::_instantiate() 
{
  instantiate_zonotope<R>();
}





} // namespace Ariadne
                                                            
                                                            
                                                            

namespace {

using Geometry::verbosity;

template<class R, class Tag>
void
adjoin_subdivision(ListSet< Zonotope<R,Tag> >& ls, const Zonotope<R,Tag>& z) 
{
  Zonotope<R,Tag> z1,z2;
  make_lpair(z1,z2)=subdivide_pair(z);
  ls.adjoin(z1);
  ls.adjoin(z2);
}









/* Test vertices individually. Highly inefficient!! */
tribool 
superset(const Zonotope<Rational>& z, const Rectangle<Rational>& r)
{
  tribool result=true;
  for(Rectangle<Rational>::vertices_const_iterator rv_iter=r.vertices_begin(); 
      rv_iter!=r.vertices_end(); ++rv_iter) 
  {
    const Point<Rational>& pt=*rv_iter;
    result=result && Geometry::contains(z,pt);
    if(result==false) {
      break;
    }
  }
  return result;
}



/* Set up linear program to solve 
 *   \f[x=c+Ge;\ l\leq x\leq u;\ -1\leq e\leq1\f].
 *
 * Change variables to normalize \f$x\f$ and \f$e\f$
 *   \f[x'=x-l,\ e'=e+1;   x'-Ge' = c-G1-l;  0\leq x\leq u-l; \ 0\leq e\leq 2.\f] 
 * 
 * Introduce slack variables sx and se, and artificial variables ax. Problem in standard form
 *   \f[ \begin{matrix}I&0\\0&I\\\pm I&\mp G\end{matrix} \begin{matrix}x'\\e'\end{matrix}
 *        + \begin{matrix}I&&\\&I&\\&&I\end{matrix}\begin{matrix}sx\\se\\ax\end{matrix}
 *             = \begin{matrix}u-l\\2\\\pm(c-G1-l)\end{matrix} \f]
 * 
 */
tribool
disjoint(const Zonotope<Rational,ExactTag>& z, const Rectangle<Rational>& r)
{
  ARIADNE_LOG(8,"disjoint(Zonotope<Rational> q, Rectangle<Rational> r)\n");
  ARIADNE_LOG(9,"z="<<z<<", r="<<r<<"\n");
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(z,r,"tribool disjoint(Zonotope<Rational> z, Rectangle<Rational> r)");
  dimension_type d=z.dimension();
  size_type m=z.number_of_generators();
  
  // Construct tableau for testing intersection of zonotope and rectangle
  // Rectangle  l<=x<=u
  // Zonotope  x==c+Ge,  -1<=e<=1
  // 
  // Translate x'=x-l,  e'=e+1
  //   0<=x'<=u-l      ->  x' +     + sx'               == u-l
  //   0<=e'<=2        ->     +  e' +     + se'         == 2
  //   x'+l==c+G(e'-1) ->  x' + Ge'             +/- ax' == c-l-G1
  //  
  // Change sign of RHS of first equality if necessary
  // Introduce slack variables for last two inequalities
  typedef Numeric::Rational Q;
  LinearAlgebra::Matrix<Q> T(2*d+m+1,d+m+1);
  
  const Geometry::Point<Q>& l=r.lower_corner();
  const Geometry::Point<Q>& u=r.upper_corner();
  const Geometry::Point<Q>& c=z.centre();
  const LinearAlgebra::Matrix<Q>& G=z.generators();
  
  const LinearAlgebra::Vector<Q> qo(m,Q(1));
  const LinearAlgebra::Vector<Q> ql=l.position_vector();
  const LinearAlgebra::Vector<Q> qu=u.position_vector();
  const LinearAlgebra::Vector<Q> qd=qu-ql;
  const LinearAlgebra::Vector<Q> qc=c.position_vector();
  const LinearAlgebra::Matrix<Q> qG=G;
  const LinearAlgebra::Vector<Q> qrhs=qc-ql-qG*qo;
  
  if(Geometry::verbosity>8) { std::clog << "ql=" << ql << ", qd=" << qd <<", qc=" << qc << ", qrhs=" << qrhs << std::endl; }
  
  // Set up constraints x+sx=u-l
  for(size_type i=0; i!=d; ++i) {
    T(i,i)=1;
    T(i,d+m)=qu(i)-ql(i);
  }
  
  // Set up constraints e+se=2
  for(size_type j=0; j!=m; ++j) {
    T(d+j,d+j)=1;
    T(d+j,d+m)=2;
  }
  
  // Set up constraints x-Ge \pm ax=c-l-G1 
  for(size_type i=0; i!=d; ++i) {
    if(qrhs(i)>=Q(0)) {
      T(i+d+m,i)=1;
      for(size_type j=0; j!=m; ++j) {
        T(i+d+m,d+j)=-qG(i,j);
      }
      T(i+d+m,d+m)=qrhs(i);
    }
    else {
      T(i+d+m,i)=-1;
      for(size_type j=0; j!=m; ++j) {
        T(i+d+m,d+j)=qG(i,j);
      }
      T(i+d+m,d+m)=-qrhs(i);
    }
  } 
  
  // Set up cost function ax^T1
  for(size_type i=0; i!=d; ++i) {
    T(2*d+m,i) -= T(i+d+m,i);
    for(size_type j=0; j!=m; ++j) {
      T(2*d+m,d+j) -= T(i+d+m,d+j);
    }
    T(2*d+m,d+m) -= T(i+d+m,d+m);
  }
  
  LinearProgramming::LinearProgram<Q> lp(T);
  tribool result=!lp.is_feasible();
  
  return result;
}



tribool
disjoint(const Zonotope<Rational>& z1, const Zonotope<Rational>& z2)
{
  typedef Numeric::Rational Q;
  ARIADNE_CHECK_EQUAL_DIMENSIONS(z1,z2,"tribool disjoint(Zonotope<Rational> z1, Zonotope<Rational> z2)");
  
  dimension_type d=z1.dimension();
  Q one=1;
  size_type m1=z1.number_of_generators();
  size_type m2=z2.number_of_generators();
  
  LinearAlgebra::Matrix<Q> T(m1+m2+d+1,m1+m2+1);
  
  const LinearAlgebra::Vector<Q> qo1(m1,one);
  const LinearAlgebra::Vector<Q> qo2(m2,one);
  
  const Geometry::Point<Q>& qc1=z1.centre();
  const LinearAlgebra::Matrix<Q>& qG1=z1.generators();
  const Geometry::Point<Q>& qc2=z2.centre();
  const LinearAlgebra::Matrix<Q>& qG2=z2.generators();
  LinearAlgebra::Vector<Q> qrhs = qG1*qo1 - qG2*qo2 + (qc2 - qc1);
  
  // Set up constraints e1 + se1 = 2
  for(size_type j1=0; j1!=m1; ++j1) {
    T(j1,j1)=1;
    T(j1,m1+m2)=2;
  }
  
  // Set up constraints e2 + se2 = 2
  for(size_type j2=0; j2!=m2; ++j2) {
    T(m1+j2,m1+j2)=1;
    T(m1+j2,m1+m2)=2;
  }
  
  // Set up constraints G1*e1 - G2*e2 = (c2 - G2*1) - (c1 - G1*1)
  for(size_type i=0; i!=d; ++i) {
    if(qrhs(i)>=Q(0)) {
      for(size_type j1=0; j1!=m1; ++j1) {
        T(m1+m2+i,j1)=qG1(i,j1);
      }
      for(size_type j2=0; j2!=m2; ++j2) {
        T(m1+m2+i,m1+j2)=qG2(i,j2);
      }
      T(m1+m2+i,m1+m2)=qrhs(i);
    }
    else {
      for(size_type j1=0; j1!=m1; ++j1) {
        T(m1+m2+i,j1)=-qG1(i,j1);
      }
      for(size_type j2=0; j2!=m2; ++j2) {
        T(m1+m2+i,m1+j2)=-qG2(i,j2);
      }
      T(m1+m2+i,m1+m2)=-qrhs(i);
    }
  } 
  
  // Set up cost function ax^T1
  for(size_type i=0; i!=d; ++i) {
    for(size_type j1=0; j1!=m1; ++j1) {
      T(m1+m2+d,j1) -= T(m1+m2+i,j1);
    }
    for(size_type j2=0; j2!=m2; ++j2) {
      T(m1+m2+d,m1+j2) -= T(m1+m2+i,m1+j2);
    }
    T(m1+m2+d,m1+m2) -= T(m1+m2+i,m1+m2);
  }
  
  LinearProgramming::LinearProgram<Q> lp(T);
  
  tribool result=!lp.is_feasible();
  
  //std::clog << "disjoint(" << z1 << "," << z2 << ")=" << result << std::endl;
  return result;
}

/* Set up LP problem to solve \f$c+Ge=p\f$; \f$-1<=e<=1\f$.
 * Change variables so that the problem becomes \f$Ge=p-c-G1;\ 0\leq e\leq2\f$.
 * Change sign of \f$ Ge=p-c-G1\f$ to make right-hand side positive.
 */
tribool 
contains(const Zonotope<Rational>& z, const Point<Rational>& pt)
{ 
  //std::clog << "Zonotope<Rational>::contains(const Point<R>& )" << std::endl;
  typedef Rational Q;
  ARIADNE_CHECK_EQUAL_DIMENSIONS(z,pt,"tribool contains(Zonotope<Rational> z, Point<Rational> pt)");
  dimension_type d=z.dimension();
  dimension_type m=z.number_of_generators();
  
  Q zero=0;
  Q one=1;
  Q two=2;
  
  const Geometry::Point<Q>& qc=z.centre();
  const Geometry::Point<Q>& qp=pt;
  const LinearAlgebra::Matrix<Q>& qG=z.generators();
  const LinearAlgebra::Vector<Q> qo(m,one);
  const LinearAlgebra::Vector<Q> zv(m,zero);
  const LinearAlgebra::Vector<Q> tv(m,two);
  
  LinearAlgebra::Vector<Q> qrhs=qp-qc+qG*qo;
  
  LinearAlgebra::Matrix<Q> T(d+m+1,m+1);
  
  // Set up constraints e+se=2
  for(dimension_type j=0; j!=m; ++j) {
    T(j,j)=1;
    T(j,m)=2;
  }
  
  // Set up constraints Ge \pm ax = p-c+G1
  for(dimension_type i=0; i!=d; ++i) {
    if(qrhs(i)>=Q(0)) {
      for(dimension_type j=0; j!=m; ++j) {
        T(m+i,j)=qG(i,j); 
      }
      T(m+i,m)=qrhs(i);
    } else {
      for(dimension_type j=0; j!=m; ++j) {
        T(m+i,j)=-qG(i,j); 
      }
      T(m+i,m)=-qrhs(i);
    }
  }
  
  // Set up cost function ax^T 1
  for(dimension_type i=0; i!=d; ++i) {
    for(dimension_type j=0; j!=m; ++j) {
      T(m+d,j)-=T(m+i,j);
    }
    T(m+d,m)-=T(m+i,m);
  }
  
  LinearProgramming::LinearProgram<Q> lp(T);
  //std::clog << lp.tableau() << std::endl;
  tribool result=lp.is_feasible();
  //std::clog << lp.tableau() << std::endl;
  return result;
}



ListSet< Zonotope<Rational,ExactTag> > inline
subdivide(const Zonotope<Rational,ExactTag>& z)
{
  typedef Numeric::Rational Q;
  using namespace LinearAlgebra;
  
  ListSet< Zonotope<Rational,ExactTag> > result;

  size_type d=z.dimension();
  size_type m=z.number_of_generators();
  Point<Q> c=z.centre();
  
  Q max_norm=0;
  size_type max_column=0;
  for(size_type j=0; j<m; ++j) {
    Q norm = LinearAlgebra::norm(Vector<Q>(z.generators().column(j)));
    if(norm>max_norm) {
      max_norm=norm;
      max_column=j;
    }
  }
  
  Matrix<Q> new_generators=z.generators();
  size_type j=max_column;
  for(size_type i=0; i!=d; ++i) {
    new_generators(i,j)/=2;
  }
  Vector<Q> v=new_generators.column(j);
  Point<Q> new_centre=c-v;
  result.adjoin(Zonotope<Q>(new_centre,new_generators));
  new_centre=c+v;
  result.adjoin(Zonotope<Q>(new_centre,new_generators));

  return result;

}




template<class R> inline
ListSet< Zonotope<R,UniformErrorTag> >
subdivide(const Zonotope<R,UniformErrorTag>& z)
{
  typedef Numeric::Interval<R> I;
  using namespace LinearAlgebra;
  
  ListSet< Zonotope<R,UniformErrorTag> > result;
  
  size_type d=z.dimension();
  size_type m=z.number_of_generators();
  Point<I> c=z.centre();
  
  I max_radius=0;
  size_type max_value=0;
  for(size_type j=0; j<d; ++j) {
    R radius = c[j].radius();
    if(radius>max_radius) {
      max_radius=radius;
      max_value=j;
    }
  }
  
  I max_norm=0;
  size_type max_column=0;
  for(size_type j=0; j<m; ++j) {
    I norm = LinearAlgebra::norm(Vector<R>(z.generators().column(j)));
    if(norm>max_norm) {
      max_norm=norm;
      max_column=j;
    }
  }
  
  if(max_norm>max_radius) {
    LinearAlgebra::Matrix<R> new_generators=z.generators();
    size_type j=max_column;
    for(size_type i=0; i!=d; ++i) {
      new_generators(i,j)=div_up(new_generators(i,j),2);
    }
    
    Vector<R> v=new_generators.column(j);
    Point<I> new_centre=c-v;
    result.adjoin(Zonotope<R,UniformErrorTag>(new_centre,new_generators));
    new_centre=c+v;
    result.adjoin(Zonotope<R,UniformErrorTag>(new_centre,new_generators));
 } else {
    const I& cmv=z.centre()[max_value];
    Point<I> new_centre = z.centre();
    const Matrix<R>& new_generators = z.generators();
    new_centre[max_value]=I(cmv.lower(),cmv.midpoint());
    result.adjoin(Zonotope<R,UniformErrorTag>(new_centre,new_generators));
    new_centre[max_value]=I(cmv.midpoint(),cmv.upper());
    result.adjoin(Zonotope<R,UniformErrorTag>(new_centre,new_generators));
  }
  return result;
}




template<class R> inline
void
instantiate_zonotope()
{
  tribool tb;
  Geometry::Point<R> pt;
  Geometry::Rectangle<R> r;
  Geometry::Polyhedron<R> p;
  Geometry::Zonotope<R,ExactTag> z;
  Geometry::Zonotope<R,UniformErrorTag> ez;
  Geometry::Zonotope<R,IntervalTag> iz;
  std::ostream* os=0;
  std::istream* is=0;
  
  tb=Geometry::contains(z,pt);
  tb=Geometry::contains(ez,pt);
  tb=Geometry::disjoint(z,r);
  tb=Geometry::disjoint(ez,r);
  tb=Geometry::superset(z,r);
  tb=Geometry::superset(ez,r);
  tb=Geometry::subset(z,r);
  tb=Geometry::subset(ez,r);
  tb=Geometry::subset(iz,r);
  tb=Geometry::subset(z,p);
  tb=Geometry::subset(ez,p);
  tb=Geometry::subset(iz,p);
  
  //Geometry::subdivide(z);
  Geometry::subdivide(ez);
  
  Geometry::orthogonal_over_approximate(ez,iz);
  Geometry::orthogonal_over_approximate(ez,ez);

  Geometry::over_approximate(z,z);
  Geometry::over_approximate(z,ez);
  Geometry::over_approximate(z,iz);
  Geometry::over_approximate(ez,ez);
  Geometry::over_approximate(ez,iz);

  Geometry::approximate(z,iz);
  Geometry::approximate(z,ez);
  Geometry::approximate(z,z);

  Geometry::operator<<(*os,z);
  Geometry::operator<<(*os,ez);
  Geometry::operator<<(*os,iz);

  Geometry::operator>>(*is,z);
  Geometry::operator>>(*is,ez);
  Geometry::operator>>(*is,iz);

}

template<> inline
void
instantiate_zonotope<Rational>()
{
  typedef Rational R;
  tribool* tb=0;
  std::ostream* os=0;
  std::istream* is=0;
  Geometry::Point<R>* pt=0;
  Geometry::Rectangle<R>* r=0;
  Geometry::Zonotope<R>* z=0;
  Geometry::ListSet< Zonotope<R> >* zls=0;
  
  *tb=Geometry::contains(*z,*pt);
  *tb=Geometry::disjoint(*z,*r);
  *tb=Geometry::superset(*z,*r);
  *tb=Geometry::subset(*z,*r);

  *zls=Geometry::subdivide(*z);
  
  Geometry::approximate(*z,*z);

  Geometry::operator<<(*os,*z);
  Geometry::operator>>(*is,*z);

}


} // namespace 



