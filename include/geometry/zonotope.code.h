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

#include "base/tuple.h"
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
#include "function/affine_model.h"

#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/box.h"
#include "geometry/polyhedron.h"
#include "geometry/constraint_set.h"
#include "geometry/list_set.h"

#include "output/logging.h"



namespace {
  
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;

inline Rational med_approx(const Rational& ql, const Rational& qu) {
  return (ql+qu)/2;
}

inline Rational rad_up(const Rational& ql, const Rational& qu) {
  return (ql+qu)/2;
}

inline Rational add_up(const Rational& q1, const Rational& q2) {
  return q1+q2;
}

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
LinearAlgebra::Vector<R>
row_norms(const Matrix< Interval<R> >& A)
{
  size_type const& m=A.number_of_rows();
  size_type const& n=A.number_of_columns();
  Vector<R> e(m);
  for(size_type i=0; i!=m; ++i) {
    for(size_type j=0; j!=n; ++j) {
      e[i]=add_up(e[i],mag(A(i,j)));
    }
  }
  return e;
}

template<class R>
LinearAlgebra::Vector<R>
row_errors(const Matrix< Interval<R> >& A)
{
  size_type const& m=A.number_of_rows();
  size_type const& n=A.number_of_columns();
  Vector<R> e(m);
  for(size_type i=0; i!=m; ++i) {
    for(size_type j=0; j!=n; ++j) {
      e[i]=add_up(e[i],A(i,j).radius());
    }
  }
  return e;
}

template<class R> inline
LinearAlgebra::Vector<R>
errors(const Point< Interval<R> >& pt)
{  
  Vector<R> result(pt.dimension());
  for(size_type i=0; i!=pt.dimension(); ++i) {
    result[i]=pt[i].radius();
  }
  return result;
}


template<class R> inline
LinearAlgebra::Vector<R>
row_errors(const Point< Interval<R> >& pt, const Matrix< Interval<R> >& A)
{
  assert(pt.dimension()==A.number_of_rows());
  Vector<R> result(pt.dimension());
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    result[i]=pt[i].width();
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      result[i]=add_up(result[i],A(i,j).width());
    }
    result[i]=div_up(result[i],2);
  }
  return result;
}
  
template<class R> inline
Vector<R>
add_up(const Vector<R>& v1, const Vector<R>& v2) 
{
  Vector<R> result;
  for(size_type i=0; i!=v1.size(); ++i) {
    result[i]=add_up(v1[i],v2[i]);
  }
  return result;
}

template<class R> inline
tribool norm_grtr(const LinearAlgebra::Vector<R>& v1, const LinearAlgebra::Vector<R>& v2) 
{
  return LinearAlgebra::norm(v1)>LinearAlgebra::norm(v2);
}





tribool contains(const Zonotope<Rational>& z, const Point<Rational>& r);
tribool disjoint(const Zonotope<Rational>& z, const Box<Rational>& r);
tribool superset(const Zonotope<Rational>& z, const Box<Rational>& r);
tribool subset(const Zonotope<Rational>& z, const Box<Rational>& r);

template<class R> void instantiate_zonotope();
template<> void instantiate_zonotope<Numeric::Rational>();

} // namespace 






namespace Ariadne {

extern int Geometry::verbosity;


template<class R> 
Geometry::Zonotope<R>::Zonotope()
  : _centre(), _generators(), _error()
{
}

template<class R> 
Geometry::Zonotope<R>::Zonotope(dimension_type d)
  : _centre(d), _generators(d,0), _error(d)
{
}

template<class R> 
Geometry::Zonotope<R>::Zonotope(dimension_type d, size_type m)
  : _centre(d), _generators(d,m), _error(d)
{
}

template<class R> 
Geometry::Zonotope<R>::Zonotope(const Point<R>& c, const LinearAlgebra::Matrix<R>& G, const LinearAlgebra::Vector<R>& e)
  : _centre(c), _generators(G), _error(e)
{
  ARIADNE_ASSERT(c.dimension()==G.number_of_rows());
  ARIADNE_ASSERT(c.dimension()==e.size());
}

template<class R> 
Geometry::Zonotope<R>::Zonotope(const Point<R>& c, const LinearAlgebra::Matrix<R>& G)
  : _centre(c), _generators(G), _error(c.dimension())
{
  ARIADNE_ASSERT(c.dimension()==G.number_of_rows());
}

template<class R> 
Geometry::Zonotope<R>::Zonotope(const Point<I>& c, const LinearAlgebra::Matrix<R>& G)
  : _centre(midpoint(c)), _generators(G), _error(errors(c))
{
  ARIADNE_ASSERT(c.dimension()==G.number_of_rows());
}

template<class R> 
Geometry::Zonotope<R>::Zonotope(const Point<R>& c, const LinearAlgebra::Matrix<I>& G)
  : _centre(c), _generators(midpoint(G)), _error(row_errors(G))
{
  ARIADNE_ASSERT(c.dimension()==G.number_of_rows());
}

template<class R> 
Geometry::Zonotope<R>::Zonotope(const Point<I>& c, const LinearAlgebra::Matrix<I>& G)
  : _centre(midpoint(c)), _generators(midpoint(G)), _error(row_errors(c,G))
{
  ARIADNE_ASSERT(c.dimension()==G.number_of_rows());
}


template<class R>       
Geometry::Zonotope<R>::Zonotope(const Zonotope<R>& z) 
  : _centre(z._centre), _generators(z._generators), _error(z._error)
{
}

template<class R>       
Geometry::Zonotope<R>&
Geometry::Zonotope<R>::operator=(const Zonotope<R>& z) 
{ 
  if(this!=&z) {
    this->_centre=z._centre;
    this->_generators=z._generators;
    this->_error=z._error;
  }
  return *this;
}

template<class R>       
dimension_type
Geometry::Zonotope<R>::dimension() const
{
  return this->_centre.dimension();
}

template<class R>       
size_type
Geometry::Zonotope<R>::number_of_generators() const
{
  return this->_generators.number_of_columns();
}

template<class R>       
const Geometry::Point<R>&
Geometry::Zonotope<R>::centre() const
{
  return this->_centre;
}

template<class R>       
const LinearAlgebra::Matrix<R>&
Geometry::Zonotope<R>::generators() const
{
  return this->_generators;
}

template<class R>       
const LinearAlgebra::Vector<R>&
Geometry::Zonotope<R>::error() const
{
  return this->_error;
}

template<class R>       
Geometry::Box<R>
Geometry::Zonotope<R>::domain() const
{
  return Box<R>::unit_box(this->number_of_generators());
}

template<class R>       
Geometry::Box<R>
Geometry::Zonotope<R>::bounding_box() const
{
  const Zonotope<R>& z=*this;
  LinearAlgebra::Vector<I> b=z.centre().position_vector()+z.generators()*z.domain().position_vectors()+z.error()*I(-1,1);
  return Box<R>(b);
}

template<class R>       
R
Geometry::Zonotope<R>::radius() const
{
  return this->bounding_box().radius();
}

template<class R>       
tribool
Geometry::Zonotope<R>::contains(const Point<R>& pt) const
{
  return Geometry::contains(*this,pt);
}




template<class R>       
tribool
Geometry::empty(const Zonotope<R>& z) 
{
  return false;
}

template<class R>       
tribool
Geometry::bounded(const Zonotope<R>& z) 
{
  return true;
}

template<class R>       
tribool
Geometry::subset(const Box<R>& r, const Zonotope<R>& z) 
{
  return superset(z,r);
}



template<class R>       
R 
Geometry::radius(const Zonotope<R>& z) 
{
  typedef Numeric::Interval<R> I;
  return radius(z.centre()+z.generators()*z.domain().position_vectors()+z.error()*I(-1,1));
}


template<class R>       
tribool
Geometry::disjoint(const Box<R>& r, Zonotope<R>& z) 
{
  return disjoint(z,r);
}

template<class R>       
tribool
Geometry::subset(const Box<R>& r, Zonotope<R>& z) 
{
  return superset(z,r);
}





template<class R>
Geometry::Box<R>
Geometry::bounding_box(const Zonotope<R>& z)
{
  return z.bounding_box();
}


template<class R>
Geometry::ListSet< Zonotope<R> >
Geometry::split(const Zonotope<R>& z)
{
  // FIXME: Not quite guarenteed to give an over-approximation
  typedef Numeric::Interval<R> I;
  using namespace LinearAlgebra;
  
  ListSet< Zonotope<R>  > result;
  
  size_type d=z.dimension();
  size_type m=z.number_of_generators();
  Point<R> const& c=z.centre();
  Matrix<R> const& G=z.generators();
  Vector<R> const& e=z.error();
  
  array<R> norms(m,0);
  for(size_type j=0; j!=m; ++j) {
    norms[j]=norm(Vector<R>(G.column(j)));
  }

  R max_norm=0;
  size_type longest_generator=0;
  for(size_type j=0; j<m; ++j) {
    if(norms[j]>max_norm) {
      max_norm=norms[j];
      longest_generator=j;
    }
  }
  for(size_type k=0; k<d; ++k) {
    if(e[k]>max_norm) {
      max_norm=e[k];
      longest_generator=m+k;
    }
  }
  
  if(longest_generator<m) {
    LinearAlgebra::Matrix<R> new_generators=z.generators();
    size_type j=longest_generator;
    for(size_type i=0; i!=d; ++i) {
      new_generators(i,j)=div_up(new_generators(i,j),2);
    }
    
    Vector<R> v=new_generators.column(j);
    Point<R> new_centre=sub_approx(c,v);
    result.adjoin(Zonotope<R>(new_centre,new_generators,e));
    new_centre=add_approx(c,v);
    result.adjoin(Zonotope<R>(new_centre,new_generators,e));
 } else {
    dimension_type k=longest_generator-m;
    Point<R> new_centre = z.centre();
    const Matrix<R>& new_generators = z.generators();
    Vector<R> new_error=e;
    new_error[k]=div_up(new_error[k],2);
    new_centre[k]=add_approx(z.centre()[k],new_error[k]);
    result.adjoin(Zonotope<R>(new_centre,new_generators,new_error));
    new_centre[k]=sub_approx(z.centre()[k],new_error[k]);
    result.adjoin(Zonotope<R>(new_centre,new_generators,new_error));
  }
  return result;
} 




template<class R>       
Geometry::Zonotope<R>::Zonotope(const Box<R>& r) 
  : _centre(r.dimension()), _generators(r.dimension(),r.dimension()), _error(r.dimension())
{
  dimension_type d=r.dimension();
  Point<R>& c=this->_centre;
  LinearAlgebra::Matrix<R>& G=this->_generators;
  LinearAlgebra::Vector<R>& e=this->_error;
  for(size_type i=0; i!=d; ++i) {
    c[i]=med_approx(r.lower_bound(i),r.upper_bound(i));
    for(size_type j=0; j!=d; ++j) {
      G(i,j)=0;
    }
    G(i,i)=rad_up(r.lower_bound(i),r.upper_bound(i));
    e(i)=0;
  }
}


template<class R>
Zonotope<R>
Geometry::apply(const Function::AffineModel<R>& am,
                const Zonotope<R>& z)
{
  using namespace LinearAlgebra;
  typedef Interval<R> I;
  
  ARIADNE_ASSERT(z.centre()==am.centre());
  if(!subset(z,am.domain())) {
    std::cerr<<"z="<<z<<"\nz.bounding_box()="<<z.bounding_box()<<"\nam.domain()="<<am.domain()<<std::endl;
  }
  ARIADNE_ASSERT(possibly(subset(z,am.domain())));
  
  dimension_type d=z.dimension();
  size_type m=z.number_of_generators();
  dimension_type nd=am.result_size();

  //Point<R> const& c=z.centre();
  Matrix<R> const& G=z.generators();
  Vector<R> const& e=z.error();

  Point<I> const& nic=am.value();
  Matrix<I> const& iDf=am.jacobian();
  Matrix<I> niG=iDf*G;
  
  Point<R> nc=midpoint(nic);
  Matrix<R> nG=midpoint(niG);
  Vector<R> ne(nd);

  for(size_type i=0; i!=nd; ++i) {
    R& err=ne[i];
    err=add_up(err,nic[i].radius());
    for(size_type j=0; j!=m; ++j) {
      err=add_up(err,niG(i,j).radius());
    }
    for(size_type k=0; k!=d; ++k) {
      err=add_up(err,mul_up(mag(iDf(i,k)),e[k]));
    }
  }
  
  return Zonotope<R>(nc,nG,ne);

}



template<class R>
tribool 
Geometry::contains(const Zonotope<R>& z, const Point<R>& pt)  
{
  return ::contains(Zonotope<Rational>(z),Point<Rational>(pt));
}


template<class R>
tribool
Geometry::disjoint(const Zonotope<R>& z, const Box<R>& r)
{
  return ::disjoint(Zonotope<Rational>(z),Box<Rational>(r));
}

template<class R>
tribool
Geometry::intersects(const Zonotope<R>& z, const Box<R>& r)
{
  return not ::disjoint(Zonotope<Rational>(z),Box<Rational>(r));
}

template<class R>
tribool
Geometry::superset(const Zonotope<R>& z, const Box<R>& r)
{
  return ::superset(Zonotope<Rational>(z),Box<Rational>(r));
}

template<class R>
tribool
Geometry::subset(const Zonotope<R>& z, const Box<R>& r)
{
  return Geometry::subset(bounding_box(z),r);
}

template<class R>
tribool
Geometry::subset(const Zonotope<R>& z, const Polyhedron<R>& p)
{
  typedef Numeric::Interval<R> I;
  if(z.error()==0) {
    LinearAlgebra::Vector<I> im=(p.A()*z.generators())*z.domain().position_vectors()+(p.A()*z.centre().position_vector()+p.b());
    return im>=0;
  } else {
    return subset(error_free_over_approximation(z),p);
  }
}





template<class R> 
Geometry::Zonotope<R> 
Geometry::approximation(const Zonotope<R>& z)
{
  return Zonotope<R>(z.centre(),z.generators());
}

template<class R> 
Zonotope<R>
Geometry::over_approximation(const Zonotope<R>& z)
{
  return z;
}

template<class R> 
Zonotope<R>
Geometry::error_free_over_approximation(const Zonotope<R>& z)
{
  if(z.error()==0) {
    return z;
  }
  dimension_type d=z.dimension();
  size_type m=z.number_of_generators();
  Matrix<R> nG(d,m+d);
  nG(slice(0,d),slice(0,m))=z.generators();
  for(size_type i=0; i!=d; ++i) {
    nG(i,i+m)=z.error()[i];
  }
  return Zonotope<R>(z.centre(),nG);
}

template<class R> 
Zonotope<R>
Geometry::nonsingular_over_approximation(const Zonotope<R>& z)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}  



template<class R> 
Geometry::Zonotope<R>
Geometry::orthogonal_over_approximation(const Zonotope<R>& z)
{
  //assert(iz.dimension()==iz.number_of_generators());
  typedef Numeric::Interval<R> I;
  Zonotope<R> ez=error_free_over_approximation(z);

  const Point<R>& c=ez.centre();
  const Matrix<R>& G=ez.generators();
  
  Matrix<R> aQ,aR;
  make_lpair(aQ,aR)=qr_approx(G);

  Matrix<I> aQinv=inverse(aQ);
  Matrix<I> iR=aQinv*G;
  DiagonalMatrix<R> aD(row_norms(iR));

  Matrix<I> niG=aQ*aD;

  return Zonotope<R>(c,niG);
}  

template<class R>
Geometry::Zonotope<R>
Geometry::cascade_over_approximation(const Zonotope<R>& z, size_type cs)
{
  using namespace std;
  using namespace LinearAlgebra;
  if(z.number_of_generators()<=z.dimension()*cs) { return z; }  

  assert(z.number_of_generators()%z.dimension()==0);

  dimension_type d=z.dimension();
  size_type nb=z.number_of_generators()/z.dimension(); // number of generator blocks
   

  const Matrix<R>& G=z.generators();
  array<R> norms(nb);
  for(size_type i=0; i!=nb; ++i) {
    norms[i]=norm(G(slice(0,d),slice(i*d,d))).upper();
  }
  
  // Compute the new number of blocks
  size_type nnb=cs;
  R sum=0;
  for(size_type i=nb-1; i!=0; --i) {
    sum=add_approx(sum,norms[i]);
    if(sum>norms[i-1]) {
      nnb=i;
    }
  }
  nnb=min(nnb,cs);
  // Reduce generators
  Matrix<R> rG(d,d*nnb);
  rG(slice(0,d),slice(0,d*(nnb-1)))=G(slice(0,d),slice(0,d*(nnb-1)));
  for(size_type i=0; i!=d; ++i) {
    R& err=rG(i,d*(nnb-1)+i);
    for(size_type j=d*(nnb-1); j!=G.number_of_columns(); ++j) {
      err=add_up(err,abs(G(i,j)));
    }
  }
  return Zonotope<R>(z.centre(),rG);
}


template<class R>
tribool
Geometry::disjoint(const Zonotope<R>& z, const ConstraintSet<R>& cs)
{
  //TODO: Change disjoint(Zonotope,Box) to accept unbounded boxes.
  ARIADNE_CHECK_EQUAL_DIMENSIONS(z,cs,"disjoint(Zonotope,ConstraintSet)");
  Zonotope<R> fz=apply(Function::AffineModel<R>(z.bounding_box(),z.centre(),cs.function()),z);
  Box<R> bcd=closed_intersection(cs.codomain(),fz.bounding_box());
  return disjoint(fz,bcd);
}

template<class R>
tribool
Geometry::subset(const Zonotope<R>& z, const ConstraintSet<R>& cs)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(z,cs,"disjoint(Zonotope,ConstraintSet)");
  Zonotope<R> fz=apply(Function::AffineModel<R>(z.bounding_box(),z.centre(),cs.function()),z);
  return subset(fz.bounding_box(),cs.codomain());
}

template<class R>
tribool
Geometry::intersects(const Zonotope<R>& z, const ConstraintSet<R>& cs)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(z,cs,"disjoint(Zonotope,ConstraintSet)");
  Zonotope<R> fz=apply(Function::AffineModel<R>(z.bounding_box(),z.centre(),cs.function()),z);
  Box<R> bcd=closed_intersection(cs.codomain(),fz.bounding_box());
  return intersects(fz,bcd);
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





template<class R>
std::ostream&
Geometry::operator<<(std::ostream& os, const Zonotope<R>& z) 
{
  typedef Numeric::Interval<R> I;
  os << "["<<z.centre();
  for(size_type j=0; j!=z.number_of_generators(); ++j) {
    os << ";";
    os << z.generators().column(j);
  }
  os << "]";
  return os;
}



template<class R>
std::istream& 
Geometry::operator>>(std::istream& is, Zonotope<R>& z)
{
  Point<R> centre;
  LinearAlgebra::Matrix<R> generators;
  char c0,c1,c2;
  is >> c0 >> centre >> c1 >> generators >> c2;
  z = Zonotope<R>(centre,generators);
  return is;
}







template<class R>
void
Geometry::Zonotope<R>::_instantiate() 
{
  ::instantiate_zonotope<R>();
}





} // namespace Ariadne
                                                            
                                                            
                                                            

namespace {

using Geometry::verbosity;

template<class R>
void
adjoin_subdivision(ListSet< Zonotope<R> >& ls, const Zonotope<R>& z) 
{
  ls.adjoin(Geometry::subdivide(z));
}









/* Test vertices individually. Highly inefficient!! */
tribool 
superset(const Zonotope<Rational>& z, const Box<Rational>& r)
{
  tribool result=true;
  for(Box<Rational>::vertices_const_iterator rv_iter=r.vertices_begin(); 
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
disjoint(const Zonotope<Rational>& z, const Box<Rational>& r)
{
  ARIADNE_LOG(8,"disjoint(Zonotope<Rational> q, Box<Rational> r)\n");
  ARIADNE_LOG(9,"z="<<z<<", r="<<r<<"\n");
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(z,r,"tribool disjoint(Zonotope<Rational> z, Box<Rational> r)");
  dimension_type d=z.dimension();
  size_type m=z.number_of_generators();
  
  // Construct tableau for testing intersection of zonotope and rectangle
  // Box  l<=x<=u
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
  
  const Geometry::Point<Q> l=r.lower_corner()-z.error();
  const Geometry::Point<Q> u=r.upper_corner()+z.error();
  const Geometry::Point<Q>& c=z.centre();
  const LinearAlgebra::Matrix<Q>& G=z.generators();
  //const LinearAlgebra::Vector<Q>& e=z.error();
  
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










template<class R> 
void
instantiate_zonotope()
{
  tribool* tb=0;
  R* r=0;
  Geometry::Point<R>* pt=0;
  Geometry::Box<R>* bx=0;
  Geometry::Polyhedron<R>* p=0;
  Geometry::Zonotope<R>* z=0;
  Geometry::ConstraintSet<R>* cs=0;
  Function::AffineModel<R>* am=0;
  std::ostream* os=0;
  std::istream* is=0;

  *r=Geometry::radius(*z);
  *bx=Geometry::bounding_box(*z);
  
  *tb=Geometry::contains(*z,*pt);
  *tb=Geometry::disjoint(*z,*bx);
  *tb=Geometry::superset(*z,*bx);
  *tb=Geometry::subset(*bx,*z);
  *tb=Geometry::subset(*z,*bx);
  *tb=Geometry::subset(*z,*p);
  
  *tb=Geometry::disjoint(*z,*cs);
  *tb=Geometry::intersects(*z,*cs);
  *tb=Geometry::subset(*z,*cs);

  Geometry::split(*z);

  Geometry::approximation(*z);
  Geometry::over_approximation(*z);
  Geometry::error_free_over_approximation(*z);
  Geometry::orthogonal_over_approximation(*z);
  Geometry::nonsingular_over_approximation(*z);
  Geometry::cascade_over_approximation(*z,1);

  Geometry::apply(*am,*z);

  Geometry::operator<<(*os,*z);
  Geometry::operator>>(*is,*z);
}

template<> 
void
instantiate_zonotope<Rational>()
{
  typedef Rational R;
  tribool* tb=0;
  R* r=0;
  std::ostream* os=0;
  std::istream* is=0;
  Geometry::Point<R>* pt=0;
  Geometry::Box<R>* bx=0;
  Geometry::Zonotope<R>* z=0;
  //Geometry::ListSet< Zonotope<R> >* zls=0;
  
  *r=Geometry::radius(*z);
  *bx=Geometry::bounding_box(*z);
  
  *tb=Geometry::contains(*z,*pt);
  *tb=Geometry::contains(*z,*pt);
  *tb=Geometry::disjoint(*z,*bx);
  *tb=Geometry::superset(*z,*bx);
  *tb=Geometry::subset(*bx,*z);
  *tb=Geometry::subset(*z,*bx);

  *z=Geometry::approximation(*z);

  Geometry::operator<<(*os,*z);
  Geometry::operator>>(*is,*z);

}


} // namespace 



