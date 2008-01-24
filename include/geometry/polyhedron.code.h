/***************************************************************************
 *            polyhedron.code.h
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

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "polyhedron.h"

#include "macros/assert.h"

#include "base/sequence_io.h"

#include "numeric/interval.h"
#include "numeric/arithmetic.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_programming/linear_program.h"

#include "geometry/ddconv.h"
#include "geometry/ddconv.code.h"
#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/box.h"
#include "geometry/rectangle.h"
#include "geometry/polytope.h"

#include "output/logging.h"


// Specializations of polyhedral operations
namespace {

using namespace Ariadne;
using Geometry::verbosity;

inline
void
assign(Geometry::Polyhedron<Numeric::Rational>& plhd, const Geometry::PointList<Numeric::Rational>& vl)
{
  plhd=Geometry::Polyhedron<Numeric::Rational>(Geometry::Polytope<Numeric::Rational>(vl));
}

template<class X> inline
void
assign(Geometry::Polyhedron<X>& plhd, const Geometry::PointList<X>& pltp)
{
  throw NotImplemented(__PRETTY_FUNCTION__); 
}



inline
tribool 
empty(const Geometry::Polyhedron<Numeric::Rational>& plhd) 
{
  try {
    Geometry::Polytope<Numeric::Rational> pltp(plhd);
    ARIADNE_LOG(7,"empty(plhd): plhd="<<plhd<<" pltp="<<pltp);
    return pltp.empty();
  }
  catch(Geometry::UnboundedSet& e) {
    return false;
  }
}

template<class R> inline
tribool 
empty(const Geometry::Polyhedron<R>& plhd) 
{
  return empty(Geometry::Polyhedron<Numeric::Rational>(plhd));
}

template<class R> inline
tribool 
empty(const Geometry::Polyhedron< Numeric::Interval<R> >& plhd) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



inline
tribool 
bounded(const Geometry::Polyhedron<Numeric::Rational>& plhd)
{
  try {
    Geometry::Polytope<Numeric::Rational> pltp(plhd);
    return pltp.bounded();
  }
  catch(Geometry::UnboundedSet& e) {
    return false;
  }
}

template<class R> inline
tribool 
bounded(const Geometry::Polyhedron<R>& plhd) 
{
  return bounded(Geometry::Polyhedron<Numeric::Rational>(plhd));
}


template<class R> inline
tribool 
bounded(const Geometry::Polyhedron< Numeric::Interval<R> >& plhd) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


tribool 
disjoint(const Geometry::Polyhedron<Numeric::Rational>& plhd, const Geometry::Box<Numeric::Rational>& r)
{
  typedef Numeric::Rational R;
  typedef Numeric::Rational F;
  
  if(verbosity>7) { std::clog << "tribool disjoint(Polyhedron<Rational>,Polyhedron<Rational>)" << std::endl; }
  if(verbosity>8) { std::clog << plhd << " " << r << std::endl; }
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(plhd,r,"tribool disjoint(Polyhedron<Rational> plhd, Box<Rational> r)");
  dimension_type d=plhd.dimension();
  size_type nc=plhd.number_of_constraints();
  size_type nnc=0; // number of negative constraints Ax<=b
  
  // Translate x'=x-l;  Ax<=b A(x'+l)<=b Ax'<= b-Al b'=b-Al 
  LinearAlgebra::Vector<F> u=r.upper_corner()-r.lower_corner();
  LinearAlgebra::Matrix<R> A=plhd.A();
  LinearAlgebra::Vector<F> b=plhd.b()-A*r.lower_corner().position_vector();
  
  if(verbosity>8) { std::clog << "u'=" << u << " A'=" << A << " b'=" << b << std::endl; }
  
  // count number of negative constraints
  for(size_type i=0; i!=nc; ++i) {
    if(b(i)<0) {
      ++nnc;
    }
  }
  
  LinearAlgebra::Matrix<F> T(nc+d+1,d+nnc+1);
  // set up constraints Ax<=b
  size_type k=0; // negative constraint number
  for(size_type i=0; i!=nc; ++i) {
    if(b(i)<0) {
      for(size_type j=0; j!=d; ++j) {
        T(i,j)=-A(i,j);
      }
      T(i,d+k)=-1;
      T(i,d+nnc)=-b(i);
      ++k;
    } else {
      for(size_type j=0; j!=d; ++j) {
        T(i,j)=A(i,j);
      }
      T(i,d+nnc)=b(i);
    }          
  }
  //
  for(size_type i=0; i!=d; ++i) {
    T(nc+i,i)=1;
    T(nc+i,d+nnc)=u(i);
  }
  
  // SetInterface value function for feasibility problem
  k=0;
  for(size_type i=0; i!=nc; ++i) {
    if(b(i)<0) {
      for(size_type j=0; j!=d+nnc; ++j) {
        T(nc+d,j)-=T(i,j);
      }
      T(nc+d,d+nnc)-=T(i,d+nnc);
    }
  }
  
  if(verbosity>8) { std::clog << T << std::endl; }
  
  LinearProgramming::LinearProgram<F> lp(T);
  tribool result=!lp.is_feasible();
  return result;
}

template<class T>
tribool 
disjoint(const Geometry::Polyhedron< Numeric::Float<T> >& plhd, 
         const Geometry::Box< Numeric::Float<T> >& r)
{
  return ::disjoint(Geometry::Polyhedron<Numeric::Rational>(plhd),
                  Geometry::Box<Numeric::Rational>(r));
}
                  
template<class R>
tribool 
disjoint(const Geometry::Polyhedron< Numeric::Interval<R> >& plhd, const Geometry::Box<R>& bx)
{
  throw NotImplemented("disjoint(Polyhedron<Fuzzy>,Box<Fuzzy>");
}
                  

tribool 
disjoint(const Geometry::Polyhedron<Numeric::Rational>& plhd, const Geometry::Polytope<Numeric::Rational>& pltp)
{
  return disjoint(plhd,Geometry::Polyhedron<Numeric::Rational>(pltp));
}

template<class T>
tribool 
disjoint(const Geometry::Polyhedron< Numeric::Float<T> >& plhd, 
         const Geometry::Polytope< Numeric::Float<T> >& pltp)
{
  return ::disjoint(Geometry::Polyhedron<Numeric::Rational>(plhd),
                    Geometry::Polytope<Numeric::Rational>(pltp));
}
                  
template<class R>
tribool 
disjoint(const Geometry::Polyhedron< Numeric::Interval<R> >& plhd, const Geometry::Polytope< Numeric::Interval<R> >& pltp)
{
  throw NotImplemented("disjoint(Polyhedron<Fuzzy>,Polytope<Fuzzy>");
}
                  




template<class R, class BS> inline
tribool 
subset(const BS& A, const Geometry::Polyhedron<R>& B)
{
  tribool result=true;
  for(typename BS::vertices_const_iterator v=A.vertices_begin();
      v!=A.vertices_end(); ++v)
  {
    for(typename Geometry::Polyhedron<R>::constraints_const_iterator c=B.constraints_begin();
        c!=B.constraints_end(); ++c)
    {
      result=result && (satisfies(*v,*c));
      if(!result) { return result; }
    }
  }
  return result;
}



inline
tribool 
subset(const Geometry::Polyhedron<Numeric::Rational>& plhd, const Geometry::Box<Numeric::Rational>& r)
{
  ARIADNE_LOG(3,"Geometry::subset(plhd,r): ""plhd="<<plhd<<", r="<<r<<"\n");
  typedef Numeric::Rational Q;
  tribool result=true;
  dimension_type d=plhd.dimension();
  array<Q> data(d+1u);
  for(size_type i=0; i!=d+1u; ++i) {
    data[i]=0;
  }
  Geometry::Polyhedron<Q> halfspace;
  ARIADNE_LOG(9,"data="<<data<<"\n");
  
  for(size_type j=0; j!=d; ++j) {
    // Check if disjoint from lower halfspace
    data[j]=-1;
    data[d]=r.lower_bound(j);
    ARIADNE_LOG(9,"data="<<data<<"\n");
    halfspace=Geometry::Polyhedron<Q>(d,1u,data.begin());
    ARIADNE_LOG(8,"halfspace="<<halfspace<<"\n");
    result=result && disjoint(plhd,halfspace);
    // Check if disjoint from upper halfspace
    data[j]=1;
    data[d]=-r.upper_bound(j);
    halfspace=Geometry::Polyhedron<Q>(d,1u,data.begin());
    ARIADNE_LOG(8,"halfspace="<<halfspace<<"\n");
    result=result && disjoint(plhd,halfspace);
    // Early return
    if(result==false) {
      return result;
    }
    data[j]=0;
  }
  return result;
}

template<class T>
inline
tribool 
subset(const Geometry::Polyhedron< Numeric::Float<T> >& plhd, 
       const Geometry::Box< Numeric::Float<T> >& r)
{
  return subset(Geometry::Polyhedron<Numeric::Rational>(plhd),Geometry::Box<Numeric::Rational>(r));
}

template<class R>
inline
tribool 
subset(const Geometry::Polyhedron< Numeric::Interval<R> >& plhd, const Geometry::Box<R>& r)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


inline
tribool 
subset(const Geometry::Polyhedron<Numeric::Rational>& plhd1, const Geometry::Polyhedron<Numeric::Rational>& plhd2)
{
  return Geometry::subset(Geometry::Polytope<Numeric::Rational>(plhd1),plhd2);
}

template<class R>
inline
tribool 
subset(const Geometry::Polyhedron<R>& plhd1, const Geometry::Polyhedron<R>& plhd2)
{
  return ::subset(Geometry::Polyhedron<Numeric::Rational>(plhd1),Geometry::Polyhedron<Numeric::Rational>(plhd2));
}

template<class R>
inline
tribool 
subset(const Geometry::Polyhedron< Numeric::Interval<R> >& plhd1, const Geometry::Polyhedron< Numeric::Interval<R> >& plhd2)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}






inline
Geometry::Box<Numeric::Rational> 
bounding_box(const Geometry::Polyhedron<Numeric::Rational>& plhd)
{
  return Geometry::Polytope<Numeric::Rational>(plhd).bounding_box();
}

template<class R> inline
Geometry::Box<R> 
bounding_box(const Geometry::Polyhedron<R>& plhd)
{
  using Numeric::Rational;
  Geometry::Box<Numeric::Rational> qbb=bounding_box(Geometry::Polyhedron<Numeric::Rational>(plhd));
  Geometry::Box<R> bb(plhd.dimension());
  for(dimension_type i=0; i!=bb.dimension(); ++i) {
    bb.set_lower_bound(i,R(qbb.lower_bound(i),Numeric::round_down));
    bb.set_upper_bound(i,R(qbb.upper_bound(i),Numeric::round_up));
  }
  return bb;
}

template<class R> inline
Geometry::Box<R> 
bounding_box(const Geometry::Polyhedron< Numeric::Interval<R> >& plhd)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}




} // namespace

  

namespace Ariadne {

template<class X> 
tribool 
Geometry::empty(const Polyhedron<X>& plhd) 
{
  return ::empty(plhd);
}

template<class X> 
tribool 
Geometry::bounded(const Polyhedron<X>& plhd) 
{
  return ::bounded(plhd);
}







template<class X> 
Geometry::Polyhedron<X>::Polyhedron(const Box<R>& bx)
  : _dimension(bx.dimension()), 
    _number_of_constraints(bx.dimension()*2u),
    _data((bx.dimension()+1u)*bx.dimension()*2u,static_cast<X>(0))
{
  dimension_type d=bx.dimension();
  LinearAlgebra::MatrixSlice<X> constraints=this->_constraints();
  for(size_type i=0; i!=d; ++i) {
    constraints(i,i)=static_cast<X>(1);
    constraints(i,d)=-bx.lower_bound(i);
    constraints(i+d,i)=static_cast<X>(-1); 
    constraints(i+d,d)=bx.upper_bound(i);
  }
}





template<class X>
Geometry::Polyhedron<X>::Polyhedron(const std::string& str)
{
  std::stringstream ss(str);
  ss >> *this;
}


template<class X>
Geometry::Polyhedron<X>::Polyhedron(dimension_type d)
  : _dimension(d), 
    _number_of_constraints(0), 
    _data()
{
}


template<class X>
Geometry::Polyhedron<X>::Polyhedron(dimension_type d, size_type nc, const X* data)
  : _dimension(d), 
    _number_of_constraints(nc), 
    _data(data,data+(d+1)*nc)
{ 
}


template<class X>
Geometry::Polyhedron<X>::Polyhedron(const LinearAlgebra::Matrix<X>& A,
                                    const LinearAlgebra::Vector<X>& b) 
  : _dimension(A.number_of_columns()), 
    _number_of_constraints(A.number_of_rows()), 
    _data((A.number_of_columns()+1)*A.number_of_rows())
{
  using namespace LinearAlgebra;
  ARIADNE_CHECK_SIZE(b,A.number_of_rows(),"Polyhedron::Polyhedron(Matrix A, Vector b)");
  dimension_type d=this->dimension();
  dimension_type nc=this->number_of_constraints();
  LinearAlgebra::MatrixSlice<X>(nc,d,this->begin(),d+1u,1u)=-A;
  LinearAlgebra::VectorSlice<X>(nc,this->begin()+d,d+1u)=b;
}


template<class X>
Geometry::Polyhedron<X>::Polyhedron(const LinearAlgebra::Matrix<X>& C) 
  : _dimension(C.number_of_columns()-1), 
    _number_of_constraints(C.number_of_rows()), 
    _data(C.data())
{
}


template<class X>
Geometry::Polyhedron<X>::Polyhedron(const PointList<X>& pts)
  : _dimension(pts.dimension()), _number_of_constraints(0), _data()
{
  assign(*this,pts);
}



template<class X>
void
Geometry::Polyhedron<X>::new_constraint(const Halfspace<X>& c)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,c,"void Polyhedron::new_constraint(PolyhedralConstraint& c)");
  dimension_type d=this->_dimension;
  size_type sz=this->_data.size();
  this->_data.resize(sz+d+1u);
  for(dimension_type i=0; i!=d; ++i) {
    this->_data[sz+i]=c.data()[i];
    this->_data[sz+d]=c.data()[d];
  }
}


template<class X>
Geometry::Halfspace<X>
Geometry::Polyhedron<X>::constraint(size_type i) const
{
  ARIADNE_ASSERT(i<this->_number_of_constraints);
  return Halfspace<X>(this->_dimension,this->_data.begin()+i*(this->_dimension+1u));
}


template<class X>
LinearAlgebra::MatrixSlice<X>
Geometry::Polyhedron<X>::_constraints()
{
  return LinearAlgebra::MatrixSlice<X>(this->_number_of_constraints,
                                       this->_dimension+1u,
                                       this->begin(),
                                       this->_dimension+1u,
                                       1u);
}




template<class X>
typename Geometry::Polyhedron<X>::constraints_const_iterator
Geometry::Polyhedron<X>::constraints_begin() const
{
  return constraints_const_iterator(*this,0u);
}

template<class X>
typename Geometry::Polyhedron<X>::constraints_const_iterator
Geometry::Polyhedron<X>::constraints_end() const
{
  return constraints_const_iterator(*this,this->number_of_constraints());
}

template<class X>
dimension_type
Geometry::Polyhedron<X>::dimension() const
{
  return this->_dimension;
}



template<class X> 
Geometry::Box<typename Geometry::Polyhedron<X>::real_type> 
Geometry::Polyhedron<X>::bounding_box() const
{
  return ::bounding_box(*this);
}

template<class X>
tribool 
Geometry::Polyhedron<X>::empty() const
{
  return ::empty(*this);
}

template<class X>
tribool 
Geometry::Polyhedron<X>::bounded() const
{
  return ::bounded(*this);
}




template<class X, class R>
tribool 
Geometry::disjoint(const Polyhedron<X>& plhd, const Box<R>& bx)
{
  return ::disjoint(plhd,bx);
}


template<class X, class R>
tribool 
Geometry::disjoint(const Box<R>& bx, const Polyhedron<X>& plhd)
{
  return ::disjoint(plhd,bx);
}



template<class X>
tribool 
Geometry::disjoint(const Polyhedron<X>& plhd, const Polytope<X>& pltp)
{
  return ::disjoint(plhd,pltp);
}


template<class X>
tribool 
Geometry::disjoint(const Polytope<X>& pltp, const Polyhedron<X>& plhd)
{
  return ::disjoint(plhd,pltp);
}



template<class X>
tribool 
Geometry::disjoint(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2)
{
  return closed_intersection(plhd1,plhd2).empty();
}

template<class X>
tribool 
Geometry::subset(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2)
{
  return ::subset(plhd1,plhd2);
}

template<class X,class R>
tribool 
Geometry::subset(const Polyhedron<X>& plhd, const Box<R>& r)
{
  return ::subset(plhd,r);
}


template<class X,class R>
tribool 
Geometry::subset(const Box<R>& r, const Polyhedron<X>& plhd)
{
  return ::subset(r,plhd);
}


template<class X1,class X2>
tribool 
Geometry::subset(const Polytope<X1>& pltp, const Polyhedron<X2>& plhd)
{
  return ::subset(pltp,plhd);
}



template<class X>
tribool 
Geometry::equal(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2)
{
  return Geometry::subset(plhd1,plhd2) && Geometry::subset(plhd2,plhd1); 
}


template<class X>
Geometry::Polyhedron<X> 
Geometry::open_intersection(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(plhd1,plhd2,"Polyhedron open_intersection(Polyhedron plhd1, Polyhedron plhd2)");
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class X>
Geometry::Polyhedron<X> 
Geometry::closed_intersection(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(plhd1,plhd2,"Polyhedron closed_intersection(Polyhedron plhd1, Polyhedron plhd2)");
  dimension_type d=plhd1.dimension();
  size_type nc1=plhd1.number_of_constraints();
  size_type nc2=plhd2.number_of_constraints();
  LinearAlgebra::Matrix<X> A(nc1+nc2,d);
  LinearAlgebra::MatrixSlice<X>(nc1,d,A.begin(),d,1)=plhd1.A();
  LinearAlgebra::MatrixSlice<X>(nc2,d,A.begin()+nc1*d,d,1)=plhd2.A();
  LinearAlgebra::Vector<X> b=direct_sum(plhd1.b(),plhd2.b());
  return Polyhedron<X>(A,b);
}


template<class X>
Geometry::Polyhedron<X> 
Geometry::closed_intersection(const Rectangle<X>& r, const Polyhedron<X>& plhd)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r,plhd,"Polyhedron closed_intersection(Rectangle r, Polyhedron plhd)");
  return closed_intersection(Polyhedron<X>(r),plhd);
}


template<class X>
Geometry::Polyhedron<X> 
Geometry::closed_intersection(const Polyhedron<X>& plhd, const Rectangle<X>& r)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(plhd,r,"Polyhedron closed_intersection(Polyhedron plhd, Rectangle r)");
  return closed_intersection(plhd,Polyhedron<X>(r));
}


template<class X>
Geometry::Polyhedron<X> 
Geometry::polyhedron(const Rectangle<X>& r)
{
  return Polyhedron<X>(r);
}


template<class X>
Geometry::Polyhedron<typename Numeric::traits<X>::arithmetic_type>
Geometry::polyhedron(const Polytope<X>& pltp)
{
  typedef typename Numeric::traits<X>::arithmetic_type F;
  
  dimension_type d=pltp.dimension();
  size_type nv=pltp.number_of_vertices();
  
  if(nv==0) {
    // empty polytope; return empty polyhedron
    array<F> data(d+1,F(0));
    data[d]=-1;
    return Polyhedron<F>(d,1,data.begin());
  }
  
  const LinearAlgebra::Matrix<X> G=pltp.generators();
  
  std::vector< LinearAlgebra::Vector<F> > result;
  std::vector< LinearAlgebra::Vector<F> > argument;
  
  LinearAlgebra::Vector<F> tmp(d+1);
  
  for(size_type j=0; j!=nv; ++j) {
    for(size_type i=0; i!=d+1u; ++i) {
      tmp(i)=G(i,j);
    }
    argument.push_back(tmp);
  }
  
  ddconv(result,argument);     
  
  size_type nc=result.size();
  array<F> data((d+1)*nc);
  for(size_type i=0; i!=nc; ++i) {
    for(size_type j=0; j!=d+1u; ++j) {
      data[i*(d+1u)+j]=result[i](j);
    }
  }
  
  return Polyhedron<F>(d,nc,data.begin());
}



template<class X>
std::string
Geometry::Polyhedron<X>::name()
{
  return std::string("Polyhedron")+"<"+Numeric::name<X>()+">";
}

template<class X>  
std::ostream& 
Geometry::Polyhedron<X>::write(std::ostream& os) const
{
  //return os << "Polyhedron( A=" << this->A() << ", b=" << this->b() << " )";
  os << "Polyhedron( constraints=";
  dimension_type d=this->dimension();
  size_type nc=this->number_of_constraints();
  for(size_type i=0; i!=nc; ++i) {
    os << ( i==0 ? "[" : "," );
    for(size_type j=0; j!=d; ++j) {
      os << ( j==0 ? "(" : ",");
      os << this->_data[i*(d+1)+j]; 
    }
    os << ":" << this->_data[i*(d+1)+d] << ")";
  }
  os << "] )";
  return os;
}

template<class X>  
std::istream& 
Geometry::Polyhedron<X>::read(std::istream& is) 
{
  std::vector< std::vector<X> > Alst;
  std::vector< X > Blst;
  
  std::vector<X> a;
  R b;
  
  char c;
  is >> c;
  assert(c=='[');
  
  c=is.peek();
  while(c=='[') {
    // Read constraint ax<=b in form [a_1,a_2,...,a_n;b];
    Base::read_sequence(is,a,'[',';',',');
    is >> b;
    is >> c;
    assert(c==']');
    Alst.push_back(a);
    Blst.push_back(b);
  }
  
  size_type m=Alst.size();
  size_type n=Alst[0].size();
  LinearAlgebra::Matrix<X> A(m,n);
  LinearAlgebra::Vector<X> B(m);
  for(uint i=0; i!=m; ++i) {
    for(size_type j=0; j!=n; ++j) {
      A(i,j)=Alst[i][j];
    }
    B(i)=Blst[i];
  }
  
  *this=Polyhedron<X>(A,B);
  
  return is;
}


template<class X>
void
Geometry::Polyhedron<X>::_instantiate() 
{
  typedef typename Numeric::traits<X>::number_type R;
  typedef typename Numeric::traits<X>::arithmetic_type F;
  tribool tb;
  Box<R> bx;
  Rectangle<X> r;
  Polytope<X> c;
  Polyhedron<X> p;
  Polyhedron<F> ip;
  
  p=Polyhedron<X>(r);
  p=Polyhedron<X>(p);
  ip=Polyhedron<F>(c);
  
  tb=equal(p,p);
  tb=disjoint(bx,p);
  tb=disjoint(p,bx);
  tb=disjoint(c,p);
  tb=disjoint(p,c);
  tb=disjoint(p,p);
  tb=subset(bx,p);
  tb=subset(c,p);
  tb=subset(p,p);
  tb=subset(p,bx);
  
  closed_intersection(p,p);
  closed_intersection(r,p);
  closed_intersection(p,r);
  open_intersection(p,p);
  
  p=polyhedron(r);
  ip=polyhedron(c);
}



} // namespace Ariadne

