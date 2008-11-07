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

#include "numeric/rational.h"
#include "numeric/float.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_programming/linear_program.h"

#include "geometry/ddconv.h"
#include "geometry/ddconv.code.h"
#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/box.h"
#include "geometry/polytope.h"

#include "output/logging.h"


// Specializations of polyhedral operations
namespace {

using namespace Ariadne;

inline
void
assign(Polyhedron<Rational>& plhd, const PointList<Rational>& vl)
{
  plhd=Polyhedron<Rational>(Polytope<Rational>(vl));
}

template<class X> inline
void
assign(Polyhedron<X>& plhd, const PointList<X>& pltp)
{
  throw NotImplemented(__PRETTY_FUNCTION__); 
}



inline
tribool 
empty(const Polyhedron<Rational>& plhd) 
{
  try {
    Polytope<Rational> pltp(plhd);
    ARIADNE_LOG(7,"empty(plhd): plhd="<<plhd<<" pltp="<<pltp);
    return pltp.empty();
  }
  catch(UnboundedSet& e) {
    return false;
  }
}

template<class R> inline
tribool 
empty(const Polyhedron<R>& plhd) 
{
  return empty(Polyhedron<Rational>(plhd));
}

template<class R> inline
tribool 
empty(const Polyhedron< Interval<R> >& plhd) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



inline
tribool 
bounded(const Polyhedron<Rational>& plhd)
{
  try {
    Polytope<Rational> pltp(plhd);
    return pltp.bounded();
  }
  catch(UnboundedSet& e) {
    return false;
  }
}

template<class R> inline
tribool 
bounded(const Polyhedron<R>& plhd) 
{
  return bounded(Polyhedron<Rational>(plhd));
}


template<class R> inline
tribool 
bounded(const Polyhedron< Interval<R> >& plhd) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


tribool 
disjoint(const Polyhedron<Rational>& plhd, const Box<Rational>& r)
{
  typedef Rational R;
  typedef Rational F;
  
  if(verbosity>7) { std::clog << "tribool disjoint(Polyhedron<Rational>,Polyhedron<Rational>)" << std::endl; }
  if(verbosity>8) { std::clog << plhd << " " << r << std::endl; }
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(plhd,r,"tribool disjoint(Polyhedron<Rational> plhd, Box<Rational> r)");
  dimension_type d=plhd.dimension();
  size_type nc=plhd.number_of_constraints();
  size_type nnc=0; // number of negative constraints Ax<=b
  
  // Translate x'=x-l;  Ax<=b A(x'+l)<=b Ax'<= b-Al b'=b-Al 
  Vector<F> u=r.upper_corner()-r.lower_corner();
  Matrix<R> A=plhd.A();
  Vector<F> b=plhd.b()-A*r.lower_corner().position_vector();
  
  if(verbosity>8) { std::clog << "u'=" << u << " A'=" << A << " b'=" << b << std::endl; }
  
  // count number of negative constraints
  for(size_type i=0; i!=nc; ++i) {
    if(b(i)<0) {
      ++nnc;
    }
  }
  
  Matrix<F> T(nc+d+1,d+nnc+1);
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
  
  LinearProgram<F> lp(T);
  tribool result=!lp.is_feasible();
  return result;
}

template<class T>
tribool 
disjoint(const Polyhedron< Float<T> >& plhd, 
         const Box< Float<T> >& r)
{
  return ::disjoint(Polyhedron<Rational>(plhd),
                    Box<Rational>(r));
}
                  
template<class R>
tribool 
disjoint(const Polyhedron< Interval<R> >& plhd, const Box<R>& bx)
{
  throw NotImplemented("disjoint(Polyhedron<Fuzzy>,Box<Fuzzy>");
}
                  

tribool 
disjoint(const Polyhedron<Rational>& plhd, const Polytope<Rational>& pltp)
{
  return disjoint(plhd,Polyhedron<Rational>(pltp));
}

template<class T>
tribool 
disjoint(const Polyhedron< Float<T> >& plhd, 
         const Polytope< Float<T> >& pltp)
{
  return ::disjoint(Polyhedron<Rational>(plhd),
                    Polytope<Rational>(pltp));
}
                  
template<class R>
tribool 
disjoint(const Polyhedron< Interval<R> >& plhd, const Polytope< Interval<R> >& pltp)
{
  throw NotImplemented("disjoint(Polyhedron<Fuzzy>,Polytope<Fuzzy>");
}
                  




template<class R, class BS> inline
tribool 
subset(const BS& A, const Polyhedron<R>& B)
{
  tribool result=true;
  for(typename BS::vertices_const_iterator v=A.vertices_begin();
      v!=A.vertices_end(); ++v)
  {
    for(typename Polyhedron<R>::constraints_const_iterator c=B.constraints_begin();
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
subset(const Polyhedron<Rational>& plhd, const Box<Rational>& r)
{
  ARIADNE_LOG(3,"subset(plhd,r): ""plhd="<<plhd<<", r="<<r<<"\n");
  typedef Rational Q;
  tribool result=true;
  dimension_type d=plhd.dimension();
  array<Q> data(d+1u);
  for(size_type i=0; i!=d+1u; ++i) {
    data[i]=0;
  }
  Polyhedron<Q> halfspace;
  ARIADNE_LOG(9,"data="<<data<<"\n");
  
  for(size_type j=0; j!=d; ++j) {
    // Check if disjoint from lower halfspace
    data[j]=-1;
    data[d]=r.lower_bound(j);
    ARIADNE_LOG(9,"data="<<data<<"\n");
    halfspace=Polyhedron<Q>(d,1u,data.begin());
    ARIADNE_LOG(8,"halfspace="<<halfspace<<"\n");
    result=result && disjoint(plhd,halfspace);
    // Check if disjoint from upper halfspace
    data[j]=1;
    data[d]=-r.upper_bound(j);
    halfspace=Polyhedron<Q>(d,1u,data.begin());
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
subset(const Polyhedron< Float<T> >& plhd, 
       const Box< Float<T> >& r)
{
  return subset(Polyhedron<Rational>(plhd),Box<Rational>(r));
}

template<class R>
inline
tribool 
subset(const Polyhedron< Interval<R> >& plhd, const Box<R>& r)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


inline
tribool 
subset(const Polyhedron<Rational>& plhd1, const Polyhedron<Rational>& plhd2)
{
  return subset(Polytope<Rational>(plhd1),plhd2);
}

template<class R>
inline
tribool 
subset(const Polyhedron<R>& plhd1, const Polyhedron<R>& plhd2)
{
  return ::subset(Polyhedron<Rational>(plhd1),Polyhedron<Rational>(plhd2));
}

template<class R>
inline
tribool 
subset(const Polyhedron< Interval<R> >& plhd1, const Polyhedron< Interval<R> >& plhd2)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}






inline
Box<Rational> 
bounding_box(const Polyhedron<Rational>& plhd)
{
  return Polytope<Rational>(plhd).bounding_box();
}

template<class R> inline
Box<R> 
bounding_box(const Polyhedron<R>& plhd)
{
  Box<Rational> qbb=bounding_box(Polyhedron<Rational>(plhd));
  Box<R> bb(plhd.dimension());
  for(dimension_type i=0; i!=bb.dimension(); ++i) {
    bb.set_lower_bound(i,R(qbb.lower_bound(i),round_down));
    bb.set_upper_bound(i,R(qbb.upper_bound(i),round_up));
  }
  return bb;
}

template<class R> inline
Box<R> 
bounding_box(const Polyhedron< Interval<R> >& plhd)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}




} // namespace

  

namespace Ariadne {

template<class X> 
tribool 
empty(const Polyhedron<X>& plhd) 
{
  return ::empty(plhd);
}

template<class X> 
tribool 
bounded(const Polyhedron<X>& plhd) 
{
  return ::bounded(plhd);
}







template<class X> 
Polyhedron<X>::Polyhedron(const Box<R>& bx)
  : _dimension(bx.dimension()), 
    _number_of_constraints(bx.dimension()*2u),
    _data((bx.dimension()+1u)*bx.dimension()*2u,static_cast<X>(0))
{
  dimension_type d=bx.dimension();
  MatrixSlice<X> constraints=this->_constraints();
  for(size_type i=0; i!=d; ++i) {
    constraints(2*i,i)=static_cast<X>(1);
    constraints(2*i,d)=-bx.lower_bound(i);
    constraints(2*i+1,i)=static_cast<X>(-1); 
    constraints(2*i+1,d)=bx.upper_bound(i);
  }
}





template<class X>
Polyhedron<X>::Polyhedron(const std::string& str)
{
  std::stringstream ss(str);
  ss >> *this;
}


template<class X>
Polyhedron<X>::Polyhedron(dimension_type d)
  : _dimension(d), 
    _number_of_constraints(0), 
    _data()
{
}


template<class X>
Polyhedron<X>::Polyhedron(dimension_type d, size_type nc, const X* data)
  : _dimension(d), 
    _number_of_constraints(nc), 
    _data(data,data+(d+1)*nc)
{ 
}


template<class X>
Polyhedron<X>::Polyhedron(const Matrix<X>& A,
                          const Vector<X>& b) 
  : _dimension(A.number_of_columns()), 
    _number_of_constraints(A.number_of_rows()), 
    _data((A.number_of_columns()+1)*A.number_of_rows())
{
  
  ARIADNE_CHECK_SIZE(b,A.number_of_rows(),"Polyhedron::Polyhedron(Matrix A, Vector b)");
  dimension_type d=this->dimension();
  dimension_type nc=this->number_of_constraints();
  MatrixSlice<X>(nc,d,this->data().begin(),d+1u,1u)=-A;
  VectorSlice<X>(nc,this->data().begin()+d,d+1u)=b;
}


template<class X>
Polyhedron<X>::Polyhedron(const Matrix<X>& C) 
  : _dimension(C.number_of_columns()-1), 
    _number_of_constraints(C.number_of_rows()), 
    _data(C.data())
{
}


template<class X>
Polyhedron<X>::Polyhedron(const PointList<X>& pts)
  : _dimension(pts.dimension()), _number_of_constraints(0), _data()
{
  assign(*this,pts);
}


template<class X>
Polyhedron<X>::Polyhedron(const std::vector< Halfspace<X> >& hsv)
  : _dimension(hsv[0].dimension()), _number_of_constraints(0), _data()
{
  for(size_type i=0; i!=hsv.size(); ++i) {
    this->new_constraint(hsv[i]);
  }
}


template<class X>
Polyhedron<X>::Polyhedron(const Halfspace<X>& hs)
  : _dimension(hs.dimension()), _number_of_constraints(0), _data()
{
  this->new_constraint(hs);
}



template<class X>
void
Polyhedron<X>::new_constraint(const Halfspace<X>& c)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,c,"void Polyhedron::new_constraint(PolyhedralConstraint& c)");
  dimension_type d=this->_dimension;
  size_type sz=this->_data.size();
  this->_data.resize(sz+d+1u);
  ++this->_number_of_constraints;
  for(dimension_type i=0; i!=d; ++i) {
    this->_data[sz+i]=c.data()[i];
    this->_data[sz+d]=c.data()[d];
  }
}


template<class X>
Halfspace<X>
Polyhedron<X>::constraint(size_type i) const
{
  ARIADNE_ASSERT(i<this->_number_of_constraints);
  return Halfspace<X>(this->_dimension,this->_data.begin()+i*(this->_dimension+1u));
}


template<class X>
MatrixSlice<X>
Polyhedron<X>::_constraints()
{
  return MatrixSlice<X>(this->_number_of_constraints,
                                       this->_dimension+1u,
                                       this->data().begin(),
                                       this->_dimension+1u,
                                       1u);
}




template<class X>
typename Polyhedron<X>::constraints_const_iterator
Polyhedron<X>::constraints_begin() const
{
  return constraints_const_iterator(*this,0u);
}

template<class X>
typename Polyhedron<X>::constraints_const_iterator
Polyhedron<X>::constraints_end() const
{
  return constraints_const_iterator(*this,this->number_of_constraints());
}

template<class X>
dimension_type
Polyhedron<X>::dimension() const
{
  return this->_dimension;
}



template<class X> 
Box<typename Polyhedron<X>::real_type> 
Polyhedron<X>::bounding_box() const
{
  return ::bounding_box(*this);
}

template<class X>
tribool 
Polyhedron<X>::empty() const
{
  return ::empty(*this);
}

template<class X>
tribool 
Polyhedron<X>::bounded() const
{
  return ::bounded(*this);
}




template<class X>
tribool 
contains(const Polyhedron<X>& plhd, const Point<typename Polyhedron<X>::real_type>& pt)
{
  return plhd.contains(pt);
}

template<class X>
tribool 
disjoint(const Polyhedron<X>& plhd, const Box<typename Polyhedron<X>::real_type>& bx)
{
  return ::disjoint(plhd,bx);
}


template<class X>
tribool 
disjoint(const Box<typename Polyhedron<X>::real_type>& bx, const Polyhedron<X>& plhd)
{
  return ::disjoint(plhd,bx);
}



template<class X>
tribool 
disjoint(const Polyhedron<X>& plhd, const Polytope<X>& pltp)
{
  return ::disjoint(plhd,pltp);
}


template<class X>
tribool 
disjoint(const Polytope<X>& pltp, const Polyhedron<X>& plhd)
{
  return ::disjoint(plhd,pltp);
}



template<class X>
tribool 
disjoint(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2)
{
  return closed_intersection(plhd1,plhd2).empty();
}

template<class X>
tribool 
subset(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2)
{
  return ::subset(plhd1,plhd2);
}

template<class X>
tribool 
subset(const Polyhedron<X>& plhd, const Box<typename Polyhedron<X>::real_type>& r)
{
  return ::subset(plhd,r);
}


template<class X>
tribool 
subset(const Box<typename Polyhedron<X>::real_type>& r, const Polyhedron<X>& plhd)
{
  return ::subset(r,plhd);
}


template<class X1,class X2>
tribool 
subset(const Polytope<X1>& pltp, const Polyhedron<X2>& plhd)
{
  return ::subset(pltp,plhd);
}




template<class X>
Polyhedron<X> 
open_intersection(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(plhd1,plhd2,"Polyhedron open_intersection(Polyhedron plhd1, Polyhedron plhd2)");
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class X>
Polyhedron<X> 
closed_intersection(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(plhd1,plhd2,"Polyhedron closed_intersection(Polyhedron plhd1, Polyhedron plhd2)");
  dimension_type d=plhd1.dimension();
  size_type nc1=plhd1.number_of_constraints();
  size_type nc2=plhd2.number_of_constraints();
  Matrix<X> A(nc1+nc2,d);
  MatrixSlice<X>(nc1,d,A.begin(),d,1)=plhd1.A();
  MatrixSlice<X>(nc2,d,A.begin()+nc1*d,d,1)=plhd2.A();
  Vector<X> b=direct_sum(plhd1.b(),plhd2.b());
  return Polyhedron<X>(A,b);
}


template<class X>
Polyhedron<X> 
closed_intersection(const Box<typename Polyhedron<X>::real_type>& r, const Polyhedron<X>& plhd)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r,plhd,"Polyhedron closed_intersection(Rectangle r, Polyhedron plhd)");
  return closed_intersection(Polyhedron<X>(r),plhd);
}

 
template<class X>
Polyhedron<X> 
closed_intersection(const Polyhedron<X>& plhd, const Box<typename Polyhedron<X>::real_type>& r)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(plhd,r,"Polyhedron closed_intersection(Polyhedron plhd, Rectangle r)");
  return closed_intersection(plhd,Polyhedron<X>(r));
}


template<class X>
Polyhedron<X> 
polyhedron(const Halfspace<X>& hs)
{
  return Polyhedron<X>(hs);
}


template<class R>
Polyhedron<R> 
polyhedron(const Box<R>& r)
{
  return Polyhedron<R>(r);
}


template<class X>
Polyhedron<typename traits<X>::arithmetic_type>
polyhedron(const Polytope<X>& pltp)
{
  typedef typename traits<X>::arithmetic_type F;
  
  dimension_type d=pltp.dimension();
  size_type nv=pltp.number_of_vertices();
  
  if(nv==0) {
    // empty polytope; return empty polyhedron
    array<F> data(d+1,F(0));
    data[d]=-1;
    return Polyhedron<F>(d,1,data.begin());
  }
  
  const Matrix<X> G=pltp.generators();
  
  std::vector< Vector<F> > result;
  std::vector< Vector<F> > argument;
  
  Vector<F> tmp(d+1);
  
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


template<class R>
Polyhedron<R>
approx_polyhedron(const Polytope<R>& pltp)
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  
  dimension_type d=pltp.dimension();
  size_type nv=pltp.number_of_vertices();
  
  if(nv==0) {
    // empty polytope; return empty polyhedron
    array<R> data(d+1,R(0));
    data[d]=-1;
    return Polyhedron<R>(d,1,data.begin());
  }
  
  const Matrix<A> G=pltp.generators();
  
  std::vector< Vector<A> > result;
  std::vector< Vector<A> > argument;
  
  Vector<A> tmp(d+1);
  
  for(size_type j=0; j!=nv; ++j) {
    for(size_type i=0; i!=d+1u; ++i) {
      tmp(i)=G(i,j);
    }
    argument.push_back(tmp);
  }
  
  ddconv(result,argument);     
  
  size_type nc=result.size();
  array<R> data((d+1)*nc);
  for(size_type i=0; i!=nc; ++i) {
    for(size_type j=0; j!=d+1u; ++j) {
      data[i*(d+1u)+j]=R(result[i](j));
    }
  }
  
  return Polyhedron<R>(d,nc,data.begin());
}



template<class X>
std::string
Polyhedron<X>::name()
{
  return std::string("Polyhedron")+"<"+Ariadne::name<X>()+">";
}

template<class X>  
std::ostream& 
Polyhedron<X>::write(std::ostream& os) const
{
  //return os << "Polyhedron( A=" << this->A() << ", b=" << this->b() << " )";
  os << "Polyhedron( constraints=";
  dimension_type d=this->dimension();
  size_type nc=this->number_of_constraints();
  for(size_type i=0; i!=nc; ++i) {
    os << ( i==0 ? "[" : "," );
    for(size_type j=0; j!=d; ++j) {
      os << ( j==0 ? "(" : ",");
      os << X(-this->_data[i*(d+1)+j]); 
    }
    os << ";" << this->_data[i*(d+1)+d] << ")";
  }
  os << "] )";
  return os;
}

template<class X>  
std::istream& 
Polyhedron<X>::read(std::istream& is) 
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
    read_sequence(is,a,'[',';',',');
    is >> b;
    is >> c;
    assert(c==']');
    Alst.push_back(a);
    Blst.push_back(b);
  }
  
  size_type m=Alst.size();
  size_type n=Alst[0].size();
  Matrix<X> A(m,n);
  Vector<X> B(m);
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
Polyhedron<X>::_instantiate() 
{
  typedef typename traits<X>::number_type R;
  typedef typename traits<X>::arithmetic_type F;
  tribool tb;
  Point<R> pt;
  Box<R> bx;
  Polytope<R> pl;
  Halfspace<X> hs;
  Polytope<X> c;
  Polyhedron<X> p;
  Polyhedron<F> ip;
  
  p=Polyhedron<X>(bx);
  p=Polyhedron<X>(p);
  ip=Polyhedron<F>(c);
  
  tb=Ariadne::contains(p,pt);
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
  closed_intersection(bx,p);
  closed_intersection(p,bx);
  open_intersection(p,p);
  
  p=polyhedron(hs);
  p=polyhedron(bx);
  ip=polyhedron(c);

  p=approx_polyhedron(pl);
}



} // namespace Ariadne

