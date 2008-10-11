/***************************************************************************
 *            vector.inline.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
#include "linear_algebra/slice.h"

namespace Ariadne {
    

template<class R>
Vector<R>::Vector()
  : _array(0) 
{ 
}

template<class R> inline
Vector<R>::Vector(const size_type& n)
  : _array(n,static_cast<R>(0)) 
{ 
}

template<class R> inline
Vector<R>::Vector(const size_type& n, const R& x)
  : _array(n) 
{ 
  for(size_type i=0; i!=n; ++i) { (*this)(i)=x; }
}


template<class R> template<class RR> inline
Vector<R>::Vector(const array<RR>& ary)
  : _array(ary) 
{ 
}

template<class R> template<class RR> inline
Vector<R>::Vector(const size_type& n, const RR* ptr, const size_type& inc)
  : _array(n) 
{ 
  for(size_type i=0; i!=n; ++i) { 
    (*this)(i)=ptr[i*inc]; 
  }
}


template<class R> template<class E> inline
Vector<R>::Vector(const VectorExpression<E>& ve)
  : _array(ve().size()) 
{ 
  const E& v=ve(); 
  for(size_type i=0; i!=this->size(); ++i) {
    (*this)(i)=R(v(i)); 
  }
}

template<class R> inline
Vector<R>::Vector(const Vector<R>& v) 
  : _array(v._array) 
{ 
}



template<class R> inline
Vector<R>& 
Vector<R>::operator=(const Vector<R>& v) 
{
  if(this!=&v) { this->_array=v._array; } return *this; }


template<class R1, class R2> inline
bool 
operator==(const Vector<R1>& v1, const Vector<R2>& v2)  
{
  if(v1.size()!=v2.size()) { return false; }
  for(size_type i=0; i!=v1.size(); ++i) { 
    if(v1(i)!=v2(i)) { return false; } 
  }
  return true; 
}


template<class R1, class R2> inline
bool 
operator!=(const Vector<R1>& v1, const Vector<R2>& v2)  
{
  return !(v1==v2); 
}


template<class R> inline
Vector<R> 
Vector<R>::zero(const size_type& n) 
{
  Vector<R> v(n); for (dimension_type i=0; i<n; ++i) { v(i)=R(0); } return v; 
}


template<class R> inline
Vector<R> 
Vector<R>::one(const size_type& n) 
{
  Vector<R> v(n); for (dimension_type i=0; i<n; ++i) { v(i)=R(1); } return v; 
}


template<class R> inline
Vector<R> 
Vector<R>::unit(dimension_type n, dimension_type i) 
{
  Vector<R> v(n); for (dimension_type k=0; k<n; ++k) { v(k)=R(0); } v(i)=R(1); return v; }


template<class R> inline
array<R>& 
Vector<R>::data() 
{
  return this->_array; 
}


template<class R> inline
const array<R>&
Vector<R>::data() const 
{
  return this->_array; 
}


template<class R> inline
void 
Vector<R>::resize(const size_type& n) 
{
  this->_array.resize(n); 
}


template<class R> inline
size_type 
Vector<R>::size() const 
{
  return this->_array.size(); 
}


template<class R> inline
R* 
Vector<R>::begin() 
{
  return this->data().begin(); 
}

template<class R> inline
const R* 
Vector<R>::begin() const 
{
  return this->data().begin(); 
}

template<class R> inline
R* 
Vector<R>::end() 
{ 
  return this->data().end(); 
}

template<class R> inline
const R* 
Vector<R>::end() const 
{
  return this->data().end(); 
}


template<class R> inline
size_type 
Vector<R>::increment() const 
{
  return 1; 
}


template<class R> inline
VectorSlice<const R> 
Vector<R>::operator() (const Slice& i) const
{
  return VectorSlice<const R>(i.size(),this->_array.begin()+i.start(),i.stride());
}

template<class R> inline
VectorSlice<R> 
Vector<R>::operator() (const Slice& i) 
{
  return VectorSlice<R>(i.size(),this->_array.begin()+i.start(),i.stride());
}




template<class R> inline
R& 
Vector<R>::operator() (const size_type& i) 
{
  return this->_array[i]; 
}

template<class R> inline
const R& 
Vector<R>::operator() (const size_type& i) const 
{
  return this->_array[i]; 
}


template<class R> inline
R& 
Vector<R>::operator[] (const size_type& i) 
{ 
  return this->_array[i]; 
}

template<class R> inline
const R& 
Vector<R>::operator[] (const size_type& i) const 
{
  return this->_array[i]; 
}










template<class R1, class E2> inline
Vector<R1>&
operator+=(Vector<R1>& v1, const VectorExpression<E2>& e2) {
  const E2& v2=e2();
  ARIADNE_CHECK_EQUAL_SIZES(v1,v2,"Vector& operator+=(Vector,VectorExpression)");
  for(size_type i=0; i!=v1.size(); ++i) { v1(i)+=v2(i); } return v1; 
}

template<class R1, class E2> inline
Vector<R1>&
operator-=(Vector<R1>& v1, const VectorExpression<E2>& e2) {
  const E2& v2=e2();
  ARIADNE_CHECK_SIZE(v1,v2.size(),"Vector& operator-=(Vector,VectorExpression)");
  //for(size_type i=0; i!=v1.size(); ++i) { v1(i)-=v2(i); } return v1; 
  for(size_type i=0; i!=v1.size(); ++i) { v1(i)=v1(i)-v2(i); } return v1; 
}

template<class R1, class R2> inline
Vector<R1>&
operator*=(Vector<R1>& v1, const R2& s2) {
  for(size_type i=0; i!=v1.size(); ++i) { v1(i)*=s2; } return v1; 
}

template<class R1, class R2> inline
Vector<R1>&
operator/=(Vector<R1>& v1, const R2& s2) {
  for(size_type i=0; i!=v1.size(); ++i) { v1(i)/=s2; } return v1; 
}



template<class E1, class E2> inline
BinaryVectorVectorExpression<Add,E1,E2> 
operator+(const VectorExpression<E1>& e1, const VectorExpression<E2>& e2) {
  const E1& v1=e1(); const E2& v2=e2();
  if(v1.size()!=v2.size()) {
    ARIADNE_THROW(IncompatibleSizes,"VectorExpression operator+(VectorExpression ve1, VectorExpression ve2)","ve1.size()="<<v1.size()<<", ve2.size()="<<v2.size());
  }
  return BinaryVectorVectorExpression<Add,E1,E2>(Add(),v1,v2);
}


template<class E1, class E2> inline
BinaryVectorVectorExpression<Sub,E1,E2> 
operator-(const VectorExpression<E1>& e1, const VectorExpression<E2>& e2) {
  const E1& v1=e1(); const E2& v2=e2();
  if(v1.size()!=v2.size()) {
    ARIADNE_THROW(IncompatibleSizes,"VectorExpression operator-(VectorExpression ve1, VectorExpression ve2)","ve1.size()="<<v1.size()<<", ve2.size()="<<v2.size());
  }
  return BinaryVectorVectorExpression<Sub,E1,E2>(Sub(),v1,v2);
}


template<class E1, class E2> inline
BinaryVectorScalarExpression<Mul,E1,E2> 
operator*(const E2& e2, const VectorExpression<E1>& e1) {
  const E1& v1=e1(); const E2& s2=e2;
  return BinaryVectorScalarExpression<Mul,E1,E2>(Mul(),v1,s2);
}


template<class E1, class E2> inline
BinaryVectorScalarExpression<Mul,E1,E2> 
operator*(const VectorExpression<E1>& e1, const E2& e2) {
  const E1& v1=e1(); const E2& s2=e2;
  return BinaryVectorScalarExpression<Mul,E1,E2>(Mul(),v1,s2);
}


template<class E1, class E2> inline
BinaryVectorScalarExpression<Div,E1,E2> 
operator/(const VectorExpression<E1>& e1, const E2& e2) {
  const E1& v1=e1(); const E2& s2=e2;
  return BinaryVectorScalarExpression<Div,E1,E2>(Div(),v1,s2);
}




template<class R> inline 
Vector<R>
zero_vector(size_type n)
{
  return Vector<R>(n);
}


template<class R> inline 
Vector<R>
unit_vector(size_type n, size_type i)
{
  Vector<R> result(n);
  result[i]=1;
  return result;
}


template<class R> inline 
Vector<R>
midpoint(const Vector< Interval<R> >& iv) 
{
  Vector<R> result(iv.size());
  for(size_type i=0; i!=iv.size(); ++i) {
    result(i) = midpoint(iv(i));
  }
  return result;
}

template<class R> inline
Vector< Interval<R> >
intersection(const Vector< Interval<R> >& iv1, const Vector< Interval<R> >& iv2) 
{
  ARIADNE_ASSERT(iv1.size()==iv2.size());
  Vector< Interval<R> > ivr(iv1.size());
  for(uint i=0; i!=iv1.size(); ++i) {
    ivr[i]=intersection(iv1[i],iv2[i]);
  }
  return ivr;
}


template<class T> inline 
Vector< Float<T> >
radius(const Vector< Interval< Float<T> > >& iv) 
{
  Vector< Float<T> > result(iv.size());
  for(size_type i=0; i!=iv.size(); ++i) {
    result(i) = iv(i).radius();
  }
  return result;
}



template<class R> inline
bool
disjoint(const Vector< Interval<R> >& iv1, const Vector< Interval<R> >& iv2) 
{
  ARIADNE_ASSERT(iv1.size()==iv2.size());
  for(uint i=0; i!=iv1.size(); ++i) {
    if(disjoint(iv1[i],iv2[i])) {
      return true;
    }
  }
  return false;
}


template<class R> inline
bool
encloses(const Vector< Interval<R> >& iv, const Vector<R>& v) 
{
  ARIADNE_CHECK_EQUAL_SIZES(iv,v,"bool contains_value(Vector<Interval>,Vector<Float>)");
  for(size_type i=0; i!=v.size(); ++i) {
    if(!encloses(iv(i),v(i))) {
      return false;
    }
  }
  return true;
}


template<class R> inline
bool
subset(const Vector< Interval<R> >& iv1, const  Vector< Interval<R> >& iv2) 
{
  ARIADNE_CHECK_EQUAL_SIZES(iv1,iv2,"bool subset(Vector<Interval>,Vector<Interval>)");
  for(size_type i=0; i!=iv1.size(); ++i) {
    if(!refines(iv1(i),iv2(i))) {
      return false;
    }
  }
  return true;
}

template<class R> inline
bool
refines(const Vector< Interval<R> >& iv1, const  Vector< Interval<R> >& iv2) 
{
  ARIADNE_CHECK_EQUAL_SIZES(iv1,iv2,"bool refines(Vector<Interval>,Vector<Interval>)");
  return subset(iv1,iv2);
}


template<class R> inline 
Vector<R>
approximation(const Vector<R>& v) 
{
  return v;
}

template<class R> inline 
Vector<R>
approximation(const Vector< Interval<R> >& iv) 
{
  return midpoint(iv);
}

template<class R> inline 
bool
operator==(const Vector<R>& v, int n) 
{
  assert(n==0);
  for(size_type i=0; i!=v.size(); ++i) {
    if(v[i]!=n) {
      return false;
    }
  }
  return true;
}

template<class R> inline 
tribool
operator>=(const Vector< Interval<R> >& v, int n) 
{
  tribool result=true;
  for(size_type i=0; i!=v.size(); ++i) {
    if(v[i].upper()<n) {
      return false;
    }
    if(v[i].lower()<=n) {
      result=indeterminate;
    }
  }
  return result;
}


template<class R> inline 
tribool
operator<=(const Vector< Interval<R> >& v, int n) 
{
  tribool result=true;
  for(size_type i=0; i!=v.size(); ++i) {
    if(v[i].upper()>n) {
      return false;
    }
    if(v[i].lower()>=n) {
      result=indeterminate;
    }
  }
  return result;
}


template<class R> inline
std::ostream&
operator<<(std::ostream& os, const Vector<R>& v)
{
  return v.write(os);
}

template<class R> inline
std::istream&
operator>>(std::istream& is, Vector<R>& v)
{
  return v.read(is);
}





template<class R> inline
Vector<R> 
operator-(const Vector<R>& v) 
{
  Vector<R> result(v.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=-v(i);
  }
  return result;
}


template<class R1, class R2> inline
Vector<typename traits<R1,R2>::arithmetic_type> 
operator+(const Vector<R1>& v1, const Vector<R2>& v2) 
{
  ARIADNE_CHECK_EQUAL_SIZES(v1,v2,"Vector operator+(Vector,Vector)");
  Vector<typename traits<R1,R2>::arithmetic_type> result(v1.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v1(i)+v2(i);
  }
  return result;
}


template<class R1, class R2> inline
Vector<class traits<R1,R2>::arithmetic_type> 
operator-(const Vector<R1>& v1, const Vector<R2>& v2) 
{
  ARIADNE_CHECK_EQUAL_SIZES(v1,v2,"Vector operator-(Vector,Vector)");
  Vector<typename traits<R1,R2>::arithmetic_type> result(v1.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v1(i)-v2(i);
  }
  return result;
}


template<class R1, class R2> inline
Vector<typename traits<R1,R2>::arithmetic_type> 
operator*(const R1& s, const Vector<R2>& v) 
{
  return v*s;
}


template<class R1, class R2> inline
Vector<typename traits<R1,R2>::arithmetic_type> 
operator*(const Vector<R1>& v, const R2& s) 
{
  Vector<typename traits<R1,R2>::arithmetic_type> result(v.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v(i)*s;
  }
  return result;
}


template<class R1, class R2> inline
Vector<typename traits<R1,R2>::arithmetic_type> 
operator/(const Vector<R1>& v, const R2& s) 
{
  Vector<typename traits<R1,R2>::arithmetic_type> result(v.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v(i)/s;
  }
  return result;
}



template<class R1, class R2> inline
Vector<typename traits<R1,R2>::arithmetic_type> 
operator*(const R1& s, const VectorSlice<R2>& v) 
{
  return v*s;
}

template<class R1, class R2> inline
Vector<typename traits<R1,R2>::arithmetic_type> 
operator*(const VectorSlice<R1>& v, const R2& s) 
{
  Vector<typename traits<R1,R2>::arithmetic_type> result(v.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v(i)*s;
  }
  return result;
}

template<class R1, class R2> inline
Vector<typename traits<R1,R2>::arithmetic_type> 
operator/(const VectorSlice<R1>& v, const R2& s) 
{
  Vector<typename traits<R1,R2>::arithmetic_type> result(v.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v(i)/s;
  }
  return result;
}






template<class R> inline
Vector<R> 
add_approx(const Vector<R>& u, const Vector<R>& v) 
{
  ARIADNE_CHECK_EQUAL_SIZES(u,v,"Vector add_approx(Vector,Vector)");
  Vector<R> result(u.size());
  for(size_type i=0; i!=u.size(); ++i) {
    result(i)=add<RoundApprox>(u(i),v(i));
  }
  return result;
}


template<class R> inline
Vector<R> 
sub_approx(const Vector<R>& u, const Vector<R>& v)
{
  ARIADNE_CHECK_EQUAL_SIZES(u,v,"Vector sub_approx(Vector,Vector)");
  Vector<R> result(u.size());
  for(size_type i=0; i!=u.size(); ++i) {
    result(i)=sub<RoundApprox>(u(i),v(i));
  }
  return result;
}


template<class R> inline
Vector<R> 
mul_approx(const R& s, const Vector<R>& v)
{
  Vector<R> result(v.size());
  for(size_type i=0; i!=v.size(); ++i) {
    result(i)=mul<RoundApprox>(v(i),s);
  }
  return result;
}


template<class R> inline
Vector<R> mul_approx(const Vector<R>& v, const R& s)
{
  Vector<R> result(v.size());
  for(size_type i=0; i!=v.size(); ++i) {
    result(i)=mul<RoundApprox>(v(i),s);
  }
  return result;
}


template<class R> inline
Vector<R> 
div_approx(const Vector<R>& v, const R& s)
{
  Vector<R> result(v.size());
  for(size_type i=0; i!=v.size(); ++i) {
    result(i)=div<RoundApprox>(v(i),s);
  }
  return result;
}



template<class R> inline 
R 
inner_product(const Vector<R>& u, const Vector<R>& v) 
{
  ARIADNE_CHECK_EQUAL_SIZES(u,v,"Scalar inner_product(Vector,Vector)");
  R result=0;
  for(size_type i=0; i!=u.size(); ++i) {
    result+=u(i)*v(i);
  }
  return result;
}


template<class R> inline 
Vector<R> 
direct_sum(const Vector<R>& v1, const Vector<R>& v2) 
{
  return concatenate(v1,v2);
}


template<class R> inline 
Vector<R> 
concatenate(const Vector<R>& v1, const Vector<R>& v2) 
{
  Vector<R> result(v1.size()+v2.size());
  for(size_type i=0; i!=v1.size(); ++i) {
    result(i)=v1(i);
  }
  for(size_type i=0; i!=v2.size(); ++i) {
    result(i+v1.size())=v2(i);
  }
  return result;
}

template<class R> inline 
Vector<R> 
concatenate(const Vector<R>& v, const R& s) 
{
  Vector<R> result(v.size()+1);
  for(size_type i=0; i!=v.size(); ++i) {
    result(i)=v(i);
  }
  result(v.size())=s;
  return result;
}


template<class R> inline 
Vector<R> 
join(const Vector<R>& v1, const Vector<R>& v2) 
{
  return concatenate(v1,v2);
}


template<class R> inline 
Vector<R> 
join(const Vector<R>& v, const R& s) 
{
  return concatenate(v,s);
}


template<class R> inline
R 
sup_norm(const Vector<R>& v) 
{
  R result=static_cast<R>(0);
  for(size_type i=0; i!= v.size(); ++i) {
    result=max(R(abs(v(i))),result); }
  return result; 
}


template<class R> inline 
R 
norm(const Vector<R>& v) 
{
  return sup_norm(v);
}


} //namespace Ariadne
