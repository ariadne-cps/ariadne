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
 

namespace Ariadne {
    

template<class R>
LinearAlgebra::Vector<R>::Vector()
  : _array(0) 
{ 
}

template<class R> inline
LinearAlgebra::Vector<R>::Vector(const size_type& n)
  : _array(n,static_cast<R>(0)) 
{ 
}

template<class R> inline
LinearAlgebra::Vector<R>::Vector(const size_type& n, const R& x)
  : _array(n) 
{ 
  for(size_type i=0; i!=n; ++i) { (*this)(i)=x; }
}

template<class R> inline
LinearAlgebra::Vector<R>::Vector(const array<R>& ary)
  : _array(ary) 
{ 
}

template<class R> inline
LinearAlgebra::Vector<R>::Vector(const size_type& n, const R* ptr, const size_type& inc)
  : _array(n) 
{ 
  for(size_type i=0; i!=n; ++i) { 
    (*this)(i)=ptr[i*inc]; 
  }
}


template<class R> template<class E> inline
LinearAlgebra::Vector<R>::Vector(const VectorExpression<E>& ve)
  : _array(ve().size()) 
{ 
  const E& v=ve(); 
  for(size_type i=0; i!=this->size(); ++i) {
    (*this)(i)=v(i); 
  }
}

template<class R> inline
LinearAlgebra::Vector<R>::Vector(const Vector<R>& v) 
  : _array(v._array) 
{ 
}



template<class R> inline
LinearAlgebra::Vector<R>& 
LinearAlgebra::Vector<R>::operator=(const Vector<R>& v) 
{
  if(this!=&v) { this->_array=v._array; } return *this; }


template<class R> inline
bool 
LinearAlgebra::Vector<R>::operator==(const Vector<R>& v) const 
{
  if(this->size()!=v.size()) { return false; }
  for(size_type i=0; i!=this->size(); ++i) { 
    if((*this)(i)!=v(i)) { return false; } 
  }
  return true; 
}


template<class R> inline
bool 
LinearAlgebra::Vector<R>::operator!=(const Vector<R>& v) const 
{
  return !(*this==v); 
}


template<class R> inline
LinearAlgebra::Vector<R> 
LinearAlgebra::Vector<R>::zero(const size_type& n) 
{
  Vector<R> v(n); for (dimension_type i=0; i<n; ++i) { v(i)=R(0); } return v; 
}


template<class R> inline
LinearAlgebra::Vector<R> 
LinearAlgebra::Vector<R>::one(const size_type& n) 
{
  Vector<R> v(n); for (dimension_type i=0; i<n; ++i) { v(i)=R(1); } return v; 
}


template<class R> inline
LinearAlgebra::Vector<R> 
LinearAlgebra::Vector<R>::unit(dimension_type n, dimension_type i) 
{
  Vector<R> v(n); for (dimension_type k=0; k<n; ++k) { v(k)=R(0); } v(i)=R(1); return v; }


template<class R> inline
array<R>& 
LinearAlgebra::Vector<R>::data() 
{
  return this->_array; 
}


template<class R> inline
const array<R>&
LinearAlgebra::Vector<R>::data() const 
{
  return this->_array; 
}


template<class R> inline
void 
LinearAlgebra::Vector<R>::resize(const size_type& n) 
{
  this->_array.resize(n); 
}


template<class R> inline
size_type 
LinearAlgebra::Vector<R>::size() const 
{
  return this->_array.size(); 
}


template<class R> inline
R* 
LinearAlgebra::Vector<R>::begin() 
{
  return this->data().begin(); 
}

template<class R> inline
const R* 
LinearAlgebra::Vector<R>::begin() const 
{
  return this->data().begin(); 
}

template<class R> inline
R* 
LinearAlgebra::Vector<R>::end() 
{ 
  return this->data().end(); 
}

template<class R> inline
const R* 
LinearAlgebra::Vector<R>::end() const 
{
  return this->data().end(); 
}


template<class R> inline
size_type 
LinearAlgebra::Vector<R>::increment() const 
{
  return 1; 
}


template<class R> inline
R& 
LinearAlgebra::Vector<R>::operator() (const size_type& i) 
{
  return this->_array[i]; 
}

template<class R> inline
const R& 
LinearAlgebra::Vector<R>::operator() (const size_type& i) const 
{
  return this->_array[i]; 
}


template<class R> inline
R& 
LinearAlgebra::Vector<R>::operator[] (const size_type& i) 
{ 
  return this->_array[i]; 
}

template<class R> inline
const R& 
LinearAlgebra::Vector<R>::operator[] (const size_type& i) const 
{
  return this->_array[i]; 
}






template<class R> inline
LinearAlgebra::VectorSlice<R>::VectorSlice(const size_type& size, R* begin, const size_type& increment)
  : _size(size), _begin(begin), _increment(increment) 
{
}


template<class R> inline
LinearAlgebra::VectorSlice<R>::VectorSlice(const Vector<R>& v)
  : _size(v.size()), _begin(v.begin()), _increment(v.increment())
{
}


template<class R> inline
LinearAlgebra::VectorSlice<R>::VectorSlice(const array<R>& a)
  : _size(a.size()), _begin(a.begin()), _increment(1u)
{
}


template<class R> inline
size_type 
LinearAlgebra::VectorSlice<R>::size() const 
{
  return this->_size; 
}


template<class R> inline
const R*
LinearAlgebra::VectorSlice<R>::begin() const 
{
  return this->_begin; 
}


template<class R> inline
const R*
LinearAlgebra::VectorSlice<R>::end() const 
{
  return this->_begin+this->_size*this->_increment; 
}


template<class R> inline
size_type 
LinearAlgebra::VectorSlice<R>::increment() const 
{
  return this->_increment; 
}



template<class R> inline     
const R&
LinearAlgebra::VectorSlice<R>::operator() (const size_type& i) const 
{
  return this->_begin[i*this->_increment]; 
}


template<class R> inline    
R& 
LinearAlgebra::VectorSlice<R>::operator() (const size_type& i) 
{
  return this->_begin[i*this->_increment]; 
}


template<class R> inline     
const R&
LinearAlgebra::VectorSlice<R>::operator[] (const size_type& i) const 
{
  return this->_begin[i*this->_increment]; 
}


template<class R> inline   
R&
LinearAlgebra::VectorSlice<R>::operator[] (const size_type& i) 
{
  return this->_begin[i*this->_increment]; 
}


template<class R> template<class E> inline
LinearAlgebra::VectorSlice<R>&
LinearAlgebra::VectorSlice<R>::operator=(const VectorExpression< E >& v) 
{
  const E& e=v(); 
  ARIADNE_CHECK_EQUAL_SIZES(*this,e,"VectorSlice& VectorSlice::operator=(VectorExcpression)");
  for(size_type i=0; i!=e.size(); ++i) {
    this->_begin[i*this->_increment]=e(i); 
  }
  return *this;
}


template<class R> inline     
std::ostream& 
LinearAlgebra::VectorSlice<R>::write(std::ostream& os) const 
{ 
  return Vector<R>(*this).write(os); 
}


template<class R1, class E2> inline
LinearAlgebra::Vector<R1>&
LinearAlgebra::operator+=(Vector<R1>& v1, const VectorExpression<E2>& e2) {
  const E2& v2=e2();
  ARIADNE_CHECK_EQUAL_SIZES(v1,v2,"Vector& operator+=(Vector,VectorExpression)");
  for(size_type i=0; i!=v1.size(); ++i) { v1(i)+=v2(i); } return v1; 
}

template<class R1, class E2> inline
LinearAlgebra::Vector<R1>&
LinearAlgebra::operator-=(Vector<R1>& v1, const VectorExpression<E2>& e2) {
  const E2& v2=e2();
  ARIADNE_CHECK_SIZE(v1,v2.size(),"Vector& operator-=(Vector,VectorExpression)");
  //for(size_type i=0; i!=v1.size(); ++i) { v1(i)-=v2(i); } return v1; 
  for(size_type i=0; i!=v1.size(); ++i) { v1(i)=v1(i)-v2(i); } return v1; 
}

template<class R1, class R2> inline
LinearAlgebra::Vector<R1>&
LinearAlgebra::operator*=(Vector<R1>& v1, const R2& s2) {
  for(size_type i=0; i!=v1.size(); ++i) { v1(i)*=s2; } return v1; 
}

template<class R1, class R2> inline
LinearAlgebra::Vector<R1>&
LinearAlgebra::operator/=(Vector<R1>& v1, const R2& s2) {
  for(size_type i=0; i!=v1.size(); ++i) { v1(i)/=s2; } return v1; 
}



template<class E1, class E2> inline
LinearAlgebra::BinaryVectorVectorExpression<E1,E2,LinearAlgebra::plus> 
LinearAlgebra::operator+(const VectorExpression<E1>& e1, const VectorExpression<E2>& e2) {
  const E1& v1=e1(); const E2& v2=e2();
  if(v1.size()!=v2.size()) {
    ARIADNE_THROW(IncompatibleSizes,"VectorExpression operator+(VectorExpression ve1, VectorExpression ve2)","ve1.size()="<<v1.size()<<", ve2.size()="<<v2.size());
  }
  return BinaryVectorVectorExpression<E1,E2,plus>(v1,v2,plus());
}


template<class E1, class E2> inline
LinearAlgebra::BinaryVectorVectorExpression<E1,E2,LinearAlgebra::minus> 
LinearAlgebra::operator-(const VectorExpression<E1>& e1, const VectorExpression<E2>& e2) {
  const E1& v1=e1(); const E2& v2=e2();
  if(v1.size()!=v2.size()) {
    ARIADNE_THROW(IncompatibleSizes,"VectorExpression operator-(VectorExpression ve1, VectorExpression ve2)","ve1.size()="<<v1.size()<<", ve2.size()="<<v2.size());
  }
  return BinaryVectorVectorExpression<E1,E2,minus>(v1,v2,minus());
}


template<class E1, class E2> inline
LinearAlgebra::BinaryVectorScalarExpression<E1,E2,LinearAlgebra::times> 
LinearAlgebra::operator*(const E2& e2, const VectorExpression<E1>& e1) {
  const E1& v1=e1(); const E2& s2=e2;
  return BinaryVectorScalarExpression<E1,E2,times>(v1,s2,times());
}


template<class E1, class E2> inline
LinearAlgebra::BinaryVectorScalarExpression<E1,E2,LinearAlgebra::times> 
LinearAlgebra::operator*(const VectorExpression<E1>& e1, const E2& e2) {
  const E1& v1=e1(); const E2& s2=e2;
  return BinaryVectorScalarExpression<E1,E2,times>(v1,s2,times());
}


template<class E1, class E2> inline
LinearAlgebra::BinaryVectorScalarExpression<E1,typename Numeric::traits<E2>::closure_type,LinearAlgebra::divides> 
LinearAlgebra::operator/(const VectorExpression<E1>& e1, const E2& e2) {
  const E1& v1=e1(); const E2& s2=e2;
  return BinaryVectorScalarExpression<E1,E2,divides>(v1,s2,divides());
}




template<class R> inline 
LinearAlgebra::Vector<R>
LinearAlgebra::zero_vector(size_type n)
{
  return Vector<R>(n);
}


template<class R> inline 
LinearAlgebra::Vector<R>
LinearAlgebra::unit_vector(size_type n, size_type i)
{
  Vector<R> result(n);
  result[i]=1;
  return result;
}


template<class R> inline 
LinearAlgebra::Vector<R>
LinearAlgebra::approximate_value(const Vector< Numeric::Interval<R> >& iv) 
{
  Vector<R> result(iv.size());
  for(size_type i=0; i!=iv.size(); ++i) {
    result(i) = approximate_value(iv(i));
  }
  return result;
}


template<class R> inline
bool
LinearAlgebra::contains_value(const Vector< Numeric::Interval<R> >& iv,const Vector<R>& v) 
{
  ARIADNE_CHECK_EQUAL_SIZES(iv,v,"bool contains_value(Vector<Interval>,Vector<Float>)");
  for(size_type i=0; i!=v.size(); ++i) {
    if(!Numeric::contains_value(iv(i),v(i))) {
      return false;
    }
  }
  return true;
}


template<class R> inline
std::ostream&
LinearAlgebra::operator<<(std::ostream& os, const Vector<R>& v)
{
  return v.write(os);
}

template<class R> inline
std::istream&
LinearAlgebra::operator>>(std::istream& is, Vector<R>& v)
{
  return v.read(is);
}

template<class R> inline
std::ostream&
LinearAlgebra::operator<<(std::ostream& os, const VectorSlice<R>& vs)
{
  return vs.write(os);
}





template<class R> inline
LinearAlgebra::Vector<R> 
LinearAlgebra::operator-(const Vector<R>& v) 
{
  Vector<R> result(v.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=-v(i);
  }
  return result;
}


template<class R1, class R2> inline
LinearAlgebra::Vector<typename Numeric::traits<R1,R2>::arithmetic_type> 
LinearAlgebra::operator+(const Vector<R1>& v1, const Vector<R2>& v2) 
{
  ARIADNE_CHECK_EQUAL_SIZES(v1,v2,"Vector operator+(Vector,Vector)");
  Vector<typename Numeric::traits<R1,R2>::arithmetic_type> result(v1.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v1(i)+v2(i);
  }
  return result;
}


template<class R1, class R2> inline
LinearAlgebra::Vector<class Numeric::traits<R1,R2>::arithmetic_type> 
LinearAlgebra::operator-(const Vector<R1>& v1, const Vector<R2>& v2) 
{
  ARIADNE_CHECK_EQUAL_SIZES(v1,v2,"Vector operator-(Vector,Vector)");
  Vector<typename Numeric::traits<R1,R2>::arithmetic_type> result(v1.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v1(i)-v2(i);
  }
  return result;
}


template<class R1, class R2> inline
LinearAlgebra::Vector<typename Numeric::traits<R1,R2>::arithmetic_type> 
LinearAlgebra::operator*(const R1& s, const Vector<R2>& v) 
{
  return v*s;
}


template<class R1, class R2> inline
LinearAlgebra::Vector<typename Numeric::traits<R1,R2>::arithmetic_type> 
LinearAlgebra::operator*(const Vector<R1>& v, const R2& s) 
{
  Vector<typename Numeric::traits<R1,R2>::arithmetic_type> result(v.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v(i)*s;
  }
  return result;
}


template<class R1, class R2> inline
LinearAlgebra::Vector<typename Numeric::traits<R1,R2>::arithmetic_type> 
LinearAlgebra::operator/(const Vector<R1>& v, const R2& s) 
{
  Vector<typename Numeric::traits<R1,R2>::arithmetic_type> result(v.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v(i)/s;
  }
  return result;
}



template<class R> inline
LinearAlgebra::Vector<R> 
LinearAlgebra::add_approx(const Vector<R>& u, const Vector<R>& v) 
{
  ARIADNE_CHECK_EQUAL_SIZES(u,v,"Vector add_approx(Vector,Vector)");
  Vector<R> result(u.size());
  for(size_type i=0; i!=u.size(); ++i) {
    result(i)=add_approx(u(i),v(i));
  }
  return result;
}


template<class R> inline
LinearAlgebra::Vector<R> 
LinearAlgebra::sub_approx(const Vector<R>& u, const Vector<R>& v)
{
  ARIADNE_CHECK_EQUAL_SIZES(u,v,"Vector sub_approx(Vector,Vector)");
  Vector<R> result(u.size());
  for(size_type i=0; i!=u.size(); ++i) {
    result(i)=sub_approx(u(i),v(i));
  }
  return result;
}


template<class R> inline
LinearAlgebra::Vector<R> 
LinearAlgebra::mul_approx(const R& s, const Vector<R>& v)
{
  Vector<R> result(v.size());
  for(size_type i=0; i!=v.size(); ++i) {
    result(i)=mul_approx(v(i),s);
  }
  return result;
}


template<class R> inline
LinearAlgebra::Vector<R> LinearAlgebra::mul_approx(const Vector<R>& v, const R& s)
{
  Vector<R> result(v.size());
  for(size_type i=0; i!=v.size(); ++i) {
    result(i)=mul_approx(v(i),s);
  }
  return result;
}


template<class R> inline
LinearAlgebra::Vector<R> 
LinearAlgebra::div_approx(const Vector<R>& v, const R& s)
{
  Vector<R> result(v.size());
  for(size_type i=0; i!=v.size(); ++i) {
    result(i)=div_approx(v(i),s);
  }
  return result;
}



template<class R> inline 
R 
LinearAlgebra::inner_product(const Vector<R>& u, const Vector<R>& v) 
{
  ARIADNE_CHECK_EQUAL_SIZES(u,v,"Scalar inner_product(Vector,Vector)");
  R result=0;
  for(size_type i=0; i!=u.size(); ++i) {
    result+=u(i)*v(i);
  }
  return result;
}


template<class R> inline 
LinearAlgebra::Vector<R> 
LinearAlgebra::direct_sum(const Vector<R>& v1, const Vector<R>& v2) 
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
R 
LinearAlgebra::sup_norm(const Vector<R>& v) 
{
  R result=static_cast<R>(0);
  for(size_type i=0; i!= v.size(); ++i) {
    result=Numeric::max(Numeric::abs(v(i)),result); }
  return result; 
}


template<class R> inline 
R 
LinearAlgebra::norm(const Vector<R>& v) 
{
  return sup_norm(v);
}


} //namespace Ariadne
