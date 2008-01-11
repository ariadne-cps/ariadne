/***************************************************************************
 *            covector.inline.h
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
LinearAlgebra::Covector<R>::Covector()
  : _data(0) 
{ 
}

template<class R> inline
LinearAlgebra::Covector<R>::Covector(const size_type& n)
  : _data(n,static_cast<R>(0)) 
{ 
}

template<class R> inline
LinearAlgebra::Covector<R>::Covector(const size_type& n, const R& x)
  : _data(n) 
{ 
  for(size_type i=0; i!=n; ++i) { (*this)(i)=x; }
}


template<class R> template<class RR> inline
LinearAlgebra::Covector<R>::Covector(const array<RR>& ary)
  : _data(ary) 
{ 
}

template<class R> template<class RR> inline
LinearAlgebra::Covector<R>::Covector(const size_type& n, const RR* ptr, const size_type& inc)
  : _data(n) 
{ 
  for(size_type i=0; i!=n; ++i) { 
    (*this)(i)=ptr[i*inc]; 
  }
}


template<class R> template<class RR> inline
LinearAlgebra::Covector<R>::Covector(const Covector<RR>& v) 
  : _data(v.data()) 
{ 
}



template<class R> inline
LinearAlgebra::Covector<R>::Covector(const Covector<R>& v) 
  : _data(v._data) 
{ 
}



template<class R> inline
LinearAlgebra::Covector<R>& 
LinearAlgebra::Covector<R>::operator=(const Covector<R>& v) 
{
  if(this!=&v) { this->_data=v._data; } return *this; 
}


template<class R> template<class RR> inline
LinearAlgebra::Covector<R>& 
LinearAlgebra::Covector<R>::operator=(const Covector<RR>& v) 
{
  this->_data=v.data(); return *this; 
}



template<class R> inline
bool 
LinearAlgebra::Covector<R>::operator==(const Covector<R>& v) const 
{
  if(this->size()!=v.size()) { return false; }
  for(size_type i=0; i!=this->size(); ++i) { 
    if((*this)(i)!=v(i)) { return false; } 
  }
  return true; 
}


template<class R> inline
bool 
LinearAlgebra::Covector<R>::operator!=(const Covector<R>& v) const 
{
  return !(*this==v); 
}


template<class R> inline
LinearAlgebra::Covector<R> 
LinearAlgebra::Covector<R>::zero(const size_type& n) 
{
  Covector<R> v(n); for (dimension_type i=0; i<n; ++i) { v(i)=R(0); } return v; 
}


template<class R> inline
LinearAlgebra::Covector<R> 
LinearAlgebra::Covector<R>::one(const size_type& n) 
{
  Covector<R> v(n); for (dimension_type i=0; i<n; ++i) { v(i)=R(1); } return v; 
}


template<class R> inline
LinearAlgebra::Covector<R> 
LinearAlgebra::Covector<R>::unit(dimension_type n, dimension_type i) 
{
  Covector<R> v(n); for (dimension_type k=0; k<n; ++k) { v(k)=R(0); } v(i)=R(1); return v; }


template<class R> inline
array<R>& 
LinearAlgebra::Covector<R>::data() 
{
  return this->_data; 
}


template<class R> inline
const array<R>&
LinearAlgebra::Covector<R>::data() const 
{
  return this->_data; 
}


template<class R> inline
void 
LinearAlgebra::Covector<R>::resize(const size_type& n) 
{
  this->_data.resize(n); 
}


template<class R> inline
size_type 
LinearAlgebra::Covector<R>::size() const 
{
  return this->_data.size(); 
}


template<class R> inline
R* 
LinearAlgebra::Covector<R>::begin() 
{
  return this->data().begin(); 
}

template<class R> inline
const R* 
LinearAlgebra::Covector<R>::begin() const 
{
  return this->data().begin(); 
}

template<class R> inline
R* 
LinearAlgebra::Covector<R>::end() 
{ 
  return this->data().end(); 
}

template<class R> inline
const R* 
LinearAlgebra::Covector<R>::end() const 
{
  return this->data().end(); 
}


template<class R> inline
size_type 
LinearAlgebra::Covector<R>::increment() const 
{
  return 1; 
}


template<class R> inline
R& 
LinearAlgebra::Covector<R>::operator() (const size_type& i) 
{
  return this->_data[i]; 
}

template<class R> inline
const R& 
LinearAlgebra::Covector<R>::operator() (const size_type& i) const 
{
  return this->_data[i]; 
}


template<class R> inline
R& 
LinearAlgebra::Covector<R>::operator[] (const size_type& i) 
{ 
  return this->_data[i]; 
}

template<class R> inline
const R& 
LinearAlgebra::Covector<R>::operator[] (const size_type& i) const 
{
  return this->_data[i]; 
}







template<class R1, class R2> inline
LinearAlgebra::Covector<R1>&
LinearAlgebra::operator+=(Covector<R1>& v1, const Covector<R2>& v2) 
{
  ARIADNE_CHECK_EQUAL_SIZES(v1,v2,"Covector& operator+=(Covector,Covector)");
  for(size_type i=0; i!=v1.size(); ++i) { v1(i)+=v2(i); } return v1; 
}

template<class R1, class R2> inline
LinearAlgebra::Covector<R1>&
LinearAlgebra::operator-=(Covector<R1>& v1, const Covector<R2>& v2) 
{
  ARIADNE_CHECK_SIZE(v1,v2.size(),"Covector& operator-=(Covector,Covector)");
  for(size_type i=0; i!=v1.size(); ++i) { v1(i)=v1(i)-=v2(i); } return v1; 
}

template<class R1, class R2> inline
LinearAlgebra::Covector<R1>&
LinearAlgebra::operator*=(Covector<R1>& v1, const R2& s2) 
{
  for(size_type i=0; i!=v1.size(); ++i) { v1(i)*=s2; } return v1; 
}

template<class R1, class R2> inline
LinearAlgebra::Covector<R1>&
LinearAlgebra::operator/=(Covector<R1>& v1, const R2& s2) {
  for(size_type i=0; i!=v1.size(); ++i) { v1(i)/=s2; } return v1; 
}





template<class R> inline 
LinearAlgebra::Covector<R>
LinearAlgebra::zero_covector(size_type n)
{
  return Covector<R>(n);
}


template<class R> inline 
LinearAlgebra::Covector<R>
LinearAlgebra::unit_covector(size_type n, size_type i)
{
  Covector<R> result(n);
  result[i]=1;
  return result;
}


template<class R> inline 
LinearAlgebra::Covector<R>
LinearAlgebra::midpoint(const Covector< Numeric::Interval<R> >& iv) 
{
  Covector<R> result(iv.size());
  for(size_type i=0; i!=iv.size(); ++i) {
    result(i) = midpoint(iv(i));
  }
  return result;
}



template<class R> inline
bool
LinearAlgebra::encloses(const Covector< Numeric::Interval<R> >& iv, const Covector<R>& v) 
{
  ARIADNE_CHECK_EQUAL_SIZES(iv,v,"bool contains_value(Covector<Interval>,Covector<Float>)");
  for(size_type i=0; i!=v.size(); ++i) {
    if(!Numeric::encloses(iv(i),v(i))) {
      return false;
    }
  }
  return true;
}


template<class R> inline
bool
LinearAlgebra::refines(const Covector< Numeric::Interval<R> >& iv1, const  Covector< Numeric::Interval<R> >& iv2) 
{
  ARIADNE_CHECK_EQUAL_SIZES(iv1,iv2,"bool refines(Covector<Interval>,Covector<Interval>)");
  for(size_type i=0; i!=iv1.size(); ++i) {
    if(!Numeric::refines(iv1(i),iv2(i))) {
      return false;
    }
  }
  return true;
}


template<class R1,class R2> inline 
LinearAlgebra::Covector<R1>
LinearAlgebra::approximation(const Covector<R2>& iv) 
{
  Covector<R1> result(iv.size());
  for(size_type i=0; i!=iv.size(); ++i) {
    set_(result(i),iv(i),Numeric::round_approx);
  }
  return result;
}








template<class R> inline
LinearAlgebra::Covector<R> 
LinearAlgebra::operator-(const Covector<R>& v) 
{
  Covector<R> result(v.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=-v(i);
  }
  return result;
}


template<class R1, class R2> inline
LinearAlgebra::Covector<typename Numeric::traits<R1,R2>::arithmetic_type> 
LinearAlgebra::operator+(const Covector<R1>& v1, const Covector<R2>& v2) 
{
  ARIADNE_CHECK_EQUAL_SIZES(v1,v2,"Covector operator+(Covector,Covector)");
  Covector<typename Numeric::traits<R1,R2>::arithmetic_type> result(v1.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v1(i)+v2(i);
  }
  return result;
}


template<class R1, class R2> inline
LinearAlgebra::Covector<class Numeric::traits<R1,R2>::arithmetic_type> 
LinearAlgebra::operator-(const Covector<R1>& v1, const Covector<R2>& v2) 
{
  ARIADNE_CHECK_EQUAL_SIZES(v1,v2,"Covector operator-(Covector,Covector)");
  Covector<typename Numeric::traits<R1,R2>::arithmetic_type> result(v1.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v1(i)-v2(i);
  }
  return result;
}


template<class R1, class R2> inline
LinearAlgebra::Covector<typename Numeric::traits<R1,R2>::arithmetic_type> 
LinearAlgebra::operator*(const R1& s, const Covector<R2>& v) 
{
  return v*s;
}


template<class R1, class R2> inline
LinearAlgebra::Covector<typename Numeric::traits<R1,R2>::arithmetic_type> 
LinearAlgebra::operator*(const Covector<R1>& v, const R2& s) 
{
  Covector<typename Numeric::traits<R1,R2>::arithmetic_type> result(v.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v(i)*s;
  }
  return result;
}


template<class R1, class R2> inline
LinearAlgebra::Covector<typename Numeric::traits<R1,R2>::arithmetic_type> 
LinearAlgebra::operator/(const Covector<R1>& v, const R2& s) 
{
  Covector<typename Numeric::traits<R1,R2>::arithmetic_type> result(v.size());
  for(size_type i=0; i!=result.size(); ++i) {
    result(i)=v(i)/s;
  }
  return result;
}



template<class R> inline
LinearAlgebra::Covector<R> 
LinearAlgebra::add_approx(const Covector<R>& u, const Covector<R>& v) 
{
  ARIADNE_CHECK_EQUAL_SIZES(u,v,"Covector add_approx(Covector,Covector)");
  Covector<R> result(u.size());
  for(size_type i=0; i!=u.size(); ++i) {
    result(i)=add_approx(u(i),v(i));
  }
  return result;
}


template<class R> inline
LinearAlgebra::Covector<R> 
LinearAlgebra::sub_approx(const Covector<R>& u, const Covector<R>& v)
{
  ARIADNE_CHECK_EQUAL_SIZES(u,v,"Covector sub_approx(Covector,Covector)");
  Covector<R> result(u.size());
  for(size_type i=0; i!=u.size(); ++i) {
    result(i)=sub_approx(u(i),v(i));
  }
  return result;
}


template<class R> inline
LinearAlgebra::Covector<R> 
LinearAlgebra::mul_approx(const R& s, const Covector<R>& v)
{
  Covector<R> result(v.size());
  for(size_type i=0; i!=v.size(); ++i) {
    result(i)=mul_approx(v(i),s);
  }
  return result;
}


template<class R> inline
LinearAlgebra::Covector<R> LinearAlgebra::mul_approx(const Covector<R>& v, const R& s)
{
  Covector<R> result(v.size());
  for(size_type i=0; i!=v.size(); ++i) {
    result(i)=mul_approx(v(i),s);
  }
  return result;
}


template<class R> inline
LinearAlgebra::Covector<R> 
LinearAlgebra::div_approx(const Covector<R>& v, const R& s)
{
  Covector<R> result(v.size());
  for(size_type i=0; i!=v.size(); ++i) {
    result(i)=div_approx(v(i),s);
  }
  return result;
}





template<class R> inline 
LinearAlgebra::Covector<R> 
LinearAlgebra::direct_sum(const Covector<R>& v1, const Covector<R>& v2) 
{
  return concatenate(v1,v2);
}


template<class R> inline 
LinearAlgebra::Covector<R> 
LinearAlgebra::concatenate(const Covector<R>& v1, const Covector<R>& v2) 
{
  Covector<R> result(v1.size()+v2.size());
  for(size_type i=0; i!=v1.size(); ++i) {
    result(i)=v1(i);
  }
  for(size_type i=0; i!=v2.size(); ++i) {
    result(i+v1.size())=v2(i);
  }
  return result;
}

template<class R> inline 
LinearAlgebra::Covector<R> 
LinearAlgebra::concatenate(const Covector<R>& v, const R& s) 
{
  Covector<R> result(v.size()+1);
  for(size_type i=0; i!=v.size(); ++i) {
    result(i)=v(i);
  }
  result(v.size())=s;
  return result;
}


template<class R> inline
R 
LinearAlgebra::sup_norm(const Covector<R>& v) 
{
  R result=static_cast<R>(0);
  for(size_type i=0; i!= v.size(); ++i) {
    result=Numeric::max(Numeric::abs(v(i)),result); }
  return result; 
}


template<class R> inline 
R 
LinearAlgebra::norm(const Covector<R>& v) 
{
  return sup_norm(v);
}


} //namespace Ariadne
