/***************************************************************************
 *            sparse_differential.h
 *
 *  Copyright 2008  Pieter Collins
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

#ifndef ARIADNE_SPARSE_DIFFERENTIAL_H
#define ARIADNE_SPARSE_DIFFERENTIAL_H

#include <map>

#include "macros/throw.h"
#include "base/array.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "differentiation/multi_index.h"
#include "differentiation/power_series.h"
#include "differentiation/function_series.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class PowerSeries;

template<class X> struct AffineTransformation;

template<class X> class SparseDifferential;
template<class X> class SparseDifferentialVector;

template<class X> SparseDifferential<X> operator+(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> operator-(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> operator+(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
template<class X> SparseDifferential<X> operator-(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
template<class X> SparseDifferential<X> operator*(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
template<class X> SparseDifferential<X> operator/(const SparseDifferential<X>& x, const SparseDifferential<X>& y);

template<class X, class R> SparseDifferential<X> operator+(const SparseDifferential<X>& x, const R& c);
template<class X, class R> SparseDifferential<X> operator-(const SparseDifferential<X>& x, const R& c);
template<class X, class R> SparseDifferential<X> operator*(const SparseDifferential<X>& x, const R& c);
template<class X, class R> SparseDifferential<X> operator/(const SparseDifferential<X>& x, const R& c);
template<class X, class R> SparseDifferential<X> operator+(const R& c, const SparseDifferential<X>& x);
template<class X, class R> SparseDifferential<X> operator-(const R& c, const SparseDifferential<X>& x);
template<class X, class R> SparseDifferential<X> operator*(const R& c, const SparseDifferential<X>& x);
template<class X, class R> SparseDifferential<X> operator/(const R& c, const SparseDifferential<X>& x);

template<class X> SparseDifferential<X> operator+(const SparseDifferential<X>& x, const X& c);
template<class X> SparseDifferential<X> operator-(const SparseDifferential<X>& x, const X& c);
template<class X> SparseDifferential<X> operator*(const SparseDifferential<X>& x, const X& c);
template<class X> SparseDifferential<X> operator/(const SparseDifferential<X>& x, const X& c);
template<class X> SparseDifferential<X> operator+(const X& c, const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> operator-(const X& c, const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> operator*(const X& c, const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> operator/(const X& c, const SparseDifferential<X>& x);

template<class X> SparseDifferential<X> pow(const SparseDifferential<X>& x, int n);
template<class X> SparseDifferential<X> sqr(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> sqrt(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> exp(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> log(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> sin(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> cos(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> tan(const SparseDifferential<X>& x);

template<class X, class Y> Y evaluate(const SparseDifferential<X>& y, const Vector<Y>& z);
template<class X> SparseDifferential<X> embed(const SparseDifferential<X>& x, uint size, uint start);
template<class X> SparseDifferential<X> compose(const PowerSeries<X>& x, const SparseDifferential<X>& y);
template<class X> SparseDifferential<X> translate(const SparseDifferential<X>& x, const Vector<X>& c);
template<class X> SparseDifferential<X> derivative(const SparseDifferential<X>& x, uint i);
template<class X> SparseDifferential<X> antiderivative(const SparseDifferential<X>& x, uint i);


/*! \brief A class representing the derivatives of a scalar variable depending on multiple arguments. */
template<class X>
class SparseDifferential
{
 public:
  typedef MultiIndex index_type;
  typedef X value_type;
  typedef SparseDifferentialVector<X> vector_type;
  typedef typename std::map<MultiIndex,X>::iterator iterator;
  typedef typename std::map<MultiIndex,X>::const_iterator const_iterator;

  SparseDifferential() : _as(1), _deg(0), _data() { }
  //explicit SparseDifferential(int c) : _as(0), _deg(0), _data() { this->set_value(c); }
  SparseDifferential(uint as, ushort deg) : _as(as), _deg(deg), _data() { _data[MultiIndex(as)]=0; }
  template<class XX> SparseDifferential(uint as, ushort deg, const XX* ptr) : _as(as), _deg(deg), _data() { 
    for(MultiIndex j(as); j.degree()<=deg; ++j) { if(*ptr!=0) { _data[j]=*ptr; } ++ptr; } }
  
  SparseDifferential<X>& operator=(const X& c) { this->_data.clear(); this->_data[MultiIndex(this->argument_size())]=c; return *this; }

  SparseDifferential<X>& operator+=(const SparseDifferential<X>& x);
  SparseDifferential<X>& operator-=(const SparseDifferential<X>& x);
  template<class R> SparseDifferential<X>& operator+=(const R& c);
  template<class R> SparseDifferential<X>& operator-=(const R& c);
  template<class R> SparseDifferential<X>& operator*=(const R& c);
  template<class R> SparseDifferential<X>& operator/=(const R& c);

  void set_degree(ushort d) { this->_deg = d; }

  X& operator[](const uint& j) { return this->_data[MultiIndex::unit(this->_as,j)]; }
  X& operator[](const MultiIndex& a) { ARIADNE_ASSERT(a.number_of_variables()==this->argument_size()); return this->_data[a]; }
  X& value() { return this->operator[](MultiIndex(this->_as)); }
  X& gradient(uint j) { return this->operator[](MultiIndex::unit(this->_as,j)); }

  void set_value(const X& c) { this->value()=c; }
  void set_gradient(uint j, const X& d) { this->gradient(j)=d; }

  const_iterator begin() const { return this->_data.begin(); }
  const_iterator end() const { return this->_data.end(); }
  uint argument_size() const { return this->_as; }
  ushort degree() const { return this->_deg; }
  X operator[](const MultiIndex& a) const { 
    ARIADNE_ASSERT(a.number_of_variables()==this->argument_size()); 
    const_iterator iter=this->_data.find(a); 
    if(iter==this->_data.end()) { return X(0); } else { return iter->second; } }
  X value() const { return this->operator[](MultiIndex(this->_as)); }
  X gradient(uint j) const { return this->operator[](MultiIndex::unit(this->_as,j)); }

  static SparseDifferential<X> zero(uint as, ushort d) {
    return SparseDifferential<X>(as,d); }
  static SparseDifferential<X> constant(uint as, ushort d, const X& c) {
    SparseDifferential<X> r(as,d); r.value()=c; return r; }
  static SparseDifferential<X> variable(uint as, ushort d, const X& x, uint i) {
    SparseDifferential<X> r(as,d); r._data[MultiIndex::zero(as)]=x; r._data[MultiIndex::unit(as,i)]=1.0; return r; }

  bool operator==(const SparseDifferential<X>& sd) const {
    if(this->argument_size()!=sd.argument_size()) { return false; }
    for(MultiIndex j(this->argument_size()); j.degree()<=std::max(this->degree(),sd.degree()); ++j) {
      if((*this)[j]!=sd[j]) { return false; } }
    return true;
  }
  bool operator!=(const SparseDifferential<X>& sd) const { return !(*this==sd); }

  friend SparseDifferential<X> operator+<>(const SparseDifferential<X>& x);
  friend SparseDifferential<X> operator-<>(const SparseDifferential<X>& x);
  friend SparseDifferential<X> operator+<>(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
  friend SparseDifferential<X> operator-<>(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
  friend SparseDifferential<X> operator*<>(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
  friend SparseDifferential<X> operator/<>(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
  friend SparseDifferential<X> compose<>(const PowerSeries<X>& x, const SparseDifferential<X>& y);
  friend SparseDifferential<X> derivative<>(const SparseDifferential<X>& x, uint j);
  friend SparseDifferential<X> antiderivative<>(const SparseDifferential<X>& x, uint j);
 public:
  void cleanup() { 
    std::map<MultiIndex,X> _new_data;
    // Important to make sure not empty for some algorithms
    _new_data[MultiIndex::zero(this->argument_size())]=_data[MultiIndex::zero(this->argument_size())];
    for(typename std::map<MultiIndex,X>::iterator iter=this->_data.begin(); iter!=this->_data.end(); ++iter) {
      if(iter->second!=0) { _new_data[iter->first]=iter->second; } }
    std::swap(_new_data,_data); }
 private:
  uint _as;
  ushort _deg;
  std::map<MultiIndex,X> _data;
 private:
  //BOOST_CONCEPT_ASSERT((DifferentialConcept< SparseDifferential<X> >));
};

template<class X>
SparseDifferential<X>& SparseDifferential<X>::operator+=(const SparseDifferential<X>& x)
{
  for(const_iterator iter=x._data.begin(); iter!=x._data.end(); ++iter) {
    this->_data[iter->first]+=iter->second;
  }
  return *this;
}

template<class X>
SparseDifferential<X>& SparseDifferential<X>::operator-=(const SparseDifferential<X>& x)
{
  for(const_iterator iter=x._data.begin(); iter!=x._data.end(); ++iter) {
    this->_data[iter->first]-=iter->second;
  }
  return *this;
}

template<class X> template<class R>
SparseDifferential<X>& SparseDifferential<X>::operator+=(const R& c)
{
  this->_data[MultiIndex(this->_as)]+=c; return *this;
}

template<class X> template<class R>
SparseDifferential<X>& SparseDifferential<X>::operator-=(const R& c)
{
  this->_data[MultiIndex(this->_as)]-=c; return *this;
}

template<class X> template<class R>
SparseDifferential<X>& SparseDifferential<X>::operator*=(const R& c)
{
  if(c==0) {
    X zero=this->_data.begin()->second; zero*=0;
    this->_data.clear();
    this->_data[MultiIndex(this->_as)]=zero;
  } else {
    for(iterator iter=this->_data.begin(); iter!=this->_data.end(); ++iter) {
      iter->second*=c;
    }
  }
  return *this;
}


template<class X> template<class R>
SparseDifferential<X>& SparseDifferential<X>::operator/=(const R& c)
{
  for(iterator iter=this->_data.begin(); iter!=this->_data.end(); ++iter) {
    iter->second/=c;
  }
  return *this;
}

template<class X>
SparseDifferential<X> operator-(const SparseDifferential<X>& x)
{
  SparseDifferential<X> r(x.argument_size(),x.degree()); r-=x; return r; 
}


template<class X, class R>
SparseDifferential<X> operator+(const SparseDifferential<X>& x, const R& c)
{
  SparseDifferential<X> r(x); r+=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator+(const R& c, const SparseDifferential<X>& x)
{
  SparseDifferential<X> r(x); r+=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator-(const SparseDifferential<X>& x, const R& c)
{
  SparseDifferential<X> r(x); r-=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator-(const R& c, const SparseDifferential<X>& x)
{
  SparseDifferential<X> r(-x); r+=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator*(const SparseDifferential<X>& x, const R& c)
{
  SparseDifferential<X> r(x); r*=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator*(const R& c, const SparseDifferential<X>& x)
{
  SparseDifferential<X> r(x); r*=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator/(const SparseDifferential<X>& x, const R& c)
{
  SparseDifferential<X> r(x); r/=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator/(const R& c, const SparseDifferential<X>& x)
{
  SparseDifferential<X> r=rec(x); r*=c; return r; 
}


template<class X>
SparseDifferential<X> operator+(const SparseDifferential<X>& x, const X& c)
{
  SparseDifferential<X> r(x); r+=c; return r; 
}

template<class X>
SparseDifferential<X> operator+(const X& c, const SparseDifferential<X>& x)
{
  SparseDifferential<X> r(x); r+=c; return r; 
}

template<class X>
SparseDifferential<X> operator-(const SparseDifferential<X>& x, const X& c)
{
  SparseDifferential<X> r(x); r-=c; return r; 
}

template<class X>
SparseDifferential<X> operator-(const X& c, const SparseDifferential<X>& x)
{
  SparseDifferential<X> r(-x); r+=c; return r; 
}

template<class X>
SparseDifferential<X> operator*(const SparseDifferential<X>& x, const X& c)
{
  SparseDifferential<X> r(x); r*=c; return r; 
}

template<class X>
SparseDifferential<X> operator*(const X& c, const SparseDifferential<X>& x)
{
  SparseDifferential<X> r(x); r*=c; return r; 
}

template<class X>
SparseDifferential<X> operator/(const SparseDifferential<X>& x, const X& c)
{
  SparseDifferential<X> r(x); r/=c; return r; 
}

template<class X>
SparseDifferential<X> operator/(const X& c, const SparseDifferential<X>& x)
{
  SparseDifferential<X> r=rec(x); r*=c; return r; 
}






template<class X>
SparseDifferential<X> operator+(const SparseDifferential<X>& x, const SparseDifferential<X>& y)
{
  SparseDifferential<X> r(x); r+=y; r.cleanup(); return r; 
}

template<class X>
SparseDifferential<X> operator-(const SparseDifferential<X>& x, const SparseDifferential<X>& y)
{
  SparseDifferential<X> r(x); r-=y; r.cleanup(); return r; 
}

template<class X>
SparseDifferential<X> operator*(const SparseDifferential<X>& x, const SparseDifferential<X>& y)
{
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  typedef typename SparseDifferential<X>::const_iterator const_iterator;
  assert(x._as==y._as);
  SparseDifferential<X> r(x._as,std::min(x._deg,y._deg));
  for(const_iterator xiter=x._data.begin(); xiter!=x._data.end() && xiter->first.degree()<=r.degree(); ++xiter) {
    for(const_iterator yiter=y._data.begin(); yiter!=y._data.end() && xiter->first.degree()+yiter->first.degree()<=r.degree(); ++yiter) {
      r._data[xiter->first+yiter->first]+=(xiter->second*yiter->second);
    }
  }
  r.cleanup();
  return r;
}

template<class X>
SparseDifferential<X> operator/(const SparseDifferential<X>& x, const SparseDifferential<X>& y)
{
  return x*rec(y);
}

template<class X>
SparseDifferential<X> rec(const SparseDifferential<X>& x)
{
  return compose(FunctionSeries<X>::rec(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> sqr(const SparseDifferential<X>& x)
{
  return pow(x,2);
}

template<class X>
SparseDifferential<X> pow(const SparseDifferential<X>& x, int n)
{
  return compose(FunctionSeries<X>::pow(x.degree(),x.value(),n),x);
}

template<class X>
SparseDifferential<X> sqrt(const SparseDifferential<X>& x)
{
  return compose(FunctionSeries<X>::sqrt(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> exp(const SparseDifferential<X>& x)
{
  return compose(FunctionSeries<X>::exp(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> log(const SparseDifferential<X>& x)
{
  return compose(FunctionSeries<X>::log(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> sin(const SparseDifferential<X>& x)
{
  return compose(FunctionSeries<X>::sin(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> cos(const SparseDifferential<X>& x)
{
  return compose(FunctionSeries<X>::cos(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> tan(const SparseDifferential<X>& x)
{
  return compose(FunctionSeries<X>::tan(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> compose(const PowerSeries<X>& x, const SparseDifferential<X>& y)
{
  uint as=y.argument_size();
  ushort d=std::min(x.degree(),y.degree());

  SparseDifferential<X> w=y;
  w[MultiIndex(as)]=0;
  SparseDifferential<X> r(as,d);
  r[MultiIndex(as)]=x[d];
  for(ushort n=1; n<=d; ++n) {
    r=r*w; r+=x[d-n];
  }
  return r;
}


template<class X>
SparseDifferential<X> 
embed(const SparseDifferential<X>& x, 
      uint size, uint start)
{  
  assert(start+x.argument_size()<=size);
  SparseDifferential<X> r(size,x.degree());
  MultiIndex jr(size);
  for(typename SparseDifferential<X>::const_iterator iter=x.begin();
      iter!=x.end(); ++iter)
  {
    const MultiIndex& jx=iter->first;
    for(uint k=0; k!=x.argument_size(); ++k) {
      jr.set(start+k,jx[k]);
    }
    r[jr]=iter->second;
  }
  return r;
}


template<class X>
SparseDifferential<X> derivative(const SparseDifferential<X>& x, uint j)
{
  if(x.degree()==0) { return SparseDifferential<X>(x.argument_size(),0u); }
  SparseDifferential<X> r(x.argument_size(), x.degree()-1); 
  MultiIndex da=MultiIndex::unit(x.argument_size(),j); 
  MultiIndex ej=MultiIndex::unit(x.argument_size(),j);
  for(typename SparseDifferential<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) 
  {
    const MultiIndex& a=iter->first;
    if(a[j]!=0) { 
      da=a-ej;
      int aj=a[j];
      const X& xa=x[a];
      r[da]=xa*aj;
    }
  }
  return r;
}

template<class X>
SparseDifferential<X> antiderivative(const SparseDifferential<X>& x, uint j)
{
  SparseDifferential<X> r(x.argument_size(), x.degree()+1); 
  MultiIndex da=MultiIndex::zero(x.argument_size()); 
  MultiIndex aj=MultiIndex::unit(x.argument_size(),j);
  for(typename SparseDifferential<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) 
  {
    const MultiIndex& a=iter->first;
    const X& xj=x[a];
    da=a+aj;
    uint daj=da[j]; r[da]=xj/daj;
    //r[da]=x[a]/da[i];
  }
  return r;
}


//! Translate the polynomial given by \a x to one with centre \a v.
template<class X> 
SparseDifferential<X>  
translate(const SparseDifferential<X>& x, const Vector<X>& v)
{
  ARIADNE_ASSERT(x.argument_size()==v.size());
  uint as=x.argument_size();
  ushort d=x.degree();
  SparseDifferentialVector<X> t=SparseDifferentialVector<X>::variable(as,as,d,v);
  return evaluate(x,t);
}

template<class X>
std::ostream& operator<<(std::ostream& os, const SparseDifferential<X>& x)
{
  for(typename SparseDifferential<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) {
    if(iter==x.begin()) { os << "D{"; } else { os << ","; }
    os << iter->first << ":" << iter->second ;
  }
  return os << "}";
}




} //namespace Ariadne

#endif /* ARIADNE_SPARSE_DIFFERENTIAL_H */
