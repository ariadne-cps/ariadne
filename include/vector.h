/***************************************************************************
 *            vector.h
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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
 
/*! \file vector.h
 *  \brief Vectors in Euclidean space.
 */

#ifndef ARIADNE_VECTOR_H
#define ARIADNE_VECTOR_H 

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <string>
#include <sstream>

#include "macros.h"
#include "numeric.h"
#include "stlio.h"

using namespace boost::numeric;

namespace Ariadne {

/// A vector over a field.
template<class X>
class Vector
  : public boost::numeric::ublas::vector<X>
{
 public:
  Vector()
    : ublas::vector<X>() { }
  Vector(size_t n)
    : ublas::vector<X>(n) { for(size_t i=0; i!=this->size(); ++i) { (*this)[i]=0; } }
  Vector(size_t n, const X& t)
    : ublas::vector<X>(n) { for(size_t i=0; i!=this->size(); ++i) { (*this)[i]=t; } }
  template<class XX> Vector(size_t n, const XX* ptr)
    : ublas::vector<X>(n) { for(size_t i=0; i!=this->size(); ++i) { (*this)[i]=ptr[i]; } }
  template<class XX> Vector(const Vector<XX>& v)
    : ublas::vector<X>(v) { }
  template<class E> Vector(const ublas::vector_expression<E> &ve)
    : ublas::vector<X>(ve) { }
  template<class E> Vector<X>& operator=(const ublas::vector_expression<E> &ve) { 
    this->ublas::vector<X>::operator=(ve); return *this; }
  static Vector<X> zero(size_t n) { return Vector<Float>(n,0.0); }
  static Vector<X> one(size_t n) { return Vector<Float>(n,1.0); }
  static Vector<X> unit(size_t n,size_t i) { 
    ARIADNE_ASSERT(i<n); Vector<Float> result(n,0.0); result[i]=1.0; return result; }
  const X& get(size_t i) const { return (*this)[i]; }
  template<class T> void set(size_t i, const T& x) { (*this)[i] = x; }
};

template<class X> Vector<X> make_vector(const std::string&);

typedef ublas::slice Slice;
typedef ublas::range Range;
using ublas::range;
using ublas::slice;
using ublas::identity_matrix;


template<class X>
X sup_norm(const Vector<X>& v)
{
  X r=0;
  for(size_t i=0; i!=v.size(); ++i) {
    r=max(r,v[i]);
  }
  return r;
}

template<class X>
X norm(const Vector<X>& v)
{
  return Ariadne::sup_norm(v);
}

template<class X>
Vector<X> join(const Vector<X>& v1, const Vector<X>& v2) 
{
  size_t n1=v1.size();
  size_t n2=v2.size();
  Vector<X> r(n1+n2);
  ublas::project(r,range(0,n1))=v1;
  ublas::project(r,range(n1,n1+n2))=v2;
  return r;
}

template<class X>
Vector<X> join(const Vector<X>& v1, const X& s2) 
{
  size_t n1=v1.size();
  Vector<X> r(n1+1);
  ublas::project(r,range(0,n1))=v1;
  r[n1]=s2;
  return r;
}

  

template<class X1, class X2>
bool operator==(const Vector<X1>& v1, const Vector<X2>& v2)
{
  if(v1.size()!=v2.size()) { return false; }
  for(size_t i=0; i!=v1.size(); ++i) {
    if(v1[i]!=v2[i]) { return false; }
  }
  return true;
}


template<class X> std::ostream& operator<<(std::ostream& os, const Vector<X>& v) {
  if(v.size()==0) { os << '['; }
  for(size_t i=0; i!=v.size(); ++i) { 
    os << (i==0 ? '[' : ',') << v[i]; }
  return os << ']';
}

template<class X> std::istream& operator>>(std::istream& is, Vector<X>& v) {
  std::vector<X> vec;
  is >> vec;
  X* ptr=&vec[0];
  v=Vector<X>(vec.size(),ptr);
  return is;
}

bool subset(const Vector<Float>& v1, const Vector<Interval>& v2);
bool subset(const Vector<Interval>& v1, const Vector<Interval>& v2);
bool disjoint(const Vector<Interval>& v1, const Vector<Interval>& v2);

Vector<Interval> hull(const Vector<Interval>& v1, const Vector<Interval>& v2);
Vector<Interval> intersection(const Vector<Interval>& v1, const Vector<Interval>& v2);
Vector<Float> midpoint(const Vector<Interval>& v);
Vector<Float> lower(const Vector<Interval>& v);
Vector<Float> upper(const Vector<Interval>& v);
Float radius(const Vector<Interval>& z);
Float volume(const Vector<Interval>& z);

inline Vector<Float> add_approx(const Vector<Float>& v1, const Vector<Float>& v2) { return v1+v2; }
inline Vector<Float> sub_approx(const Vector<Float>& v1, const Vector<Float>& v2) { return v1-v2; }



template<class X, class XX> inline Vector<X> vector(size_t d, const XX* ptr) {
  return Vector<Float>(d,ptr); }
inline Vector<Float> point(size_t d, Float* ptr) { 
  return Vector<Float>(d,ptr); }
inline Vector<Interval> box(size_t d, Float* ptr) { 
  Vector<Interval> bx(d); 
  for(size_t i=0; i!=d; ++i) { 
    bx[i]=Interval(ptr[2*i],ptr[2*i+1]); }
  return bx;
}

template<class X> Vector<X> 
make_vector(const std::string& str)
{
  Vector<X> res;
  std::stringstream ss(str);
  ss >> res;
  return res;
}



} // namespace Ariadne

#endif // ARIADNE_VECTOR_H

