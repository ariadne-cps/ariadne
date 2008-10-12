#ifndef ARIADNE_VECTOR_H
#define ARIADNE_VECTOR_H 

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "numeric.h"

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
  Vector(int n)
    : ublas::vector<X>(n) { for(uint i=0; i!=this->size(); ++i) { (*this)[i]=0; } }
  Vector(uint n)
    : ublas::vector<X>(n) { for(uint i=0; i!=this->size(); ++i) { (*this)[i]=0; } }
  Vector(uint n, const X& t)
    : ublas::vector<X>(n) { for(uint i=0; i!=this->size(); ++i) { (*this)[i]=t; } }
  template<class XX> Vector(uint n, const XX* ptr)
    : ublas::vector<X>(n) { for(uint i=0; i!=this->size(); ++i) { (*this)[i]=ptr[i]; } }
  template<class XX> Vector(const Vector<XX>& v)
    : ublas::vector<X>(v) { }
  template<class E> Vector(const E& ve)
    : ublas::vector<X>(ve) { }
  uint dimension() const { return this->size(); }
  const X& get(uint i) const { return (*this)[i]; }
  template<class T> void set(uint i, const T& x) { (*this)[i] = x; }
};

typedef ublas::slice Slice;
typedef ublas::range Range;
using ublas::range;
using ublas::slice;
using ublas::project;
using ublas::identity_matrix;


template<class X>
X sup_norm(const Vector<X>& v)
{
  X r=0;
  for(uint i=0; i!=v.size(); ++i) {
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
  uint n1=v1.size();
  uint n2=v2.size();
  Vector<X> r(n1+n2);
  project(r,range(0,n1))=v1;
  project(r,range(n1,n1+n2))=v2;
  return r;
}

template<class X>
Vector<X> join(const Vector<X>& v1, const X& s2) 
{
  uint n1=v1.size();
  Vector<X> r(n1+1);
  project(r,range(0,n1))=v1;
  r[n1]=s2;
  return r;
}



template<class X1, class X2>
bool operator==(const Vector<X1>& v1, const Vector<X2>& v2)
{
  if(v1.size()!=v2.size()) { return false; }
  for(uint i=0; i!=v1.size(); ++i) {
    if(v1[i]!=v2[i]) { return false; }
  }
  return true;
}


template<class X> std::ostream& operator<<(std::ostream& os, const Vector<X>& v) {
  if(v.size()==0) { os << '['; }
  for(uint i=0; i!=v.size(); ++i) { 
    os << (i==0 ? '[' : ',') << v[i]; }
  return os << ']';
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

template<class X, class XX> inline Vector<X> operator*(const XX& s, const Vector<X>& v) { return v*s; }


template<class X, class XX> inline Vector<X> vector(uint d, const XX* ptr) {
  return Vector<Float>(d,ptr); }
inline Vector<Float> point(uint d, Float* ptr) { 
  return Vector<Float>(d,ptr); }
inline Vector<Interval> box(uint d, Float* ptr) { 
  Vector<Interval> bx(d); 
  for(uint i=0; i!=d; ++i) { 
    bx[i]=Interval(ptr[2*i],ptr[2*i+1]); }
  return bx;
}

} // namespace Ariadne

#endif
