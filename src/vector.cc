#include "macros.h"
#include "numeric.h"
#include "vector.h"

template class boost::numeric::ublas::vector<Ariadne::Float>;
template class boost::numeric::ublas::vector<Ariadne::Interval>;

namespace Ariadne {

bool subset(const Vector<Float>& v1, const Vector<Interval>& v2)
{
  ARIADNE_ASSERT(v1.size()==v2.size());
  for(uint i=0; i!=v1.size(); ++i) {
    if(!subset(v1[i],v2[i])) { return false; }
  }
  return true;
}

bool subset(const Vector<Interval>& v1, const Vector<Interval>& v2) 
{
  ARIADNE_ASSERT(v1.size()==v2.size());
  for(uint i=0; i!=v1.size(); ++i) {
    if(!subset(v1[i],v2[i])) { return false; }
  }
  return true;
}

Vector<Float> midpoint(const Vector<Interval>& v) 
{
  Vector<Float> r(v.size());
  for(uint i=0; i!=v.size(); ++i) {
    r[i]=v[i].midpoint();
  }
  return r;
}

Vector<Float> lower(const Vector<Interval>& v) 
{
  Vector<Float> r(v.size());
  for(uint i=0; i!=v.size(); ++i) {
    r[i]=v[i].lower();
  }
  return r;
}

Vector<Float> upper(const Vector<Interval>& v) 
{
  Vector<Float> r(v.size());
  for(uint i=0; i!=v.size(); ++i) {
    r[i]=v[i].upper();
  }
  return r;
}

Vector<Interval> intersection(const Vector<Interval>& v1, const Vector<Interval>& v2) 
{
  ARIADNE_ASSERT(v1.size()==v2.size());
  Vector<Interval> r(v1.size());
  for(uint i=0; i!=v1.size(); ++i) {
    r[i]=intersection(v1[i],v2[i]);
  }
  return r;
}

Float radius(const Vector<Interval>& v) 
{
  Float r=0;
  for(uint i=0; i!=v.size(); ++i) {
    r=Ariadne::max(r,v[i].radius());
  }
  return r;
}

bool disjoint(const Vector<Interval>& v1, const Vector<Interval>& v2) 
{
  ARIADNE_ASSERT(v1.size()==v2.size());
  for(uint i=0; i!=v1.size(); ++i) {
    if(v1[i].u<v2[i].l || v1[i].l>v2[i].u) {
      return true;
    }
  }
  return false;
}

Vector<Interval> intersection(const Vector<Interval>& v1, const Vector<Interval>& v2);
Vector<Float> midpoint(const Vector<Interval>& v);

} // namespace Ariadne
