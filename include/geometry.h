#ifndef ARIADNE_GEOMETRY_H
#define ARIANDE_GEOMETRY_H

#include <vector>

#include "box.h"
#include "function_interface.h"

namespace Ariadne {

template<class BS>
class ListSet
  : public std::vector<BS>
{
  void adjoin(const BS& bs) { 
    this->push_back(bs); }
  void adjoin(const std::vector<BS>& bs) {
    this->insert(this->end(),bs.begin(),bs.end()); }
};

inline
uint irmax(const Box& bx) {
  Float dmax=0.0;
  uint imax=0;
  for(uint i=0; i!=bx.size(); ++i) {
    Float d=bx[i].width();
    if(d>dmax) {
      imax=i;
      dmax=d;
    }
  }
  return imax;
}
 
inline
Box split(const Box& bx, uint i, bool lr) {
  Box result(bx);
  Float c=med_approx(bx[i].lower(),bx[i].upper());
  if(lr==false) { result[i].u=c; }
  else { result[i].l=c; }
  return result;
}

inline
std::pair<Box,Box> split(const Box& bx, uint i)
{
  std::pair<Box,Box> result(bx,bx);
  Float c=med_approx(bx[i].lower(),bx[i].upper());
  result.first[i].u=c;
  result.second[i].l=c;
  return result;
}

inline
std::pair<Box,Box> split(const Box& bx) {
  return split(irmax(bx));
}


template<class F>
tribool 
disjoint(const Box& d, const F& f, const Box& b, const Float& eps)
{
  Box fd=f(d);
  if(disjoint(fd,b)) { 
    return true;
  } else if(subset(fd,b)) {
    return false;
  } else if(d.radius()<eps) {
    return indeterminate;
  } else {
    uint i=irmax(d);
    return disjoint(split(d,i,0),f,b,eps) || disjoint(split(d,i,0),f,b,eps);
  }
}

template<class F>
tribool 
subset(const Box& d, const F& f, const Box& b, const Float& eps)
{
  Box fd=f(d);
  if(subset(fd,b)) { 
    return true;
  } else if(disjoint(fd,b)) {
    return false;
  } else if(d.radius()<eps) {
    return indeterminate;
  } else {
    uint i=irmax(d);
    return subset(split(d,i,0),f,b,eps) && subset(split(d,i,0),f,b,eps);
  }
}

template<class DS>
DS remove_subsets(const DS& ls)
{
  DS result;
  for(uint i=0; i!=ls.size(); ++i) {
    for(uint j=0; j!=ls.size(); ++j) {
      if(subset(ls[i],ls[j])) {
        break; 
      }
    }
    result.adjoin(ls[i]);
  }
}

template<class DS>
DS remove_supersets(const DS& ls)
{
  DS result;
  for(uint i=0; i!=ls.size(); ++i) {
    for(uint j=0; j!=ls.size(); ++j) {
      if(subset(ls[j],ls[i])) {
        break; 
      }
    }
    result.adjoin(ls[i]);
  }
}

} // namespace Ariadne

#endif // ARIADNE_GEOMETRY_H
