#ifndef ARIADNE_POINT_H
#define ARIADNE_POINT_H

#include "numeric.h"
#include "vector.h"

namespace Ariadne {

class Point : public Vector<Float>
{
 public:
  typedef Float real_type;
  Point() : Vector<Float>() { }
  template<class T> Point(const T& t) : Vector<Float>(t) { }
  template<class T1, class T2> Point(const T1& t1, const T2& t2) : Vector<Float>(t1,t2) { }
  uint dimension() const { return this->size(); }
  Vector<Float> centre() const { return *this; }
};

} // namespace Ariadne

#endif /* ARIADNE_POINT_H */
