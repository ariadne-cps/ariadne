#ifndef ARIADNE_BOX_H
#define ARIADNE_BOX_H

#include "numeric.h"
#include "vector.h"

namespace Ariadne {

class Box : public Vector<Interval>
{
 public:
  typedef Float real_type;
  Box() : Vector<Interval>() { }
  template<class T> Box(const T& t) : Vector<Interval>(t) { }
  template<class T1, class T2> Box(const T1& t1, const T2& t2) : Vector<Interval>(t1,t2) { }
  uint dimension() const { return this->size(); }
  Vector<Float> centre() const { return midpoint(*this); }
};

} // namespace Ariadne

#endif /* ARIADNE_BOX_H */
