#ifndef ARIADNE_FUNCTION_INTERFACE_H
#define ARIADNE_FUNCTION_INTERFACE_H

#include <iosfwd>
#include <iostream>
#include "numeric.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class SparseDifferentialVector;

//! \brief Interface for functions whose derivatives can be computed.
class FunctionInterface {
 public:
  virtual ~FunctionInterface() { };
  virtual FunctionInterface* clone() const = 0;
     
  virtual ushort smoothness() const = 0;
  virtual uint argument_size() const = 0;
  virtual uint result_size() const = 0;

  virtual Vector<Interval> evaluate(const Vector<Interval>& x) const = 0;
  virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const = 0;
  virtual SparseDifferentialVector<Float> expansion(const Vector<Float>& x, const ushort& s) const = 0;
  virtual SparseDifferentialVector<Interval> expansion(const Vector<Interval>& x, const ushort& s) const = 0;
  
  virtual std::ostream& write(std::ostream& os) const = 0;
};

inline std::ostream& operator<<(std::ostream& os, const FunctionInterface& f) {
  return f.write(os); 
}


} // namespace Ariadne

#endif
