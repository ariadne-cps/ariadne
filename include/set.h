#ifndef ARIADNE_SET_H
#define ARIADNE_SET_H

#include <iosfwd>

#include "boost/shared_ptr.hpp"
#include "macros.h"
#include "numeric.h"
#include "vector.h"
#include "set_interface.h"
#include "function_interface.h"


namespace Ariadne {

class ImageSet
  : public LocatedSetInterface
{
  Vector<Interval> _domain;
  boost::shared_ptr<FunctionInterface> _function;
 public:
  typedef Vector<Interval> BoxType;
  ImageSet(const BoxType& dom, const FunctionInterface& fn);
  const BoxType& domain() { return this->_domain; }
  const FunctionInterface& function() { return *this->_function; }
};

class ConstraintSet
  : public RegularSetInterface
{
  Vector<Interval> _codomain;
  boost::shared_ptr<FunctionInterface> _function;
 public:
  typedef Vector<Interval> BoxType;
  ConstraintSet(const BoxType& codom, const FunctionInterface& fn);
  const BoxType& codomain() { return this->_codomain; }
  const FunctionInterface& function() { return *this->_function; };
};


}

#endif
