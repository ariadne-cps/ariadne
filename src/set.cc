#include "set.h"

namespace Ariadne {

 
ImageSet::ImageSet(const BoxType& dom, const FunctionInterface& fn)
  : _domain(dom), _function(fn.clone()) 
{ 
  ARIADNE_ASSERT(dom.size()==fn.argument_size()); 
}


ConstraintSet::ConstraintSet(const BoxType& codom, const FunctionInterface& fn) 
  : _codomain(codom), _function(fn.clone()) 
{ 
  ARIADNE_ASSERT(codom.size()==fn.result_size());
}


} // namespace Ariadne
