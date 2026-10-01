/***************************************************************************
 *            function/declarations.hpp
 *
 *  Copyright  2011-26  Pieter Collins
 *
 ****************************************************************************/

#ifndef ARIADNE_FUNCTION_DECLARATIONS_HPP
#define ARIADNE_FUNCTION_DECLARATIONS_HPP

#include <iosfwd>

#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"

#include "paradigm/paradigm.hpp"
#include "paradigm/logical.decl.hpp"

#include "function/function.decl.hpp"

namespace Ariadne {


template<class P, class F> class AffineModel;
template<class P, class F> class TaylorModel;
template<class X> class Differential;
template<class X> class ElementaryAlgebra;
template<class X> class Formula;

} // namespace Ariadne

#endif /* ARIADNE_FUNCTION_DECLARATIONS_HPP */
