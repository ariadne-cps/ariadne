/***************************************************************************
 *            algebra/declarations.hpp
 *
 *  Copyright  2011-26  Pieter Collins
 *
 ****************************************************************************/

#ifndef ARIADNE_ALGEBRA_DECLARATIONS_HPP
#define ARIADNE_ALGEBRA_DECLARATIONS_HPP

#include <iosfwd>

#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"

#include "paradigm/paradigm.hpp"
#include "paradigm/logical.decl.hpp"

#include "algebra/linear_algebra.decl.hpp"
#include "algebra/differential.decl.hpp"

namespace Ariadne {


template<class X> class Algebra;
template<class X> class ElementaryAlgebra;
template<class X> class Series;

} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_DECLARATIONS_HPP */
