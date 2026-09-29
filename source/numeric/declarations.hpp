/***************************************************************************
 *            numeric/declarations.hpp
 *
 *  Copyright  2011-26  Pieter Collins
 *
 ****************************************************************************/

#ifndef ARIADNE_NUMERIC_DECLARATIONS_HPP
#define ARIADNE_NUMERIC_DECLARATIONS_HPP

#include <iosfwd>

#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"

#include "foundation/paradigm.hpp"
#include "foundation/logical.decl.hpp"

#include "numeric/number.decl.hpp"
#include "numeric/float.decl.hpp"

namespace Ariadne {


template<class X> struct InformationTypedef;
template<> struct InformationTypedef<Real> { typedef EffectiveTag Type; };
template<class P> struct InformationTypedef<Number<P>> { typedef P Type; };
template<class X> using InformationTag = typename InformationTypedef<X>::Type;

} // namespace Ariadne

#endif /* ARIADNE_NUMERIC_DECLARATIONS_HPP */
