/***************************************************************************
 *            numeric/logical.decl.hpp
 *
 *  Copyright  2013-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file numeric/logical.decl.hpp
 *  \brief
 */



#ifndef ARIADNE_LOGICAL_DECL_HPP
#define ARIADNE_LOGICAL_DECL_HPP

#include "../numeric/paradigm.hpp"

namespace Ariadne {

typedef void Void;
typedef bool Bool;

class Indeterminate;
extern const Indeterminate indeterminate;

class Boolean;
class Sierpinskian;
class Kleenean;

template<class L> class Negate;
template<class L> class Lower;
template<class L> class Upper;
template<class L> class Validated;
template<class L> class Approximate;

class NegatedSierpinskian;
class LowerKleenean;
class UpperKleenean;
class ValidatedSierpinskian;
class ValidatedNegatedSierpinskian;
class ValidatedKleenean;
class ValidatedLowerKleenean;
class ValidatedUpperKleenean;
class ApproximateKleenean;

using Fuzzy = ApproximateKleenean;

template<class P> struct LogicalTypedef;
template<class P> using LogicalType = typename LogicalTypedef<P>::Type;
template<> struct LogicalTypedef<ExactTag> { typedef Boolean Type; };
template<> struct LogicalTypedef<EffectiveTag> { typedef Kleenean Type; };
template<> struct LogicalTypedef<EffectiveLowerTag> { typedef LowerKleenean Type; };
template<> struct LogicalTypedef<EffectiveUpperTag> { typedef UpperKleenean Type; };
template<> struct LogicalTypedef<ValidatedTag> { typedef ValidatedKleenean Type; };
template<> struct LogicalTypedef<ValidatedLowerTag> { typedef ValidatedLowerKleenean Type; };
template<> struct LogicalTypedef<ValidatedUpperTag> { typedef ValidatedUpperKleenean Type; };
template<> struct LogicalTypedef<ApproximateTag> { typedef ApproximateKleenean Type; };

using Decidable = Boolean;
using Quasidecidable = Kleenean;
using Verifyable = Sierpinskian;
using Falsifyable = NegatedSierpinskian;

Boolean operator!(Boolean const&);
NegatedSierpinskian operator!(Sierpinskian const&);
Kleenean operator!(Kleenean const&);

namespace Detail {
enum class LogicalValue : char;
}
using Detail::LogicalValue;

}

#endif
