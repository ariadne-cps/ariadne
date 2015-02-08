/***************************************************************************
 *            numeric/paradigm.h
 *
 *  Copyright 2013-14  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file numeric/paradigm.h
 *  \brief
 */

#ifndef ARIADNE_PARADIGM_H
#define ARIADNE_PARADIGM_H

#include "utility/metaprogramming.h"

namespace Ariadne {

typedef Void Void;

class ParadigmError { };

enum class ParadigmCode : char {
    ERROR, EXACT,
    EFFECTIVE, EFFECTIVE_UPPER, EFFECTIVE_LOWER,
    VALIDATED, METRIC, BOUNDED, UPPER, LOWER,
    APPROXIMATE, RAW
};

//! \defgroup ParadigmSubModule Computational Paradigms
//!   \ingroup LogicModule
//! \brief Tag classes describing the computational paradigm.
//! \details In order to indicate the kind of guarantees on the approximation provided by a concrete object,
//!   every type in %Ariadne has an associated Paradigm tag.
//!
//! %Ariadne also needs to handle built-in C++ types, notably \c Bool, \c Int and \c double (floating-point) classes.
//! Since C++ internally allows unsafe and inexact conversions e.g
//!     \code Int n=1.5; // n is set to 1! \endcode
//! and
//!     \code double x=1.3; // x not exactly equal to 1.3! \endcode
//! the Builtin tag is used to describe these objects.
//! %Builtin objects should be immediately converted to %Ariadne objects in user code.


//! \ingroup ParadigmSubModule
//! \brief The <em>computational paradigm</em> supported by the object.
//! User paradigms are Exact, Effective, Validated,ValidatedBounded, ValidatedUpper, ValidatedLower or Approximate.
//! Internal paradigms are Builtin and Raw.
template<class T> using Paradigm = typename T::Paradigm;

//! \ingroup ParadigmSubModule
//! \brief The <em>computational paradigm</em> supported by the object. Equivalent to Paradigm<T>.
template<class T> using ParadigmTag = typename T::Paradigm;

struct Exact;
struct Effective;
struct EffectiveUpper;
struct EffectiveLower;
struct Validated;
struct ValidatedMetric;
struct ValidatedBounded;
struct ValidatedUpper;
struct ValidatedLower;
struct Approximate;
struct PositiveExact;
struct PositiveValidated;
struct PositiveValidatedUpper;
struct PositiveValidatedLower;
struct PositiveApproximate;
using PositiveUpper=PositiveValidatedUpper;
using PositiveLower=PositiveValidatedLower;

using Error=PositiveUpper;
using Metric=ValidatedMetric;
using Bounded=ValidatedBounded;
using Upper=ValidatedUpper;
using Lower=ValidatedLower;

using Valid=Validated;
using Metrc=ValidatedMetric;
using Bound=ValidatedBounded;
using Apprx=Approximate;

using ExactTag = Exact;
using ErrorTag = Error;
using EffectiveTag = Effective;
using EffectiveUpperTag = EffectiveUpper;
using EffectiveLowerTag = EffectiveLower;
using ValidatedTag = Validated;
using ValidatedMetricTag = ValidatedMetric;
using ValidatedBoundedTag = ValidatedBounded;
using ValidatedUpperTag = ValidatedUpper;
using ValidatedLowerTag = ValidatedLower;
using MetricTag = Metric;
using BoundedTag = Bounded;
using UpperTag = Upper;
using LowerTag = Lower;
using ApproximateTag = Approximate;

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object is of a builtin type. Such objects should be converted to %Ariadne internal types before use.
struct Builtin { Builtin() { } };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object decribes raw data. Such objects should not be used in high-level code, as they probably do not provide safe guarantees on their values.
struct Raw { Raw() { } static const ParadigmCode code = ParadigmCode::RAW; };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, and equality is decidable.
//! Only available for discrete types i.e. elements of countable spaces, or for computational types such as floating-point numbers.
struct Exact {
    static const ParadigmCode code = ParadigmCode::EXACT;
    Exact() { }
    typedef Effective NextWeakerParadigm;
};

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, but equality is undecidable.
struct Effective {
    static const ParadigmCode code = ParadigmCode::EFFECTIVE;
    Effective() { } Effective(Exact) { }
    typedef Validated NextWeakerParadigm;
};

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents upper bounds to a quantity arbitrarily accurately.
struct EffectiveUpper {
    static const ParadigmCode code = ParadigmCode::EFFECTIVE_UPPER;
    EffectiveUpper() { } EffectiveUpper(Exact) { } EffectiveUpper(Effective) { }
    typedef ValidatedUpper NextWeakerParadigm;
};

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, but equality is undecidable.
struct EffectiveLower {
    static const ParadigmCode code = ParadigmCode::EFFECTIVE_LOWER;
    EffectiveLower() { } EffectiveLower(Exact) { } EffectiveLower(Effective) { }
    typedef ValidatedLower NextWeakerParadigm;
};

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents an upper bound for a positive quantity.
struct PositiveValidatedUpper {
    static const ParadigmCode code = ParadigmCode::ERROR;
    PositiveValidatedUpper() { }
};

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents an approximation to a quantity with a bound on the error in some metric, or lower and upper bounds on a quantity.
struct Validated {
    static const ParadigmCode code = ParadigmCode::VALIDATED;
    Validated() { } Validated(Exact) { } Validated(Effective) { } Validated(ValidatedMetric); Validated(ValidatedBounded);
    typedef Approximate NextWeakerParadigm;
};

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides an approximation with a bound on the metric error.
//! For numbers, a specialisation of Validated.
struct ValidatedMetric {
    static const ParadigmCode code = ParadigmCode::METRIC;
    ValidatedMetric() { } ValidatedMetric(Exact) { } ValidatedMetric(Effective) { }
    ValidatedMetric(Validated) { } ValidatedMetric(ValidatedBounded);
    typedef Approximate NextWeakerParadigm;
};

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides lower and upper bounds for a quantity.
//! For numbers, is a specialisation of Validated.
struct ValidatedBounded {
    static const ParadigmCode code = ParadigmCode::BOUNDED;
    ValidatedBounded() { } ValidatedBounded(Exact) { } ValidatedBounded(Effective) { }
    ValidatedBounded(Validated) { } ValidatedBounded(ValidatedMetric);
    typedef Approximate NextWeakerParadigm;
};

inline Validated::Validated(ValidatedMetric) { }
inline Validated::Validated(ValidatedBounded) { }
inline ValidatedMetric::ValidatedMetric(ValidatedBounded) { }
inline ValidatedBounded::ValidatedBounded(ValidatedMetric) { }

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides an upper bound for a quantity.
struct ValidatedUpper {
    static const ParadigmCode code = ParadigmCode::UPPER;
    ValidatedUpper() { } ValidatedUpper(Exact) { } ValidatedUpper(Effective) { } ValidatedUpper(EffectiveUpper) { }
    ValidatedUpper(Validated) { } ValidatedUpper(ValidatedMetric) { } ValidatedUpper(ValidatedBounded) { }
};

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides a lower bound for a quantity.
struct ValidatedLower {
    static const ParadigmCode code = ParadigmCode::LOWER;
    ValidatedLower() { } ValidatedLower(Exact) { } ValidatedLower(Effective) { } ValidatedLower(EffectiveLower) { }
    ValidatedLower(Validated) { } ValidatedLower(ValidatedMetric) { } ValidatedLower(ValidatedBounded) { }
};

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides an approximation to a quantity with no guarantees on the error.
struct Approximate {
    static const ParadigmCode code = ParadigmCode::APPROXIMATE;
    Approximate() { } Approximate(Exact) { }
    Approximate(Effective) { } Approximate(EffectiveUpper) { } Approximate(EffectiveLower) { }
    Approximate(Validated) { }Approximate(ValidatedBounded) { } Approximate(ValidatedMetric) { } Approximate(ValidatedUpper) { } Approximate(ValidatedLower) { }
    typedef Void NextWeakerParadigm;
};

namespace Detail {
Exact weaker_paradigm(Exact,Exact);
Effective weaker_paradigm(Effective,Effective);
EffectiveUpper weaker_paradigm(EffectiveUpper,EffectiveUpper);
EffectiveLower weaker_paradigm(EffectiveLower,EffectiveLower);
Validated weaker_paradigm(Validated,Validated);
Validated weaker_paradigm(ValidatedMetric,ValidatedMetric);
Validated weaker_paradigm(ValidatedBounded,ValidatedBounded);
ValidatedUpper weaker_paradigm(ValidatedUpper,ValidatedUpper);
ValidatedLower weaker_paradigm(ValidatedLower,ValidatedLower);
Approximate weaker_paradigm(Approximate,Approximate);

Approximate equality_paradigm(Approximate,Approximate);
ValidatedLower equality_paradigm(ValidatedLower,ValidatedUpper);
ValidatedLower equality_paradigm(ValidatedUpper,ValidatedLower);
ValidatedLower equality_paradigm(Validated,Validated);
ValidatedLower equality_paradigm(ValidatedMetric,ValidatedMetric);
ValidatedLower equality_paradigm(Effective,Effective);
Exact equality_paradigm(Exact,Exact);

Exact negate_paradigm(Exact);
Effective negate_paradigm(Effective);
Validated negate_paradigm(Validated);
ValidatedBounded negate_paradigm(ValidatedBounded);
ValidatedMetric negate_paradigm(ValidatedMetric);
ValidatedLower negate_paradigm(ValidatedUpper);
ValidatedUpper negate_paradigm(ValidatedLower);
Approximate negate_paradigm(Approximate);

ValidatedUpper error_paradigm(Exact);
ValidatedUpper error_paradigm(Effective);
ValidatedUpper error_paradigm(Validated);
ValidatedUpper error_paradigm(ValidatedMetric);
ValidatedUpper error_paradigm(ValidatedBounded);
ValidatedUpper error_paradigm(ValidatedUpper);
Approximate error_paradigm(ValidatedLower);
Approximate error_paradigm(Approximate);

template<class T> T widen_paradigm(T);
ValidatedBounded widen_paradigm(Exact);

}

//! \ingroup ParadigmSubModule
//! \brief Inherits from \c TrueType if \a P is a paradigm tag class.
template<class P> struct IsParadigm : IsConvertible<P,Approximate> { };
//! \brief Inherits from \c TrueType if paradigm \a P1 is weaker than \a P2.
template<class P1, class P2> struct IsWeaker : IsConvertible<P2,P1> { };
//! \brief Inherits from \c TrueType if paradigm \a P1 is stronger than \a P2.
template<class P1, class P2> struct IsStronger : IsConvertible<P1,P2> { };

//! \ingroup ParadigmSubModule
//! \brief The strongest paradigm which is weaker than both paradigms \a P1 and \a P2.
template<class P1, class P2> using Weaker = decltype(Detail::weaker_paradigm(declval<P1>(),declval<P2>()));

//! \ingroup ParadigmSubModule
//! \brief The strongest paradigm which is strictly weaker than \a P. If \a P is \c ValidatedBounded, is defined as \c Approximate.
template<class P> using Weaken = typename P::NextWeakerParadigm;



//! \ingroup ParadigmSubModule
//! \brief The paradigm associated with the negation of an object.
template<class P> using Negated = decltype(Detail::negate_paradigm(declval<P>()));
template<class P> using Opposite = decltype(Detail::negate_paradigm(declval<P>()));

//! \ingroup ParadigmSubModule
//! \brief The paradigm obtained by widening a type due to roundoff error.
template<class P> using Widen = decltype(Detail::widen_paradigm(declval<P>()));

template<class P1, class P2=Negated<P1>> using Equality = decltype(Detail::equality_paradigm(declval<P1>(),declval<P2>()));

template<class P> using Errors = decltype(Detail::error_paradigm(declval<P>()));

template<class P1, class P2=Negated<P1>> using LessThan = Weaker<P1,Negated<P2>>;


} // namespace Ariadne

#endif /* ARIADNE_PARADIGM_H */
