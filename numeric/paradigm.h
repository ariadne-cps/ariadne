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

enum class ParadigmCode : short {
    APPROXIMATE_FLAG=1,
    LOWER_FLAG=2,
    UPPER_FLAG=4,
    METRIC_FLAG=8,
    EFFECTIVE_FLAG=16,
    EXACT_FLAG=32,
    POSITIVE_FLAG=64,

    RAW=0,
    APPROXIMATE=ParadigmCode::APPROXIMATE_FLAG,
    LOWER=ParadigmCode::LOWER_FLAG|ParadigmCode::APPROXIMATE,
    UPPER=ParadigmCode::UPPER_FLAG|ParadigmCode::APPROXIMATE,
    BOUNDED=ParadigmCode::LOWER|ParadigmCode::UPPER,
    METRIC=ParadigmCode::METRIC_FLAG|ParadigmCode::APPROXIMATE,
    VALIDATED=ParadigmCode::BOUNDED|ParadigmCode::METRIC,
    EFFECTIVE_UNINFORMATIVE=ParadigmCode::EFFECTIVE_FLAG,
    EFFECTIVE_LOWER=ParadigmCode::EFFECTIVE_FLAG|ParadigmCode::LOWER,
    EFFECTIVE_UPPER=ParadigmCode::EFFECTIVE_FLAG|ParadigmCode::UPPER,
    EFFECTIVE=ParadigmCode::EFFECTIVE_FLAG|ParadigmCode::BOUNDED|ParadigmCode::METRIC,
    EXACT=ParadigmCode::EXACT_FLAG|ParadigmCode::EFFECTIVE,

    POSITIVE_APPROXIMATE=ParadigmCode::POSITIVE_FLAG|ParadigmCode::APPROXIMATE,
    POSITIVE_LOWER=ParadigmCode::POSITIVE_FLAG|ParadigmCode::LOWER,
    POSITIVE_UPPER=ParadigmCode::POSITIVE_FLAG|ParadigmCode::UPPER,
    POSITIVE_EFFECTIVE_UPPER=ParadigmCode::POSITIVE_FLAG|ParadigmCode::EFFECTIVE_UPPER,
    POSITIVE_EXACT=ParadigmCode::POSITIVE_FLAG|ParadigmCode::EXACT,
};
inline constexpr ParadigmCode operator^(ParadigmCode p1, ParadigmCode p2) {
    return ParadigmCode(static_cast<char>(p1) ^ static_cast<char>(p2)); }
inline constexpr ParadigmCode operator&(ParadigmCode p1, ParadigmCode p2) {
    return ParadigmCode(static_cast<char>(p1) & static_cast<char>(p2)); }
inline constexpr ParadigmCode operator|(ParadigmCode p1, ParadigmCode p2) {
    return ParadigmCode(static_cast<char>(p1) | static_cast<char>(p2)); }
inline constexpr ParadigmCode operator!(ParadigmCode p) {
    return ParadigmCode(!static_cast<char>(p)); }
inline constexpr bool operator>=(ParadigmCode p1, ParadigmCode p2) {
    return p2 == (p1&p2); }
inline constexpr ParadigmCode strengthen(ParadigmCode p) {
    return (p >= ParadigmCode::BOUNDED) or (p >= ParadigmCode::METRIC) ? (p | ParadigmCode::VALIDATED) : p; }
constexpr bool is_weaker(ParadigmCode p1, ParadigmCode p2) {
    return (p1 & strengthen(p2)) == p1; }
constexpr ParadigmCode next_weaker(ParadigmCode p) {
    return is_weaker(ParadigmCode::EFFECTIVE_FLAG,p) ? (p & ParadigmCode::VALIDATED) : ParadigmCode::APPROXIMATE; }
constexpr ParadigmCode widen(ParadigmCode p) {
    return p & ParadigmCode::VALIDATED; }
constexpr ParadigmCode positive(ParadigmCode p) {
    return p | ParadigmCode::POSITIVE_FLAG; }



//! \defgroup ParadigmSubModule Computational Paradigms
//!   \ingroup LogicModule
//! \brief Tag classes describing the information provided by a type.
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

template<ParadigmCode CODE> struct InformationLevel { static constexpr ParadigmCode code() { return CODE; } };

//! \ingroup ParadigmSubModule
//! \brief The <em>computational paradigm</em> supported by the object.
//! User paradigms are Exact, Effective, Validated,ValidatedBounded, ValidatedUpper, ValidatedLower or Approximate.
//! Internal paradigms are Builtin and Raw.
template<class T> using Paradigm = typename T::Paradigm;

//! \ingroup ParadigmSubModule
//! \brief The <em>computational paradigm</em> supported by the object. Equivalent to Paradigm<T>.
template<class T> using ParadigmTag = typename T::Paradigm;


//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object is of a builtin type. Such objects should be converted to %Ariadne internal types before use.
struct Builtin : InformationLevel<ParadigmCode::RAW>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object decribes raw data. Such objects should not be used in high-level code, as they probably do not provide safe guarantees on their values.
struct Raw : InformationLevel<ParadigmCode::RAW>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, and equality is decidable.
//! Only available for discrete types i.e. elements of countable spaces, or for computational types such as floating-point numbers.
struct Exact : InformationLevel<ParadigmCode::EXACT>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, but equality is undecidable.
struct Effective : InformationLevel<ParadigmCode::EFFECTIVE>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, but only convergent upper bounds can be computed.
struct EffectiveUpper : InformationLevel<ParadigmCode::EFFECTIVE_UPPER>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, but only convergent lower bounds can be computed.
struct EffectiveLower : InformationLevel<ParadigmCode::EFFECTIVE_LOWER>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, but no information can be extracted in finite time.
struct EffectiveUninformative : InformationLevel<ParadigmCode::EFFECTIVE_UNINFORMATIVE>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents an approximation to a quantity with a bound on the error in some metric, or lower and upper bounds on a quantity.
struct Validated : InformationLevel<ParadigmCode::VALIDATED>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides an approximation with a bound on the metric error.
//! For numbers, a specialisation of Validated.
struct Metric : InformationLevel<ParadigmCode::METRIC>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides lower and upper bounds for a quantity.
//! For numbers, is a specialisation of Validated.
struct Bounded : InformationLevel<ParadigmCode::BOUNDED>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides an upper bound for a quantity.
struct Upper : InformationLevel<ParadigmCode::UPPER>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides a lower bound for a quantity.
struct Lower : InformationLevel<ParadigmCode::LOWER>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides an approximation to a quantity with no guarantees on the error.
struct Approximate : InformationLevel<ParadigmCode::APPROXIMATE>  { };


//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents an upper bound for a positive quantity.
struct PositiveUpper;

struct PositiveApproximate : InformationLevel<ParadigmCode::POSITIVE_APPROXIMATE>  { };
struct PositiveLower : InformationLevel<ParadigmCode::POSITIVE_LOWER>  { };
struct PositiveUpper : InformationLevel<ParadigmCode::POSITIVE_UPPER>  { };
struct PositiveEffectiveUpper : InformationLevel<ParadigmCode::POSITIVE_EFFECTIVE_UPPER>  { };
struct PositiveExact : InformationLevel<ParadigmCode::POSITIVE_EXACT>  { };

using ValidatedLower = Lower;
using ValidatedUpper = Upper;
using ValidatedBounded = Bounded;
using ValidatedMetric = Metric;
using PositiveValidatedLower = PositiveLower;
using PositiveValidatedUpper = PositiveUpper;
using Error = PositiveUpper;

using ApproximateTag = Approximate;
using LowerTag = Lower;
using UpperTag = Upper;
using BoundedTag = Bounded;
using MetricTag = Metric;
using ValidatedTag = Validated;
using EffectiveLowerTag = EffectiveLower;
using EffectiveUpperTag = EffectiveUpper;
using EffectiveTag = Effective;
using ExactTag = Exact;
using ErrorTag = Error;

template<ParadigmCode PC> struct InformationTraits { };
template<> struct InformationTraits<ParadigmCode::RAW> { typedef Void Paradigm; };
template<> struct InformationTraits<ParadigmCode::APPROXIMATE> { typedef Approximate Paradigm; };
template<> struct InformationTraits<ParadigmCode::LOWER> { typedef Lower Paradigm; };
template<> struct InformationTraits<ParadigmCode::UPPER> { typedef Upper Paradigm; };
template<> struct InformationTraits<ParadigmCode::BOUNDED> { typedef Bounded Paradigm; };
template<> struct InformationTraits<ParadigmCode::METRIC> { typedef Metric Paradigm; };
template<> struct InformationTraits<ParadigmCode::VALIDATED> { typedef Validated Paradigm; };
template<> struct InformationTraits<ParadigmCode::EFFECTIVE_LOWER> { typedef EffectiveLower Paradigm; };
template<> struct InformationTraits<ParadigmCode::EFFECTIVE_UPPER> { typedef EffectiveUpper Paradigm; };
template<> struct InformationTraits<ParadigmCode::EFFECTIVE> { typedef Effective Paradigm; };
template<> struct InformationTraits<ParadigmCode::EXACT> { typedef Exact Paradigm; };
template<> struct InformationTraits<ParadigmCode::POSITIVE_APPROXIMATE> { typedef PositiveApproximate Paradigm; };
template<> struct InformationTraits<ParadigmCode::POSITIVE_UPPER> { typedef PositiveUpper Paradigm; };
template<> struct InformationTraits<ParadigmCode::POSITIVE_EXACT> { typedef PositiveExact Paradigm; };
template<ParadigmCode PC> using ParadigmClass = typename InformationTraits<PC>::Paradigm;

namespace Detail {
Approximate equality_paradigm(Approximate,Approximate);
ValidatedLower equality_paradigm(ValidatedLower,ValidatedUpper);
ValidatedLower equality_paradigm(ValidatedUpper,ValidatedLower);
ValidatedLower equality_paradigm(Validated,Validated);
ValidatedLower equality_paradigm(ValidatedBounded,ValidatedBounded);
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

Approximate negate_paradigm(PositiveApproximate);
ValidatedLower negate_paradigm(PositiveValidatedUpper);
ValidatedUpper negate_paradigm(PositiveValidatedLower);
Exact negate_paradigm(PositiveExact);

Exact invert_paradigm(Exact);
Effective invert_paradigm(Effective);
Validated invert_paradigm(Validated);
ValidatedBounded invert_paradigm(ValidatedBounded);
ValidatedMetric invert_paradigm(ValidatedMetric);
ValidatedLower invert_paradigm(ValidatedUpper);
ValidatedUpper invert_paradigm(ValidatedLower);
Approximate invert_paradigm(Approximate);

PositiveApproximate invert_paradigm(PositiveApproximate);
PositiveValidatedLower invert_paradigm(PositiveValidatedUpper);
PositiveValidatedUpper invert_paradigm(PositiveValidatedLower);
PositiveExact invert_paradigm(PositiveExact);

template<class T> T strengthen_paradigm(T);
Validated strengthen_paradigm(ValidatedBounded);
Validated strengthen_paradigm(ValidatedMetric);


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
ValidatedBounded widen_paradigm(Effective);
ValidatedBounded widen_paradigm(Validated);

template<class T> T unsigned_paradigm(T);
PositiveExact unsigned_paradigm(Exact);
PositiveEffectiveUpper unsigned_paradigm(EffectiveUpper);
PositiveValidatedUpper unsigned_paradigm(ValidatedUpper);
PositiveValidatedLower unsigned_paradigm(ValidatedLower);
PositiveApproximate unsigned_paradigm(Approximate);

template<class T> T signed_paradigm(T);
Exact signed_paradigm(PositiveExact);
EffectiveUpper signed_paradigm(PositiveEffectiveUpper);
ValidatedUpper signed_paradigm(PositiveValidatedUpper);
ValidatedLower signed_paradigm(PositiveValidatedLower);
Approximate signed_paradigm(PositiveApproximate);

}

template<bool b> using BooleanConstant = std::integral_constant<bool,b>;

//! \ingroup ParadigmSubModule
//! \brief Inherits from \c TrueType if \a P is a paradigm tag class.
template<class P> struct IsParadigm : IsSame<decltype(P::code()),ParadigmCode> { };
//! \brief Inherits from \c TrueType if paradigm \a P1 is weaker than \a P2.
template<class P1, class P2> struct IsWeaker : BooleanConstant<is_weaker(P1::code(),P2::code())> { };

//! \brief Inherits from \c TrueType if paradigm \a P1 is stronger than \a P2.
template<class P1, class P2> struct IsStronger : IsWeaker<P2,P1> { };

template<class P1, class P2=P1> struct ParadigmTraits {
    using Weaker = typename InformationTraits<P1::code()&P2::code()>::Paradigm;
    using Stronger = typename InformationTraits<P1::code()&P2::code()>::Paradigm;
};

template<> struct ParadigmTraits<Metric,Bounded> {
    using Weaker = typename InformationTraits<ParadigmCode::BOUNDED>::Paradigm;
    using Stronger = typename InformationTraits<ParadigmCode::METRIC>::Paradigm;
};

template<> struct ParadigmTraits<Bounded,Metric> {
    using Weaker = ParadigmClass<ParadigmCode::BOUNDED>;
    using Stronger = ParadigmClass<ParadigmCode::METRIC>;
};

template<> struct ParadigmTraits<Metric,Upper> {
    using Weaker = typename InformationTraits<ParadigmCode::UPPER>::Paradigm;
    using Stronger = typename InformationTraits<ParadigmCode::METRIC>::Paradigm;
};

template<> struct ParadigmTraits<Upper,Metric> {
    using Weaker = ParadigmClass<ParadigmCode::UPPER>;
    using Stronger = ParadigmClass<ParadigmCode::METRIC>;
};

template<> struct ParadigmTraits<Metric,Lower> {
    using Weaker = typename InformationTraits<ParadigmCode::LOWER>::Paradigm;
    using Stronger = typename InformationTraits<ParadigmCode::METRIC>::Paradigm;
};

template<> struct ParadigmTraits<Lower,Metric> {
    using Weaker = ParadigmClass<ParadigmCode::LOWER>;
    using Stronger = ParadigmClass<ParadigmCode::METRIC>;
};

template<class P> struct ParadigmTraits<P,P> {
    using Weaker = typename InformationTraits<P::code()>::Paradigm;
    using Stronger = typename InformationTraits<P::code()>::Paradigm;
    using NextWeaker = typename InformationTraits<next_weaker(P::code())>::Paradigm;
};

template<> struct ParadigmTraits<Approximate,Approximate> {
  private:
    typedef Approximate P;
  public:
    using Weaker = typename InformationTraits<P::code()>::Paradigm;
    using Stronger = typename InformationTraits<P::code()>::Paradigm;
    using NextWeaker = Void;
};

//! \ingroup ParadigmSubModule
//! \brief The strongest paradigm which is weaker than both paradigms \a P1 and \a P2.
template<class P1, class P2> using Weaker = typename ParadigmTraits<P1,P2>::Weaker;
//! \ingroup ParadigmSubModule
//! \brief The weakest paradigm which is stronger than both paradigms \a P1 and \a P2.
template<class P1, class P2> using Stronger = typename ParadigmTraits<P1,P2>::Stronger;
//! \ingroup ParadigmSubModule
//! \brief The strongest paradigm which is strictly weaker than \a P. If \a P is \c ValidatedBounded, is defined as \c Approximate.
template<class P> using NextWeaker = typename ParadigmTraits<P>::NextWeaker;
template<class P> using Weaken = typename ParadigmTraits<P>::NextWeaker;



//! \ingroup ParadigmSubModule
//! \brief The paradigm associated with the negation of an object.
template<class P> using Negated = decltype(Detail::negate_paradigm(declval<P>()));
template<class P> using Inverted = decltype(Detail::invert_paradigm(declval<P>()));
template<class P> using Opposite = decltype(Detail::negate_paradigm(declval<P>()));
template<class P> using Generic = decltype(Detail::strengthen_paradigm(declval<P>()));
template<class P> using Unsigned = decltype(Detail::unsigned_paradigm(declval<P>()));
template<class P> using Signed = decltype(Detail::signed_paradigm(declval<P>()));
template<class P> using Unorder = Weaker<P,Opposite<P>>;

//! \ingroup ParadigmSubModule
//! \brief The paradigm obtained by widening a type due to roundoff error.
template<class P> using Widen = decltype(Detail::widen_paradigm(declval<P>()));

template<class P1, class P2=Negated<P1>> using Equality = decltype(Detail::equality_paradigm(declval<P1>(),declval<P2>()));

template<class P> using Errors = decltype(Detail::error_paradigm(declval<P>()));

template<class P1, class P2=Negated<P1>> using LessThan = Weaker<P1,Negated<P2>>;


} // namespace Ariadne

#endif /* ARIADNE_PARADIGM_H */
