/***************************************************************************
 *            numeric/paradigm.hpp
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

/*! \file numeric/paradigm.hpp
 *  \brief
 */

#ifndef ARIADNE_PARADIGM_HPP
#define ARIADNE_PARADIGM_HPP

#include <cstdint>
#include "../utility/metaprogramming.hpp"

namespace Ariadne {

typedef Void Void;

class ParadigmError { };

typedef std::uint16_t  ParadigmCodeType;

enum class ParadigmCode : ParadigmCodeType {
    APPROXIMATE_FLAG=1,
    LOWER_FLAG=2,
    UPPER_FLAG=4,
    BOUNDS_FLAGS=LOWER_FLAG|UPPER_FLAG,
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
    POSITIVE_BOUNDED=ParadigmCode::BOUNDED|ParadigmCode::POSITIVE_FLAG,
    POSITIVE_METRIC=ParadigmCode::METRIC|ParadigmCode::POSITIVE_FLAG,
    POSITIVE_VALIDATED=ParadigmCode::BOUNDED|ParadigmCode::METRIC|ParadigmCode::POSITIVE_FLAG,
    POSITIVE_EFFECTIVE_UPPER=ParadigmCode::POSITIVE_FLAG|ParadigmCode::EFFECTIVE_UPPER,
    POSITIVE_EXACT=ParadigmCode::POSITIVE_FLAG|ParadigmCode::EXACT,
};
inline constexpr ParadigmCode operator^(ParadigmCode p1, ParadigmCode p2) {
    return ParadigmCode(static_cast<ParadigmCodeType>(p1) ^ static_cast<ParadigmCodeType>(p2)); }
inline constexpr ParadigmCode operator&(ParadigmCode p1, ParadigmCode p2) {
    return ParadigmCode(static_cast<ParadigmCodeType>(p1) & static_cast<ParadigmCodeType>(p2)); }
inline constexpr ParadigmCode operator|(ParadigmCode p1, ParadigmCode p2) {
    return ParadigmCode(static_cast<ParadigmCodeType>(p1) | static_cast<ParadigmCodeType>(p2)); }
inline constexpr ParadigmCode operator!(ParadigmCode p) {
    return ParadigmCode(~static_cast<ParadigmCodeType>(p)); }
inline constexpr bool operator>=(ParadigmCode p1, ParadigmCode p2) {
    return p2 == (p1&p2); }

inline constexpr ParadigmCode strengthen(ParadigmCode p) {
    return (p >= ParadigmCode::BOUNDED) or (p >= ParadigmCode::METRIC) ? (p | ParadigmCode::VALIDATED) : p; }
constexpr bool is_weaker(ParadigmCode p1, ParadigmCode p2) {
    return (p1 & strengthen(p2)) == p1; }
constexpr ParadigmCode next_weaker(ParadigmCode p) {
    return is_weaker(ParadigmCode::EFFECTIVE_FLAG,p) ? (p & ParadigmCode::VALIDATED) : ParadigmCode::APPROXIMATE; }

constexpr Bool is_exact(ParadigmCode p) {
    return (p >= ParadigmCode::EXACT_FLAG) ; }
constexpr Bool is_effective(ParadigmCode p) {
    return (p >= ParadigmCode::EFFECTIVE_FLAG) ; }
constexpr Bool is_validated(ParadigmCode p) {
    return (p >= ParadigmCode::BOUNDS_FLAGS) xor (p >= ParadigmCode::METRIC_FLAG); }
constexpr Bool is_directed(ParadigmCode p) {
    return Bool(p & ParadigmCode::LOWER_FLAG) xor Bool(p & ParadigmCode::UPPER_FLAG); }
constexpr Bool is_unsigned(ParadigmCode p) {
    return Bool(p | ParadigmCode::POSITIVE_FLAG); }
constexpr ParadigmCode unsign(ParadigmCode p) {
    return p | ParadigmCode::POSITIVE_FLAG; }
constexpr ParadigmCode sign(ParadigmCode p) {
    return p & !ParadigmCode::POSITIVE_FLAG; }
constexpr ParadigmCode opposite(ParadigmCode p) {
    return is_directed(p) ? p ^ ParadigmCode::BOUNDS_FLAGS : p; }
constexpr ParadigmCode undirect(ParadigmCode p) {
    return is_directed(p) ? p & !ParadigmCode::BOUNDS_FLAGS : p; }
constexpr ParadigmCode negate(ParadigmCode p) {
    return sign(is_directed(p) ? p ^ ParadigmCode::BOUNDS_FLAGS : p); }
constexpr ParadigmCode invert(ParadigmCode p) {
    return is_directed(p) ? (is_unsigned(p) ? opposite(p) : undirect(p)) : p; }
constexpr ParadigmCode widen(ParadigmCode p) {
    return Bool(p & ParadigmCode::EFFECTIVE_FLAG) ? p & ParadigmCode::POSITIVE_BOUNDED : p; }
//    return Bool(p & ParadigmCode::EFFECTIVE_FLAG) ? p & ParadigmCode::POSITIVE_METRIC : p; }
constexpr ParadigmCode null(ParadigmCode p) {
    return is_exact(p) ? sign(p) : is_validated(p) ? p & ParadigmCode::LOWER : p & ParadigmCode::APPROXIMATE; }



//! \defgroup ParadigmSubModule Computational paradigms
//!   \ingroup LogicalModule
//! \brief Tag classes describing the information provided by a type.
//! \details In order to indicate the kind of guarantees on the approximation provided by a concrete object,
//!   every type in %Ariadne has an associated Paradigm tag.
//!
//! %Ariadne also needs to handle built-in C++ types, notably \c Bool, \c Int and \c double (floating-point) classes.
//! Since C++ internally allows unsafe and inexact conversions e.g
//!     \code Int n=1.5; // n is set to 1! \endcode
//! and
//!     \code double x=1.3; // x not exactly equal to 1.3! \endcode
//! the BuiltinTag tag is used to describe these objects.
//! %BuiltinTag objects should be immediately converted to %Ariadne objects in user code.

template<ParadigmCode CODE> struct InformationLevel { static constexpr ParadigmCode code() { return CODE; } };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object is of a builtin type. Such objects should be converted to %Ariadne internal types before use.
struct BuiltinTag : InformationLevel<ParadigmCode::RAW>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object decribes raw data. Such objects should not be used in high-level code, as they probably do not provide safe guarantees on their values.
struct RawTag : InformationLevel<ParadigmCode::RAW>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, and equality is decidable.
//! Only available for discrete types i.e. elements of countable spaces, or for computational types such as floating-point numbers.
struct ExactTag : InformationLevel<ParadigmCode::EXACT>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, but equality is undecidable.
struct EffectiveTag : InformationLevel<ParadigmCode::EFFECTIVE>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, but only convergent upper bounds can be computed.
struct EffectiveUpperTag : InformationLevel<ParadigmCode::EFFECTIVE_UPPER>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, but only convergent lower bounds can be computed.
struct EffectiveLowerTag : InformationLevel<ParadigmCode::EFFECTIVE_LOWER>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, but no information can be extracted in finite time.
struct EffectiveUninformativeTag : InformationLevel<ParadigmCode::EFFECTIVE_UNINFORMATIVE>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents an approximation to a quantity with a bound on the error in some metric, or lower and upper bounds on a quantity.
struct ValidatedTag : InformationLevel<ParadigmCode::VALIDATED>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides an approximation with a bound on the metric error.
//! For numbers, a specialisation of ValidatedTag.
struct MetricTag : InformationLevel<ParadigmCode::METRIC>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides lower and upper bounds for a quantity.
//! For numbers, is a specialisation of ValidatedTag.
struct BoundedTag : InformationLevel<ParadigmCode::BOUNDED>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides an upper bound for a quantity.
struct UpperTag : InformationLevel<ParadigmCode::UPPER>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides a lower bound for a quantity.
struct LowerTag : InformationLevel<ParadigmCode::LOWER>  { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides an approximation to a quantity with no guarantees on the error.
struct ApproximateTag : InformationLevel<ParadigmCode::APPROXIMATE>  { };


//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents an upper bound for a positive quantity.
struct PositiveUpperTag;

struct PositiveApproximateTag : InformationLevel<ParadigmCode::POSITIVE_APPROXIMATE>  { };
struct PositiveLowerTag : InformationLevel<ParadigmCode::POSITIVE_LOWER>  { };
struct PositiveUpperTag : InformationLevel<ParadigmCode::POSITIVE_UPPER>  { };
struct PositiveBoundedTag : InformationLevel<ParadigmCode::POSITIVE_BOUNDED>  { };
struct PositiveMetricTag : InformationLevel<ParadigmCode::POSITIVE_METRIC>  { };
struct PositiveEffectiveUpperTag : InformationLevel<ParadigmCode::POSITIVE_EFFECTIVE_UPPER>  { };
struct PositiveExactTag : InformationLevel<ParadigmCode::POSITIVE_EXACT>  { };

using ValidatedLowerTag = LowerTag;
using ValidatedUpperTag = UpperTag;
using ValidatedBoundedTag = BoundedTag;
using ValidatedMetricTag = MetricTag;
using PositiveValidatedLowerTag = PositiveLowerTag;
using PositiveValidatedUpperTag = PositiveUpperTag;
using ErrorTag = PositiveUpperTag;

template<ParadigmCode PC> struct InformationTraits { };
template<> struct InformationTraits<ParadigmCode::RAW> { typedef Void Paradigm; };
template<> struct InformationTraits<ParadigmCode::APPROXIMATE> { typedef ApproximateTag Paradigm; };
template<> struct InformationTraits<ParadigmCode::LOWER> { typedef LowerTag Paradigm; };
template<> struct InformationTraits<ParadigmCode::UPPER> { typedef UpperTag Paradigm; };
template<> struct InformationTraits<ParadigmCode::BOUNDED> { typedef BoundedTag Paradigm; };
template<> struct InformationTraits<ParadigmCode::METRIC> { typedef MetricTag Paradigm; };
template<> struct InformationTraits<ParadigmCode::VALIDATED> { typedef ValidatedTag Paradigm; };
template<> struct InformationTraits<ParadigmCode::EFFECTIVE_LOWER> { typedef EffectiveLowerTag Paradigm; };
template<> struct InformationTraits<ParadigmCode::EFFECTIVE_UPPER> { typedef EffectiveUpperTag Paradigm; };
template<> struct InformationTraits<ParadigmCode::EFFECTIVE> { typedef EffectiveTag Paradigm; };
template<> struct InformationTraits<ParadigmCode::EXACT> { typedef ExactTag Paradigm; };
template<> struct InformationTraits<ParadigmCode::POSITIVE_APPROXIMATE> { typedef PositiveApproximateTag Paradigm; };
template<> struct InformationTraits<ParadigmCode::POSITIVE_LOWER> { typedef PositiveLowerTag Paradigm; };
template<> struct InformationTraits<ParadigmCode::POSITIVE_UPPER> { typedef PositiveUpperTag Paradigm; };
template<> struct InformationTraits<ParadigmCode::POSITIVE_EXACT> { typedef PositiveExactTag Paradigm; };
template<ParadigmCode PC> using ParadigmClass = typename InformationTraits<PC>::Paradigm;

template<bool b> using BooleanConstant = std::integral_constant<bool,b>;

//! \ingroup ParadigmSubModule
//! \brief Inherits from \c TrueType if \a P is a paradigm tag class.
template<class P> struct IsParadigm : IsSame<decltype(P::code()),ParadigmCode> { };
//! \ingroup ParadigmSubModule
//! \brief Inherits from \c TrueType if paradigm \a P1 is weaker than \a P2.
template<class P1, class P2> struct IsWeaker : BooleanConstant<is_weaker(P1::code(),P2::code())> { };
//! \ingroup ParadigmSubModule
//! \brief Inherits from \c TrueType if paradigm \a P1 is stronger than \a P2.
template<class P1, class P2> struct IsStronger : IsWeaker<P2,P1> { };

template<class P1, class P2=P1> struct ParadigmTraits;

//@{
//! \ingroup ParadigmSubModule
//! \name Information traits
//! \nosubgrouping

#ifdef DOXYGEN
//! \brief Traits class describing information levels
//!   \param P May be ExactTag, EffectiveTag, ValidatedTag, ApproximateTag etc.
//! \ingroup ParadigmSubModule
template<class P1, class P2=P1> struct ParadigmTraits {
  public:
    typedef typename Weaker; //!< The strongest paradigm (most information) which is weaker than (can be converted from) both \a P1 and \a P2.
    typedef typename Stronger; //!< The weakest paradigm (least information) which stronger than (can be converted to) both \a P1 and \a P2.
};
//! \brief Traits class describing information levels
//! \ingroup ParadigmSubModule
template<class P> struct ParadigmTraits<P,P> {
  public:
    typedef P Weaker; //!< The weaker of \a P and itself is \a P.
    typedef P Stronger; //!< The stronger of \a P and itself is \a P.
    typedef typename NextWeaker; //!< The paradigm (information level) directly below (weaker) in the information hierarchy to \a P.
};
#endif

//! \brief The <em>computational paradigm</em> supported by the object.
//! User paradigms are ExactTag, EffectiveTag, ValidatedTag,ValidatedBoundedTag, ValidatedUpperTag, ValidatedLowerTag or ApproximateTag.
//! Internal paradigms are BuiltinTag and RawTag.
//! \ingroup ParadigmSubModule
template<class T> using Paradigm = typename T::Paradigm;

//! \brief The <em>computational paradigm</em> supported by the object. Equivalent to Paradigm<T>.
//! \ingroup ParadigmSubModule
template<class T> using ParadigmTag = typename T::Paradigm;

//! \brief The strongest paradigm (most information) which is weaker than (can be converted from) both \a P1 and \a P2.
//! \ingroup ParadigmSubModule
template<class P1, class P2> using Weaker = typename ParadigmTraits<P1,P2>::Weaker;

//! \brief The weakest paradigm (least information) which stronger than (can be converted to) both \a P1 and \a P2.
//! \ingroup ParadigmSubModule
template<class P1, class P2> using Stronger = typename ParadigmTraits<P1,P2>::Stronger;
//! \brief The strongest paradigm which is strictly weaker than \a P. e.g. If \a P is \c ValidatedTag, is defined as \c ApproximateTag.
//! \ingroup ParadigmSubModule
template<class P> using NextWeaker = typename ParadigmTraits<P>::NextWeaker;
//! \brief Synonym for NextWeaker.
//! \ingroup ParadigmSubModule
template<class P> using Weaken = typename ParadigmTraits<P>::NextWeaker;
//@}


template<class P1, class P2> struct ParadigmTraits {
    using Weaker = typename InformationTraits<P1::code()&P2::code()>::Paradigm;
    using Stronger = typename InformationTraits<P1::code()&P2::code()>::Paradigm;
};

template<> struct ParadigmTraits<MetricTag,BoundedTag> {
    using Weaker = typename InformationTraits<ParadigmCode::BOUNDED>::Paradigm;
    using Stronger = typename InformationTraits<ParadigmCode::METRIC>::Paradigm;
};

template<> struct ParadigmTraits<BoundedTag,MetricTag> {
    using Weaker = ParadigmClass<ParadigmCode::BOUNDED>;
    using Stronger = ParadigmClass<ParadigmCode::METRIC>;
};

template<> struct ParadigmTraits<MetricTag,UpperTag> {
    using Weaker = typename InformationTraits<ParadigmCode::UPPER>::Paradigm;
    using Stronger = typename InformationTraits<ParadigmCode::METRIC>::Paradigm;
};

template<> struct ParadigmTraits<UpperTag,MetricTag> {
    using Weaker = ParadigmClass<ParadigmCode::UPPER>;
    using Stronger = ParadigmClass<ParadigmCode::METRIC>;
};

template<> struct ParadigmTraits<MetricTag,LowerTag> {
    using Weaker = typename InformationTraits<ParadigmCode::LOWER>::Paradigm;
    using Stronger = typename InformationTraits<ParadigmCode::METRIC>::Paradigm;
};

template<> struct ParadigmTraits<LowerTag,MetricTag> {
    using Weaker = ParadigmClass<ParadigmCode::LOWER>;
    using Stronger = ParadigmClass<ParadigmCode::METRIC>;
};

template<class P> struct ParadigmTraits<P,P> {
    using Weaker = typename InformationTraits<P::code()>::Paradigm;
    using Stronger = typename InformationTraits<P::code()>::Paradigm;
    using NextWeaker = typename InformationTraits<next_weaker(P::code())>::Paradigm;
};

template<> struct ParadigmTraits<ApproximateTag,ApproximateTag> {
  private:
    typedef ApproximateTag P;
  public:
    using Weaker = typename InformationTraits<P::code()>::Paradigm;
    using Stronger = typename InformationTraits<P::code()>::Paradigm;
    using NextWeaker = Void;
};




template<class P> using Opposite = ParadigmClass<opposite(P::code())>;
template<class P> using Negated = ParadigmClass<negate(P::code())>;
template<class P> using Inverted = ParadigmClass<invert(P::code())>;
template<class P> using Generic = ParadigmClass<strengthen(P::code())>;
template<class P> using Unsigned = ParadigmClass<unsign(P::code())>;
template<class P> using Signed = ParadigmClass<sign(P::code())>;
template<class P> using Widen = ParadigmClass<widen(P::code())>;
template<class P> using Unorder = ParadigmClass<undirect(P::code())>;
template<class P> using Undirect = ParadigmClass<undirect(P::code())>;
template<class P> using Null = ParadigmClass<null(P::code())>;
template<class P1, class P2=Negated<P1>> using LessThan = Weaker<Negated<P1>,P2>;
template<class P1, class P2=Negated<P1>> using Equality = Null<Weaker<P1,Negated<P2>>>;

} // namespace Ariadne

#endif /* ARIADNE_PARADIGM_HPP */
