/***************************************************************************
 *            numeric/logical.h
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

//! \file logical.h
//! \brief %Logical classes

#ifndef ARIADNE_LOGICAL_H
#define ARIADNE_LOGICAL_H

#include "utility/stdlib.h"
#include "utility/typedefs.h"
#include "numeric/paradigm.h"

#include "logical.decl.h"

namespace Ariadne {

class Effort {
    Nat _m;
  public:
    static Effort get_default() { return Effort(0u); }
    explicit Effort(Nat m) : _m(m) { }
    operator Nat() const { return _m; }
};

template<class P> class Logical;
template<> class Logical<Effective>;



//!  \ingroup LogicalTypes
//! \brief An enumeration containing the possible values of a logical variable.
enum class LogicalValue : char {
    FALSE=-2, //!< Definitely not true.
    UNLIKELY=-1, //!< Considered unlikely to be true.
    INDETERMINATE= 0, //!< Truth is unknown, possibly undecidable.
    LIKELY=+1, //!< Considered likely to be true.
    TRUE=+2 //!< Definitely true.
};

inline Bool definitely(const LogicalValue& l) { return l==LogicalValue::TRUE; }
inline Bool possibly(const LogicalValue& l) { return l!=LogicalValue::FALSE; }
inline Bool decide(const LogicalValue& l) { return l==LogicalValue::TRUE || l==LogicalValue::LIKELY; };

LogicalValue equal(LogicalValue l1, LogicalValue l2);
inline LogicalValue disjunction(LogicalValue l1, LogicalValue l2) { return (l1<l2 ? l1 : l2); }
inline LogicalValue conjunction(LogicalValue l1, LogicalValue l2) { return (l1>l2 ? l1 : l2); }
inline LogicalValue negation(LogicalValue l) { return static_cast<LogicalValue>(-static_cast<char>(l)); }

OutputStream& operator<<(OutputStream& os, LogicalValue b);


//!  \ingroup LogicalTypes
//!  \brief A logical variable for the paradigm \a P, which must be %Exact, %Validated, %Upper, %Lower or %Approximate.
//!  Used as a base of named logical types Boolean, Tribool, Sierpinski and Fuzzy. Implemented in terms of LogicalValue.
template<class P> class Logical
{
    template<class PP> friend class Logical;
    LogicalValue _v;
  public:
    explicit Logical(LogicalValue v);
    explicit operator LogicalValue () const { return _v; }
  public:
    constexpr Logical() : Logical(LogicalValue::FALSE) { }
    template<class PP, EnableIf<IsWeaker<P,PP>> =dummy> constexpr Logical(Logical<PP> l) : Logical(l._v) { }
    //! \brief Convert from a builtin boolean value.
    constexpr Logical(Bool b) : Logical(b?LogicalValue::TRUE:LogicalValue::FALSE) { }

    //! \brief Convert to a builtin boolean value. Calls the decide() function.
    // TODO: This should be explicit (except for Boolean)

    // explicit operator Bool () const { return decide(this->_v); }

    //! \brief Equality of two logical values.
    friend inline Logical<P> operator==(Logical<P> l1, Logical<P> l2) { return Logical<P>(equal(l1._v,l2._v)); }
    friend inline Logical<P> operator!=(Logical<P> l1, Logical<P> l2) { return Logical<P>(negation(equal(l1._v,l2._v))); }
    //! \brief %Logical disjunction.
    friend inline Logical<P> operator&&(Logical<P> l1, Logical<P> l2) { return Logical<P>(disjunction(l1._v,l2._v)); }
    //! \brief %Logical conjunction.
    friend inline Logical<P> operator||(Logical<P> l1, Logical<P> l2) { return Logical<P>(conjunction(l1._v,l2._v)); }
    //! \brief %Logical exclusive or.
    friend inline Logical<P> operator^(Logical<P> l1, Logical<P> l2) { return Logical<P>(negation(equal(l1._v,l2._v))); }
    //! \brief %Logical negation.
    friend inline Logical<Negated<P>> operator!(Logical<P> l) { return Logical<Negated<P>>(negation(l._v)); }
    //! \brief Returns \c true only if \a l represents the result of a logical predicate which is definitely true.
    //! Returns \c false for values other than LogicalValue::TRUE.
    friend inline Bool definitely(Logical<P> l) { return l._v == LogicalValue::TRUE; }
    //! \brief Returns \c true only if \a l represents the result of a logical predicate which may be true.
    //!  Returns \c false  only for the value LogicalValue::FALSE.
    friend inline Bool possibly(Logical<P> l) { return l._v != LogicalValue::FALSE; }
    //! \brief Converts the logical value into a true/false boolean context.
    //! Returns \c true for LogicalValue::TRUE or LogicalValue::LIKELY.
    //! Note that decide(indeterminate) is false.
    friend inline Bool decide(Logical<P> l) { return l._v >= LogicalValue::LIKELY; }
    //! \brief Returns \c true if the value is definitely TRUE or FALSE.
    friend inline Bool is_determinate(Logical<P> l) { return l._v == LogicalValue::TRUE || l._v == LogicalValue::FALSE; }
    friend inline Bool is_indeterminate(Logical<P> l) { return l._v != LogicalValue::TRUE && l._v != LogicalValue::FALSE; }
    //! \brief Returns \c true if the values have the same code; does not mean they represent equal results.
    friend inline Bool same(Logical<P> l1, Logical<P> l2) { return l1._v == l2._v; }
    //! \brief Write to an output stream.
    friend inline OutputStream& operator<<(OutputStream& os, Logical<P> l) { return os << l._v; }
  private:
    friend class Tribool;
};


template<> inline Logical<Exact>::Logical(LogicalValue v)
     : _v(v) { assert(v==LogicalValue::FALSE || v==LogicalValue::TRUE); }
template<> inline Logical<Validated>::Logical(LogicalValue v)
    : _v(v) { }
template<> inline Logical<Upper>::Logical(LogicalValue v)
    : _v(v==LogicalValue::FALSE?LogicalValue::UNLIKELY:v) { }
template<> inline Logical<Lower>::Logical(LogicalValue v)
    : _v(v==LogicalValue::TRUE?LogicalValue::LIKELY:v) { }
template<> inline Logical<Approximate>::Logical(LogicalValue v)
    : _v(v==LogicalValue::TRUE?LogicalValue::LIKELY:v==LogicalValue::FALSE?LogicalValue::UNLIKELY:v) { }
class LogicalHandle;
LogicalValue check(LogicalHandle const& l, Effort e);

class LogicalInterface {
    friend LogicalValue check(LogicalHandle const& l, Effort e);
    friend OutputStream& operator<<(OutputStream& os, LogicalHandle const& l);
  public:
    virtual ~LogicalInterface() = default;
  private:
    virtual LogicalValue _check(Effort) const = 0;
    virtual OutputStream& _write(OutputStream&) const = 0;
};

class LogicalHandle {
    SharedPointer<const LogicalInterface> _ptr;
  public:
    explicit LogicalHandle(SharedPointer<const LogicalInterface> p) : _ptr(p) { }
    explicit LogicalHandle(LogicalValue v);
    friend LogicalValue check(LogicalHandle const& l, Effort e) { return l._ptr->_check(e); }
    friend OutputStream& operator<<(OutputStream& os, LogicalHandle const& l) { return l._ptr->_write(os); }
};

template<> class Logical<Effective>
{
    LogicalHandle _v;
    template<class P> friend class Logical;
  public:
    explicit Logical<Effective>(SharedPointer<const LogicalInterface> p) : _v(p) { }
    template<class B, EnableIf<IsSame<B,Bool>> =dummy> Logical(B b) : Logical(Logical<Exact>(b)) { }
    Logical<Effective>(Logical<Exact> l);
    Logical<Validated> check(Effort e) const { return Logical<Validated>(Ariadne::check(_v,e)); }
    friend Bool decide(Logical<Effective> l, Effort e=Effort::get_default()) { return decide(l.check(e)); }
    friend Bool definitely(Logical<Effective> l, Effort e=Effort::get_default()) { return definitely(l.check(e)); }
    friend Bool possibly(Logical<Effective> l, Effort e=Effort::get_default()) { return possibly(l.check(e)); }
};

template<> class Logical<EffectiveUpper>
{
    LogicalHandle _v;
  public:
    explicit Logical<EffectiveUpper>(SharedPointer<const LogicalInterface> p) : _v(p) { }
    Logical<EffectiveUpper>(Logical<Effective> l) : _v(l._v) { }
    Logical<ValidatedUpper> check(Effort e) const { return Logical<ValidatedUpper>(Ariadne::check(_v,e)); }
    friend Bool decide(Logical<EffectiveUpper> l, Effort e=Effort::get_default()) { return decide(l.check(e)); }
};

template<> class Logical<EffectiveLower>
{
    LogicalHandle _v;
  public:
    Logical<EffectiveLower>(SharedPointer<const LogicalInterface> p) : _v(p) { }
    Logical<EffectiveLower>(Logical<Effective> l) : _v(l._v) { }
    Logical<ValidatedLower> check(Effort e) const { return Logical<ValidatedLower>(Ariadne::check(_v,e)); }
    friend Bool decide(Logical<EffectiveLower> l, Effort e=Effort::get_default()) { return decide(l.check(e)); }
};

typedef Logical<Effective> Quasidecidable;
typedef Logical<EffectiveUpper> Verifyable;
typedef Logical<EffectiveLower> Falsifyable;

inline Logical<EffectiveLower> operator&&(Logical<EffectiveLower> l1, Logical<Exact> l2) {
    if(decide(l2)) { return l1; } else { return Logical<Effective>(false); } }
inline Logical<EffectiveUpper> operator||(Logical<EffectiveUpper> l1, Logical<Exact> l2) {
    if(decide(l2)) { return Logical<Effective>(true); } else { return l1; } }

//! \ingroup LogicalTypes
//! \brief The logical constant representing an unknown value.
static const Logical<Validated> indeterminate = Logical<Validated>(LogicalValue::INDETERMINATE);

//! \ingroup LogicalTypes
//! \brief The logical constant representing an value which is deemed likely to be true, but for which truth has not been confirmed.
static const Logical<Approximate> likely = Logical<Approximate>(LogicalValue::LIKELY);

//! \ingroup LogicalTypes
//! \brief The logical constant representing an value which is deemed unlikely to be true, but for which truth has not been ruled out.
static const Logical<Approximate> unlikely = Logical<Approximate>(LogicalValue::UNLIKELY);

inline Bool definitely(Bool b) { return b; }
inline Bool possibly(Bool b) { return b; }
inline Bool decide(Bool b) { return b; }

//! \ingroup LogicalTypes
//! \brief A logical variable representing the result of a decidable proposition.
class Boolean : public Logical<Exact> {
  public:
    Boolean(Bool b=false) : Logical<Exact>(b) { }
    Boolean(Logical<Exact> l) : Logical<Exact>(l) { }
    friend Boolean operator&&(Boolean l1, Boolean l2) { return Logical<Exact>(l1) && Logical<Exact>(l2); }
    friend Boolean operator&&(Boolean l1, Bool l2) { return l1 && Boolean(l2); }
    friend Boolean operator&&(Bool l1, Boolean l2) { return Boolean(l1) && l2; }
    friend Boolean operator||(Boolean l1, Boolean l2) { return Logical<Exact>(l1) || Logical<Exact>(l2); }
    friend Boolean operator||(Boolean l1, Bool l2) { return l1 || Boolean(l2); }
    friend Boolean operator||(Bool l1, Boolean l2) { return Boolean(l1) || l2; }
    friend Boolean operator!(Boolean l) { return !Logical<Exact>(l); }
    operator Bool () const { return decide(*this); }
};

//! \ingroup LogicalTypes
//! \brief A logical variable representing the result of a proposition with some undecidable instances.
//! Takes value \c INDETERMINATE for undecidable instances, for instances for which the information provided is insufficient to obtain a definite result, or for algorithms for which obtaining a result would be deemed to take unacceptably long.
//! The concrete type of the constant indeterminate.
class Tribool : public Logical<Validated> {
 public:
    using Logical<Validated>::Logical;
    Tribool() :  Logical<Validated>(LogicalValue::INDETERMINATE) { }
    Tribool(Logical<Effective> l) : Logical<Validated>(l.check(Effort::get_default())) { }
    // FIXME: Currently needed for Real<Real comparison; Is there a better name?
    explicit Tribool(Logical<Lower> l) : Logical<Validated>(l._v) { }
    explicit Tribool(Logical<Upper> l) : Logical<Validated>(l._v) { }
  public:
    operator Logical<Effective>() const { return Logical<Effective>(*this); }
};

//! \ingroup LogicalTypes
//! \brief A logical variable representing the result of a verifyable proposition. May not take the value FALSE.
class Sierpinski : public Logical<Upper> {
  public:
    using Logical<Upper>::Logical;
    Sierpinski(Logical<EffectiveUpper> l) : Logical<ValidatedUpper>(l.check(Effort::get_default())) { }
};

// TODO: Should this be a user class?
class NegSierpinski : public Logical<Lower> {
  public:
    using Logical<Lower>::Logical;
    NegSierpinski(Logical<EffectiveLower> l) : Logical<ValidatedLower>(l.check(Effort::get_default())) { }
};

//! \ingroup LogicalTypes
//! \brief A logical variable representing the result of a proposition
//! for which not enough information is provided to yield a definite answer.
//! May not take the values \c TRUE or \c FALSE.
//! Takes value \c LIKELY the information suggests that the result is \c true,
//! and \c UNLIKELY the information suggests that the result is \c false.
//! May take value \c INDETERMINATE if the information provided does not strongly suggest either result.
class Fuzzy : public Logical<Approximate> {
    using Logical<Approximate>::Logical;
};


}

#endif
