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
    explicit Effort() : _m(0u) { }
    explicit Effort(Nat m) : _m(m) { }
    operator Nat() const { return _m; }
    friend OutputStream& operator<<(OutputStream& os, Effort eff) { return os << "Effort(" << eff._m << ")"; }
};

template<class P> class Logical;

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
inline Bool probably(const LogicalValue& l) { return l==LogicalValue::TRUE || l==LogicalValue::LIKELY; };
inline Bool decide(const LogicalValue& l) { return l==LogicalValue::TRUE || l==LogicalValue::LIKELY; };

LogicalValue equality(LogicalValue l1, LogicalValue l2);
inline LogicalValue negation(LogicalValue l) { return static_cast<LogicalValue>(-static_cast<char>(l)); }
inline LogicalValue conjunction(LogicalValue l1, LogicalValue l2) { return (l1<l2 ? l1 : l2); }
inline LogicalValue disjunction(LogicalValue l1, LogicalValue l2) { return (l1>l2 ? l1 : l2); }
inline LogicalValue exclusive(LogicalValue l1, LogicalValue l2) { return negation(equality(l1,l2)); }
inline LogicalValue check(LogicalValue l) { return l; }

OutputStream& operator<<(OutputStream& os, LogicalValue b);

class LogicalHandle;
LogicalValue check(LogicalHandle const& l, Effort e);
LogicalValue check(LogicalHandle const& l); //!< DEPRECATED

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
    friend LogicalHandle equality(LogicalHandle l1, LogicalHandle l2);
    friend LogicalHandle conjunction(LogicalHandle l1, LogicalHandle l2);
    friend LogicalHandle disjunction(LogicalHandle l1, LogicalHandle l2);
    friend LogicalHandle exclusive(LogicalHandle l1, LogicalHandle l2);
    friend LogicalHandle negation(LogicalHandle l);
    friend LogicalValue check(LogicalHandle const& l, Effort e) { return l._ptr->_check(e); }
    friend LogicalValue check(LogicalHandle const& l) { return check(l, Effort::get_default()); } //!< DEPRECATED
    friend OutputStream& operator<<(OutputStream& os, LogicalHandle const& l) { return l._ptr->_write(os); }
};


template<class P> class LogicalFacade {
};

template<> class LogicalFacade<ExactTag> {
    typedef ExactTag P;
  public:
    operator Bool () const;
  public:
    friend Logical<ExactTag> operator||(Bool b1, Logical<ExactTag> l2);
    friend Logical<ExactTag> operator||(Logical<ExactTag> l1, Bool b2);
    friend Logical<ExactTag> operator&&(Bool b1, Logical<ExactTag> l2);
    friend Logical<ExactTag> operator&&(Logical<ExactTag> l1, Bool b2);
};


template<> class Logical<EffectiveTag>;

//!  \ingroup LogicalTypes
//!  \brief A logical variable for the paradigm \a P, which must be %ExactTag, %ValidatedTag, %UpperTag, %LowerTag or %ApproximateTag.
//!  Used as a base of exact, validated and approximate logical types. Implemented in terms of LogicalValue.
template<class P> class Logical
    : public LogicalFacade<P>
{
    template<class PP> friend class Logical;
    LogicalValue _v;
  public:
    explicit Logical(LogicalValue v);
    explicit operator LogicalValue () const { return _v; }
  public:
    constexpr Logical() : Logical(LogicalValue::FALSE) { }
    template<class PP, EnableIf<IsWeaker<P,PP>> =dummy> constexpr Logical(Logical<PP> l) : Logical(check(l._v)) { }
    //! \brief Convert from a builtin boolean value.
    constexpr Logical(Bool b) : Logical(b?LogicalValue::TRUE:LogicalValue::FALSE) { }

    //! \brief Convert to a builtin boolean value. Calls the decide() function.
    // TODO: This should be explicit (except for Boolean)

    //! \brief Equality of two logical values.
    friend inline Logical<P> operator==(Logical<P> l1, Logical<P> l2) { return Logical<P>(equality(l1._v,l2._v)); }
    friend inline Logical<P> operator!=(Logical<P> l1, Logical<P> l2) { return Logical<P>(negation(equality(l1._v,l2._v))); }
    //! \brief %Logical conjunction [and].
    friend inline Logical<P> operator&&(Logical<P> l1, Logical<P> l2) { return Logical<P>(conjunction(l1._v,l2._v)); }
    //! \brief %Logical disjunction [or].
    friend inline Logical<P> operator||(Logical<P> l1, Logical<P> l2) { return Logical<P>(disjunction(l1._v,l2._v)); }
    //! \brief %Logical exclusive or.
    friend inline Logical<P> operator^(Logical<P> l1, Logical<P> l2) { return Logical<P>(negation(equality(l1._v,l2._v))); }
    //! \brief %Logical negation.
    friend inline Logical<Negated<P>> operator!(Logical<P> const& l) { return Logical<Negated<P>>(negation(l._v)); }
    //! \brief Returns \c true only if \a l represents the result of a logical predicate which is definitely true.
    //! Returns \c false for values other than LogicalValue::TRUE.
    friend inline Bool definitely(Logical<P> l) { return l._v == LogicalValue::TRUE; }
    //! \brief Returns \c true only if \a l represents the result of a logical predicate which may be true.
    //!  Returns \c false  only for the value LogicalValue::FALSE.
    friend inline Bool possibly(Logical<P> l) { return l._v != LogicalValue::FALSE; }
    //! \brief Converts the logical value into a true/false boolean context.
    //! Returns \c true for LogicalValue::TRUE or LogicalValue::LIKELY.
    friend inline Bool probably(Logical<P> l) { return l._v >= LogicalValue::LIKELY; }
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
};

inline LogicalFacade<ExactTag>::operator Bool () const { return decide(static_cast<Logical<ExactTag>const&>(*this)); }

template<> inline Logical<ExactTag>::Logical(LogicalValue v)
     : _v(v) { assert(v==LogicalValue::FALSE || v==LogicalValue::TRUE); }
template<> inline Logical<ValidatedTag>::Logical(LogicalValue v)
    : _v(v) { }
template<> inline Logical<UpperTag>::Logical(LogicalValue v)
    : _v(v==LogicalValue::FALSE?LogicalValue::UNLIKELY:v) { }
template<> inline Logical<LowerTag>::Logical(LogicalValue v)
    : _v(v==LogicalValue::TRUE?LogicalValue::LIKELY:v) { }
template<> inline Logical<ApproximateTag>::Logical(LogicalValue v)
    : _v(v==LogicalValue::TRUE?LogicalValue::LIKELY:v==LogicalValue::FALSE?LogicalValue::UNLIKELY:v) { }

template<> class Logical<EffectiveTag>
{
    LogicalHandle _v;
    template<class P> friend class Logical;
  public:
    explicit Logical<EffectiveTag>() : _v(LogicalValue::INDETERMINATE) { }
    explicit Logical<EffectiveTag>(SharedPointer<const LogicalInterface> p) : _v(p) { }
    explicit Logical<EffectiveTag>(LogicalValue v) : _v(LogicalHandle(v)) { }
    explicit Logical<EffectiveTag>(LogicalHandle h) : _v(h) { }
    explicit operator LogicalHandle () const { return _v; }
    template<class B, EnableIf<IsSame<B,Bool>> =dummy> Logical(B b) : Logical(Logical<ExactTag>(b)) { }
    Logical<EffectiveTag>(Logical<ExactTag> l) : _v(static_cast<LogicalValue>(l)) { };
    Logical<ValidatedTag> check(Effort e) const { return Logical<ValidatedTag>(Ariadne::check(_v,e)); }
    friend Logical<ValidatedTag> check(Logical<EffectiveTag> l, Effort e) { return Logical<ValidatedTag>(Ariadne::check(l._v,e)); }
    friend Logical<EffectiveTag> operator==(Logical<EffectiveTag> l1, Logical<EffectiveTag> l2) {
        return Logical<EffectiveTag>(equality(l1._v,l2._v)); }
    friend Logical<EffectiveTag> operator&&(Logical<EffectiveTag> l1, Logical<EffectiveTag> l2) {
        return Logical<EffectiveTag>(conjunction(l1._v,l2._v)); }
    friend Logical<EffectiveTag> operator||(Logical<EffectiveTag> l1, Logical<EffectiveTag> l2) {
        return Logical<EffectiveTag>(disjunction(l1._v,l2._v)); }
    friend Logical<EffectiveTag> operator^(Logical<EffectiveTag> l1, Logical<EffectiveTag> l2) {
        return Logical<EffectiveTag>(exclusive(l1._v,l2._v)); }
    friend Logical<EffectiveTag> operator!(Logical<EffectiveTag> const& l) {
        return Logical<EffectiveTag>(negation(l._v)); }
    friend Bool decide(Logical<EffectiveTag> l, Effort e) { return decide(l.check(e)); }
    friend Bool definitely(Logical<EffectiveTag> l, Effort e) { return definitely(l.check(e)); }
    friend Bool possibly(Logical<EffectiveTag> l, Effort e) { return possibly(l.check(e)); }
    friend Bool decide(Logical<EffectiveTag> l) { return decide(l.check(Effort::get_default())); }  //!< DEPRECATED
    friend Bool definitely(Logical<EffectiveTag> l) { return definitely(l.check(Effort::get_default())); }  //!< DEPRECATED
    friend Bool possibly(Logical<EffectiveTag> l) { return possibly(l.check(Effort::get_default())); }  //!< DEPRECATED
    friend inline OutputStream& operator<<(OutputStream& os, Logical<EffectiveTag> l) { return os << l._v; }
};

template<> class Logical<EffectiveUpperTag>
{
    LogicalHandle _v;
    template<class P> friend class Logical;
  public:
    explicit Logical<EffectiveUpperTag>(SharedPointer<const LogicalInterface> p) : _v(p) { }
    explicit Logical<EffectiveUpperTag>(LogicalHandle h) : _v(h) { }
    explicit operator LogicalHandle () const { return _v; }
    Logical<EffectiveUpperTag>(Logical<EffectiveTag> l) : _v(l._v) { }
    Logical<ValidatedUpperTag> check(Effort e) const { return Logical<ValidatedUpperTag>(Ariadne::check(_v,e)); }
    friend Logical<EffectiveUpperTag> operator&&(Logical<EffectiveUpperTag> l1, Logical<EffectiveUpperTag> l2) {
        return Logical<EffectiveUpperTag>(conjunction(l1._v,l2._v)); }
    friend Logical<EffectiveUpperTag> operator||(Logical<EffectiveUpperTag> l1, Logical<EffectiveUpperTag> l2) {
        return Logical<EffectiveUpperTag>(disjunction(l1._v,l2._v)); }
    friend Logical<EffectiveLowerTag> operator!(Logical<EffectiveUpperTag> const& l);
    friend Logical<EffectiveUpperTag> operator!(Logical<EffectiveLowerTag> const& l);
    friend Logical<ValidatedUpperTag> check(Logical<EffectiveUpperTag> l, Effort e) { return Logical<ValidatedUpperTag>(Ariadne::check(l._v,e)); }
    friend Bool decide(Logical<EffectiveUpperTag> l, Effort e) { return decide(l.check(e)); }
    friend Bool decide(Logical<EffectiveUpperTag> l) { return decide(l.check(Effort::get_default())); } //!< DEPRECATED
    friend inline OutputStream& operator<<(OutputStream& os, Logical<EffectiveUpperTag> l) { return os << l._v; }
};

template<> class Logical<EffectiveLowerTag>
{
    LogicalHandle _v;
    template<class P> friend class Logical;
  public:
    Logical<EffectiveLowerTag>(SharedPointer<const LogicalInterface> p) : _v(p) { }
    explicit Logical<EffectiveLowerTag>(LogicalHandle h) : _v(h) { }
    explicit operator LogicalHandle () const { return _v; }
    Logical<EffectiveLowerTag>(Logical<EffectiveTag> l) : _v(l._v) { }
    Logical<ValidatedLowerTag> check(Effort e) const { return Logical<ValidatedLowerTag>(Ariadne::check(_v,e)); }
    friend Logical<ValidatedLowerTag> check(Logical<EffectiveLowerTag> l, Effort e) { return Logical<ValidatedLowerTag>(Ariadne::check(l._v,e)); }
    friend Logical<EffectiveLowerTag> operator&&(Logical<EffectiveLowerTag> l1, Logical<EffectiveLowerTag> l2) {
        return Logical<EffectiveLowerTag>(conjunction(l1._v,l2._v)); }
    friend Logical<EffectiveLowerTag> operator||(Logical<EffectiveLowerTag> l1, Logical<EffectiveLowerTag> l2) {
        return Logical<EffectiveLowerTag>(disjunction(l1._v,l2._v)); }
    friend Logical<EffectiveUpperTag> operator!(Logical<EffectiveLowerTag> const& l) {
        return Logical<EffectiveUpperTag>(negation(l._v)); }
    friend Logical<EffectiveLowerTag> operator!(Logical<EffectiveUpperTag> const& l) {
        return Logical<EffectiveLowerTag>(negation(l._v)); }
    friend Bool decide(Logical<EffectiveLowerTag> l, Effort e) { return decide(l.check(e)); }
    friend Bool decide(Logical<EffectiveLowerTag> l) { return decide(l.check(Effort::get_default())); } //!< DEPRECATED
    friend inline OutputStream& operator<<(OutputStream& os, Logical<EffectiveLowerTag> l) { return os << l._v; }
};

typedef Logical<EffectiveTag> Quasidecidable;
typedef Logical<EffectiveUpperTag> Verifyable;
typedef Logical<EffectiveLowerTag> Falsifyable;

inline Logical<EffectiveLowerTag> operator&&(Logical<EffectiveLowerTag> l1, Logical<ExactTag> l2) {
    if(decide(l2)) { return l1; } else { return Logical<EffectiveTag>(false); } }
inline Logical<EffectiveUpperTag> operator||(Logical<EffectiveUpperTag> l1, Logical<ExactTag> l2) {
    if(decide(l2)) { return Logical<EffectiveTag>(true); } else { return l1; } }

//! \ingroup LogicalTypes
//! \brief The logical constant representing an unknown value.
extern const Logical<EffectiveTag> indeterminate;

//! \ingroup LogicalTypes
//! \brief The logical constant representing an value which is deemed likely to be true, but for which truth has not been confirmed.
static const Logical<ApproximateTag> likely = Logical<ApproximateTag>(LogicalValue::LIKELY);

//! \ingroup LogicalTypes
//! \brief The logical constant representing an value which is deemed unlikely to be true, but for which truth has not been ruled out.
static const Logical<ApproximateTag> unlikely = Logical<ApproximateTag>(LogicalValue::UNLIKELY);

inline Bool definitely(Bool b) { return b; }
inline Bool possibly(Bool b) { return b; }
inline Bool decide(Bool b) { return b; }

}

#endif
