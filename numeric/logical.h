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

/*! \file numeric/logical.h
 *  \brief
 */



//! \file logical.h
//! \brief %Logical classes

#ifndef ARIADNE_LOGICAL_H
#define ARIADNE_LOGICAL_H

#include "utility/stdlib.h"
#include "numeric/paradigm.h"

#include "logical.decl.h"

//! \brief Main %Ariadne namespace.
namespace Ariadne {

using OutputStream = std::ostream;

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

inline bool definitely(const LogicalValue& l) { return l==LogicalValue::TRUE; }
inline bool possibly(const LogicalValue& l) { return l!=LogicalValue::FALSE; }
inline bool decide(const LogicalValue& l) { return l==LogicalValue::TRUE || l==LogicalValue::LIKELY; };

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
    constexpr Logical(bool b) : Logical(b?LogicalValue::TRUE:LogicalValue::FALSE) { }

    //! \brief Convert to a builtin boolean value. Calls the decide() function.
    // TODO: This should be explicit (except for Boolean)
    explicit operator bool () const { return decide(this->_v); }

    //! \brief Equality of two logical values.
    friend inline Logical<P> operator==(Logical<P> l1, Logical<P> l2) { return Logical<P>(equal(l1._v,l2._v)); }
    //! \brief %Logical disjunction.
    friend inline Logical<P> operator&&(Logical<P> l1, Logical<P> l2) { return Logical<P>(disjunction(l1._v,l2._v)); }
    //! \brief %Logical conjunction.
    friend inline Logical<P> operator||(Logical<P> l1, Logical<P> l2) { return Logical<P>(conjunction(l1._v,l2._v)); }
    //! \brief %Logical negation.
    friend inline Logical<Negated<P>> operator!(Logical<P> l) { return Logical<Negated<P>>(negation(l._v)); }
    //! \brief Returns \c true only if \a l represents the result of a logical predicate which is definitely true.
    //! Returns \c false for values other than LogicalValue::TRUE.
    friend inline bool definitely(Logical<P> l) { return l._v == LogicalValue::TRUE; }
    //! \brief Returns \c true only if \a l represents the result of a logical predicate which may be true.
    //!  Returns \c false  only for the value LogicalValue::FALSE.
    friend inline bool possibly(Logical<P> l) { return l._v != LogicalValue::FALSE; }
    //! \brief Converts the logical value into a true/false boolean context.
    //! Returns \c true for LogicalValue::TRUE or LogicalValue::LIKELY.
    //! Note that decide(indeterminate) is false.
    friend inline bool decide(Logical<P> l) { return l._v >= LogicalValue::LIKELY; }
    //! \brief Write to an output stream.
    friend inline OutputStream& operator<<(OutputStream& os, Logical<P> l) { return os << l._v; }
};

template<> inline Logical<Exact>::Logical(LogicalValue v)
     : _v(v) { assert(v==LogicalValue::FALSE || v==LogicalValue::TRUE); }
template<> inline Logical<Effective>::Logical(LogicalValue v)
    : _v(v) { }
template<> inline Logical<Validated>::Logical(LogicalValue v)
    : _v(v) { }
template<> inline Logical<Upper>::Logical(LogicalValue v)
    : _v(v==LogicalValue::FALSE?LogicalValue::UNLIKELY:v) { }
template<> inline Logical<Lower>::Logical(LogicalValue v)
    : _v(v==LogicalValue::TRUE?LogicalValue::LIKELY:v) { }
template<> inline Logical<Approximate>::Logical(LogicalValue v)
    : _v(v==LogicalValue::TRUE?LogicalValue::LIKELY:v==LogicalValue::FALSE?LogicalValue::UNLIKELY:v) { }

//! \ingroup LogicalTypes
//! \brief The logical constant representing an unknown value.
static const Logical<Validated> indeterminate = Logical<Validated>(LogicalValue::INDETERMINATE);

//! \ingroup LogicalTypes
//! \brief The logical constant representing an value which is deemed likely to be true, but for which truth has not been confirmed.
static const Logical<Approximate> likely = Logical<Approximate>(LogicalValue::LIKELY);

//! \ingroup LogicalTypes
//! \brief The logical constant representing an value which is deemed unlikely to be true, but for which truth has not been ruled out.
static const Logical<Approximate> unlikely = Logical<Approximate>(LogicalValue::UNLIKELY);

inline bool definitely(bool b) { return b; }
inline bool possibly(bool b) { return b; }
inline bool decide(bool b) { return b; }

//! \ingroup LogicalTypes
//! \brief A logical variable representing the result of a decidable proposition.
class Boolean : public Logical<Exact> {
  public:
    Boolean(bool b=false) : Logical<Exact>(b) { }
    Boolean(Logical<Exact> l) : Logical<Exact>(l) { }
    operator bool () const { return this->Logical<Exact>::operator bool(); }
};

//! \ingroup LogicalTypes
//! \brief A logical variable representing the result of a proposition with some undecidable instances.
//! Takes value \c INDETERMINATE for undecidable instances, for instances for which the information provided is insufficient to obtain a definite result, or for algorithms for which obtaining a result would be deemed to take unacceptably long.
//! The concrete type of the constant indeterminate.
class Tribool : public Logical<Validated> {
    using Logical<Validated>::Logical;
    // FIXME: Currently needed for Real<Real comparison; Is there a better name?
    public: operator Logical<Effective>() const { return Logical<Effective>(*this); }
};

//! \ingroup LogicalTypes
//! \brief A logical variable representing the result of a verifyable proposition. May not take the value FALSE.
class Sierpinski : public Logical<Upper> {
    using Logical<Upper>::Logical;
};

// TODO: Should this be a user class?
class NegSierpinski : public Logical<Lower> { using Logical<Lower>::Logical; };

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
