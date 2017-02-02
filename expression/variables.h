/***************************************************************************
 *            variables.h
 *
 *  Copyright 2008-16  Pieter Collins
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


/*! \file variables.h
 *  \brief Internal variables
 */

#ifndef ARIADNE_VARIABLES_H
#define ARIADNE_VARIABLES_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "utility/macros.h"
#include "utility/pointer.h"
#include "utility/container.h"
#include "utility/tribool.h"
#include "utility/string.h"

#include "numeric/logical.decl.h"
#include "numeric/number.decl.h"
#include "expression/identifier.h"
#include "expression/operations.h"

namespace Ariadne {

class UntypedVariable;
class ExtendedUntypedVariable;

template<class T> class Variable;
template<class T> class LetVariable;
template<class T> class DottedVariable;
template<class T> class PrimedVariable;

template<class T> class Variables;


template<class T> class Constant;
template<class T> class Expression;
template<class LHS,class RHS> class Assignment;

// Simplifying typedefs
typedef Variable<Boolean> BooleanVariable;
typedef Variable<Kleenean> KleeneanVariable;
typedef Variable<String> StringVariable;
typedef Variable<Integer> IntegerVariable;
typedef Variable<Real> RealVariable;
typedef Variables<Real> RealVariables;

typedef PrimedVariable<String> PrimedStringVariable;
typedef LetVariable<Integer> LetIntegerVariable;
typedef PrimedVariable<Integer> PrimedIntegerVariable;
typedef LetVariable<Real> LetRealVariable;
typedef PrimedVariable<Real> PrimedRealVariable;
typedef DottedVariable<Real> DottedRealVariable;

template<class UB> class Interval;
template<class UB> class VariableInterval;
template<class IVL> class VariablesBox;
typedef Interval<Real> RealInterval;
typedef VariableInterval<Real> RealVariableInterval;
typedef VariablesBox<RealInterval> RealVariablesBox;


enum class VariableType : char { BOOLEAN, KLEENEAN, ENUMERATED, STRING, INTEGER, REAL };

template<class T> inline VariableType variable_type() { ARIADNE_FAIL_MSG("Unknown variable type"); }
template<> inline constexpr VariableType variable_type<Boolean>() { return VariableType::BOOLEAN; }
template<> inline constexpr VariableType variable_type<Kleenean>() { return VariableType::KLEENEAN; }
template<> inline constexpr VariableType variable_type<String>() { return VariableType::STRING; }
template<> inline constexpr VariableType variable_type<Integer>() { return VariableType::INTEGER; }
template<> inline constexpr VariableType variable_type<Real>() { return VariableType::REAL; }

inline String class_name(const VariableType& tp) {
    switch(tp) {
        case VariableType::BOOLEAN: return "Boolean";
        case VariableType::KLEENEAN: return "Kleenean";
        case VariableType::ENUMERATED: return "Enumerated";
        case VariableType::STRING: return "String";
        case VariableType::INTEGER: return "Integer";
        case VariableType::REAL: return "Real";
    }
    assert(false); // To prevent warnings on reaching end of function
}

//! A named variable of unknown type.
class UntypedVariable {
  public:
    //! \brief The name of the variable.
    const Identifier& name() const { return this->_name; }
    const VariableType& type() const { return this->_type; }
    Bool operator==(const UntypedVariable& other) const {
        return (this->name()==other.name()) && (this->type()==other.type()); }
    Bool operator!=(const UntypedVariable& other) const { return !(*this==other); }
    Bool operator<(const UntypedVariable& other) const {
        return this->name()<other.name() || (this->name()==other.name() && this->type() < other.type()); }
    friend OutputStream& operator<<(OutputStream& os, const UntypedVariable& var) {
        return os << var.name(); }
  protected:
    explicit UntypedVariable(const Identifier& nm, VariableType tp)
        : _name(nm), _type(tp) { }
  private:
    Identifier _name;
    VariableType _type;
};


//! \ingroup ExpressionModule
//! \brief A named variable of type \a T.
//! \sa Expression \sa Assignment
template<class T> class Variable
    : public UntypedVariable
    , public DeclareExpressionOperations<T>
{
  public:
    typedef T Type;
    typedef Variable<T> BaseType;
    //! \brief Construct a variable with name \a nm.
    explicit Variable(const Identifier& nm) : UntypedVariable(nm,variable_type<T>()) { }
    Variable<T> const& base() const { return *this; }
    inline Variable<T>& operator=(const Variable<T>& e) = default;
    inline Assignment<Variable<T>,T> operator=(const T& c) const;
    template<class XL, class XU> inline VariableInterval<XU> in(const XL& l, const XU& u);
    template<class XU> inline VariableInterval<XU> in(const Interval<XU>& ivl);
    Expression<T> create_zero() const { return Expression<T>::constant(0); }
  public:
    friend LetVariable<T> let(const Variable<T>&);
    friend PrimedVariable<T> prime(const Variable<T>&);
    friend DottedVariable<Real> dot(const Variable<Real>&);
};

class TimeVariable : public Variable<Real> {
  public:
    TimeVariable() : Variable<Real>(" t ") { }
};

//! \ingroup ExpressionModule
//! \brief A list of variables of type \a T.
//! \sa Variable
template<class T> class Variables : public List<Variable<T>> {
  public:
    Variables(Identifier name, SizeType num) : List<Variable<T>>() {
        this->reserve(num); for(SizeType i=0; i!=num; ++i) { this->append(Variable<T>(name+to_str(i))); } }
    Variables<T>& operator=(Variables<T> const&) = default;
    inline List<Assignment<Variable<T>,T>> operator=(const List<T>& c) const;
    template<class IVL> inline VariablesBox<IVL> in(const List<IVL>& bx) const;
};



enum class VariableCategory : char { SIMPLE, DOTTED, PRIMED };

//! A named variable of type \a T, possibly decorated by a "let", "dot" or "prime"
//! representing a time derivative or updated value.
class ExtendedUntypedVariable
    : public UntypedVariable
{
    VariableCategory _category;
  public:
    VariableCategory category() const { return this->_category; }
  protected:
    explicit ExtendedUntypedVariable(const Identifier& name, VariableType type, VariableCategory category)
        : UntypedVariable(name, type), _category(category) { }
};

inline OutputStream& operator<<(OutputStream& os, ExtendedUntypedVariable const& var) {
    switch(var.category()) {
        case VariableCategory::SIMPLE: os << var.name(); break;
        case VariableCategory::PRIMED: os << "prime("<<var.name()<<")"; break;
        case VariableCategory::DOTTED: os << "dot("<<var.name()<<")"; break;
    }
    return os;
}



template<class T> class ExtendedVariable
    : public ExtendedUntypedVariable
{
  public:
    typedef Variable<T> BaseType;
    Variable<T> base() const { return Variable<T>(this->name()); }
  protected:
    explicit ExtendedVariable(const Variable<T>& var, VariableCategory category) : ExtendedUntypedVariable(var.name(),variable_type<T>(),category) { }
};

//! A named variable of type \a T to be used on the left-hand-side of an assignment denoting an algebraic equation.
template<class T> class LetVariable
    : public ExtendedVariable<T>
{
  public:
    //! \brief Decorate a simple variable to denote the left-hand-side of an algebraic equation.
    friend LetVariable<T> let(const Variable<T>& var) { return LetVariable<T>(var); }
    //! \brief Construct an assignment statement representing the algebraic equation \a var := \a expr.
    inline Assignment<Variable<T>,Expression<T>> operator=(const Expression<T>& expr) const;
    inline Assignment<Variable<T>,Expression<T>> operator=(const Variable<T>& expr) const;
    inline Assignment<Variable<T>,Expression<T>> operator=(const T& c) const;
  private:
    explicit LetVariable(const Variable<T>& var) : ExtendedVariable<T>(var,VariableCategory::SIMPLE) { }
};
template<class T> inline LetVariable<T> set(const Variable<T>& var) { return let(var); }


//! \brief A named variable of type \a T decorated by a prime representing a value after a discrete jump.
template<class T> class PrimedVariable
    : public ExtendedVariable<T>
{
  public:
    //! \brief Decorate a simple variable with a prime.
    friend PrimedVariable<T> prime(const Variable<T>& var) { return PrimedVariable<T>(var); }
    //! \brief Construct an assignment statement representing the differential equation \a var' := \a expr.
    inline Assignment<PrimedVariable<T>,Expression<T>> operator=(const Expression<T>& e) const;
    inline Assignment<PrimedVariable<T>,Expression<T>> operator=(const T& c) const;
  private:
    explicit PrimedVariable(const Variable<T>& var) : ExtendedVariable<T>(var,VariableCategory::PRIMED) { }
};
template<class T> inline PrimedVariable<T> next(const Variable<T>& var) { return prime(var); }


//! \brief A named variable of type \a T decorated by a dot representing differentiation with respect to time.
template<class T> class DottedVariable
    : public ExtendedVariable<T>
{
  public:
    //! \brief Decorate a simple variable with a dot.
    friend DottedVariable<Real> dot(const Variable<Real>& var);
    //! \brief Construct an assignment statement representing the differential equation \a dot(var) := \a expr.
    inline Assignment<DottedVariable<T>,Expression<T>> operator=(const Expression<T>& e) const;
    inline Assignment<DottedVariable<T>,Expression<T>> operator=(const T& c) const;
  private:
    explicit DottedVariable(const Variable<T>& var) : ExtendedVariable<T>(var,VariableCategory::DOTTED) { }
};
DottedVariable<Real> inline dot(const Variable<Real>& var) { return DottedVariable<Real>(var); }



template<class T> struct LetVariables {
    const List<Variable<T>> _lhs;
    LetVariables(const List<Variable<T>>& lhs) : _lhs(lhs) { }
    List<Assignment<Variable<T>,Expression<T>>> operator=(const List<Expression<T>>&);
};
template<class T> inline LetVariables<T> let(const List<Variable<T>>& lhs) { return LetVariables<T>(lhs); }
template<class T> inline LetVariables<T> set(const List<Variable<T>>& lhs) { return LetVariables<T>(lhs); }
inline LetVariables<Real> let(const InitializerList<Variable<Real>>& lhs) { return let(List<Variable<Real>>(lhs)); }
inline LetVariables<Real> set(const InitializerList<Variable<Real>>& lhs) { return set(List<Variable<Real>>(lhs)); }

template<class T> struct PrimedVariables {
    const List<Variable<T>> _lhs;
    PrimedVariables(const List<Variable<T>>& lhs) : _lhs(lhs) { }
    List<Assignment<PrimedVariable<T>,Expression<T>>> operator=(const List<Expression<T>>&);
    friend OutputStream& operator<<(OutputStream& os, PrimedVariables<T> const& dv) { return os <<"prime("<<dv._lhs<<")"; }
};
template<class T> inline PrimedVariables<T> prime(const List<Variable<T>>& lhs) { return PrimedVariables<T>(lhs); }
template<class T> inline PrimedVariables<T> next(const List<Variable<T>>& lhs) { return PrimedVariables<T>(lhs); }
inline PrimedVariables<Real> prime(const InitializerList<Variable<Real>>& lhs) { return prime(List<Variable<Real>>(lhs)); }
inline PrimedVariables<Real> next(const InitializerList<Variable<Real>>& lhs) { return next(List<Variable<Real>>(lhs)); }

template<class T> struct DottedVariables {
    const List<Variable<T>> _lhs;
    DottedVariables(const List<Variable<T>>& lhs) : _lhs(lhs) { }
    List<Assignment<DottedVariable<T>,Expression<T>>> operator=(const List<Expression<T>>&);
    friend OutputStream& operator<<(OutputStream& os, DottedVariables<T> const& dv) { return os <<"dot("<<dv._lhs<<")"; }
};
inline DottedVariables<Real> dot(const List<Variable<Real>>& lhs) { return DottedVariables<Real>(lhs); }
inline DottedVariables<Real> dot(const InitializerList<Variable<Real>>& lhs) { return DottedVariables<Real>(List<Variable<Real>>(lhs)); }

} // namespace Ariadne

#endif /* ARIADNE_VARIABLES_H */
