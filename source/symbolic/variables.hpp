/***************************************************************************
 *            symbolic/variables.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file symbolic/variables.hpp
 *  \brief Internal variables
 */

#ifndef ARIADNE_VARIABLES_HPP
#define ARIADNE_VARIABLES_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "../utility/macros.hpp"
#include "../utility/pointer.hpp"
#include "../utility/container.hpp"
#include "../utility/tribool.hpp"
#include "../utility/string.hpp"

#include "../numeric/logical.decl.hpp"
#include "../numeric/number.decl.hpp"
#include "../symbolic/expression.decl.hpp"
#include "../symbolic/identifier.hpp"
#include "../symbolic/operations.hpp"

namespace Ariadne {

class UntypedVariable;
class ExtendedUntypedVariable;


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
        default:
            ARIADNE_FAIL_MSG("Unhandled VariableType for output stream\n");
    }
}

//! \ingroup SymbolicModule
//! \brief A named variable of unknown type.
//! \see Variable
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


//! \ingroup SymbolicModule
//! \brief A named variable of type \a T.
//! %Ariadne supports variables of type Boolean, Kleenean, String, Integer and Real.
//! \see TimeVariable \see Constant, Expression, Assignment, Space
template<class T> class Variable
    : public UntypedVariable
    , public DeclareExpressionOperations<T>
{
  public:
    typedef T Type;
    typedef Variable<T> BaseType;
    //! \brief Construct a variable with name \a name.
    explicit Variable(const Identifier& name) : UntypedVariable(name,variable_type<T>()) { }
    Variable<T> const& base() const { return *this; }
    inline Variable<T>& operator=(const Variable<T>& v) = default;
    inline Assignment<Variable<T>,T> operator=(const T& c) const;
    template<class XL, class XU> inline VariableInterval<XU> in(const XL& l, const XU& u);
    template<class XU> inline VariableInterval<XU> in(const Interval<XU>& ivl);
    Expression<T> create_zero() const { return Expression<T>::constant(0); }
  public:
    friend LetVariable<T> let(const Variable<T>&);
    friend PrimedVariable<T> prime(const Variable<T>&);
    friend DottedVariable<Real> dot(const Variable<Real>&);
};

//! \ingroup SymbolicModule
//! \brief A special variable representing time.
class TimeVariable : public Variable<Real> {
  public:
    TimeVariable() : Variable<Real>(" t ") { }
};

//@{
//! \related Variable \name Type synonyms.
using BooleanVariable = Variable<Boolean>; //!< .
using KleeneanVariable = Variable<Kleenean>; //!< .
using StringVariable = Variable<String>; //!< .
using IntegerVariable = Variable<Integer>; //!< .
using RealVariable = Variable<Real>; //!< .
//@}

//! \ingroup SymbolicModule
//! \brief A list of variables of type \a T.
//! \sa Variable
template<class T> class Variables : public List<Variable<T>> {
  public:
    //! \brief Construct \a n variables with name \a name.
    Variables(Identifier name, SizeType n) : List<Variable<T>>() {
        this->reserve(n); for(SizeType i=0; i!=n; ++i) { this->append(Variable<T>(name+to_str(i))); } }
    Variables<T>& operator=(Variables<T> const&) = default;
    inline List<Assignment<Variable<T>,T>> operator=(const List<T>& c) const;
    template<class IVL> inline VariablesBox<IVL> in(const List<IVL>& bx) const;
};



enum class VariableCategory : char { SIMPLE, DOTTED, PRIMED };

inline OutputStream& operator<<(OutputStream& os, VariableCategory const& cat) {
    switch(cat) {
        case VariableCategory::SIMPLE: os << "SIMPLE"; break;
        case VariableCategory::PRIMED: os << "DOTTED"; break;
        case VariableCategory::DOTTED: os << "PRIMED"; break;
        default:
            ARIADNE_FAIL_MSG("Unhandled VariableCategory for output streaming\n");
    }
    return os;
}

// A named variable of type \a T, possibly decorated by a "let", "dot" or "prime"
// representing a time derivative or updated value.
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
        default:
            ARIADNE_FAIL_MSG("Unhandled VariableCategory "<<var.category()<<"\n");
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

//! \ingroup SymbolicModule
//! A named variable of type \a T to be used on the left-hand-side of an assignment denoting an algebraic equation.
//! \relates Variable
template<class T> class LetVariable
    : public ExtendedVariable<T>
{
  public:
    //! \brief Decorate a simple variable to denote the left-hand-side of an algebraic equation.
    friend LetVariable<T> let(const Variable<T>& var) { return LetVariable<T>(var); }
    //! \brief Construct an assignment statement representing the algebraic equation \a var := \a expr.
    inline Assignment<Variable<T>,Expression<T>> operator=(const Expression<T>& expr) const;
  private:
    explicit LetVariable(const Variable<T>& var) : ExtendedVariable<T>(var,VariableCategory::SIMPLE) { }
};
template<class T> inline LetVariable<T> set(const Variable<T>& var) { return let(var); }


//! \ingroup SymbolicModule
//! \brief A named variable of type \a T decorated by a prime representing a value after a discrete jump.
template<class T> class PrimedVariable
    : public ExtendedVariable<T>
{
  public:
    //! \brief Decorate a simple variable with a prime.
    friend PrimedVariable<T> prime(const Variable<T>& var) { return PrimedVariable<T>(var); }
    //! \brief Construct an assignment statement representing the differential equation \a var' := \a expr.
    inline Assignment<PrimedVariable<T>,Expression<T>> operator=(const Expression<T>& e) const;
  private:
    explicit PrimedVariable(const Variable<T>& var) : ExtendedVariable<T>(var,VariableCategory::PRIMED) { }
};
template<class T> inline PrimedVariable<T> next(const Variable<T>& var) { return prime(var); }


//! \ingroup SymbolicModule
//! \brief A named variable of type \a T decorated by a dot representing differentiation with respect to time.
template<class T> class DottedVariable
    : public ExtendedVariable<T>
{
  public:
    //! \brief Decorate a simple variable with a dot.
    friend DottedVariable<Real> dot(const Variable<Real>& var);
    //! \brief Construct an assignment statement representing the differential equation \a dot(var) := \a expr.
    inline Assignment<DottedVariable<T>,Expression<T>> operator=(const Expression<T>& e) const;
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

#endif /* ARIADNE_VARIABLES_HPP */
