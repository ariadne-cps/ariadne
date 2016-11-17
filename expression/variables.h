/***************************************************************************
 *            variables.h
 *
 *  Copyright 2008-9  Pieter Collins
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
#include "expression/operations.h"

namespace Ariadne {

class Identifier;

template<class T> class Set;

class String;

//! \ingroup ExpressionModule
//! \brief A class representing the name of a variable.
//! \details A proxy for a standard string; used to distinguish a string used as a variable name from a value.
//! \sa Variable
class Identifier : public String
{
  public:
    Identifier() : String() { }
    Identifier(const char* cstr) : String(cstr) { }
    //! \brief Construct an identifier from a standard string.
    Identifier(const StringType& str) : String(str) { }
};

class UntypedVariable;
template<class T> class ExtendedVariable;
template<class T> class Variable;
template<class T> class LetVariable;
template<class T> class DottedVariable;
template<class T> class PrimedVariable;

template<class T> class Variables;

template<class T> class Constant;
template<class T> class Expression;
template<class LHS,class RHS> class Assignment;

// Simplifying typedefs
typedef Variable<String> StringVariable;
typedef Constant<String> StringConstant;
typedef PrimedVariable<String> PrimedStringVariable;
typedef Variable<Integer> IntegerVariable;
typedef PrimedVariable<Integer> PrimedIntegerVariable;
typedef ExtendedVariable<Real> ExtendedRealVariable;
typedef LetVariable<Real> LetRealVariable;
typedef DottedVariable<Real> DottedRealVariable;
typedef PrimedVariable<Real> PrimedRealVariable;
typedef Variable<Real> RealVariable;
typedef Constant<Real> RealConstant;
typedef Variables<Real> RealVariables;

class RealVariableInterval;
class RealVariablesBox;

//! \ingroup ExpressionModule
//! A named constant of type \a T.
template<class T> class Constant
    : public T
{
  public:
    explicit Constant(const String& str, const T& value)
        : T(value), _name(str) { }
    const Identifier& name() const { return _name; }
    const T& value() const { return *this; }
  private:
    Identifier _name;
};

template<> class Constant<String>
    : public String
{
  public:
    explicit Constant(const String& value) : String(value) { }
    const Identifier& name() const { return static_cast<const Identifier&>(static_cast<const String&>(*this)); }
    const String& value() const { return *this; }
};

enum VariableType { type_bool, type_tribool, type_enumerated, type_string, type_integer, type_real };
enum VariableCategory { simple, dotted, primed };

//! A named variable of unknown type.
class UntypedVariable {
  public:
    //! \brief The name of the variable.
    const Identifier& name() const { return *_name_ptr; }
    const VariableType& type() const { return this->_type; }
    Bool operator==(const UntypedVariable& other) const {
        return (this->name()==other.name()) && (this->_category==other._category); }
    Bool operator!=(const UntypedVariable& other) const { return !(*this==other); }
    Bool operator<(const UntypedVariable& other) const {
        return this->name()<other.name() || (this->name()==other.name() && this->_type < other._type); }
    virtual OutputStream& write(OutputStream&) const;
  public:
    static StringType name(const VariableType& tp) {
        switch(tp) {
            case type_bool: return "Bool";
            case type_tribool: return "Kleenean";
            case type_enumerated: return "Enumerated";
            case type_string: return "String";
            case type_integer: return "Integer";
            case type_real: return "Real";
        }
        return "Unknown";
    }
  protected:
    explicit UntypedVariable(const String& nm, VariableType tp, VariableCategory cat=simple)
        : _name_ptr(new Identifier(nm)), _type(tp), _category(cat) { }
  private:
    std::shared_ptr<Identifier> _name_ptr;
    VariableType _type;
    VariableCategory _category;
};

template<class T> inline VariableType variable_type() { ARIADNE_FAIL_MSG("Unknown variable type"); }
template<> inline VariableType variable_type<Boolean>() { return type_bool; }
template<> inline VariableType variable_type<Kleenean>() { return type_tribool; }
template<> inline VariableType variable_type<String>() { return type_string; }
template<> inline VariableType variable_type<Integer>() { return type_integer; }
template<> inline VariableType variable_type<Real>() { return type_real; }

inline OutputStream& UntypedVariable::write(OutputStream& os) const {
    switch(this->_category) {
        case simple: os << this->name(); break;
        case dotted: os << "dot("<<this->name()<<")"; break;
        case primed: os << "prime("<<this->name()<<")"; break;
    }
    //os << ":" << name(this->_type);
    return os;
}

inline OutputStream& operator<<(OutputStream& os, const UntypedVariable& var) {
    return var.write(os); }



//! A named variable of type \a T, possibly decorated by a "dot" or "prime"
//! representing a time derivative or updated value.
template<class T> class ExtendedVariable
    : public UntypedVariable
    , public DeclareExpressionOperations<T>
{
  public:
    //! \brief The type (class) of data held by the variable.
    typedef T Type;
    Assignment< ExtendedVariable<T>,Expression<T> > operator=(const Expression<T>& e) const;
  protected:
    explicit ExtendedVariable(const String& nm, VariableCategory cat=simple)
        : UntypedVariable(nm, variable_type<T>(), cat) { }
};

template<class R, class D> struct IsRealBuiltin : False { };
template<> struct IsRealBuiltin<Real,int> : True { };
template<> struct IsRealBuiltin<Real,double> : True { };
template<class R, class D, class T=Dummy> using EnableIfRealBuiltin = EnableIf<IsRealBuiltin<R,D>,T>;


//! \ingroup ExpressionModule
//! \brief A named variable of type \a T.
//! \sa Expression \sa Assignment
template<class T> class Variable
    : public ExtendedVariable<T>
{
    typedef Assignment< Variable<T>, Expression<T> > AssignmentType;
    typedef Assignment< Variable<T>, T > ConstantAssignmentType;
  public:
    typedef Variable<T> BaseType;
    typedef T ValueType;
    //! \brief Construct a variable with name \a nm.
    explicit Variable(const Identifier& nm) : ExtendedVariable<T>(nm) { }
    Variable<T> const& base() const { return *this; }
    inline ConstantAssignmentType operator=(const T& e) const;
    inline AssignmentType operator=(const Constant<T>& cnst) const;
    inline AssignmentType operator=(const Variable<T>& e) const;
    inline AssignmentType operator=(const Expression<T>& e) const;
    template<class XL, class XU> inline RealVariableInterval in(const XL& l, const XU& u);
    template<class D> inline EnableIfRealBuiltin<T,D,ConstantAssignmentType> operator=(D e) const;
    Expression<T> create_zero() const { return Expression<T>::constant(0); }
};

class TimeVariable : public Variable<Real> {
  public:
    TimeVariable() : Variable<Real>(" t ") { }
};

//! \ingroup ExpressionModule
//! \brief A list of variables of type \a T.
//! \sa Variable
template<class T> class Variables : public List<Variable<T>> {
    typedef Assignment< Variable<T>, Expression<T> > AssignmentType;
  public:
    Variables(Identifier name, SizeType num) : List<Variable<T>>() {
        this->reserve(num); for(SizeType i=0; i!=num; ++i) { this->append(Variable<T>(name+to_str(i))); } }
    Variables<T>& operator=(Variables<T> const&) = default;
    inline List<AssignmentType> operator=(const List<Expression<T>>& e) const;
    template<class IVL> inline RealVariablesBox in(const List<IVL>& bx) const;
};

template<class T> LetVariable<T> let(const Variable<T>&);

//! A named variable of type \a T to be used on the left-hand-side of an assignment denoting an algebraic equation.
template<class T> class LetVariable
    : public ExtendedVariable<T>
{
    typedef Assignment<Variable<T>,Expression<T>> AssignmentType;
  public:
    typedef Variable<T> BaseType;
    //! \brief Decorate a simple variable to denote the left-hand-side of an algebraic equation.
    friend LetVariable<T> let<>(const Variable<T>&);
    Variable<T> base() const { return Variable<T>(this->name()); }
    inline AssignmentType operator=(const T& val) const;
    inline AssignmentType operator=(const Constant<T>& cnst) const;
    inline AssignmentType operator=(const Variable<T>& var) const;
    //! \brief Construct an assignment statement representing the algebraic equation \a var := \a expr.
    inline Assignment<Variable<T>,Expression<T>> operator=(const Expression<T>& expr) const;
    template<class D> inline EnableIfRealBuiltin<T,D,AssignmentType> operator=(D e) const;
  private:
    explicit LetVariable(const Variable<T>& var) : ExtendedVariable<T>(var.name(),simple) { }
};

//! \relates Variable \brief Decorate a simple variable to denote the left-hand-side of an algebraic equation.
template<class T> inline LetVariable<T> let(const Variable<T>& var) {
    return LetVariable<T>(var); }



DottedVariable<Real> dot(const Variable<Real>&);

template<class T> class DottedVariable;

//! \brief A named variable of type \a T decorated by a dot representing differentiation with respect to time.
template<class T> class DottedVariable
    : public ExtendedVariable<T>
{
    typedef Assignment<DottedVariable<T>,Expression<T>> AssignmentType;
  public:
    typedef Variable<T> BaseType;
    //! \brief Decorate a simple variable with a dot.
    friend DottedVariable<Real> dot(const Variable<Real>&);
    Variable<T> base() const { return Variable<T>(this->name()); }
    inline AssignmentType operator=(const T& e) const;
    inline AssignmentType operator=(const Constant<T>& cnst) const;
    inline AssignmentType operator=(const Variable<T>& e) const;
    //! \brief Construct an assignment statement representing the differential equation \a dot(var) := \a expr.
    inline Assignment<DottedVariable<T>,Expression<T>> operator=(const Expression<T>& e) const;
    template<class D> inline EnableIfRealBuiltin<T,D,AssignmentType> operator=(D e) const;
  private:
    explicit DottedVariable(const Variable<T>& var) : ExtendedVariable<T>(var.name(),dotted) { }
};

//! \relates Variable  \brief Decorate a simple variable with a dot.
inline DottedVariable<Real> dot(const Variable<Real>& var) {
    return DottedVariable<Real>(var); }


template<class T> PrimedVariable<T> prime(const Variable<T>&);
template<class T> PrimedVariable<T> next(const Variable<T>&);

//! \brief A named variable of type \a T decorated by a prime representing a value after a discrete jump.
template<class T> class PrimedVariable
    : public ExtendedVariable<T>
{
    typedef Assignment<PrimedVariable<T>,Expression<T>> AssignmentType;
  public:
    typedef Variable<T> BaseType;
    //! \brief Decorate a simple variable with a prime.
    friend PrimedVariable<T> prime<>(const Variable<T>&);
    Variable<T> base() const { return Variable<T>(this->name()); }
    inline AssignmentType operator=(const T& val) const;
    inline AssignmentType operator=(const Constant<T>& cnst) const;
    inline AssignmentType operator=(const Variable<T>& var) const;
    //! \brief Construct an assignment statement representing the differential equation \a var' := \a expr.
    inline Assignment<PrimedVariable<T>,Expression<T>> operator=(const Expression<T>& expr) const;
    template<class D> inline EnableIfRealBuiltin<T,D,AssignmentType> operator=(D e) const;
  private:
    explicit PrimedVariable(const Variable<T>& var) : ExtendedVariable<T>(var.name(),primed) { }
};

//! \relates Variable \brief Decorate a simple variable with a prime.
template<class T> inline PrimedVariable<T> prime(const Variable<T>& var) {
    return PrimedVariable<T>(var); }
template<class T> inline PrimedVariable<T> next(const Variable<T>& var) {
    return prime(var); }

} // namespace Ariadne

#endif /* ARIADNE_VARIABLES_H */
