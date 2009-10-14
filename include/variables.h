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

#include "macros.h"
#include "pointer.h"
#include "container.h"

namespace Ariadne {


template<class T> class Set;

typedef bool Boolean;
typedef tribool Tribool;
typedef std::string String;
class Integer;
class Real;
class EnumeratedValue;

typedef String Identifier;

class UntypedVariable;
template<class T> class ExtendedVariable;
template<class T> class Variable;
template<class T> class DottedVariable;
template<class T> class PrimedVariable;

template<class R> class Expression;
template<class LHS,class RHS> class Assignment;


template<class T> class Constant
{
  public:
    template<class X> explicit Constant(const String& str, const X& value)
        : _name_ptr(new String(str)), _value_ptr(new T(value)) { }
    //explicit Constant(const String& str, const T& value)
    //    : _name_ptr(new String(str)), _value_ptr(new T(value)) { }
    const String& name() const { return *_name_ptr; }
    const T& value() const { return *_value_ptr; }
    bool operator==(const Constant<T>& other) const {
        if(this->name()==other.name()) { assert(this->value()==other.value()); return true; } else { return false; } }
  private:
    shared_ptr<String> _name_ptr;
    shared_ptr<T> _value_ptr;
};

enum VariableType { type_bool, type_tribool, type_enumerated, type_string, type_integer, type_real };
enum VariableCategory { simple, dotted, primed };

class UntypedVariable {
  public:
    const String& name() const { return *_name_ptr; }
    bool operator==(const UntypedVariable& other) const {
        return (this->name()==other.name()) && (this->_category==other._category); }
    bool operator!=(const UntypedVariable& other) const { return !(*this==other); }
    bool operator<(const UntypedVariable& other) const {
        return this->name()<other.name() || (this->name()==other.name() && this->_type < other._type); }
    virtual std::ostream& write(std::ostream&) const;
  public:
    static std::string name(const VariableType& tp) {
        switch(tp) {
            case type_bool: return "Bool";
            case type_tribool: return "Tribool";
            case type_enumerated: return "Enumerated";
            case type_string: return "String";
            case type_integer: return "Integer";
            case type_real: return "Real";
        }
        return "Unknown";
    }
  protected:
    explicit UntypedVariable(const std::string& nm, VariableType tp, VariableCategory cat=simple) 
        : _name_ptr(new std::string(nm)), _type(tp), _category(cat) { }
  private:
    shared_ptr<std::string> _name_ptr;
    VariableType _type;
    VariableCategory _category;
};

template<class T> inline VariableType variable_type() { ARIADNE_FAIL_MSG("Unknown variable type"); }
template<> inline VariableType variable_type<bool>() { return type_bool; }
template<> inline VariableType variable_type<tribool>() { return type_tribool; }
template<> inline VariableType variable_type<String>() { return type_string; }
template<> inline VariableType variable_type<Integer>() { return type_integer; }
template<> inline VariableType variable_type<Real>() { return type_real; }

inline std::ostream& UntypedVariable::write(std::ostream& os) const {
    switch(this->_category) {
        case simple: os << this->name(); break;
        case dotted: os << "dot("<<this->name()<<")"; break;
        case primed: os << "next("<<this->name()<<")"; break;
    }
    os << ":" << name(this->_type);
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const UntypedVariable& var) {
    return var.write(os); }



template<class T> class ExtendedVariable
    : public UntypedVariable
{
  public:
    Assignment< ExtendedVariable<T>,Expression<T> > operator=(const Expression<T>& e) const;
  protected:
    explicit ExtendedVariable(const String& nm, VariableCategory cat=simple)
        : UntypedVariable(nm, variable_type<T>(), cat) { }
};

template<class T> class Variable
    : public ExtendedVariable<T>
{
  public:
    explicit Variable(const String& nm) : ExtendedVariable<T>(nm) { }
    Assignment< Variable<T>, Expression<T> > operator=(const double& e) const;
    Assignment< Variable<T>, Expression<T> > operator=(const T& e) const;
    Assignment< Variable<T>, Expression<T> > operator=(const Variable<T>& e) const;
    Assignment< Variable<T>, Expression<T> > operator=(const Expression<T>& e) const;
};


DottedVariable<Real> dot(const Variable<Real>&);

template<class T> class DottedVariable
    : public ExtendedVariable<T>
{
  public:
    friend DottedVariable<Real> dot(const Variable<Real>&);
    Variable<T> base() const { return Variable<T>(this->name()); }
    Assignment< DottedVariable<T>, Expression<T> > operator=(const double& e) const;
    Assignment< DottedVariable<T>, Expression<T> > operator=(const T& e) const;
    Assignment< DottedVariable<T>, Expression<T> > operator=(const Variable<T>& e) const;
    Assignment< DottedVariable<T>, Expression<T> > operator=(const Expression<T>& e) const;
  private:
    explicit DottedVariable(const Variable<T>& var) : ExtendedVariable<T>(var.name(),dotted) { }
};

inline DottedVariable<Real> dot(const Variable<Real>& var) {
    return DottedVariable<Real>(var); }


template<class T> PrimedVariable<T> next(const Variable<T>&);

template<class T> class PrimedVariable
    : public ExtendedVariable<T>
{
    typedef Assignment< PrimedVariable<T>, Expression<T> > AssignmentType;
  public:
    friend PrimedVariable<T> next<>(const Variable<T>&);
    Variable<T> base() const { return Variable<T>(this->name()); }
    AssignmentType operator=(const double& val) const;
    AssignmentType operator=(const T& val) const;
    AssignmentType operator=(const Variable<T>& var) const;
    AssignmentType operator=(const Expression<T>& expr) const;
  private:
    explicit PrimedVariable(const Variable<T>& var) : ExtendedVariable<T>(var.name(),primed) { }
};

template<class T> inline PrimedVariable<T> next(const Variable<T>& var) {
    return PrimedVariable<T>(var); }

inline std::ostream& operator<<(std::ostream& os, const Variable<EnumeratedValue>& v) { return os << v.name() << ":Enumerated"; }
inline std::ostream& operator<<(std::ostream& os, const Variable<String>& v) { return os << v.name() << ":String"; }
inline std::ostream& operator<<(std::ostream& os, const Variable<Integer>& v) { return os << v.name() << ":Integer"; }
inline std::ostream& operator<<(std::ostream& os, const Variable<Real>& v) { return os << v.name() << ":Real"; }
inline std::ostream& operator<<(std::ostream& os, const Variable<Boolean>& v) { return os << v.name() << ":Boolean"; }
inline std::ostream& operator<<(std::ostream& os, const Variable<Tribool>& v) { return os << v.name() << ":Tribool"; }


} // namespace Ariadne

#endif /* ARIADNE_VARIABLES_H */
