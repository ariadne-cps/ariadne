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

#include "operators.h"

namespace Ariadne {


template<class T> class Set;

typedef bool Boolean;
typedef tribool Tribool;
typedef std::string String;
class Integer;
class Real;
class EnumeratedValue;

typedef String Identifier;

template<class T> class Variable;
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

template<class T> class Variable
{
  public:
    explicit Variable(const String& str) : _name_ptr(new String(str)) { }
    const String& name() const { return *_name_ptr; }
    bool operator==(const Variable<T>& other) const { return this->name()==other.name(); }
    template<class X> bool operator==(const Variable<X>& other) const { ARIADNE_ASSERT(this->name()!=other.name()); return false; }
    template<class X> bool operator!=(const Variable<X>& other) const { return !(*this==other); }
    Assignment< Variable<T>,Expression<T> > operator=(const Expression<T>& e) const;
  private:
    shared_ptr<String> _name_ptr;
};

template<class T> inline bool operator<(const Variable<T> v1, const Variable<T>& v2) { return v1.name()<v2.name(); }

inline std::ostream& operator<<(std::ostream& os, const Variable<EnumeratedValue>& v) { return os << v.name() << ":Enumerated"; }
inline std::ostream& operator<<(std::ostream& os, const Variable<String>& v) { return os << v.name() << ":String"; }
inline std::ostream& operator<<(std::ostream& os, const Variable<Integer>& v) { return os << v.name() << ":Integer"; }
inline std::ostream& operator<<(std::ostream& os, const Variable<Real>& v) { return os << v.name() << ":Real"; }

template<class T> class Space;
template<class T> std::ostream& operator<<(std::ostream& os, const Space<T>& spc);

/*! \brief A space defined as a list of named variables of type \a T
 *  \details \sa Variable
 */
template<class T> struct Space
{
    typedef unsigned int SizeType;
    typedef Variable<T> VariableType;
  public:
    //! \brief The trivial space \f$\R^0\f$.
    Space() : _variables() { }
    Space(const Set<String>& vs) : _variables(vs.begin(),vs.end()) { }
    Space(const List<String>& vl) { for(uint i=0; i!=vl.size(); ++i) { this->append(VariableType(vl[i])); } }
    Space(const List<VariableType>& vl) { for(uint i=0; i!=vl.size(); ++i) { this->append(vl[i]); } }

    SizeType size() const { return _variables.size(); }
    //! \brief The dimension of the space.
    SizeType dimension() const { return _variables.size(); }
    //! \brief The \a i<sup>th</sup> named variable.
    const VariableType& operator[](SizeType i) const { return _variables.at(i); }
    const VariableType& variable(SizeType i) const { return _variables.at(i); }

    //! \brief A list giving ordered variables.
    List<VariableType> variables() const { return this->_variables; }
    //! \brief A map giving the index of a given variable.
    Map<VariableType,SizeType> indices() const {
        Map<VariableType,SizeType> indices;
        for(uint i=0; i!=this->_variables.size(); ++i) {
            ARIADNE_ASSERT_MSG(!indices.has_key(_variables[i]),"Repeated variable "<<_variables[i]<<" in space "<<_variables)
            indices.insert(this->_variables[i],i);
        }
        return indices; }

    //! \brief The index of the named variable \a v.
    SizeType index(const VariableType& v) const {
        for(uint i=0; i!=_variables.size(); ++i) {
            if(v==_variables[i]) { return i; } }
        ARIADNE_ASSERT_MSG(false,"Variable "<<v<<" is not in the Space "<<*this);
        return _variables.size(); }
    //! \brief Append the named variable \a v to the variables defining the space; ignores if the variable is already in the space.
    Space<T>& insert(const VariableType& v) {
        for(uint i=0; i!=_variables.size(); ++i) {
            if(_variables[i]==v) { return *this; } }
        _variables.push_back(v); return *this; }
    //! \brief Adjoins the variables in \a spc.
    Space<T>& adjoin(const Space<T>& spc) {
        for(uint i=0; i!=spc._variables.size(); ++i) { this->insert(spc._variables[i]); } return *this; }
    //! \brief Append the named variable \a v to the variables defining the space.
    Space<T>& append(const VariableType& v) {
        for(uint i=0; i!=_variables.size(); ++i) {
            ARIADNE_ASSERT_MSG(_variables[i]!=v,"Variable "<<v<<" is already a variable of the StateSpace "<<*this);
        }
        _variables.push_back(v); return *this; }
    //! \brief Append the named variable \a v to the variables defining the space.
    Space<T>& operator,(const String& v) { return this->append(VariableType(v)); }
  private:
    List<VariableType> _variables;
};

template<class T> inline std::ostream& operator<<(std::ostream& os, const Space<T>& spc) { return os << spc.variables(); }

template<class T> inline Space<T> operator,(const Identifier& v1, const Identifier& v2) {
    Space<T> r; r,Variable<T>(v1),Variable<T>(v2); return r; }

template<class T> inline Space<T> join(const Space<T>& spc1, const Space<T>& spc2) {
    Space<T> r(spc1); r.adjoin(spc2); return r; }


template<class T> Variable<T> variable(const String& s) { return Variable<T>(s); }
template<class T> Space<T> variables(const List<String>& s) { return Space<T>(s); }



} // namespace Ariadne

#endif /* ARIADNE_VARIABLES_H */
