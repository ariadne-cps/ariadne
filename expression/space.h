/***************************************************************************
 *            space.h
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


/*! \file space.h
 *  \brief Spaces formed by variables.
 */

#ifndef ARIADNE_SPACE_H
#define ARIADNE_SPACE_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "utility/macros.h"
#include "utility/pointer.h"
#include "utility/container.h"

#include "expression/variables.h"

namespace Ariadne {


template<class T> class Space;
template<class T> OutputStream& operator<<(OutputStream& os, const Space<T>& spc);

class Real;
typedef Space<Real> RealSpace;

//! \ingroup ExpressionModule
//! \brief A space defined as a list of named variables of type \a T.
//!   Allows conversion between Euclidean space \f$\mathbb{R}^n\f$ and a space defined by named variables.
//!  \details \sa Variable
template<class T> class Space
{
  public:
    typedef Variable<T> VariableType;
  public:
    //! \brief The trivial space \f$\R^0\f$.
    Space() : _variables() { }
    Space(const List<VariableType>& vl) { for(Nat i=0; i!=vl.size(); ++i) { this->append(vl[i]); } }
    Space(const InitializerList<VariableType>& vl) { for(Nat i=0; i!=vl.size(); ++i) { this->append(vl.begin()[i]); } }

    Bool operator==(const Space<T>& other) const { return this->_variables==other._variables; }
    Bool operator!=(const Space<T>& other) const { return !(*this == other); }

    SizeType size() const { return _variables.size(); }
    //! \brief The dimension of the space.
    SizeType dimension() const { return _variables.size(); }
    //! \brief The \a i<sup>th</sup> named variable.
    const VariableType operator[](SizeType i) const { return VariableType(_variables.at(i)); }
    const VariableType variable(SizeType i) const { return VariableType(_variables.at(i)); }

    //! \brief A list giving ordered variables.
    List<Identifier> variable_names() const { return this->_variables; }
    //! \brief A list giving ordered variables.
    List<VariableType> variables() const { return List<VariableType>(this->_variables); }
    //! \brief A map giving the index of a given variable.
    Map<Identifier,SizeType> indices_from_names() const {
        Map<Identifier,SizeType> indices;
        for(Nat i=0; i!=this->_variables.size(); ++i) {
            ARIADNE_ASSERT_MSG(!indices.has_key(_variables[i]),"Repeated variable "<<_variables[i]<<" in space "<<_variables)
            indices.insert(this->_variables[i],i);
        }
        return indices; }
    //! \brief A map giving the index of a given variable.
    Map<VariableType,SizeType> indices() const {
        Map<VariableType,SizeType> indices;
        for(Nat i=0; i!=this->_variables.size(); ++i) {
            ARIADNE_ASSERT_MSG(!indices.has_key(VariableType(_variables[i])),"Repeated variable "<<_variables[i]<<" in space "<<_variables)
            indices.insert(VariableType(this->_variables[i]),i);
        }
        return indices; }

    //! \brief Tests if the variable \a v is in the space.
    Bool contains(const VariableType& v) const {
        for(Nat i=0; i!=_variables.size(); ++i) {
            if(v.name()==_variables[i]) { return true; } }
        return false; }
    //! \brief The index of the named variable \a v.
    SizeType index(const VariableType& v) const {
        for(Nat i=0; i!=_variables.size(); ++i) {
            if(v.name()==_variables[i]) { return i; } }
        ARIADNE_ASSERT_MSG(false,"Variable "<<v<<" is not in the Space "<<*this);
        return _variables.size(); }
    SizeType index(const String& n) const {
        for(Nat i=0; i!=_variables.size(); ++i) {
            if(n==_variables[i]) { return i; } }
        ARIADNE_ASSERT_MSG(false,"Variable named "<<n<<" is not in the Space "<<*this);
        return _variables.size(); }
    //! \brief Append the named variable \a v to the variables defining the space; ignores if the variable is already in the space.
    Space<T>& insert(const VariableType& v) {
        for(Nat i=0; i!=_variables.size(); ++i) {
            if(_variables[i]==v.name()) { return *this; } }
        _variables.push_back(v.name()); return *this; }
    //! \brief Adjoins the variables in \a spc.
    Space<T>& adjoin(const Space<T>& spc) {
        for(Nat i=0; i!=spc._variables.size(); ++i) { this->insert(VariableType(spc._variables[i])); } return *this; }
    //! \brief Append the named variable \a v to the variables defining the space.
    Space<T>& append(const VariableType& v) {
        for(Nat i=0; i!=_variables.size(); ++i) {
            ARIADNE_ASSERT_MSG(_variables[i]!=v.name(),"Variable "<<v<<" is already a variable of the StateSpace "<<*this);
        }
        _variables.push_back(v.name()); return *this; }
  private:
    List<Identifier> _variables;
};

template class Space<Real>;

template<class T> inline OutputStream& operator<<(OutputStream& os, const Space<T>& spc) { return os << spc.variables(); }

template<class T> inline Space<T> join(const Space<T>& spc1, const Space<T>& spc2) {
    Space<T> r(spc1); r.adjoin(spc2); return r; }
template<class T> inline Space<T> join(const Space<T>& spc1, const Variable<T>& var2) {
    Space<T> r(spc1); r.append(var2); return r; }

// Compiled conversion operators to allow conversion between expression and function.
SizeType dimension(const Space<Real>& spc);
Space<Real> space(const List< Variable<Real> >& vars);

template<class T> Variable<T> variable(const String& s) { return Variable<T>(s); }
template<class T> Space<T> variables(const List<String>& s) { return Space<T>(s); }



} // namespace Ariadne

#endif /* ARIADNE_SPACE_H */
