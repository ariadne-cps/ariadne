/***************************************************************************
 *            space.hpp
 *
 *  Copyright 2008-17 Pieter Collins
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

/*! \file space.hpp
 *  \brief Spaces formed by variables.
 */

#ifndef ARIADNE_SPACE_HPP
#define ARIADNE_SPACE_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "../utility/macros.hpp"
#include "../utility/pointer.hpp"
#include "../utility/container.hpp"

#include "../symbolic/variables.hpp"

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
    Space();
    Space(const List<VariableType>& vl);
    Space(const List<Identifier>& vl);
    Space(const InitializerList<VariableType>& vl);

    Bool operator==(const Space<T>& other) const;
    Bool operator!=(const Space<T>& other) const;

    SizeType size() const;
    //! \brief The dimension of the space.
    SizeType dimension() const;
    //! \brief The \a i<sup>th</sup> named variable.
    const VariableType operator[](SizeType i) const;
    const VariableType variable(SizeType i) const;

    //! \brief A list giving ordered variables.
    List<Identifier> variable_names() const;
    //! \brief A list giving ordered variables.
    List<VariableType> variables() const;
    //! \brief A map giving the index of a given variable.
    Map<Identifier,SizeType> indices_from_names() const;
    //! \brief A map giving the index of a given variable.
    Map<VariableType,SizeType> indices() const;

    //! \brief Tests if the variable \a v is in the space.
    Bool contains(const VariableType& v) const;
    //! \brief Tests if all the variables \a vs is in the space.
    Bool contains(const Set<VariableType>& vs) const;
    //! \brief The index of the named variable \a v.
    SizeType index(const VariableType& v) const;
    SizeType index(const Identifier& n) const;
    //! \brief Append the named variable \a v to the variables defining the space; ignores if the variable is already in the space.
    Space<T>& insert(const VariableType& v);
    //! \brief Adjoins the variables in \a spc.
    Space<T>& adjoin(const Space<T>& spc);
    //! \brief Append the named variable \a v to the variables defining the space.
    Space<T>& append(const VariableType& v);
  private:
    List<Identifier> _variables;
};

template<class T> OutputStream& operator<<(OutputStream& os, const Space<T>& spc) { return os << spc.variables(); }

template<class T> Space<T> join(const Space<T>& spc1, const Space<T>& spc2) {
    Space<T> r(spc1); r.adjoin(spc2); return r; }
template<class T> Space<T> join(const Space<T>& spc1, const Variable<T>& var2) {
    Space<T> r(spc1); r.append(var2); return r; }

// Compiled conversion operators to allow conversion between expression and function.
SizeType dimension(const Space<Real>& spc);
List<Identifier> variable_names(const Space<Real>& spc);
List<Identifier> variable_names(const List<Variable<Real>>& spc);
Space<Real> real_space(const List<Identifier>& vars);

template<class T> Variable<T> variable(const Identifier& s);
template<class T> Space<T> variables(const List<Identifier>& s);


} // namespace Ariadne

#endif /* ARIADNE_SPACE_HPP */
