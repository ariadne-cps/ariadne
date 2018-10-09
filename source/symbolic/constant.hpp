/***************************************************************************
 *            constant.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file constant.hpp
 *  \brief Named constants
 */

#ifndef ARIADNE_CONSTANT_HPP
#define ARIADNE_CONSTANT_HPP

#include "../utility/string.hpp"
#include "../symbolic/identifier.hpp"

namespace Ariadne {

template<class T> class Constant;
typedef Constant<String> StringConstant;
typedef Constant<Integer> IntegerConstant;
typedef Constant<Real> RealConstant;

//! \ingroup ExpressionModule
//! A named constant of type \a T.
template<class T> class Constant
    : public T
{
  public:
    explicit Constant(const String& str, const T& value) : T(value), _name(str) { }
    const Identifier& name() const { return _name; }
    const T& value() const { return *this; }
    friend OutputStream& operator<<(OutputStream& os, Constant<T> const& c) {
        return os << c._name << "(=" << static_cast<T const&>(c) << ")"; }
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

} // namespace Ariadne

#endif /* ARIADNE_CONSTANT_HPP */
