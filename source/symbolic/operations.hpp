/***************************************************************************
 *            expression.hpp
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

/*! \file expression.hpp
 *  \brief Internal expressions
 */

#ifndef ARIADNE_EXPRESSION_OPERATIONS_HPP
#define ARIADNE_EXPRESSION_OPERATIONS_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

namespace Ariadne {

class String;

template<class T> class Set;

class Identifier;

template<class T> class Variable;
template<class T> class Space;
template<class T> class Expression;
template<class LHS,class RHS> class Assignment;

template<class X> class Vector;
template<class X> class Formula;
template<class X> class Algebra;


typedef Expression<Boolean> DiscretePredicate;
typedef Expression<Kleenean> ContinuousPredicate;
typedef Expression<String> StringExpression;
typedef Expression<Integer> IntegerExpression;
typedef Expression<Real> RealExpression;



template<class X> struct ExpressionNode;

//! \name Methods for building expressions
//! \related Expression
template<class T> struct DeclareExpressionOperations;

template<> struct DeclareExpressionOperations<Boolean>  {
    //! \related Expression \brief Logical disjunction.
    friend Expression<Boolean> operator&&(Expression<Boolean> const& e1, Expression<Boolean> const& e2);
    //! \related Expression \brief Logical conjunction.
    friend Expression<Boolean> operator||(Expression<Boolean> const& e1, Expression<Boolean> const& e2);
    //! \related Expression \brief Logical negation.
    friend Expression<Boolean> operator!(Expression<Boolean> const& e);
};

template<> struct DeclareExpressionOperations<Kleenean> {
    //! \related Expression \brief Fuzzy logical disjunction.
    friend Expression<Kleenean> operator&&(Expression<Kleenean> const& e1, Expression<Kleenean> const& e2);
    //! \related Expression \brief Fuzzy logical conjunction.
    friend Expression<Kleenean> operator||(Expression<Kleenean> const& e1, Expression<Kleenean> const& e2);
    //! \related Expression \brief Fuzzy logical negation.
    friend Expression<Kleenean> operator!(Expression<Kleenean> const& e);
};

template<> struct DeclareExpressionOperations<String>  {
    //! \related Expression \brief %String equality.
    friend Expression<Boolean> operator==(Variable<String> v1, const String& s2);
    //! \related Expression \brief %String inequality.
    friend Expression<Boolean> operator!=(Variable<String> v1, const String& s2);
};

template<> struct DeclareExpressionOperations<Integer>  {
    //! \related Expression \brief %Integer equality predicate.
    friend Expression<Boolean> operator==(Expression<Integer> const& e1, Expression<Integer> const& e2);
    //! \related Expression \brief %Integer inequality predicate.
    friend Expression<Boolean> operator!=(Expression<Integer> const& e1, Expression<Integer> const& e2);
    //! \related Expression \brief %Integer comparison predicate (greater or equal).
    friend Expression<Boolean> operator>=(Expression<Integer> const& e1, Expression<Integer> const& e2);
    //! \related Expression \brief %Integer comparison (less or equal)..
    friend Expression<Boolean> operator<=(Expression<Integer> const& e1, Expression<Integer> const& e2);
    //! \related Expression \brief %Integer comparison (greater).
    friend Expression<Boolean> operator> (Expression<Integer> const& e1, Expression<Integer> const& e2);
    //! \related Expression \brief %Integer comparison (less).
    friend Expression<Boolean> operator< (Expression<Integer> const& e1, Expression<Integer> const& e2);

    //! \related Expression \brief %Integer unary plus expression (identity).
    friend Expression<Integer> operator+(Expression<Integer> const& e);
    //! \related Expression \brief %Integer unary minus expression.
    friend Expression<Integer> operator-(Expression<Integer> const& e);
    //! \related Expression \brief %Integer addition expression.
    friend Expression<Integer> operator+(Expression<Integer> const& e1, Expression<Integer> const& e2);
    //! \related Expression \brief %Integer subtraction expression.
    friend Expression<Integer> operator-(Expression<Integer> const& e1, Expression<Integer> const& e2);
    //! \related Expression \brief %Integer multiplication expression.
    friend Expression<Integer> operator*(Expression<Integer> const& e1, Expression<Integer> const& e2);

    friend Expression<Integer>& operator+=(Expression<Integer>& e1, Expression<Integer> const& e2);
    friend Expression<Integer>& operator-=(Expression<Integer>& e1, Expression<Integer> const& e2);
    friend Expression<Integer>& operator*=(Expression<Integer>& e1, Expression<Integer> const& e2);
};

template<> struct DeclareExpressionOperations<Real> {
    //! \related Expression \brief Positivity test.
    //! Returns \c indeterminate if the value cannot be distinguished from zero.
    friend Expression<Kleenean> sgn(Expression<Real> const& e);
    //! \related Expression \brief Fuzzy inequality comparison predicate (less) of real expressions.
    friend Expression<Kleenean> operator<=(Expression<Real> const& e1, Expression<Real> const& e2);
    //! \related Expression \brief Fuzzy inequality comparison predicate (greater) of real expressions.
    friend Expression<Kleenean> operator>=(Expression<Real> const& e1, Expression<Real> const& e2);
    //! \related Expression \brief Fuzzy inequality comparison predicate (less) of real expressions.
    friend Expression<Kleenean> operator< (Expression<Real> const& e1, Expression<Real> const& e2);
    //! \related Expression \brief Fuzzy inequality comparison predicate (greater) of real expressions.
    friend Expression<Kleenean> operator> (Expression<Real> const& e1, Expression<Real> const& e2);

    //! \related Expression \brief Equality comparison predicate of real expressions.
    friend Expression<Kleenean> operator==(Expression<Real> const& e1, Expression<Real> const& e2);
    //! \related Expression \brief Negated equality comparison predicate of real expressions.
    friend Expression<Kleenean> operator!=(Expression<Real> const& e1, Expression<Real> const& e2);

    //! \related Expression \brief %Real unary plus expression.
    friend Expression<Real> operator+(Expression<Real> const& e);
    //! \related Expression \brief %Real unary minus expression.
    friend Expression<Real> operator-(Expression<Real> const& e);
    //! \related Expression \brief %Real addition expression.
    friend Expression<Real> operator+(Expression<Real> const& e1, Expression<Real> const& e2);
    //! \related Expression \brief %Real subtraction expression.
    friend Expression<Real> operator-(Expression<Real> const& e1, Expression<Real> const& e2);
    //! \related Expression \brief %Real multiplication expression.
    friend Expression<Real> operator*(Expression<Real> const& e1, Expression<Real> const& e2);
    //! \related Expression \brief %Real division expression.
    friend Expression<Real> operator/(Expression<Real> const& e1, Expression<Real> const& e2);

    friend Expression<Real>& operator+=(Expression<Real>& e1, Expression<Real> const& e2);
    friend Expression<Real>& operator-=(Expression<Real>& e1, Expression<Real> const& e2);
    friend Expression<Real>& operator*=(Expression<Real>& e1, Expression<Real> const& e2);
    friend Expression<Real>& operator/=(Expression<Real>& e1, Expression<Real> const& e2);

    //! \related Expression \brief %Real unary plus expression.
    //! Equivalent to +\a e.
    friend Expression<Real> pos(Expression<Real> const& e);
    //! \related Expression \brief %Real negation expression.
    //! Equivalent to -\a e.
    friend Expression<Real> neg(Expression<Real> const& e);
     //! \related Expression \brief %Real addition expression.
    friend Expression<Real> add(Expression<Real> const& e1, Expression<Real> const& e2);
    //! \related Expression \brief %Real subtraction expression.
    friend Expression<Real> sub(Expression<Real> const& e1, Expression<Real> const& e2);
    //! \related Expression \brief %Real multiplication expression.
    friend Expression<Real> mul(Expression<Real> const& e1, Expression<Real> const& e2);
    //! \related Expression \brief %Real division expression.
    friend Expression<Real> div(Expression<Real> const& e1, Expression<Real> const& e2);
    //! \related Expression \brief %Real integer power expression.
    friend Expression<Real> pow(Expression<Real> const& e, Int n);
   //! \related Expression \brief %Real reciprocal expression.
    //! Equivalent to 1/\a e.
    friend Expression<Real> rec(Expression<Real> const& e);
    //! \related Expression \brief %Real square expression.
    friend Expression<Real> sqr(Expression<Real> const& e);
    //! \related Expression \brief %Real square root expression.
    friend Expression<Real> sqrt(Expression<Real> const& e);
    //! \related Expression \brief %Real exponential expression.
    friend Expression<Real> exp(Expression<Real> const& e);
    //! \related Expression \brief %Real natural logarithm expression.
    friend Expression<Real> log(Expression<Real> const& e);
    //! \related Expression \brief %Real sine expression.
    friend Expression<Real> sin(Expression<Real> const& e);
    //! \related Expression \brief %Real cosine expression.
    friend Expression<Real> cos(Expression<Real> const& e);
    //! \related Expression \brief %Real tangent expression.
    friend Expression<Real> tan(Expression<Real> const& e);
    //! \related Expression \brief %Real arctangent expression.
    friend Expression<Real> atan(Expression<Real> const& e);

    //! \related Expression \brief Real maximum expression.
    friend Expression<Real> max(Expression<Real> const& e1, Expression<Real> const& e2);
    //! \related Expression \brief %Real minimum expression.
    friend Expression<Real> min(Expression<Real> const& e1, Expression<Real> const& e2);
    //! \related Expression \brief %Real absolute value expression.
    friend Expression<Real> abs(Expression<Real> const& e);
};


} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_OPERATIONS_HPP */
