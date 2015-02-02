/***************************************************************************
 *            expression.h
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


/*! \file expression.h
 *  \brief Internal expressions
 */

#ifndef ARIADNE_EXPRESSION_H
#define ARIADNE_EXPRESSION_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include <boost/pool/singleton_pool.hpp>

#include "utility/macros.h"
#include "utility/declarations.h"
#include "utility/pointer.h"
#include "utility/container.h"

#include "expression/valuation.h"
#include "expression/operators.h"
#include "expression/variables.h"

namespace Ariadne {

class Boolean;
class Tribool;

class String;
class Integer;
class Real;

template<class T> class Set;

class Identifier;

template<class T> class Variable;
template<class T> class Space;
template<class T> class Expression;
template<class LHS,class RHS> class Assignment;

template<class T, class V> class Valuation;
typedef Valuation<String,String> StringValuation;
typedef Valuation<Integer,Integer> IntegerValuation;
class DiscreteValuation;
template<class X> class ContinuousValuation;

template<class X> class ScalarFunction;
template<class X> class Vector;
template<class X> class Formula;


typedef Expression<Boolean> DiscretePredicate;
typedef Expression<Tribool> ContinuousPredicate;
typedef Expression<String> StringExpression;
typedef Expression<Integer> IntegerExpression;
typedef Expression<Real> RealExpression;



template<class X> struct ExpressionNode;

//! \ingroup ExpressionModule
//! \brief A simple expression in named variables.
//!
//! %Ariadne supports expressions of type Boolean, Tribool, String, Integer and Real.
//! The class Real is a dummy type which can be implemented in many different ways.
//!
//! The independent variables are given string names, rather than an integer index.
//! Formulae in different variables may be combined; the variables of the resulting formula
//! are all variables occuring in all formulae.
//! Formulae may be manipulated symbolically.
//!
//! \sa Variable \sa Assignment
template<class T>
class Expression {
    typedef SharedPointer<const ExpressionNode<T>> Pointer;
  public:
    typedef T ValueType;
    typedef Constant<T> ConstantType;
    typedef Variable<T> VariableType;
  public:
    explicit Expression(SharedPointer<const ExpressionNode<T>> eptr) : _root(eptr) { }
  public:
    //! \brief Default expression is a constant with value \c 0.
    Expression();
    //! \brief Construct an expression from a numerical value.
    Expression(const T& c);
    //! \brief Construct an expression from a named constant.
    Expression(const Constant<T>& c);
    //! \brief Construct an expression from a variable.
    Expression(const Variable<T>& v);
    //! \brief Construct a constant expression from a value.
    static Expression<T> constant(const ValueType& c);
    //! \brief Construct an expression from a name.
    static Expression<T> variable(const Identifier& c);
  public:
    const Operator& op() const;
    OperatorCode code() const;
    OperatorKind kind() const;
    const ValueType& val() const;
    const Identifier& var() const;
    const Expression<T>& arg() const;
    const Int& num() const;
    const Expression<T>& arg1() const;
    const Expression<T>& arg2() const;
    template<class A> const Expression<A>& cmp1(A* dummy=0) const;
    template<class A> const Expression<A>& cmp2(A* dummy=0) const;
    friend OutputStream& operator<<(OutputStream& os, Expression<T> const& e) { return e._write(os); }
  public:
    //! \brief The variables needed to compute the expression.
    Set<UntypedVariable> arguments() const;
  public:
    SharedPointer<const ExpressionNode<T>> node_ptr() const { return _root; }
    const ExpressionNode<T>* node_raw_ptr() const { return _root.operator->(); }
  private:
    OutputStream& _write(OutputStream& os) const;
  private:
    SharedPointer<const ExpressionNode<T>> _root;
};


//@{
//! \name List operations.
//! \related Expression
template<class T, class X> struct IsExpression;

template<class T, class X> struct IsExpression : public False { };
template<class T> struct IsExpression< T, T > : public True { };
template<class T> struct IsExpression< T, Constant<T> > : public True { };
template<class T> struct IsExpression< T, Variable<T> > : public True { };
template<class T> struct IsExpression< T, Expression<T> > : public True { };

template<class X1, class X2, EnableIf<And<IsExpression<Real,X1>,IsExpression<Real,X2>>> =dummy> inline
List<Expression<Real>> operator,(const X1& e1, const X2& e2) {
    List< Expression<Real> > r; r.append(e1); r.append(e2); return r; }

template<class X, EnableIf<IsExpression<Real,X>> =dummy> inline
List<Expression<Real>> operator,(Int c, const X& e) {
    List< Expression<Real> > r; r.append(Expression<Real>::constant(c)); r.append(Expression<Real>(e)); return r; }

template<class T, class X, EnableIf<IsExpression<T,X>> =dummy> inline
List<Expression<T>> operator,(List<Expression<T>> l, const X& e) {
    List< Expression<T> > r(l); r.append(Expression<T>(e)); return r; }
//@}

//@{
//! \name Evaluation and related operations.
//! \related Expression

Boolean evaluate(const Expression<Boolean>& e, const DiscreteValuation& q);
String evaluate(const Expression<String>& e, const StringValuation& q);
Integer evaluate(const Expression<Integer>& e, const IntegerValuation& q);
Real evaluate(const Expression<Integer>& e, const ContinuousValuation<Real>& q);
Tribool evaluate(const Expression<Tribool>& e, const ContinuousValuation<Real>& q);

//! \brief Evaluate expression \a e on argument \a x which is a map of variable identifiers to values of type \c A.
template<class A> typename Logic<A>::Type evaluate(const Expression<typename Logic<A>::Type>& e, const Map<Identifier,A>& x);

//! \brief Evaluate expression \a e on argument \a x which is a map of variable identifiers to values of type \c A.
template<class T> T evaluate(const Expression<T>& e, const Map<Identifier,T>& x);

//! \brief Extract the arguments of expression \a e.
template<class T> Set<Identifier> arguments(const Expression<T>& e);

//! \brief Returns \a true if the expression\a e is syntactically equal to the constant \a c.
template<class T> Bool is_constant(const Expression<T>& e, const typename Expression<T>::ValueType& c);

//! \brief Returns \a true if the expression \a e is syntactically equal to the variable with name \a vn.
template<class T> Bool is_variable(const Expression<T>& e, const Identifier& vn);

//! \brief Returns \a true if the expression \a e is syntactically equal to the variable \a v.
template<class T> Bool is_variable(const Expression<T>& e, const Variable<T>& v);

//! \brief Simplify the expression \a e.
template<class T> Expression<T> simplify(const Expression<T>& e);

//! \brief Tests whether two expressions are identical.
template<class T> Bool identical(const Expression<T>& e1, const Expression<T>& e2);

//! \brief Returns true if the expressions are mutual negations.
//!
//! Currently can only test for pairs of the form (a1<=a2; a1>=a2),  (a1<=a2; a2<=a1)
//! or (a1>=a2; a2>=a1).
Bool opposite(Expression<Tribool> p, Expression<Tribool> q);

//! \brief Given \a sign when the predicate \a p is true.
Expression<Real> indicator(Expression<Tribool> p, Sign sign=POSITIVE);

//! \brief Substitute all occurrences of variable \a v of type \c Y with constant value \a c.
template<class T, class Y> Expression<T> substitute(const Expression<T>& e, const Variable<Y>& v, const Y& c);

//! \brief Substitute all occurrences of variable \a v of type \c Y with expression value \a se.
template<class T, class Y> Expression<T> substitute(const Expression<T>& e, const Variable<Y>& v, const Expression<Y>& se);

template<class T, class Y> Expression<T> substitute(const Expression<T>& e, const List< Assignment< Variable<Y>,Expression<Y> > >& a);


//! \brief Make a function on a Euclidean domain given an ordered list including all argument variables.
ScalarFunction<EffectiveTag> make_function(const Expression<Real>& e, const Space<Real>& s);
//! \brief Make a function on coordinates given a mapping from variable names to indices.
Formula<Real> formula(const Expression<Real>& e, const Map<Identifier,Nat>& v);
Formula<Real> formula(const Expression<Real>& e, const List<Variable<Real>>& vars);
Formula<Real> formula(const Expression<Real>& res, const List<Assignment<Variable<Real>,Expression<Real>>>& aux, const Space<Real> spc);
List< Formula<Real> > formula(const List<Expression<Real>>& res, const List<Assignment<Variable<Real>,Expression<Real>>>& aux, const Space<Real> spc);

//@}


//@{
//! \name Methods for building expressions
//! \related Expression

//! \related Expression \brief Logical disjunction.
Expression<Boolean> operator&&(Expression<Boolean> e1, Expression<Boolean> e2);
//! \related Expression \brief Logical conjunction.
Expression<Boolean> operator||(Expression<Boolean> e1, Expression<Boolean> e2);
//! \related Expression \brief Logical negation.
Expression<Boolean> operator!(Expression<Boolean> e);

//! \related Expression \brief Fuzzy logical disjunction.
Expression<Tribool> operator&&(Expression<Tribool> e1, Expression<Tribool> e2);
//! \related Expression \brief Fuzzy logical conjunction.
Expression<Tribool> operator||(Expression<Tribool> e1, Expression<Tribool> e2);
//! \related Expression \brief Fuzzy logical negation.
Expression<Tribool> operator!(Expression<Tribool> e);

//! \related Expression \brief %String equality.
Expression<Boolean> operator==(Variable<String> v1, const String& s2);
//! \related Expression \brief %String inequality.
Expression<Boolean> operator!=(Variable<String> v1, const String& s2);


//! \related Expression \brief %Integer equality predicate.
Expression<Boolean> operator==(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief %Integer inequality predicate.
Expression<Boolean> operator!=(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief %Integer comparison predicate (greater or equal).
Expression<Boolean> operator>=(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief %Integer comparison (less or equal)..
Expression<Boolean> operator<=(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief %Integer comparison (greater).
Expression<Boolean> operator> (Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief %Integer comparison (less).
Expression<Boolean> operator< (Expression<Integer> e1, Expression<Integer> e2);

//! \related Expression \brief %Integer unary plus expression (identity).
Expression<Integer> operator+(Expression<Integer> e);
//! \related Expression \brief %Integer unary minus expression.
Expression<Integer> operator-(Expression<Integer> e);
//! \related Expression \brief %Integer addition expression.
Expression<Integer> operator+(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief %Integer subtraction expression.
Expression<Integer> operator-(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief %Integer multiplication expression.
Expression<Integer> operator*(Expression<Integer> e1, Expression<Integer> e2);



//! \related Expression \brief Positivity test.
//! Returns \c indeterminate if the value cannot be distinguished from zero.
Expression<Tribool> sgn(Expression<Real> e);
//! \related Expression \brief Fuzzy inequality comparison predicate (less) of real expressions.
Expression<Tribool> operator<=(Expression<Real> e1, Expression<Real> e2);
//! \related Expression \brief Fuzzy inequality comparison predicate (greater) of real expressions.
Expression<Tribool> operator>=(Expression<Real> e1, Expression<Real> e2);
//! \related Expression \brief Fuzzy inequality comparison predicate (less) of real expressions.
Expression<Tribool> operator< (Expression<Real> e1, Expression<Real> e2);
//! \related Expression \brief Fuzzy inequality comparison predicate (greater) of real expressions.
Expression<Tribool> operator> (Expression<Real> e1, Expression<Real> e2);

//! \related Expression \brief %Real unary plus expression.
Expression<Real> operator+(Expression<Real> e);
//! \related Expression \brief %Real unary minus expression.
Expression<Real> operator-(Expression<Real> e);
//! \related Expression \brief %Real addition expression.
Expression<Real> operator+(Expression<Real> e1, Expression<Real> e2);
//! \related Expression \brief %Real subtraction expression.
Expression<Real> operator-(Expression<Real> e1, Expression<Real> e2);
//! \related Expression \brief %Real multiplication expression.
Expression<Real> operator*(Expression<Real> e1, Expression<Real> e2);
//! \related Expression \brief %Real division expression.
Expression<Real> operator/(Expression<Real> e1, Expression<Real> e2);

//! \related Expression \brief %Real integer power expression.
Expression<Real> pow(Expression<Real> e, Int n);

//! \related Expression \brief %Real negation expression.
//! Equivalent to -\a e.
Expression<Real> neg(Expression<Real> e);
//! \related Expression \brief %Real reciprocal expression.
//! Equivalent to 1/\a e.
Expression<Real> rec(Expression<Real> e);
//! \related Expression \brief %Real square expression.
Expression<Real> sqr(Expression<Real> e);
//! \related Expression \brief %Real square root expression.
Expression<Real> sqrt(Expression<Real> e);
//! \related Expression \brief %Real exponential expression.
Expression<Real> exp(Expression<Real> e);
//! \related Expression \brief %Real natural logarithm expression.
Expression<Real> log(Expression<Real> e);
//! \related Expression \brief %Real sine expression.
Expression<Real> sin(Expression<Real> e);
//! \related Expression \brief %Real cosine expression.
Expression<Real> cos(Expression<Real> e);
//! \related Expression \brief %Real tangent expression.
Expression<Real> tan(Expression<Real> e);

//! \related Expression \brief Real maximum expression.
Expression<Real> max(Expression<Real> e1, Expression<Real> e2);
//! \related Expression \brief %Real minimum expression.
Expression<Real> min(Expression<Real> e1, Expression<Real> e2);
//! \related Expression \brief %Real absolute value expression.
Expression<Real> abs(Expression<Real> e);

Expression<Real> operator+(Expression<Real> e, Real c);
Expression<Real> operator-(Expression<Real> e, Real c);
Expression<Real> operator*(Expression<Real> e, Real c);
Expression<Real> operator/(Expression<Real> e, Real c);
Expression<Real> operator+(Real c, Expression<Real> e);
Expression<Real> operator-(Real c, Expression<Real> e);
Expression<Real> operator*(Real c, Expression<Real> e);
Expression<Real> operator/(Real c, Expression<Real> e);

Expression<Tribool> operator<=(Expression<Real> e, Real c);
Expression<Tribool> operator< (Expression<Real> e, Real c);
Expression<Tribool> operator>=(Expression<Real> e, Real c);
Expression<Tribool> operator> (Expression<Real> e, Real c);
Expression<Tribool> operator<=(Real c, Expression<Real> e);
Expression<Tribool> operator< (Real c, Expression<Real> e);
Expression<Tribool> operator>=(Real c, Expression<Real> e);
Expression<Tribool> operator> (Real c, Expression<Real> e);


//@}

} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_H */
