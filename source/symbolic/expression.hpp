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

#ifndef ARIADNE_EXPRESSION_HPP
#define ARIADNE_EXPRESSION_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "../utility/macros.hpp"
#include "../utility/declarations.hpp"
#include "../utility/pointer.hpp"
#include "../utility/container.hpp"

#include "../numeric/logical.decl.hpp"
#include "../numeric/number.decl.hpp"

#include "../numeric/operators.hpp"
#include "../symbolic/constant.hpp"
#include "../symbolic/variables.hpp"
#include "../symbolic/valuation.hpp"
#include "../symbolic/operations.hpp"

namespace Ariadne {

template<class T> class Set;
template<class T> struct DeclareExpressionOperations;
template<class X> struct ExpressionNode;

//! \ingroup ExpressionModule
//! \brief A simple expression in named variables.
//!
//! %Ariadne supports expressions of type Boolean, Kleenean, String, Integer and Real.
//! The class Real is a dummy type which can be implemented in many different ways.
//!
//! The independent variables are given string names, rather than an integer index.
//! Formulae in different variables may be combined; the variables of the resulting formula
//! are all variables occuring in all formulae.
//! Formulae may be manipulated symbolically.
//!
//! \sa Variable \sa Assignment
template<class T>
class Expression
    : public DeclareExpressionOperations<T>
{
    typedef SharedPointer<const ExpressionNode<T>> Pointer;
  public:
    typedef Real NumericType;
    typedef T ValueType;
    typedef Constant<T> ConstantType;
    typedef Variable<T> VariableType;
  public:
    // Use template formulation to avoid ambiguity treating Expression(0) as a pointer construction.
    template<class P, EnableIf<IsConvertible<P,SharedPointer<const ExpressionNode<T>>>> =dummy>
        explicit Expression(P const& eptr) : _root(eptr) { }
  public:
    //! \brief Default expression is a constant with value \c 0.
    Expression();
    //! \brief Construct an expression from a numerical value.
    Expression(const T& c);
    template<class U, EnableIf<IsConvertible<U,T>> =dummy> Expression(const U& c) : Expression(T(c)) { }
    //! \brief Construct an expression from a named constant.
    Expression(const Constant<T>& c);
    //! \brief Construct an expression from a variable.
    Expression(const Variable<T>& v);
    //! \brief Construct a constant expression from a value.
    static Expression<T> constant(const ValueType& c);
    //! \brief Construct a constant expression from a named variable value.
    static Expression<T> constant(const Constant<ValueType>& c);
    //! \brief Construct an expression from a name.
    static Expression<T> variable(const Identifier& c);

    //! \brief Create the zero element.
    Expression<T> create_zero() const { return Expression<T>::constant(T()); }
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
    template<class X, EnableIf<And<IsSame<T,Real>,IsSame<X,EffectiveNumericType>>> =dummy> operator Algebra<X>() const;
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
//! \name Evaluation and related operations.
//! \related Expression

template<class T> LogicType<T> evaluate(const Expression<LogicType<T>>& e, const Valuation<T>& x);
template<class T> T evaluate(const Expression<T>& e, const Valuation<T>& x);

Boolean evaluate(const Expression<Boolean>& e, const DiscreteValuation& q);
template<class X> X evaluate(const Expression<Real>& e, const ContinuousValuation<X>&);
template<class X> Kleenean evaluate(const Expression<Kleenean>&, const ContinuousValuation<X>&);

//! \brief Evaluate expression \a e on argument \a x which is a map of variable identifiers to values of type \c A.
template<class T> LogicType<T> evaluate(const Expression<LogicType<T>>& e, const Map<Identifier,T>& x);
//! \brief Evaluate expression \a e on argument \a x which is a map of variable identifiers to values of type \c A.
template<class T> T evaluate(const Expression<T>& e, const Map<Identifier,T>& x);

//! \brief Extract the arguments of expression \a e.
template<class T> Set<Identifier> arguments(const Expression<T>& e);

//! \brief Substitute all occurrences of variable \a v of type \c Y with constant value \a c.
template<class T, class Y> Expression<T> substitute(const Expression<T>& e, const Variable<Y>& v, const Y& c);
//! \brief Substitute all occurrences of variable \a v of type \c Y with expression value \a se.
template<class T, class Y> Expression<T> substitute(const Expression<T>& e, const List< Assignment< Variable<Y>,Expression<Y> > >& a);
template<class T, class Y> Vector<Expression<T>> substitute(const Vector<Expression<T>>& e, const List< Assignment< Variable<Y>,Expression<Y> > >& a);

//! \brief Returns \a true if the expression\a e is syntactically equal to the constant \a c.
template<class T> Bool is_constant(const Expression<T>& e, const SelfType<T>& c);

//! \brief Returns \a true if the expression \a e is syntactically equal to the variable with name \a vn.
template<class T> Bool is_variable(const Expression<T>& e, const Identifier& vn);

//! \brief Returns \a true if the expression \a e is syntactically equal to the variable \a v.
template<class T> Bool is_variable(const Expression<T>& e, const Variable<T>& v);

//! \brief Returns \a true if the expression \a e is syntactically constant in the variables \a vs.
Bool is_constant_in(const Expression<Real>& e, const Set<Variable<Real>>& vs);
//! \brief Returns \a true if the expression \a e is syntactically affine in the variables \a vs.
Bool is_affine_in(const Expression<Real>& e, const Set<Variable<Real>>& vs);
//! \brief Returns \a true if the vector expression \a e is syntactically affine in the variables \a vs.
Bool is_affine_in(const Vector<Expression<Real>>& e, const Set<Variable<Real>>& vs);
//! \brief Returns \a true if the vector expression \a e is syntactically additive in the variables \a vs.
Bool is_additive_in(const Vector<Expression<Real>>& e, const Set<Variable<Real>>& vs);

//! \brief Simplify the expression \a e.
template<class T> Expression<T> simplify(const Expression<T>& e);

//! \brief Tests whether two expressions are identical.
template<class T> Bool identical(const Expression<T>& e1, const Expression<T>& e2);

//! \brief Returns true if the expressions are mutual negations.
//!
//! Currently can only test for pairs of the form (a1<=a2; a1>=a2),  (a1<=a2; a2<=a1)
//! or (a1>=a2; a2>=a1).
Bool opposite(Expression<Kleenean> p, Expression<Kleenean> q);

//! \brief Given \a sign when the predicate \a p is true.
Expression<Real> indicator(Expression<Kleenean> p, Sign sign=Sign::POSITIVE);

//! \brief Simplify the expression \a e.
Expression<Real> derivative(const Expression<Real>& e, Variable<Real> v);

//! \brief Make a formula in terms of numbered coordinates from an expression in named variables.
Formula<EffectiveNumber> make_formula(const Expression<Real>& e, const Map<Identifier,Nat>& v);
Formula<EffectiveNumber> make_formula(const Expression<Real>& e, const Space<Real>& spc);
Vector<Formula<EffectiveNumber>> make_formula(const Vector<Expression<Real>>& e, const Space<Real>& spc);
Formula<EffectiveNumber> make_formula(const Expression<Real>& e, const Variable<Real>& var);
Formula<EffectiveNumber> make_formula(const Expression<Real>& e, const List<Variable<Real>>& vars);
Formula<EffectiveNumber> make_formula(const Expression<Real>& out, const List<Assignment<Variable<Real>,Expression<Real>>>& aux, const Space<Real> spc);
Vector<Formula<EffectiveNumber>> make_formula(const Vector<Expression<Real>>& out, const List<Assignment<Variable<Real>,Expression<Real>>>& aux, const Space<Real> spc);

const Formula<EffectiveNumber>& cached_make_formula(const Expression<Real>& e, const Map<Identifier,Nat>& v, Map< const Void*, Formula<EffectiveNumber> >& cache);

//! \brief Make a function on the real line given an expression in a single argument variable.
ScalarUnivariateFunction<EffectiveTag> make_function(const Variable<Real>& v, const Expression<Real>& e);
//! \brief Make a function on a Euclidean domain given an ordered list including all argument variables.
ScalarFunction<EffectiveTag> make_function(const Space<Real>& s, const Expression<Real>& e);
//! \brief Make a function on a Euclidean domain given an ordered list including all argument variables.
VectorFunction<EffectiveTag> make_function(const Space<Real>& s, const Vector<Expression<Real>>& e);

//! \brief Make a function on a Euclidean domain given an ordered list including all argument variables. // DEPRECATED
ScalarFunction<EffectiveTag> make_function(const Expression<Real>& e, const Space<Real>& s);

//! \brief Make a function on a Euclidean domain given an ordered list including all argument variables.
Expression<Real> make_expression(const ScalarFunction<EffectiveTag>& f, const Space<Real>& s);

Expression<Real> make_expression(const Formula<Real>& f, const Space<Real>& s);

//@}




} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_HPP */
