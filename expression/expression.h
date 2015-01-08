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


template<class X> struct Logic { };
template<> struct Logic<String> { typedef Boolean Type; };
template<> struct Logic<Integer> { typedef Boolean Type; };
template<> struct Logic<Real> { typedef Tribool Type; };


//@{
//! \name Concrete expression classes.
//! \related Expression

template<class X> struct ExpressionNode;

template<class D, class X, class T=Void> struct EnableIfDoubleReal { };
template<class T> struct EnableIfDoubleReal<double,Real,T> { typedef T Type; };


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
    typedef Identifier V;
  public:
    typedef T ConstantType;
  public:
    explicit Expression(const ExpressionNode<T>* eptr) : _root(eptr) { }
    explicit Expression(counted_pointer< const ExpressionNode<T> > eptr) : _root(eptr) { }
  public:
    //! \brief Construct the constant expression with the default value of \a T.
    Expression();
    //! \brief Construct a constant expression from a value.
    Expression(const T& c);
    //! \brief Construct an expression from a name.
    Expression(const Identifier& v);
    //! \brief Construct an expression from a constant.
    Expression(const Constant<T>& c);
    //! \brief Construct an expression from a variable.
    Expression(const Variable<T>& v);
    //! \brief Construct a constant expression from a value.
    static Expression<T> constant(const T& v);

  public:
    const Operator& op() const;
    OperatorCode code() const;
    OperatorKind kind() const;
    const T& val() const;
    const V& var() const;
    const Expression<T>& arg() const;
    const Int& num() const;
    const Expression<T>& arg1() const;
    const Expression<T>& arg2() const;
    template<class A> const Expression<A>& cmp1(A* dummy=0) const;
    template<class A> const Expression<A>& cmp2(A* dummy=0) const;
  public:
    //! \brief The variables needed to compute the expression.
    Set<UntypedVariable> arguments() const;
  public:
    const ExpressionNode<T>* node_ptr() const { return _root.operator->(); }
  private:
    counted_pointer< const ExpressionNode<T> > _root;
};


template<class T>
struct ExpressionNode {
    mutable uint count;
    Operator op;
    virtual ~ExpressionNode();
    explicit ExpressionNode(const Operator& o) : count(0u), op(o) { }
    explicit ExpressionNode(OperatorCode cd, OperatorKind knd) : count(0u), op(cd,knd) { }
};

template<class T> struct ConstantExpressionNode : public ExpressionNode<T> {
    T val;
    ConstantExpressionNode(const T& v) : ExpressionNode<T>(CNST,NULLARY), val(v) { }
};
template<class T> struct VariableExpressionNode : public ExpressionNode<T> {
    Identifier var;
    VariableExpressionNode(const Identifier& v) : ExpressionNode<T>(VAR,VARIABLE), var(v) { }
};
template<class T, class A=T> struct UnaryExpressionNode : public ExpressionNode<T> {
    Expression<A> arg;
    UnaryExpressionNode(const Operator& op, Expression<A> const& a) : ExpressionNode<T>(op), arg(a) { }
};
template<class T, class A1=T, class A2=A1> struct BinaryExpressionNode {
    Expression<T> arg1; Expression<T> arg2;
};
template<class T> struct BinaryExpressionNode<T> : public ExpressionNode<T> {
    Expression<T> arg1; Expression<T> arg2;
    BinaryExpressionNode(const Operator& op, Expression<T> const& a1, Expression<T> const& a2)
        : ExpressionNode<T>(op), arg1(a1), arg2(a2) { }
};
template<class T> struct BinaryExpressionNode<typename Logic<T>::Type,T,T> : public ExpressionNode<typename Logic<T>::Type> {
    typedef typename Logic<T>::Type R; typedef T A;
    Expression<A> arg1; Expression<A> arg2;
    BinaryExpressionNode(const Operator& op, Expression<A> const& a1, Expression<A> const& a2)
        : ExpressionNode<R>(op), arg1(a1), arg2(a2) { }
};
template<class R, class A=R, class N=Int> struct ScalarExpressionNode : public UnaryExpressionNode<R,A> {
    N num;
    ScalarExpressionNode(const Operator& op, Expression<R> const& a, N n)
        : UnaryExpressionNode<R,A>(op,a), num(n) { }
};

template<class T> ExpressionNode<T>::~ExpressionNode() { }

template<class T> const Operator& Expression<T>::op() const {
    return node_ptr()->op; }
template<class T> OperatorCode Expression<T>::code() const {
    return node_ptr()->op.code(); }
template<class T> OperatorKind Expression<T>::kind() const {
    return node_ptr()->op.kind(); }
template<class T> const T& Expression<T>::val() const {
    return static_cast<const ConstantExpressionNode<T>*>(node_ptr())->val; }
template<class T> const Identifier& Expression<T>::var() const {
    return static_cast<const VariableExpressionNode<T>*>(node_ptr())->var; }
template<class T> const Expression<T>& Expression<T>::arg() const {
    return static_cast<const UnaryExpressionNode<T>*>(node_ptr())->arg; }
template<class T> const Int& Expression<T>::num() const {
    return static_cast<const ScalarExpressionNode<T>*>(node_ptr())->num; }
template<class T> const Expression<T>& Expression<T>::arg1() const {
    return static_cast<const BinaryExpressionNode<T>*>(node_ptr())->arg1; }
template<class T> const Expression<T>& Expression<T>::arg2() const {
    return static_cast<const BinaryExpressionNode<T>*>(node_ptr())->arg2; }
template<class R> template<class A> const Expression<A>& Expression<R>::cmp1(A*) const {
    return static_cast<const BinaryExpressionNode<R,A>*>(node_ptr())->arg1; }
template<class R> template<class A> const Expression<A>& Expression<R>::cmp2(A*) const {
    return static_cast<const BinaryExpressionNode<R,A>*>(node_ptr())->arg2; }


template<class T> inline Expression<T>::Expression() : _root(new ConstantExpressionNode<T>(T())) { }
template<class T> inline Expression<T>::Expression(const T& c) : _root(new ConstantExpressionNode<T>(c)) { }
template<class T> inline Expression<T>::Expression(const Identifier& v) : _root(new VariableExpressionNode<T>(v)) { }
template<class T> inline Expression<T>::Expression(const Constant<T>& c): _root(new ConstantExpressionNode<T>(c.value())) { };
template<class T> inline Expression<T>::Expression(const Variable<T>& v) : _root(new VariableExpressionNode<T>(v.name())) { }
template<class T> inline Expression<T> Expression<T>::constant(const T& c) { return Expression<T>(c); }

template<class R> inline
Expression<R> make_expression(const R& c) {
    return Expression<R>(new ConstantExpressionNode<R>(c)); }
template<class R, class A> inline
Expression<R> make_expression(OperatorCode op, const Expression<A>& e) {
    return Expression<R>(new UnaryExpressionNode<R,A>(op,e)); }
template<class R, class A, class N> inline
Expression<R> make_expression(OperatorCode op, const Expression<A>& e, N n) {
    return Expression<R>(new ScalarExpressionNode<R,A,N>(op,e,n)); }
template<class R, class A1, class A2> inline
Expression<R> make_expression(OperatorCode op, const Expression<A1>& e1, Expression<A2> e2) {
    return Expression<R>(new BinaryExpressionNode<R,A1,A2>(op,e1,e2)); }

template<class R, class Op, class A> inline
Expression<R> make_expression(Op op, const Expression<A> e) {
    return make_expression<R,A>(op.code(),e); }
template<class R, class Op, class A, class N> inline
Expression<R> make_expression(Op op, const Expression<A>& e, N n) {
    return make_expression<R>(op.code(),e,n); }
template<class R, class Op, class A1, class A2> inline
Expression<R> make_expression(Op op, const Expression<A1>& e1, Expression<A2> e2) {
    return make_expression<R>(op.code(),e1,e2); }



template<class T> inline OutputStream& operator<<(OutputStream& os, const ExpressionNode<T>* e) {
    return os << (void*)(e);
}


//@}



//@{
//! \name Input / output operations.
//! \related Expression

template<class T> OutputStream& operator<<(OutputStream& os, const Expression<T>& f);
template<class T> inline OutputStream& _write_comparison(OutputStream& os, const Expression<T>& f) {
    ARIADNE_FAIL_MSG("Comparison must return a logical type."); }
template<> inline OutputStream& _write_comparison(OutputStream& os, const Expression<Tribool>& f) {
    Real* real_ptr=0; return os << "(" << f.cmp1(real_ptr) << symbol(f.op()) << f.cmp2(real_ptr) << ")"; }
template<> inline OutputStream& _write_comparison(OutputStream& os, const Expression<Boolean>& f) {
    String* string_ptr=0; return os << "(" << f.cmp1(string_ptr) << symbol(f.op()) << f.cmp2(string_ptr) << ")"; }
//FIXME: Distinguish String and Integer comparisons

//! \brief Write to an output stream
template<class T> OutputStream& operator<<(OutputStream& os, const Expression<T>& f) {
    switch(f.op()) {
        //case CNST: return os << std::fixed << std::setprecision(4) << fptr->val;
        case CNST:
            os << f.val(); return os;
            //if(f.val()==0.0) { return os << 0.0; } if(abs(f.val())<1e-4) { os << std::fixed << f.val(); } else { os << f.val(); } return os;
        case VAR:
            return os << f.var();
        case ADD:
            return os << f.arg1() << '+' << f.arg2();
        case SUB:
            os << f.arg1() << '-';
            switch(f.arg2().op()) { case ADD: case SUB: os << '(' << f.arg2() << ')'; break; default: os << f.arg2(); }
            return os;
        case MUL:
            switch(f.arg1().op()) { case ADD: case SUB: case DIV: os << '(' << f.arg1() << ')'; break; default: os << f.arg1(); }
            os << '*';
            switch(f.arg2().op()) { case ADD: case SUB: os << '(' << f.arg2() << ')'; break; default: os << f.arg2(); }
            return os;
        case DIV:
            switch(f.arg1().op()) { case ADD: case SUB: case DIV: os << '(' << f.arg1() << ')'; break; default: os << f.arg1(); }
            os << '/';
            switch(f.arg2().op()) { case ADD: case SUB: case MUL: case DIV: os << '(' << f.arg2() << ')'; break; default: os << f.arg2(); }
            return os;
        case POW:
            return os << "pow" << '(' << f.arg() << ',' << f.num() << ')';
        default:
            switch(f.kind()) {
                case UNARY: return os << f.op() << "(" << f.arg() << ")";
                case BINARY: return os << f.op() << "(" << f.arg1() << "," << f.arg2() << ")";
                // FIXME: Type-cast comparison arguments correctly
                case COMPARISON: return _write_comparison(os,f);
                default: ARIADNE_FAIL_MSG("Cannot output expression with operator "<<f.op()<<" of kind "<<f.kind()<<"\n");
            }
    }
}

//@}

//@{
//! \name Evaluation and related operations.
//! \related Expression


//! \brief Evaluate expression \a e on argument \a x which is a map of variable identifiers to values of type \c A.
template<class A>
typename Logic<A>::Type
evaluate(const Expression<typename Logic<A>::Type>& e, const Map<Identifier,A>& x)
{
    typedef typename Logic<A>::Type R;
    A* aptr=0;
    switch(e.kind()) {
        case NULLARY: return static_cast<R>(e.val());
        case UNARY: return compute(e.op(),evaluate(e.arg(),x));
        case BINARY: return compute(e.op(),evaluate(e.arg1(),x),evaluate(e.arg2(),x));
        case COMPARISON: return compare<R>(e.op(),evaluate(e.cmp1(aptr),x),evaluate(e.cmp2(aptr),x));
        default: ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" on "<<x<<"\n");
    }
}

//! \brief Evaluate expression \a e on argument \a x which is a map of variable identifiers to values of type \c A.
template<class T>
T
evaluate(const Expression<T>& e, const Map<Identifier,T>& x)
{
    switch(e.kind()) {
        case VARIABLE: return x[e.var()];
        case NULLARY: return e.val();
        case UNARY: return compute(e.op(),evaluate(e.arg(),x));
        case BINARY: return compute(e.op(),evaluate(e.arg1(),x),evaluate(e.arg2(),x));
        case SCALAR: return compute(e.op(),evaluate(e.arg(),x),e.num());
        default: ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" on "<<x<<"\n");
    }
}



template<class T> Set<UntypedVariable> Expression<T>::arguments() const {
    const Expression<T>& e=*this;
    switch(e.kind()) {
        case VARIABLE: return Set<UntypedVariable>{Variable<T>(e.var())};
        case NULLARY: return Set<UntypedVariable>();
        case UNARY: return e.arg().arguments();
        case BINARY: return join(e.arg1().arguments(),e.arg2().arguments());
        case COMPARISON: {
            const BinaryExpressionNode<T,Real>* rlp = dynamic_cast<const BinaryExpressionNode<T,Real>*>(e.node_ptr());
            if(rlp) { return join(rlp->arg1.arguments(),rlp->arg2.arguments()); }
            const BinaryExpressionNode<T,String>* strp = dynamic_cast<const BinaryExpressionNode<T,String>*>(e.node_ptr());
            if(strp) { return join(strp->arg1.arguments(),strp->arg2.arguments()); }
        }
        default: ARIADNE_FAIL_MSG("Cannot compute arguments of expression "<<e<<" of kind "<<e.kind()<<"\n");
    }
}

//! \brief Extract the arguments of expression \a e.
template<class T> Set<Identifier> arguments(const Expression<T>& e)
{
    switch(e.kind()) {
        case VARIABLE: return Set<Identifier>{e.var()};
        case NULLARY: return Set<Identifier>();
        case UNARY: return arguments(e.arg());
        case BINARY: return join(arguments(e.arg1()),arguments(e.arg2()));
        case COMPARISON: {
            const BinaryExpressionNode<T,Real>* rlp = dynamic_cast<const BinaryExpressionNode<T,Real>*>(e.node_ptr());
            if(rlp) { return join(arguments(rlp->arg1),arguments(rlp->arg2)); }
            const BinaryExpressionNode<T,String>* strp = dynamic_cast<const BinaryExpressionNode<T,String>*>(e.node_ptr());
            if(strp) { return join(arguments(strp->arg1),arguments(strp->arg2)); }
        }
        default: ARIADNE_FAIL_MSG("Cannot compute arguments of expression "<<e<<" of kind "<<e.kind()<<"\n");
    }
}


//! \brief Returns \a true if the expression\a e is syntactically equal to the constant \a c.
template<class T> Bool is_constant(const Expression<T>& e, const T& c) {
    switch(e.op()) {
        case CNST: return decide(e.val()==c);
        default: return false;
    }
}

template<class T, class X> Bool is_constant(const Expression<T>& e, const X& c) {
    return is_constant(e,static_cast<T>(c));
}

//! \brief Returns \a true if the expression \a e is syntactically equal to the variable with name \a vn.
template<class T> Bool is_variable(const Expression<T>& e, const Identifier& vn) {
    switch(e.op()) {
        case VAR: return e.var()==vn;
        default: return false;
    }
}

//! \brief Returns \a true if the expression \a e is syntactically equal to the variable \a v.
template<class T> Bool is_variable(const Expression<T>& e, const Variable<T>& v) {
    switch(e.op()) {
        case VAR: return e.var()==v.name();
        default: return false;
    }
}

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

//! \brief Convert the expression in named variables to a formula in numbered coordinates.
Formula<Real> formula(const Expression<Real>& e, const Map<Identifier,Nat>& v);


//@}


//@{
//! \name Deprecated conversions and Expression/Formula operators.
//! \related Expression

Boolean evaluate(const Expression<Boolean>& e, const DiscreteValuation& q);
String evaluate(const Expression<String>& e, const StringValuation& q);
Integer evaluate(const Expression<Integer>& e, const IntegerValuation& q);

template<class T> Tribool evaluate(const Expression<Tribool>& e, const ContinuousValuation<T>& x);
template<class T> T evaluate(const Expression<Real>& e, const ContinuousValuation<T>& x);
//template<class T> T evaluate(const Expression<Real>& e, const Map<ExtendedVariable<Real>,T>& x);

Formula<Real> formula(const Expression<Real>& e, const List< Variable<Real> >& vars);
Formula<Real> formula(const Expression<Real>& res, const List< Assignment< Variable<Real>, Expression<Real> > >& aux, const Space<Real> spc);
List< Formula<Real> > formula(const List< Expression<Real> >& res, const List< Assignment< Variable<Real>, Expression<Real> > >& aux, const Space<Real> spc);


ScalarFunction<EffectiveTag> make_function(const Expression<Real>& e, const Space<Real>& s);


//@}

//@{
//! \name Metaprogramming and sequencing operators for making lists
//! \related Expression

//! \brief Inherits from True type if class \a X is a constant, variable or expression in type \a T, otherwise inherits from False.
template<class T, class X> struct IsExpression;

template<class T, class X> struct IsExpression : public False { };
template<class T> struct IsExpression< T, T > : public True { };
template<class T> struct IsExpression< T, Constant<T> > : public True { };
template<class T> struct IsExpression< T, Variable<T> > : public True { };
template<class T> struct IsExpression< T, Expression<T> > : public True { };
template<class T, class X1, class X2, class R=Dummy> using EnableIfExpressions = EnableIf< And<IsExpression<T,X1>,IsExpression<T,X2> >, R>;

template<class X1, class X2> inline
EnableIfExpressions< Real, X1, X2, List<Expression<Real> > >
operator,(const X1& e1, const X2& e2) {
    List< Expression<Real> > r; r.append(e1); r.append(e2); return r; }

template<class X> inline
EnableIf< IsExpression<Real,X>, List<Expression<Real> > >
operator,(Int c, const X& e) {
    List< Expression<Real> > r; r.append(Expression<Real>(c)); r.append(Expression<Real>(e)); return r; }

template<class T, class X> inline
EnableIf< IsExpression<T,X>, List<Expression<T> > >
operator,(List<Expression<T> > l, const T& e) {
    List< Expression<T> > r(l); r.append(Expression<Real>(e)); return r; }

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

//@}

inline Expression<Real> operator+(Expression<Real> e, Real c) { return e + Expression<Real>::constant(c); }
inline Expression<Real> operator-(Expression<Real> e, Real c) { return e - Expression<Real>::constant(c); }
inline Expression<Real> operator*(Expression<Real> e, Real c) { return e * Expression<Real>::constant(c); }
inline Expression<Real> operator/(Expression<Real> e, Real c) { return e / Expression<Real>::constant(c); }
inline Expression<Real> operator+(Real c, Expression<Real> e) { return Expression<Real>::constant(c) + e; }
inline Expression<Real> operator-(Real c, Expression<Real> e) { return Expression<Real>::constant(c) - e; }
inline Expression<Real> operator*(Real c, Expression<Real> e) { return Expression<Real>::constant(c) * e; }
inline Expression<Real> operator/(Real c, Expression<Real> e) { return Expression<Real>::constant(c) / e; }

inline Expression<Tribool> operator<=(Expression<Real> e, Real c) { return e <= Expression<Real>::constant(c); }
inline Expression<Tribool> operator< (Expression<Real> e, Real c) { return e <  Expression<Real>::constant(c); }
inline Expression<Tribool> operator>=(Expression<Real> e, Real c) { return e >= Expression<Real>::constant(c); }
inline Expression<Tribool> operator> (Expression<Real> e, Real c) { return e >  Expression<Real>::constant(c); }
inline Expression<Tribool> operator<=(Real c, Expression<Real> e) { return Expression<Real>::constant(c) <= e; }
inline Expression<Tribool> operator< (Real c, Expression<Real> e) { return Expression<Real>::constant(c) <  e; }
inline Expression<Tribool> operator>=(Real c, Expression<Real> e) { return Expression<Real>::constant(c) >= e; }
inline Expression<Tribool> operator> (Real c, Expression<Real> e) { return Expression<Real>::constant(c) >  e; }


} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_H */
