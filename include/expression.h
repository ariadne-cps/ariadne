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

#include "macros.h"
#include "declarations.h"
#include "pointer.h"
#include "container.h"

#include "operators.h"
#include "variables.h"

namespace Ariadne {

//! Internal name for boolean objects in expressions.
typedef bool Boolean;
typedef tribool Tribool;
typedef std::string String;
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

//! The sign of a numerical value.
enum Sign { NEGATIVE=-1, ZERO=0, POSITIVE=+1 };

typedef Expression<Boolean> DiscretePredicate;
typedef Expression<Tribool> ContinuousPredicate;
typedef Expression<String> StringExpression;
typedef Expression<Integer> IntegerExpression;
typedef Expression<Real> RealExpression;


template<class X> struct Logic { };
template<> struct Logic<String> { typedef Bool Type; };
template<> struct Logic<Integer> { typedef Bool Type; };
template<> struct Logic<Real> { typedef Tribool Type; };


//@{
//! \name Concrete expression classes.
//! \related Expression

template<class X> struct ExpressionNode;

struct Index {
    Nat _i;
  public:
    explicit Index(Nat i) : _i(i) { }
    operator Nat () const { return _i; }
};

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
template<class X>
class Expression {
    typedef Index I;
    typedef Identifier V;
  public:
    typedef X ConstantType;
  public:
    explicit Expression(const ExpressionNode<X>* eptr) : _root(eptr) { }
    explicit Expression(counted_pointer< const ExpressionNode<X> > eptr) : _root(eptr) { }
  public:
    //! \brief Construct the constant expression with the default value of \a X.
    Expression();
    //! \brief Construct an expression from a value.
    Expression(const X& c);
    //! \brief Construct an expression from an index.
    Expression(const Index& i);
    //! \brief Construct an expression from a name.
    Expression(const Identifier& v);
    //! \brief Construct an expression from a constant.
    Expression(const Constant<X>& c);
    //! \brief Construct an expression from a variable.
    Expression(const Variable<X>& v);

    Expression(OperatorCode op, const Expression<X>& arg);
    Expression(OperatorCode op, const Expression<X>& arg1, int num);
    Expression(OperatorCode op, const Expression<X>& arg1, const Expression<X>& arg2);
    template<class A> Expression(OperatorCode op, const Expression<A>& arg1, const Expression<A>& arg2);
  public:
    OperatorCode op() const;
    OperatorKind kind() const;
    const X& val() const;
    const I& ind() const;
    const V& var() const;
    const Expression<X>& arg() const;
    const Int& num() const;
    const Expression<X>& arg1() const;
    const Expression<X>& arg2() const;
    template<class A> const Expression<A>& cmp1(A* dummy=0) const;
    template<class A> const Expression<A>& cmp2(A* dummy=0) const;
  public:
    Set<UntypedVariable> arguments() const;
  public:
    const ExpressionNode<X>* node_ptr() const { return _root.operator->(); }
  private:
    counted_pointer< const ExpressionNode<X> > _root;
};


template<class X>
struct ExpressionNode {
    mutable uint count;
    OperatorCode op;
    OperatorKind knd;
    virtual ~ExpressionNode();
    explicit ExpressionNode(OperatorCode o) : count(0u), op(o), knd(kind(op)) { }
    explicit ExpressionNode(OperatorCode o, OperatorKind k) : count(0u), op(o), knd(k) { }
};

template<class X> struct ConstantExpressionNode : public ExpressionNode<X> {
    X val;
    ConstantExpressionNode(const X& v) : ExpressionNode<X>(CNST,NULLARY), val(v) { }
};
template<class X> struct IndexExpressionNode : public ExpressionNode<X> {
    Index ind;
    IndexExpressionNode(const Index& i) : ExpressionNode<X>(IND,COORDINATE), ind(i) { }
};
template<class X> struct VariableExpressionNode : public ExpressionNode<X> {
    Identifier var;
    VariableExpressionNode(const Identifier& v) : ExpressionNode<X>(VAR,VARIABLE), var(v) { }
};
template<class X, class A=X> struct UnaryExpressionNode : public ExpressionNode<X> {
    Expression<A> arg;
    UnaryExpressionNode(OperatorCode op, Expression<A> const& a) : ExpressionNode<X>(op,UNARY), arg(a) { }
    UnaryExpressionNode(OperatorCode op, OperatorKind knd, Expression<A> const& a) : ExpressionNode<X>(op,knd), arg(a) { }
};
template<class X, class A1=X, class A2=A1> struct BinaryExpressionNode {
    Expression<X> arg1; Expression<X> arg2;
};
template<class X> struct BinaryExpressionNode<X> : public ExpressionNode<X> {
    Expression<X> arg1; Expression<X> arg2;
    BinaryExpressionNode(OperatorCode op, Expression<X> const& a1, Expression<X> const& a2)
        : ExpressionNode<X>(op,BINARY), arg1(a1), arg2(a2) { }
};
template<class X> struct BinaryExpressionNode<typename Logic<X>::Type,X,X> : public ExpressionNode<typename Logic<X>::Type> {
    typedef typename Logic<X>::Type R; typedef X A;
    Expression<A> arg1; Expression<A> arg2;
    BinaryExpressionNode(OperatorCode op, Expression<A> const& a1, Expression<A> const& a2)
        : ExpressionNode<R>(op,COMPARISON), arg1(a1), arg2(a2) { }
};
template<class R, class A=R, class N=Int> struct ScalarExpressionNode : public UnaryExpressionNode<R,A> {
    N num;
    ScalarExpressionNode(OperatorCode op, Expression<R> const& a, N n)
        : UnaryExpressionNode<R,A>(op,POWER,a), num(n) { }
};

template<class X> ExpressionNode<X>::~ExpressionNode() { }

template<class X> OperatorCode Expression<X>::op() const {
    return node_ptr()->op; }
template<class X> OperatorKind Expression<X>::kind() const {
    return node_ptr()->knd; }
template<class X> const X& Expression<X>::val() const {
    return static_cast<const ConstantExpressionNode<X>*>(node_ptr())->val; }
template<class X> const Index& Expression<X>::ind() const {
    return static_cast<const IndexExpressionNode<X>*>(node_ptr())->ind; }
template<class X> const Identifier& Expression<X>::var() const {
    return static_cast<const VariableExpressionNode<X>*>(node_ptr())->var; }
template<class X> const Expression<X>& Expression<X>::arg() const {
    return static_cast<const UnaryExpressionNode<X>*>(node_ptr())->arg; }
template<class X> const Int& Expression<X>::num() const {
    return static_cast<const ScalarExpressionNode<X>*>(node_ptr())->num; }
template<class X> const Expression<X>& Expression<X>::arg1() const {
    return static_cast<const BinaryExpressionNode<X>*>(node_ptr())->arg1; }
template<class X> const Expression<X>& Expression<X>::arg2() const {
    return static_cast<const BinaryExpressionNode<X>*>(node_ptr())->arg2; }
template<class R> template<class A> const Expression<A>& Expression<R>::cmp1(A*) const {
    return static_cast<const BinaryExpressionNode<R,A>*>(node_ptr())->arg1; }
template<class R> template<class A> const Expression<A>& Expression<R>::cmp2(A*) const {
    return static_cast<const BinaryExpressionNode<R,A>*>(node_ptr())->arg2; }


template<class X> inline Expression<X>::Expression() : _root(new ConstantExpressionNode<X>(X())) { }
template<class X> inline Expression<X>::Expression(const X& c) : _root(new ConstantExpressionNode<X>(c)) { }
template<class X> inline Expression<X>::Expression(const Index& i) : _root(new IndexExpressionNode<X>(i)) { }
template<class X> inline Expression<X>::Expression(const Identifier& v) : _root(new VariableExpressionNode<X>(v)) { }
template<class X> inline Expression<X>::Expression(const Constant<X>& c): _root(new ConstantExpressionNode<X>(c.value())) { };
template<class X> inline Expression<X>::Expression(const Variable<X>& v) : _root(new VariableExpressionNode<X>(v.name())) { }
template<class X> inline Expression<X>::Expression(OperatorCode op, const Expression<X>& arg)
    : _root(new UnaryExpressionNode<X>(op,arg)) { }
template<class X> inline Expression<X>::Expression(OperatorCode op, const Expression<X>& arg1, const Expression<X>& arg2)
    : _root(new BinaryExpressionNode<X>(op,arg1,arg2)) { }
template<class X> inline Expression<X>::Expression(OperatorCode op, const Expression<X>& arg, Int n)
    : _root(new ScalarExpressionNode<X,X,Int>(op,arg,n)) { }
template<class X> template<class A> inline Expression<X>::Expression(OperatorCode op, const Expression<A>& arg1, const Expression<A>& arg2)
    : _root(new BinaryExpressionNode<X,A>(op,arg1,arg2)) { }


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



template<class X> inline OutputStream& operator<<(OutputStream& os, const ExpressionNode<X>* e) {
    return os << (void*)(e);
}


//@}



//@{
//! \name Input / output operations.
//! \related Expression

//! \brief Write to an output stream
template<class X> OutputStream& operator<<(std::ostream& os, const Expression<X>& f) {
    switch(f.op()) {
        //case CNST: return os << std::fixed << std::setprecision(4) << fptr->val;
        case CNST:
            os << f.val(); return os;
            //if(f.val()==0.0) { return os << 0.0; } if(abs(f.val())<1e-4) { os << std::fixed << f.val(); } else { os << f.val(); } return os;
        case IND:
            return os << "x[" << f.ind() << "]";
        case VAR:
            return os << "[" << f.var() << "]";
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
                case COMPARISON: return os << "(" << f.arg1() << symbol(f.op()) << f.arg2() << ")";
                default: ARIADNE_FAIL_MSG("Cannot output expression with operator "<<f.op()<<" of kind "<<f.kind()<<"\n");
            }
    }
}

//@}

//@{
//! \name Evaluation and related operations.
//! \related Expression


template<class R, class A> struct Get { typedef Identifier I; R operator()(const I& v, const Map<I,A>& x) { assert(false); } };
template<class R> struct Get<R,R> { typedef Identifier I; R operator()(const I& v, const Map<I,R>& x) { return x[v]; } };

template<class X, class A> struct EvaluationResult { typedef X Type; };
template<class Y> struct EvaluationResult<Real,Y> { typedef Y Type; };
template<> struct EvaluationResult<Real,String> { typedef Real Type; };
template<> struct EvaluationResult<Real,Integer> { typedef Real Type; };

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
template<class X>
typename EvaluationResult<X,X>::Type
evaluate(const Expression<X>& e, const Map<Identifier,X>& x)
{
    typedef typename EvaluationResult<X,X>::Type R;
    switch(e.kind()) {
        case VARIABLE: return x[e.var()];
        case NULLARY: return static_cast<R>(e.val());
        case UNARY: return compute(e.op(),evaluate(e.arg(),x));
        case BINARY: return compute(e.op(),evaluate(e.arg1(),x),evaluate(e.arg2(),x));
        case POWER: return compute(e.op(),evaluate(e.arg(),x),e.num());
        default: ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" on "<<x<<"\n");
    }
}

template<class X, class R> X make_constant(const Vector<X>& v, const R& c) {
    ARIADNE_ASSERT(v.size()!=0);
    return v[0]*0+c;
}

//! \brief Evaluate an expression in numbered coordinates on a vector of a compatible type.
template<class Y>
Y evaluate(const Expression<Real>& e, const Vector<Y>& x)
{
    switch(e.kind()) {
        case COORDINATE: return x[e.ind()];
        case NULLARY: return make_constant(x,e.val());
        case UNARY: return compute(e.op(),evaluate(e.arg(),x));
        case BINARY: return compute(e.op(),evaluate(e.arg1(),x),evaluate(e.arg2(),x));
        default: ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" on "<<x<<"\n");
    }
}

template<class Y> Y& cached_evaluate(const Expression<Real>& e, const Vector<Y>& x, Map<const void*,Y>& cache) {
    ExpressionNode<Real> const* eptr=e.node_ptr();
    if(cache.has_key(eptr)) { return cache[eptr]; }
    switch(e.kind()) {
        case COORDINATE: return cache[eptr]=x[e.ind()];
        case NULLARY: return cache[eptr]=make_constant(x,e.val());
        case UNARY: return cache[eptr]=compute(e.op(),cached_evaluate(e.arg(),x,cache));
        case BINARY: return cache[eptr]=compute(e.op(),cached_evaluate(e.arg1(),x,cache),cached_evaluate(e.arg2(),x,cache));
        default: ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" on "<<x<<"\n");
    }
}


//! \brief Extract the arguments of expression \a e.
template<class R> Set<Identifier> arguments(const Expression<R>& e)
{
    typedef Identifier I;
    switch(e.kind()) {
        case VARIABLE: return Set<I>(e.var());
        case NULLARY: return Set<I>();
        case UNARY: return arguments(e.arg());
        case BINARY: return join(arguments(e.arg1()),arguments(e.arg2()));
        default: ARIADNE_FAIL_MSG("Cannot compute arguments of expression "<<e<<"\n");
    }
}

template<class R> Set<UntypedVariable> Expression<R>::arguments() const {
    const Expression<R>& e=*this;
    switch(e.kind()) {
        case VARIABLE: return Set<UntypedVariable>(Variable<R>(e.var()));
        case NULLARY: return Set<UntypedVariable>();
        case UNARY: return e.arg().arguments();
        case BINARY: return join(e.arg1().arguments(),e.arg2().arguments());
        case COMPARISON: {
            const BinaryExpressionNode<R,Real>* rlp = dynamic_cast<const BinaryExpressionNode<R,Real>*>(e.node_ptr());
            if(rlp) { return join(rlp->arg1.arguments(),rlp->arg2.arguments()); }
            const BinaryExpressionNode<R,String>* strp = dynamic_cast<const BinaryExpressionNode<R,String>*>(e.node_ptr());
            if(strp) { return join(strp->arg1.arguments(),strp->arg2.arguments()); }
        }
        default: ARIADNE_FAIL_MSG("Cannot compute arguments of expression "<<e<<"\n");
    }
}

template<class I, class X, class J> inline X& insert(Map<I,X>& m, const J& k, const X& v) {
    return m.std::map<I,X>::insert(std::make_pair(k,v)).first->second; }

//! \brief Convert the expression with index type \c I to one with variables indexed by \a J.
template<class X> const Expression<X>&
cached_convert(const Expression<X>& e, const Map<Identifier,Nat>& v, Map< const Void*, Expression<X> >& cache)
{
    const ExpressionNode<X>* eptr=e.node_ptr();
    if(cache.has_key(eptr)) { return cache.get(eptr); }
    switch(e.kind()) {
        case VARIABLE: return insert( cache, eptr, Expression<X>(v[e.var()]) );
        case NULLARY: return insert( cache, eptr, Expression<X>(e.val()) );
        case UNARY: return insert( cache, eptr, Expression<X>(e.op(),convert(e.arg(),v)));
        case BINARY: return insert( cache, eptr, Expression<X>(e.op(),convert(e.arg1(),v),convert(e.arg2(),v)) );
        default: ARIADNE_FAIL_MSG("Cannot convert expression "<<e<<" to use variables "<<v<<"\n");
    }
}

//! \brief Convert the expression with index type \c I to one with variables indexed by \a J.
template<class X> Expression<X> convert(const Expression<X>& e, const Map<Identifier,Nat>& v)
{
    Map< const Void*, Expression<X> > cache;
    return cached_convert(e,v,cache);
}


//! \brief Returns \a true if the expression\a e is syntactically equal to the constant \a c.
template<class X> Bool is_constant(const Expression<X>& e, const X& c) {
    switch(e.op()) {
        case CNST: return e.val()==c;
        default: return false;
    }
}

template<class X, class R> Bool is_constant(const Expression<X>& e, const R& c) {
    return is_constant(e,static_cast<X>(c));
}

//! \brief Returns \a true if the expression\a e is syntactically equal to the constant \a c.
template<class X> Bool is_variable(const Expression<X>& e, const Identifier& v) {
    switch(e.op()) {
        case VAR: return e.var()==v;
        default: return false;
    }
}

//! \brief Simplify the expression \a e.
template<class X> Expression<X> simplify(const Expression<X>& e);

//! \brief Tests whether two expressions are identical.
template<class X> Bool identical(const Expression<X>& e1, const Expression<X>& e2);

//! \brief Compute the derivative of expression \a e with respect to the variable \a v.
template<class R, class I> Expression<R> derivative(const Expression<R>& e, const I& v);


//@}


//@{
//! \name Deprecated conversions and Expression/Formula operators.
//! \related Expression

Boolean evaluate(const Expression<Boolean>& e, const DiscreteValuation& q);
String evaluate(const Expression<String>& e, const StringValuation& q);
Integer evaluate(const Expression<Integer>& e, const IntegerValuation& q);

template<class X> Tribool evaluate(const Expression<Tribool>& e, const ContinuousValuation<X>& x);
template<class X> X evaluate(const Expression<Real>& e, const ContinuousValuation<X>& x);
//template<class X> X evaluate(const Expression<Real>& e, const Map<ExtendedVariable<Real>,X>& x);

//! \brief Convert a real expression with index type Identifier to one with variables indexed by Nat.
Formula<Real> formula(const Expression<Real>& e, const Space<Real>& spc);

Formula<Real> formula(const Expression<Real>& e, const List< Variable<Real> >& vars);
Formula<Real> formula(const Expression<Real>& res, const List< Assignment< Variable<Real>, Expression<Real> > >& aux, const Space<Real> spc);
List< Formula<Real> > formula(const List< Expression<Real> >& res, const List< Assignment< Variable<Real>, Expression<Real> > >& aux, const Space<Real> spc);

//! \related Expression \brief Returns true if the expressions are mutual negations.
//!
//! Currently can only test for pairs of the form (a1<=a2; a1>=a2),  (a1<=a2; a2<=a1)
//! or (a1>=a2; a2>=a1).
Bool opposite(Expression<Tribool> p, Expression<Tribool> q);

//! \brief Given \a sign when the predicate \a p is true.
Expression<Real> indicator(Expression<Tribool> p, Sign sign=POSITIVE);

//! \brief Substitute all occurrences of variable \a v of type \c Y with constant value \a c.
template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Y& c);

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const List< Assignment< Variable<Y>,Expression<Y> > >& a);

//! \brief Substitute all occurrences of variable \a v of type \c Y with expression value \a se.
template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Expression<Y>& se);

template<class X> class Affine;
template<class X> class Polynomial;

template<class X> Affine<X> affine(const Formula<Real>&);
template<class X> Polynomial<X> polynomial(const Formula<Real>&);

inline Expression<Real> derivative(const Expression<Real>& e, const Variable<Real>& v) {
    return derivative(e,v.name());
}

ScalarFunction<Real> make_function(const Expression<Real>& e, const Space<Real>& s);

Bool is_variable(const Expression<Real>& e, const Identifier& v);

//@}

//@{
//! \name Metaprogramming and sequencing operators for making lists
//! \related Expression

//! \brief Inherits from True type if class \a T is a constant, variable or expression in type \a X, otherwise inherits from False.
template<class X, class T> struct IsExpression;

template<class X, class T> struct IsExpression : public False { };
template<class X> struct IsExpression< X, X > : public True { };
template<class X> struct IsExpression< X, Constant<X> > : public True { };
template<class X> struct IsExpression< X, Variable<X> > : public True { };
template<class X> struct IsExpression< X, Expression<X> > : public True { };
template<class X, class T1, class T2, class R> struct EnableIfExpressions : public EnableIf< And<IsExpression<X,T1>,IsExpression<X,T2> >, R> { };

template<class T1, class T2> inline
typename EnableIfExpressions< Real, T1, T2, List<Expression<Real> > >::Type
operator,(const T1& e1, const T2& e2) {
    List< Expression<Real> > r; r.append(e1); r.append(e2); return r; }

template<class T> inline
typename EnableIf< IsExpression<Real,T>, List<Expression<Real> > >::Type
operator,(Int c, const T& e) {
    List< Expression<Real> > r; r.append(Expression<Real>(c)); r.append(Expression<Real>(e)); return r; }

template<class X, class T> inline
typename EnableIf< IsExpression<X,T>, List<Expression<X> > >::Type
operator,(List<Expression<X> > l, const T& e) {
    List< Expression<X> > r(l); r.append(Expression<Real>(e)); return r; }

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

//! \related Expression \brief String equality.
Expression<Boolean> operator==(Variable<String> v1, const String& s2);
//! \related Expression \brief String inequality.
Expression<Boolean> operator!=(Variable<String> v1, const String& s2);


//! \related Expression \brief Integer equality predicate.
Expression<Boolean> operator==(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief Integer inequality predicate.
Expression<Boolean> operator!=(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief Integer comparison predicate (greater or equal).
Expression<Boolean> operator>=(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief Integer comparison (less or equal)..
Expression<Boolean> operator<=(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief Integer comparison (greater).
Expression<Boolean> operator> (Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief Integer comparison (less).
Expression<Boolean> operator< (Expression<Integer> e1, Expression<Integer> e2);

//! \related Expression \brief Integer unary plus expression (identity).
Expression<Integer> operator+(Expression<Integer> e);
//! \related Expression \brief Integer unary minus expression.
Expression<Integer> operator-(Expression<Integer> e);
//! \related Expression \brief Integer addition expression.
Expression<Integer> operator+(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief Integer subtraction expression.
Expression<Integer> operator-(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief Integer multiplication expression.
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

//! \related Expression \brief Real unary plus expression.
Expression<Real> operator+(Expression<Real> e);
//! \related Expression \brief Real unary minus expression.
Expression<Real> operator-(Expression<Real> e);
//! \related Expression \brief Real addition expression.
Expression<Real> operator+(Expression<Real> e1, Expression<Real> e2);
//! \related Expression \brief Real subtraction expression.
Expression<Real> operator-(Expression<Real> e1, Expression<Real> e2);
//! \related Expression \brief Real multiplication expression.
Expression<Real> operator*(Expression<Real> e1, Expression<Real> e2);
//! \related Expression \brief Real division expression.
Expression<Real> operator/(Expression<Real> e1, Expression<Real> e2);

//! \related Expression \brief Real integer power expression.
Expression<Real> pow(Expression<Real> e, Int n);

//! \related Expression \brief Real negation expression.
//! Equivalent to -\a e.
Expression<Real> neg(Expression<Real> e);
//! \related Expression \brief Real reciprocal expression.
//! Equivalent to 1/\a e.
Expression<Real> rec(Expression<Real> e);
//! \related Expression \brief Real square expression.
Expression<Real> sqr(Expression<Real> e);
//! \related Expression \brief Real square root expression.
Expression<Real> sqrt(Expression<Real> e);
//! \related Expression \brief Real exponential expression.
Expression<Real> exp(Expression<Real> e);
//! \related Expression \brief Real natural logarithm expression.
Expression<Real> log(Expression<Real> e);
//! \related Expression \brief Real sine expression.
Expression<Real> sin(Expression<Real> e);
//! \related Expression \brief Real cosine expression.
Expression<Real> cos(Expression<Real> e);
//! \related Expression \brief Real tangent expression.
Expression<Real> tan(Expression<Real> e);

//! \related Expression \brief Real maximum expression.
Expression<Real> max(Expression<Real> e1, Expression<Real> e2);
//! \related Expression \brief Real minimum expression.
Expression<Real> min(Expression<Real> e1, Expression<Real> e2);
//! \related Expression \brief Real absolute value expression.
Expression<Real> abs(Expression<Real> e);

//@}

template<class R> Expression<R> make_expression(double c) { return Expression<R>(static_cast<R>(c)); }

inline Expression<Real> operator+(Expression<Real> e, double c) { return e + make_expression<Real>(c); }
inline Expression<Real> operator-(Expression<Real> e, double c) { return e - make_expression<Real>(c); }
inline Expression<Real> operator*(Expression<Real> e, double c) { return e - make_expression<Real>(c); }
inline Expression<Real> operator/(Expression<Real> e, double c) { return e - make_expression<Real>(c); }
inline Expression<Real> operator+(double c, Expression<Real> e) { return make_expression<Real>(c) + e; }
inline Expression<Real> operator-(double c, Expression<Real> e) { return make_expression<Real>(c) - e; }
inline Expression<Real> operator*(double c, Expression<Real> e) { return make_expression<Real>(c) * e; }
inline Expression<Real> operator/(double c, Expression<Real> e) { return make_expression<Real>(c) / e; }

inline Expression<Tribool> operator<=(Expression<Real> e, double c) { return e <= make_expression<Real>(c); }
inline Expression<Tribool> operator< (Expression<Real> e, double c) { return e < make_expression<Real>(c); }
inline Expression<Tribool> operator>=(Expression<Real> e, double c) { return e >= make_expression<Real>(c); }
inline Expression<Tribool> operator> (Expression<Real> e, double c) { return e >  make_expression<Real>(c); }
inline Expression<Tribool> operator<=(double c, Expression<Real> e) { return make_expression<Real>(c) <= e; }
inline Expression<Tribool> operator< (double c, Expression<Real> e) { return make_expression<Real>(c) <  e; }
inline Expression<Tribool> operator>=(double c, Expression<Real> e) { return make_expression<Real>(c) >= e; }
inline Expression<Tribool> operator> (double c, Expression<Real> e) { return make_expression<Real>(c) >  e; }



} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_H */
