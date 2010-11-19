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
#include "pointer.h"
#include "container.h"

#include "operators.h"
#include "variables.h"
#include "vector.h"

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
template<class T> class Variable;
template<class T> class Space;
template<class T> class Expression;
template<class LHS,class RHS> class Assignment;

class DiscreteValuation;
template<class X> class ContinuousValuation;
template<class X> class Vector;

class Substitution;

typedef Expression<Boolean> DiscretePredicate;
typedef Expression<Tribool> ContinuousPredicate;
typedef Expression<String> StringExpression;
typedef Expression<Integer> IntegerExpression;
typedef Expression<Real> RealExpression;

template<class T>
class ExpressionInterface
{
    friend class Expression<T>;
  public:
    virtual ~ExpressionInterface() { }
    virtual ExpressionInterface<T>* clone() const = 0;
    virtual Set<UntypedVariable> arguments() const = 0;
    virtual Operator op() const = 0;
    virtual std::ostream& write(std::ostream& os) const = 0;
  protected:
    virtual ExpressionInterface<T>* simplify() const = 0;
};

template<class R> class Expression;
template<class R> std::ostream& operator<<(std::ostream&, const Expression<R>&);

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Y& c);
template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Expression<Y>& c);
template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const List< Assignment< Variable<Y>,Expression<Y> > >& a);

template<class X> Expression<X> simplify(const Expression<X>& e);


struct ExpressionNode
{
    void* operator new(size_t) { return boost::singleton_pool<ExpressionNode,sizeof(ExpressionNode)>::malloc(); }
    void operator delete(void* ptr) { boost::singleton_pool<ExpressionNode,sizeof(ExpressionNode)>::free(ptr); }

    unsigned int _count;
    Operator _operator;
    Float _constant;
    union {
        void* _variable;
        unsigned int _coordinate;
        ExpressionNode* _arg;
        struct { ExpressionNode* _arg1; ExpressionNode* _arg2; };
        //struct { ExpressionNode* _arg; int _power; };
    };
};

/*
template<class R, class A>
R _evaluate(ExpressionNode* e, const Vector<A>& x)
{
    switch(e->_operator) {
        case CNST: return static_cast<R>(e->_constant);
        case IND: return x[e->_coordinate];
        case ADD: return add(_evaluate<R>(e->_arg1,x),_evaluate<R>(e->_arg2,x));
        case SUB: return sub(_evaluate<R>(e->_arg1,x),_evaluate<R>(e->_arg2,x));
        case MUL: return mul(_evaluate<R>(e->_arg1,x),_evaluate<R>(e->_arg2,x));
        case DIV: return div(_evaluate<R>(e->_arg1,x),_evaluate<R>(e->_arg2,x));
        case MAX: return max(_evaluate<R>(e->_arg1,x),_evaluate<R>(e->_arg2,x));
        case MIN: return min(_evaluate<R>(e->_arg1,x),_evaluate<R>(e->_arg2,x));
        //case POW: return pow(_evaluate<R>(e->_arg,x),e->_power);
        case POS: return pos(_evaluate<R>(e->_arg,x));
        case NEG: return neg(_evaluate<R>(e->_arg,x));
        case ABS: return abs(_evaluate<R>(e->_arg,x));
        case SQR: return sqr(_evaluate<R>(e->_arg,x));
        case SQRT: return sqrt(_evaluate<R>(e->_arg,x));
        case EXP: return exp(_evaluate<R>(e->_arg,x));
        case LOG: return log(_evaluate<R>(e->_arg,x));
        case SIN: return sin(_evaluate<R>(e->_arg,x));
        case COS: return cos(_evaluate<R>(e->_arg,x));
        case TAN: return tan(_evaluate<R>(e->_arg,x));
        case ASIN: return asin(_evaluate<R>(e->_arg,x));
        case ACOS: return acos(_evaluate<R>(e->_arg,x));
        case ATAN: return atan(_evaluate<R>(e->_arg,x));
        default: assert(false);
    }
}
*/

template<class R, class A>
R _evaluate(ExpressionNode* e, const ContinuousValuation<A>& x)
{
    switch(e->_operator) {
        case CNST: return static_cast<R>(e->_constant);
        case VAR: return x[e->_variable];
        case ADD: return add(_evaluate<R>(e->_arg1,x),_evaluate<R>(e->_arg2,x));
        case SUB: return sub(_evaluate<R>(e->_arg1,x),_evaluate<R>(e->_arg2,x));
        case MUL: return mul(_evaluate<R>(e->_arg1,x),_evaluate<R>(e->_arg2,x));
        case DIV: return div(_evaluate<R>(e->_arg1,x),_evaluate<R>(e->_arg2,x));
        case MAX: return max(_evaluate<R>(e->_arg1,x),_evaluate<R>(e->_arg2,x));
        case MIN: return min(_evaluate<R>(e->_arg1,x),_evaluate<R>(e->_arg2,x));
        //case POW: return pow(_evaluate<R>(e->_arg,x),e->_power);
        case POS: return pos(_evaluate<R>(e->_arg,x));
        case NEG: return neg(_evaluate<R>(e->_arg,x));
        case ABS: return abs(_evaluate<R>(e->_arg,x));
        case SQR: return sqr(_evaluate<R>(e->_arg,x));
        case SQRT: return sqrt(_evaluate<R>(e->_arg,x));
        case EXP: return exp(_evaluate<R>(e->_arg,x));
        case LOG: return log(_evaluate<R>(e->_arg,x));
        case SIN: return sin(_evaluate<R>(e->_arg,x));
        case COS: return cos(_evaluate<R>(e->_arg,x));
        case TAN: return tan(_evaluate<R>(e->_arg,x));
        case ASIN: return asin(_evaluate<R>(e->_arg,x));
        case ACOS: return acos(_evaluate<R>(e->_arg,x));
        case ATAN: return atan(_evaluate<R>(e->_arg,x));
        default: assert(false);
    }
}

//! \ingroup ExpressionModule
//! \brief A simple expression in named variables.
//!
//! Ariadne supports expressions of type Bool, Tribool, String, Integer and Real.
//! The class Real is a dummy type which can be implemented in many different ways.
//!
//! The independent variables are given string names, rather than an integer index.
//! Formulae in different variables may be combined; the variables of the resulting formula
//! are all variables occuring in all formulae.
//! Formulae may be manipulated symbolically.
//! \sa Variable, ScalarFunction
template<class R>
class Expression {
  public:
    explicit Expression(const ExpressionInterface<R>& e) : _ptr(e.clone()) { }
    explicit Expression(const ExpressionInterface<R>* eptr) : _ptr(eptr) { }
    explicit Expression(shared_ptr< const ExpressionInterface<R> > eptr) : _ptr(eptr) { }
    Expression() { R z; *this=Expression(z); }
    Expression(const R& c);
    Expression(const Constant<R>& c);
    Expression(const Variable<R>& v);
    //! \brief The variables used in the formula.
    Set<UntypedVariable> arguments() const { return _ptr->arguments(); }
    //! \brief The operator used to construct the expression from subexpressions.
    Operator op() const { return _ptr->op(); }
    //! \brief The immediate subexpressions used in the formula.
    List< Expression<R> > subexpressions() const;
     //! \brief Substitute the constant \a c for the variable \a v.
    template<class X> Expression<R> substitute(const Variable<X>& v, const X& c) const {
        return Ariadne::substitute(*this,v,c); };
    //! \brief Simplify the expression (e.g. by evaluating constants).
    Expression<R> simplify() const {
        return Ariadne::simplify(*this); }
    //! \brief Convert to a list of expressions
    operator List< Expression<R> >() const { return List< Expression<R> >(1u,*this); }
    //! \brief Write to an output stream.
    friend std::ostream& operator<< <>(std::ostream&, const Expression<R>&);
  public:
    // Return a raw pointer to the underlying representation.
    const ExpressionInterface<R>* _raw_pointer() const { return this->_ptr.operator->(); }
    const ExpressionInterface<R>* ptr() const { return _ptr.operator->(); }
  public:
    shared_ptr<const ExpressionInterface<R> > _ptr;
};

template<>
class Expression<Real> {
    typedef Real R;
  public:
    explicit Expression(const ExpressionInterface<R>& e) : _ptr(e.clone()) { }
    explicit Expression(const ExpressionInterface<R>* eptr) : _ptr(eptr) { }
    explicit Expression(shared_ptr< const ExpressionInterface<R> > eptr) : _ptr(eptr) { }
    Expression(const double& c);
    Expression(const Interval& c);
    Expression(const Real& c);
    Expression(const Constant<R>& c);
    Expression(const Variable<R>& v);
    //! \brief Test if two expressions are identical to each other.
    friend bool identical(const Expression<R>& e1, const Expression<R>& e2);
    //! \brief The variables used in the formula.
    Set<UntypedVariable> arguments() const { return _ptr->arguments(); }
    //! \brief The operator used to construct the expression from subexpressions.
    Operator op() const { return _ptr->op(); }
    //! \brief The immediate subexpressions used in the formula.
    List< Expression<R> > subexpressions() const;
    //! \brief Substitute the constant \a c for the variable \a v.
    template<class X> Expression<R> substitute(const Variable<X>& v, const X& c) const {
        return Ariadne::substitute(*this,v,c); }
    //! \brief Simplify the expression (e.g. by evaluating constants).
    Expression<R> simplify() const {
        return Ariadne::simplify(*this); }
    //! \brief Convert to a list of expressions
    operator List< Expression<R> >() const { return List< Expression<R> >(1u,*this); }
    //! \brief Write to an output stream.
    friend std::ostream& operator<< <>(std::ostream&, const Expression<R>&);
  public:
    // Return a raw pointer to the underlying representation.
    const ExpressionInterface<R>* _raw_pointer() const { return this->_ptr.operator->(); }
    const ExpressionInterface<R>* ptr() const { return _ptr.operator->(); }
  public:
    shared_ptr<const ExpressionInterface<R> > _ptr;
};


bool identical(const Expression<Real>& e1, const Expression<Real>& e2);

template<class R> inline std::ostream& operator<<(std::ostream& os, const Expression<R>& e) { return e._ptr->write(os); }

Boolean evaluate(const Expression<Boolean>& e, const DiscreteValuation& q);
String evaluate(const Expression<String>& e, const DiscreteValuation& q);
Integer evaluate(const Expression<Integer>& e, const DiscreteValuation& q);

template<class X> Tribool evaluate(const Expression<Tribool>& e, const ContinuousValuation<X>& x);
template<class X> X evaluate(const Expression<Real>& e, const ContinuousValuation<X>& x);
template<class X> X evaluate(const Expression<Real>& e, const Map<ExtendedVariable<Real>,X>& x);

template<class X> Tribool evaluate(const Expression<Tribool>& e, const Vector<X>& x);
template<class X> X evaluate(const Expression<Real>& e, const Vector<X>& x);

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Y& c);
template<class X> Expression<X> simplify(const Expression<X>& e);

bool operator==(const Expression<Tribool>&, bool);

Expression<Real> function(const Expression<Real>& e, const Space<Real>& s);

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


//! \related Expression \brief Equality predicate.
Expression<Boolean> operator==(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief Inequality predicate.
Expression<Boolean> operator!=(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief Integer comparison predicate (greater or equal).
Expression<Boolean> operator>=(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief Integer comparison (less or equal)..
Expression<Boolean> operator<=(Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief Integer comparison (greater).
Expression<Boolean> operator> (Expression<Integer> e1, Expression<Integer> e2);
//! \related Expression \brief Integer comparison (less)..
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
Expression<Real> pow(Expression<Real> e, int n);

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

enum Sign { NEGATIVE=-1, ZERO=0, POSITIVE=+1 };

//! \related Expression \brief Try to compute a real expression which has the
//! given \a sign when the predicate \a p is true.
Expression<Real> indicator(Expression<tribool> p, Sign sign=POSITIVE);

//! \related Expression \brief Returns true if the expressions are mutual negations.
//!
//! Currently can only test for pairs of the form (a1<=a2; a1>=a2),  (a1<=a2; a2<=a1)
//! or (a1>=a2; a2>=a1).
tribool opposite(Expression<tribool> p, Expression<tribool> q);


template<class X> class Affine;
template<class X> class Polynomial;

//! \related Expression \brief Compute the derivative of expression \a e with respect to the variable \a v.
Expression<Real> derivative(const Expression<Real>& e, const Variable<Real>& v);

template<class X> Affine<X> affine(const Expression<Real>&, const Space<Real>&);
template<class X> Polynomial<X> polynomial(const Expression<Real>&, const Space<Real>&);

} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_H */
