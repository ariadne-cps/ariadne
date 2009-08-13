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

//#include "untyped_expression.h"

/*! \file expression.h
 *  \brief Internal expressions
 */

#ifndef ARIADNE_EXPRESSION_H
#define ARIADNE_EXPRESSION_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function_interface.h"

#include "macros.h"
#include "pointer.h"

#include "operators.h"
#include "container.h"

namespace Ariadne {


static const int SMOOTH=255;

typedef std::string String;
class Integer;
class Real;
class EnumeratedValue;

typedef String Identifier;

template<class T> class Variable;
template<class R> class Expression;
template<class LHS,class RHS> class Assignment;

template<class T> class Variable
{
  public:
    explicit Variable(const String& str) : _name(new String(str)) { }
    const String& name() const { return *_name; }
    bool operator==(const Variable<T>& other) { return this->name()==other.name(); }
    Assignment<Variable<T>,Expression<T> > operator=(const Expression<T>& e) const;
    template<class X> bool operator==(const Variable<X>& other) { ARIADNE_ASSERT(this->name()!=other.name()); return false; }
  private:
    shared_ptr<String> _name;
};

template<class T> inline bool operator<(const Variable<T> v1, const Variable<T>& v2) { return v1.name()<v2.name(); }

inline std::ostream& operator<<(std::ostream& os, const Variable<EnumeratedValue>& v) { return os << v.name() << ":Enumerated"; }
inline std::ostream& operator<<(std::ostream& os, const Variable<String>& v) { return os << v.name() << ":String"; }
inline std::ostream& operator<<(std::ostream& os, const Variable<Integer>& v) { return os << v.name() << ":Integer"; }
inline std::ostream& operator<<(std::ostream& os, const Variable<Real>& v) { return os << v.name() << ":Real"; }


/*! \brief A simple expression in named variables.
 *  The independent variables are given string names, rather than an integer index.
 *  Formulae in different variables may be combined; the variables of the resulting formula
 *  are all variables occuring in all formulae.
 *  Formulae may be manipulated symbolically.
 *  \sa RealVariable, ScalarFunctionInterface
 */
template<class T>
class ExpressionInterface
{
  public:
    virtual ~ExpressionInterface() { }
    virtual ExpressionInterface<T>* clone() const = 0;
    virtual Set<String> arguments() const = 0;
    virtual std::ostream& write(std::ostream& os) const = 0;
};

template<class T>
inline std::ostream& operator<<(std::ostream& os, const ExpressionInterface<T>& e) {
    return e.write(os);
}

typedef shared_ptr< ExpressionInterface<Real> > RealExpressionPointer;


//! A constant, viewed as a function \f$c:\R^n\rightarrow\R\f$.
template<class R>
class ConstantExpression
    : public ExpressionInterface<R>
{
  public:
    ConstantExpression(const R& c) : _c(c) { }
    operator R () const { return _c; }
    R value() const { return _c; }
    virtual Set<String> arguments() const { return Set<String>(); }
    virtual ConstantExpression<R>* clone() const { return new ConstantExpression<R>(*this); }
    virtual std::ostream& write(std::ostream& os) const { return os << _c; }
  private:
    R _c;
};

//! A constant, viewed as a function \f$c:\R^n\rightarrow\R\f$.
template<>
class ConstantExpression<Real>
    : public ExpressionInterface<Real>
{
  public:
    ConstantExpression(const Float& c) : _c(c) { }
    ConstantExpression(const Interval& c) :  _c(c) { }
    operator Interval() const { return _c; }
    Interval value() const { return _c; }
    virtual Set<String> arguments() const { return Set<String>(); }
    virtual ConstantExpression<Real>* clone() const { return new ConstantExpression<Real>(*this); }
    virtual std::ostream& write(std::ostream& os) const {
        if(_c.lower()==_c.upper()) { os<<midpoint(_c); } else { os<<_c; } return os; }
  private:
    Interval _c;
};


//! A projection onto a named variable.
template<class R>
class VariableExpression
    : public ExpressionInterface<R>
{
  public:
    explicit VariableExpression(const String& s) : _var(s) { }
    VariableExpression(const Variable<R>& v) : _var(v) { }
    String name() const { return this->_var.name(); }
    virtual VariableExpression<R>* clone() const { return new VariableExpression<R>(*this); }
    virtual Set<String> arguments() const { Set<String> r; r.insert(this->name()); return r; }
    virtual std::ostream& write(std::ostream& os) const { return os << this->name(); }
  private:
    Variable<R> _var;
};

template<class T> class CoordinateExpression;

//! A coordinate projection \f$\pi:\R^n\rightarrow\R\f$ given by \f$\pi(x)=x_j\f$.
template<>
class CoordinateExpression<Real>
    : public ExpressionInterface<Real>
{
    typedef unsigned int SizeType;
  public:
    explicit CoordinateExpression() : _as(0), _j(0) { }
    explicit CoordinateExpression(SizeType j) : _as(0), _j(j) { }
    explicit CoordinateExpression(unsigned int as, unsigned int j) : _as(as), _j(j) { }
    SizeType argument_size() const { return _as; }
    SizeType coordinate() const { return _j; }
    String name() const { String s("x0"); s[1]+=_j; return s; }
    virtual CoordinateExpression<Real>* clone() const { return new CoordinateExpression<Real>(*this); }
    virtual Set<String> arguments() const { Set<String> r; r.insert(this->name()); return r; }
    virtual std::ostream& write(std::ostream& os) const { return os << this->name(); }
  private:
    SizeType _as;
    SizeType _j;
};



//! A composed scalar function, based on a standard operator.
template<class R, class Op, class A=R> class UnaryExpression
    : public ExpressionInterface<R>
{
  public:
    UnaryExpression(Op op, const ExpressionInterface<A>& expr) : _op(op), _arg_ptr(expr.clone()) { }
    UnaryExpression(Op op, const ExpressionInterface<A>* expr) : _op(op), _arg_ptr(const_cast<ExpressionInterface<A>*>(expr)) { }
    UnaryExpression(Op op, shared_ptr< ExpressionInterface<A> > expr) : _op(op), _arg_ptr(expr) { }
    virtual UnaryExpression<R,Op,A>* clone() const { return new UnaryExpression<R,Op,A>(_op,_arg_ptr); }
    virtual Set<String> arguments() const { return this->_arg_ptr->arguments(); }
    virtual std::ostream& write(std::ostream& os) const;
  public:
    Op _op;
    shared_ptr< ExpressionInterface<A> > _arg_ptr;
};

template<class R, class Op, class A> inline std::ostream& UnaryExpression<R,Op,A>::write(std::ostream& os) const { return os << _op << "(" << *_arg_ptr << ")"; }


//! A composed scalar function, based on an arthmetic operator.
template<class R, class Op, class A1=R, class A2=A1> class BinaryExpression
    : public ExpressionInterface<R>
{
  public:
    BinaryExpression(Op op, const ExpressionInterface<A1>& expr1, const ExpressionInterface<A2>& expr2)
        : _op(op), _arg1_ptr(expr1.clone()), _arg2_ptr(expr2.clone()) { }
    BinaryExpression(Op op, const ExpressionInterface<A1>* expr1, const ExpressionInterface<A2>* expr2)
        : _op(op), _arg1_ptr(const_cast<ExpressionInterface<A1>*>(expr1)), _arg2_ptr(const_cast<ExpressionInterface<A2>*>(expr2)) { }
    BinaryExpression(Op op, shared_ptr< ExpressionInterface<A1> > expr1, shared_ptr< ExpressionInterface<A2> > expr2)
        : _op(op), _arg1_ptr(expr1), _arg2_ptr(expr2)  { }
    virtual BinaryExpression<R,Op,A1,A2>* clone() const { return new BinaryExpression<R,Op,A1,A2>(_op,_arg1_ptr,_arg2_ptr); }
    virtual Set<String> arguments() const { return join(this->_arg1_ptr->arguments(),this->_arg2_ptr->arguments()); }
    virtual std::ostream& write(std::ostream& os) const {
        return os << "(" << *_arg1_ptr << _op << *_arg2_ptr << ")"; }
/*
  private:
    template<class R, class A> inline
    void compute(R& r, const A& a) { r=Op()(_arg1->evaluate(a),_arg2->evaluate(a)); }
*/
  public:
    Op _op;
    shared_ptr< ExpressionInterface<A1> > _arg1_ptr;
    shared_ptr< ExpressionInterface<A2> > _arg2_ptr;
};

template<class R> class Expression;
template<class R> std::ostream& operator<<(std::ostream&, const Expression<R>&);

//! \brief A formula taking values of type \a R.
//! \tparam R The type of mathematical object the formula represents; e.g. String, Integer, Real
//! \sa Variable
template<class R>
class Expression {
  public:
    explicit Expression(const R& c) : ptr(new ConstantExpression<R>(c)) { }
    explicit Expression(const Variable<R>& v) : ptr(new VariableExpression<R>(v)) { }
    Expression(const ExpressionInterface<R>& e) : ptr(e.clone()) { }
    explicit Expression(ExpressionInterface<R>* p) : ptr(p) { }
    Expression(shared_ptr< ExpressionInterface<R> > p) : ptr(p) { }
    //! \brief The variables used in the formula.
    Set<String> arguments() const { return ptr->arguments(); }
    //! \brief Write to an output stream.
    friend std::ostream& operator<< <>(std::ostream&, const Expression<R>&);
  public:
    shared_ptr< ExpressionInterface<R> > ptr;
};

template<>
class Expression<Real> {
    typedef Real R;
  public:
    explicit Expression(const double& c) : ptr(new ConstantExpression<R>(c)) { }
    explicit Expression(const Interval& c) : ptr(new ConstantExpression<R>(c)) { }
    explicit Expression(const Variable<R>& v) : ptr(new VariableExpression<R>(v)) { }
    Expression(const ExpressionInterface<R>& e) : ptr(e.clone()) { }
    explicit Expression(ExpressionInterface<R>* p) : ptr(p) { }
    Expression(shared_ptr< ExpressionInterface<R> > p) : ptr(p) { }
    //! \brief The variables used in the formula.
    Set<String> arguments() const { return ptr->arguments(); }
    //! \brief Write to an output stream.
    friend std::ostream& operator<< <>(std::ostream&, const Expression<R>&);
  public:
    shared_ptr< ExpressionInterface<R> > ptr;
};

template<class R, class V> inline bool operator==(const Expression<R>& e, const V& v) {
    const ConstantExpression<R>* expr = dynamic_cast<const ConstantExpression<R>*>(e.ptr.operator->());
    return expr && expr->value()==v;
}

template<class R> inline std::ostream& operator<<(std::ostream& os, const Expression<R>& f) { return f.ptr->write(os); }

typedef Expression<Real> RealExpression;

//! \related Expression \brief .
template<class R> inline Expression<R> operator-(Expression<R> e) {
    return Expression<R>(new UnaryExpression<R,Neg>(Neg(),e.ptr)); }
//! \related Expression \brief .
template<class R> inline Expression<R> operator+(Expression<R> e1, Expression<R> e2) {
    return Expression<R>(new BinaryExpression<R,Add>(Add(),e1.ptr,e2.ptr)); }
//! \related Expression \brief .
template<class R> inline Expression<R> operator-(Expression<R> e1, Expression<R> e2) {
    return Expression<R>(new BinaryExpression<R,Sub>(Sub(),e1.ptr,e2.ptr)); }
//! \related Expression \brief .
template<class R> inline Expression<R> operator*(Expression<R> e1, Expression<R> e2) {
    return Expression<R>(new BinaryExpression<R,Mul>(Mul(),e1.ptr,e2.ptr)); }
//! \related Expression \brief .
template<class R> inline Expression<R> operator/(Expression<R> e1, Expression<R> e2) {
    return Expression<R>(new BinaryExpression<R,Div>(Div(),e1.ptr,e2.ptr)); }


//! \related Expression \brief .
inline Expression<Real> neg(Expression<Real> e) {
    return Expression<Real>(new UnaryExpression<Real,Neg>(Neg(),e.ptr)); }
//! \related Expression \brief .
inline Expression<Real> rec(Expression<Real> e) {
    return Expression<Real>(new UnaryExpression<Real,Rec>(Rec(),e.ptr)); }
//! \related Expression \brief .
inline Expression<Real> sqr(Expression<Real> e) {
    return Expression<Real>(new UnaryExpression<Real,Sqr>(Sqr(),e.ptr)); }
//! \related Expression \brief .
inline Expression<Real> pow(Expression<Real> e, int n) {
    return Expression<Real>(new UnaryExpression<Real,Pow>(Pow(n),e.ptr)); }
//! \related Expression \brief .
inline Expression<Real> sqrt(Expression<Real> e) {
    return Expression<Real>(new UnaryExpression<Real,Sqrt>(Sqrt(),e.ptr)); }
//! \related Expression \brief .
inline Expression<Real> exp(Expression<Real> e) {
    return Expression<Real>(new UnaryExpression<Real,Exp>(Exp(),e.ptr)); }
//! \related Expression \brief .
inline Expression<Real> log(Expression<Real> e) {
    return Expression<Real>(new UnaryExpression<Real,Log>(Log(),e.ptr)); }
//! \related Expression \brief .
inline Expression<Real> sin(Expression<Real> e) {
    return Expression<Real>(new UnaryExpression<Real,Sin>(Sin(),e.ptr)); }
//! \related Expression \brief .
inline Expression<Real> cos(Expression<Real> e) {
    return Expression<Real>(new UnaryExpression<Real,Cos>(Cos(),e.ptr)); }
//! \related Expression \brief .
inline Expression<Real> tan(Expression<Real> e) {
    return Expression<Real>(new UnaryExpression<Real,Tan>(Tan(),e.ptr)); }



//! \related Expression \brief .
template<class T> inline Expression<T> operator&&(Expression<T> e1, Expression<T> e2) {
    return Expression<T>(new BinaryExpression<T,And,T,T>(And(),e1.ptr,e2.ptr)); }
//! \related Expression \brief .
template<class T> inline Expression<T> operator||(Expression<T> e1, Expression<T> e2) {
    return Expression<T>(new BinaryExpression<T,Or>(Or(),e1.ptr,e2.ptr)); }
//! \related Expression \brief .
template<class T> inline Expression<T> operator!(Expression<T> e) {
    return Expression<T>(new UnaryExpression<T,Not,T>(Not(),e.ptr)); }

/*
//! \related EnumeratedVariable \brief .
inline
Expression<bool>
operator==(const EnumeratedVariable& lhs, const EnumeratedValue& rhs) {
    return Expression<bool>(new BinaryFormula<bool,Equal,EnumeratedValue,EnumeratedValue>(Equal(),VariableFormula<EnumeratedValue>(lhs),ConstantFormula<EnumeratedValue>(rhs)));
}

//! \related EnumeratedVariable \brief .
inline
Expression<bool>
operator==(const EnumeratedVariable& lhs, const String& rhs) {
    return lhs==EnumeratedValue(rhs);
}

//! \related Expression \brief .
inline
Expression<tribool>
operator<=(const Expression& lhs, const Expression& rhs) {
    return Expression<tribool>(new BinaryFormula<tribool,Less,Real,Real>(Less(),lhs,rhs));
}

//! \related Expression \brief .
inline
Expression<tribool>
operator>=(const Expression& lhs, const Expression& rhs) {
    return Expression<tribool>(new BinaryFormula<tribool,Gtr,Real,Real>(Gtr(),lhs,rhs));
}



inline bool operator==(const ContinuousPredicate& pred, bool truth) {
    ConstantFormula<tribool> const* f=dynamic_cast<ConstantFormula<tribool> const*>(&*pred.ptr);
    return (f && f->val==false);
}
*/



class ScalarFunction {
    Array<String> _space;
    Expression<Real> _expression;
  public:
    ScalarFunction(const Expression<Real>& e, const Array<String>& s) : _space(s.begin(),s.end()), _expression(convert(e,s)) { }

    uint argument_size() const { return _space.size(); }

    template<class X> X compute(const Vector<X>& x) {
        return evaluate(_expression,x);
    }
  protected:
    Expression<Real> convert(Expression<Real> e, const Array<String>& s) {
        Map<String,uint> _variables;
        for(uint i=0; i!=s.size(); ++i) { _variables.insert(s[i],i); }
    }

    template<class X> static X evaluate(const ExpressionInterface<Real>* e, const Vector<X>& x) {
        if(dynamic_cast<const ConstantExpression<Real>*>(e)) {
            Interval c=dynamic_cast<const ConstantExpression<Real>*>(e)->value();
            return x[0]*0+c;
        } else if(dynamic_cast<const VariableExpression<Real>*>(e)) {
            ARIADNE_ASSERT(false);
        } else if(dynamic_cast<const CoordinateExpression<Real>*>(e)) {
            uint j=dynamic_cast<const CoordinateExpression<Real>*>(e)->coordinate();
            return x[j];
        } else if(dynamic_cast<const UnaryExpression<Real,Neg>*>(e)) {
            const UnaryExpression<Real,Neg>* ue=dynamic_cast<const UnaryExpression<Real,Neg>*>(e);
            return evaluate(ue->_arg_ptr,x);
        } else if(dynamic_cast<const BinaryExpression<Real,Add>*>(e)) {
            const BinaryExpression<Real,Add>* be=dynamic_cast<const BinaryExpression<Real,Add>*>(e);
            return evaluate(be->_arg1_ptr,x)+evaluate(be->_arg2_ptr,x);
        } else if(dynamic_cast<const BinaryExpression<Real,Sub>*>(e)) {
            const BinaryExpression<Real,Sub>* be=dynamic_cast<const BinaryExpression<Real,Sub>*>(e);
            return evaluate(be->_arg1_ptr,x)-evaluate(be->_arg2_ptr,x);
        } else if(dynamic_cast<const BinaryExpression<Real,Mul>*>(e)) {
            const BinaryExpression<Real,Mul>* be=dynamic_cast<const BinaryExpression<Real,Mul>*>(e);
            return evaluate(be->_arg1_ptr,x)*evaluate(be->_arg2_ptr,x);
        } else if(dynamic_cast<const BinaryExpression<Real,Div>*>(e)) {
            const BinaryExpression<Real,Div>* be=dynamic_cast<const BinaryExpression<Real,Div>*>(e);
            return evaluate(be->_arg1_ptr,x)/evaluate(be->_arg2_ptr,x);
        } else {
            ARIADNE_ASSERT_MSG(false,"Unknown expression "<<e);
        }
    }

    std::ostream& write(std::ostream& os) const {
        return os << "Function( space="<<this->_space<<", expression="<<this->_expression<<" )";
    }

};



} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_H */
