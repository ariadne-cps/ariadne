/***************************************************************************
 *            expression.cc
 *
 *  Copyright 2009  Pieter Collins
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

#include "expression.h"
#include "space.h"
#include "valuation.h"
#include "vector.h"

#include "numeric.h"
#include "taylor_model.h"
#include "differential.h"

#include "polynomial.h"
#include "affine.h"

#include "real.h"
#include "formula.h"

namespace Ariadne {


template<class T> inline
std::ostream& operator<<(std::ostream& os, const ExpressionInterface<T>& e) { return e.write(os); }


//! A constant, viewed as a function \f$c:\R^n\rightarrow\R\f$.
template<class R>
class ConstantExpression
    : public ExpressionInterface<R>
{
  public:
    ConstantExpression(const std::string& nm, const R& c) : _name(nm), _c(c) { }
    ConstantExpression(const R& c) : _name(to_string(c)), _c(c) { }
    operator R () const { return _c; }
    R value() const { return _c; }
    virtual String name() const { return this->_name; }
    virtual Operator type() const { return CNST; }
    virtual Set<UntypedVariable> arguments() const { return Set<UntypedVariable>(); }
    virtual ConstantExpression<R>* clone() const { return new ConstantExpression<R>(*this); }
    virtual std::ostream& write(std::ostream& os) const { return os << _name; }
  protected:
    virtual ExpressionInterface<R>* simplify() const { return this->clone(); }
  private:
    String _name;
    R _c;
};

//! A constant, viewed as a function \f$c:\R^n\rightarrow\R\f$.
template<>
class ConstantExpression<Real>
    : public ExpressionInterface<Real>
{
  public:
    ConstantExpression(const std::string& nm, const Real& c) : _name(nm), _c(c) { }
    ConstantExpression(const Float& c) : _name(to_string(c)), _c(c) { }
    ConstantExpression(const Interval& c) :  _name(to_string(c)), _c(c) { }
    operator Real() const { return _c; }
    Real value() const { return _c; }
    virtual String name() const { return this->_name; }
    virtual String operator_name() const { return "const"; }
    virtual Operator type() const { return CNST; }
    virtual Set<UntypedVariable> arguments() const { return Set<UntypedVariable>(); }
    virtual ConstantExpression<Real>* clone() const { return new ConstantExpression<Real>(*this); }
    virtual std::ostream& write(std::ostream& os) const { return os << _name; }
  protected:
    virtual ExpressionInterface<Real>* simplify() const { return this->clone(); }
  private:
    String _name;
    Real _c;
};


//! A projection onto a named variable.
template<class R>
class VariableExpression
    : public ExpressionInterface<R>
{
  public:
    explicit VariableExpression(const String& s) : _var(s) { }
    VariableExpression(const Variable<R>& v) : _var(v) { }
    const Variable<R>& variable() const { return this->_var; }
    String name() const { return this->_var.name(); }
    virtual String operator_name() const { return "var"; }
    virtual Operator type() const { return VAR; }
    virtual VariableExpression<R>* clone() const { return new VariableExpression<R>(*this); }
    virtual Set<UntypedVariable> arguments() const { Set<UntypedVariable> r; r.insert(this->_var); return r; }
    virtual std::ostream& write(std::ostream& os) const { return os << this->name(); }
  protected:
    virtual ExpressionInterface<R>* simplify() const {
        return this->clone(); }
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
    SizeType index() const { return _j; }
    SizeType coordinate() const { return _j; }
    String name() const { String s("x0"); s[1]+=_j; return s; }
    virtual String operator_name() const { return "coord"; }
    virtual Operator type() const { return IND; }
    virtual CoordinateExpression<Real>* clone() const { return new CoordinateExpression<Real>(*this); }
    virtual Set<UntypedVariable> arguments() const { ARIADNE_FAIL_MSG("Cannot compute arguments of a CoordinateExpression"); }
    virtual std::ostream& write(std::ostream& os) const { return os << this->name(); }
  protected:
    virtual ExpressionInterface<Real>* simplify() const { return this->clone(); }
  private:
    SizeType _as;
    SizeType _j;
};



//! A composed scalar function, based on a standard operator.
template<class R, class Op=Operator, class A=R> class UnaryExpression
    : public ExpressionInterface<R>
{
  public:
    UnaryExpression(const Op& op, const ExpressionInterface<A>& expr)
        : _op(op), _arg(expr.clone()) { }
    UnaryExpression(const Op& op, const ExpressionInterface<A>* expr)
        : _op(op), _arg(const_cast<ExpressionInterface<A>*>(expr)) { }
    UnaryExpression(const Op& op, shared_ptr< const ExpressionInterface<A> > expr)
        : _op(op), _arg(expr) { }
    virtual String operator_name() const { return name(_op); }
    virtual Operator type() const { return static_cast<Operator>(_op); }
    virtual UnaryExpression<R,Op,A>* clone() const { return new UnaryExpression<R,Op,A>(_op,_arg._ptr); }
    virtual Set<UntypedVariable> arguments() const { return this->_arg.arguments(); }
    virtual std::ostream& write(std::ostream& os) const;
  protected:
    virtual ExpressionInterface<R>* simplify() const { return this->clone(); }
  public:
    Op _op;
    Expression<A> _arg;
};

//! A composed scalar function, based on a standard operator.
template<class R> class UnaryExpression<R,Operator,R>
    : public ExpressionInterface<R>
{
    typedef Operator Op; typedef R A;
  public:
    UnaryExpression(const Op& op, const ExpressionInterface<A>& expr)
        : _op(op), _arg(expr.clone()) { }
    UnaryExpression(const Op& op, const ExpressionInterface<A>* expr)
        : _op(op), _arg(const_cast<ExpressionInterface<A>*>(expr)) { }
    UnaryExpression(const Op& op, shared_ptr< const ExpressionInterface<A> > expr)
        : _op(op), _arg(expr) { }
    virtual String operator_name() const { return name(_op); }
    virtual Operator type() const { return static_cast<Operator>(_op); }
    virtual UnaryExpression<R,Op,A>* clone() const { return new UnaryExpression<R,Op,A>(_op,_arg._ptr); }
    virtual Set<UntypedVariable> arguments() const { return this->_arg.arguments(); }
    virtual std::ostream& write(std::ostream& os) const;
  protected:
    virtual ExpressionInterface<R>* simplify() const { return this->clone(); }
  public:
    Op _op;
    Expression<A> _arg;
};

template<class R, class Op, class A> inline std::ostream& UnaryExpression<R,Op,A>::write(std::ostream& os) const {
    switch(_op) {
        case NEG: return os << '-' << _arg;
        case NOT: return os << '!' << _arg;
        default: return os << _op << "(" << _arg << ")";
    }
}

template<class R> inline std::ostream& UnaryExpression<R,Operator,R>::write(std::ostream& os) const {
    switch(_op) {
        case NEG: return os << '-' << _arg;
        case NOT: return os << '!' << _arg;
        default: return os << _op << "(" << _arg << ")";
    }
}


//! A composed scalar function, based on an arthmetic operator.
template<class R, class Op=Operator, class A1=R, class A2=A1> class BinaryExpression
    : public ExpressionInterface<R>
{
  public:
    BinaryExpression(Op op, const ExpressionInterface<A1>& expr1, const ExpressionInterface<A2>& expr2)
        : _op(op), _arg1(expr1.clone()), _arg2(expr2.clone()) { }
    BinaryExpression(Op op, const ExpressionInterface<A1>* expr1, const ExpressionInterface<A2>* expr2)
        : _op(op), _arg1(expr1), _arg2(expr2) { }
    BinaryExpression(Op op, shared_ptr< const ExpressionInterface<A1> > expr1, shared_ptr< const ExpressionInterface<A2> > expr2)
        : _op(op), _arg1(expr1), _arg2(expr2)  { }
    virtual String operator_name() const { return name(_op); }
    virtual Operator type() const { return static_cast<Operator>(_op); }
    virtual BinaryExpression<R,Op,A1,A2>* clone() const { return new BinaryExpression<R,Op,A1,A2>(_op,_arg1._ptr,_arg2._ptr); }
    virtual Set<UntypedVariable> arguments() const { return join(this->_arg1.arguments(),this->_arg2.arguments()); }
    virtual std::ostream& write(std::ostream& os) const;
  protected:
    virtual ExpressionInterface<R>* simplify() const { return this->clone(); }
  public:
    Op _op;
    Expression<A1> _arg1;
    Expression<A2> _arg2;
};


template<class R, class Op, class A1, class A2> inline
std::ostream& BinaryExpression<R,Op,A1,A2>::write(std::ostream& os) const {
    switch(_op) {
        case AND: case OR:
        case ADD: case SUB: case MUL: case DIV:
        case EQ: case NEQ: case LT: case GT: case LEQ: case GEQ:
            return os << "(" << _arg1 << symbol(_op) << _arg2 << ")";
        default:
            return os << name(_op) << "(" << _arg1 << "," << _arg2 << ")";
    }
}




//! An expression in multiple variables of the same type
template<class R, class Op=Operator, class A=R> class MultiaryExpression
    : public ExpressionInterface<R>
{
  public:
    MultiaryExpression(Op op, const ExpressionInterface<A>& expr1, const ExpressionInterface<A>& expr2)
        : _op(op), _args() { _args.append(Expression<A>(expr1.clone())); _args.append(Expression<A>(expr2.clone()));  }
    MultiaryExpression(Op op, const ExpressionInterface<A>* expr1, const ExpressionInterface<A>* expr2)
        : _op(op), _args() { _args.append(Expression<A>(expr1)); _args.append(Expression<A>(expr2));  }
    MultiaryExpression(Op op, shared_ptr< const ExpressionInterface<A> > expr1, shared_ptr< const ExpressionInterface<A> > expr2)
        : _op(op), _args() { _args.append(expr1); _args.append(expr2);  }
    virtual unsigned int number_of_arguments() const { return _args.size(); }
    virtual String operator_name() const { return name(_op); }
    virtual Operator type() const { return static_cast<Operator>(_op); }
    virtual MultiaryExpression<R,Op,A>* clone() const { return new MultiaryExpression<R,Op,A>(_op,_args); }
    virtual Set<UntypedVariable> arguments() const {
        Set<UntypedVariable> res; for(uint i=0; i!=_args.size(); ++i) { res.adjoin(_args[i].arguments()); } return res; }
    virtual std::ostream& write(std::ostream& os) const {
        return os << _op << _args; }
/*
  private:
    template<class R, class A> inline
    void compute(R& r, const A& a) { r=Op()(_arg1->evaluate(a),_arg2->evaluate(a)); }
*/
  public:
    Op _op;
    List< Expression<A> > _args;
};



bool operator==(const Expression<Tribool>& e, bool v) {
    const ConstantExpression<Tribool>* expr = dynamic_cast<const ConstantExpression<Tribool>*>(e._raw_pointer());
    return expr && expr->value()==v;
}


template<class R> Expression<R>::Expression(const R& c)
    : _ptr(new ConstantExpression<R>(c))
{
}

template<class R> Expression<R>::Expression(const Constant<R>& c)
    : _ptr(new ConstantExpression<R>(c.name(),c.value()))
{
}

template<class R> Expression<R>::Expression(const Variable<R>& v)
    : _ptr(new VariableExpression<R>(v))
{
}

Expression<Real>::Expression(const double& c) : _ptr(new ConstantExpression<Real>(c)) { }
Expression<Real>::Expression(const Interval& c) : _ptr(new ConstantExpression<Real>(c)) { }
Expression<Real>::Expression(const Real& c) : _ptr(new ConstantExpression<Real>(c)) { }
Expression<Real>::Expression(const Constant<Real>& c) : _ptr(new ConstantExpression<Real>(c.name(),c.value())) { }
Expression<Real>::Expression(const Variable<Real>& v) : _ptr(new VariableExpression<Real>(v)) { }


template<class R> List< Expression<R> > Expression<R>::subexpressions() const {
    // Return the basic subexpressions used in forming the expression
    // Ideally, this code should be part of ExpressionInterface to avoid the
    // switching logic
    List< Expression<R> > res;
    const UnaryExpression<R>* uptr;
    const BinaryExpression<R>* bptr;
    uptr=dynamic_cast<const UnaryExpression<R>*>(this->_ptr.operator->());
    if(uptr) {
        res.append(uptr->_arg);
    } else {
        bptr=dynamic_cast<const BinaryExpression<R>*>(this->_ptr.operator->());
        if(bptr) {
            res.append(bptr->_arg1);
            res.append(bptr->_arg2);
        }
        else {
            std::cerr<<"Warning: subexpressions of "<<*this<<"\n";
        }
    }
    return res;
}

List< Expression<Real> > Expression<Real>::subexpressions() const {
    // Return the basic subexpressions used in forming the expression
    // Ideally, this code should be part of ExpressionInterface to avoid the
    // switching logic
    List< Expression<R> > res;
    const UnaryExpression<R>* uptr;
    const BinaryExpression<R>* bptr;
    uptr=dynamic_cast<const UnaryExpression<R>*>(this->_ptr.operator->());
    if(uptr) {
        res.append(uptr->_arg);
    } else {
        bptr=dynamic_cast<const BinaryExpression<R>*>(this->_ptr.operator->());
        if(bptr) {
            res.append(bptr->_arg1);
            res.append(bptr->_arg2);
        }
    }
    return res;
}


bool identical(const Expression<Real>& e1, const Expression<Real>& e2) {
   if(e1._ptr==e2._ptr) { return true; }
    if(e1.op()!=e2.op()) { return false; }
    const VariableExpression<Real> *vptr1, *vptr2;
    const ConstantExpression<Real> *cptr1, *cptr2;
    switch(e1.op()) {
        case VAR:
            vptr1 = dynamic_cast<const VariableExpression<Real>*>(e1._raw_pointer());
            vptr2 = dynamic_cast<const VariableExpression<Real>*>(e2._raw_pointer());
            assert(vptr1); assert(vptr2);
            return (vptr1->variable()==vptr2->variable());
        case CNST:
            cptr1 = dynamic_cast<const ConstantExpression<Real>*>(e1._raw_pointer());
            cptr2 = dynamic_cast<const ConstantExpression<Real>*>(e2._raw_pointer());
            assert(cptr1); assert(cptr2);
            return (cptr1->value()==cptr2->value());
        default:
            List< Expression<Real> > sub1=e1.subexpressions();
            List< Expression<Real> > sub2=e2.subexpressions();
            bool result=true;
            for(uint i=0; i!=sub1.size(); ++i) {
                result = result && identical(sub1[i],sub2[i]);
            }
            return result;
    }
}



template class Expression<Boolean>;
template class Expression<Tribool>;
template class Expression<String>;
template class Expression<Integer>;
template class Expression<Real>;





template<class R, class Op, class A> inline
Expression<R> make_expression(Op op, Expression<A> e) {
    return Expression<R>(new UnaryExpression<R,Op,A>(op,e._ptr)); }
template<class R, class Op, class A1, class A2> inline
Expression<R> make_expression(Op op, Expression<A1> e1, Expression<A2> e2) {
    return Expression<R>(new BinaryExpression<R,Op,A1,A2>(op,e1._ptr,e2._ptr)); }

template<class R, class Op, class A> inline
Expression<R> make_expression(Op op, ExpressionInterface<A>* eptr) {
    return Expression<R>(new UnaryExpression<R,Op,A>(op,eptr)); }
template<class R, class Op, class A1, class A2> inline
Expression<R> make_expression(Op op, ExpressionInterface<A1>* e1ptr, ExpressionInterface<A2>* e2ptr) {
    return Expression<R>(new BinaryExpression<R,Op,A1,A2>(op,e1ptr,e2ptr)); }


Expression<Boolean> operator&&(Expression<Boolean> e1, Expression<Boolean> e2) {
    return make_expression<Boolean>(AND,e1,e2); }
Expression<Boolean> operator||(Expression<Boolean> e1, Expression<Boolean> e2) {
    return make_expression<Boolean>(OR,e1,e2); }
Expression<Boolean> operator!(Expression<Boolean> e) {
    return make_expression<Boolean>(NOT,e); }


Expression<tribool> operator&&(Expression<tribool> e1, Expression<tribool> e2) {
    return make_expression<tribool>(AND,e1,e2); }
Expression<tribool> operator||(Expression<tribool> e1, Expression<tribool> e2) {
    return make_expression<tribool>(OR,e1,e2); }
Expression<tribool> operator!(Expression<tribool> e) {
    return make_expression<tribool>(NOT,e); }


Expression<Boolean> operator==(Variable<String> v1, const String& s2) {
    return make_expression<Boolean>(EQ,Expression<String>(v1),Expression<String>(s2)); }
Expression<Boolean> operator!=(Variable<String> v1, const String& s2) {
    return make_expression<Boolean>(NEQ,Expression<String>(v1),Expression<String>(s2)); }


Expression<Boolean> operator==(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(EQ,e1,e2); }
Expression<Boolean> operator!=(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(NEQ,e1,e2); }
Expression<Boolean> operator>=(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(GEQ,e1,e2); }
Expression<Boolean> operator<=(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(LEQ,e1,e2); }
Expression<Boolean> operator>(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(GT,e1,e2); }
Expression<Boolean> operator<(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(LT,e1,e2); }



Expression<Integer> operator+(Expression<Integer> e) {
    return make_expression<Integer>(POS,e); }
Expression<Integer> operator-(Expression<Integer> e) {
    return make_expression<Integer>(NEG,e); }
Expression<Integer> operator+(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Integer>(ADD,e1,e2); }
Expression<Integer> operator-(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Integer>(SUB,e1,e2); }
Expression<Integer> operator*(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Integer>(MUL,e1,e2); }



Expression<Tribool> sgn(Expression<Real> e) {
    return make_expression<Tribool>(SGN,e); }

Expression<Tribool> operator==(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(EQ,e1,e2); }
Expression<Tribool> operator!=(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(NEQ,e1,e2); }
Expression<Tribool> operator>=(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(GEQ,e1,e2); }
Expression<Tribool> operator<=(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(LEQ,e1,e2); }
Expression<Tribool> operator>(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(GT,e1,e2); }
Expression<Tribool> operator<(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(LT,e1,e2); }


Expression<Real> operator+(Expression<Real> e) {
    return make_expression<Real>(POS,e); }
Expression<Real> operator-(Expression<Real> e) {
    return make_expression<Real>(NEG,e); }
Expression<Real> operator+(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(ADD,e1,e2); }
Expression<Real> operator-(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(SUB,e1,e2); }
Expression<Real> operator*(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(MUL,e1,e2); }
Expression<Real> operator/(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(DIV,e1,e2); }

Expression<Real> pow(Expression<Real> e, int n) { 
    ARIADNE_NOT_IMPLEMENTED;
    //return make_expression(POW,e,n);
}

Expression<Real> neg(Expression<Real> e) {
    return make_expression<Real>(NEG,e); }
Expression<Real> rec(Expression<Real> e) {
    return make_expression<Real>(REC,e); }
Expression<Real> sqr(Expression<Real> e) {
    return make_expression<Real>(SQR,e); }
Expression<Real> sqrt(Expression<Real> e) {
    return make_expression<Real>(SQRT,e); }
Expression<Real> exp(Expression<Real> e) {
    return make_expression<Real>(EXP,e); }
Expression<Real> log(Expression<Real> e) {
    return make_expression<Real>(LOG,e); }
Expression<Real> sin(Expression<Real> e) {
    return make_expression<Real>(SIN,e); }
Expression<Real> cos(Expression<Real> e) {
    return make_expression<Real>(COS,e); }
Expression<Real> tan(Expression<Real> e) {
    return make_expression<Real>(TAN,e); }

Expression<Real> max(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(MAX,e1,e2); }
Expression<Real> min(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(MIN,e1,e2); }
Expression<Real> abs(Expression<Real> e) {
    return make_expression<Real>(ABS,e); }



inline void _set_constant(Float& r, const Interval& c) { r=midpoint(c); }
inline void _set_constant(Interval& r, const Interval& c) { r=c; }
inline void _set_constant(TaylorModel& r, const Interval& c) { r.clear(); r+=c; }
inline void _set_constant(Differential<Float>& r, const Interval& c) { r.clear(); r+=midpoint(c); }
inline void _set_constant(Differential<Interval>& r, const Interval& c) { r.clear(); r+=c; }

Boolean _compare(Operator cmp, const String& s1, const String& s2) {
    switch(cmp) {
        case EQ:  return s1==s2;
        case NEQ: return s1!=s2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate comparison "<<cmp<<" on string arguments.");
    }
}

Boolean _compare(Operator cmp, const Integer& z1, const Integer& z2) {
    switch(cmp) {
        case EQ:  return z1==z2;
        case NEQ: return z1!=z2;
        case LEQ: return z1<=z2;
        case GEQ: return z1>=z2;
        case LT:  return z1< z2;
        case GT:  return z1> z2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate comparison "<<cmp<<" on integer arguments.");
    }
}

template<class X> Tribool _compare(Operator cmp, const X& x1, const X& x2) {
    switch(cmp) {
        case GT: case GEQ: return x1>x2;
        case LT: case LEQ: return x1<x2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate comparison "<<cmp<<" on real arguments.");
    }
}

Boolean _compute(Operator op, const Boolean& b) {
    switch(op) {
        case NOT: return !b;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on a boolean arguments.");
    }
}

Boolean _compute(Operator op, const Boolean& b1, const Boolean& b2) {
    switch(op) {
        case AND: return b1 && b2;
        case OR: return b1 || b2;
        case XOR: return b1 ^ b2;
        case IMPL: return !b1 || b2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on a boolean arguments.");
    }
}

Tribool _compute(Operator op, const Tribool& b) {
    switch(op) {
        case NOT: return !b;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on a boolean arguments.");
    }
}

Tribool _compute(Operator op, const Tribool& b1, const Tribool& b2) {
    switch(op) {
        case AND: return b1 && b2;
        case OR: return b1 || b2;
        case XOR: return b1 ^ b2;
        case IMPL: return !b1 || b2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on a boolean arguments.");
    }
}

Integer _compute(Operator op, const Integer& x1, const Integer& x2) {
    switch(op) {
        case ADD: return x1+x2;
        case SUB: return x1-x2;
        case MUL: return x1*x2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on two integer arguments.");
    }
}

template<class X> X _compute(Operator op, const X& x1, const X& x2) {
    switch(op) {
        case ADD: return x1+x2;
        case SUB: return x1-x2;
        case MUL: return x1*x2;
        case DIV: return x1/x2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on two real arguments.");
    }
}

Integer _compute(Operator op, const Integer& z) {
    switch(op) {
        case POS: return +z;
        case NEG: return -z;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on one integer argument.");
    }
}

template<class X>
X _compute(Operator op, const X& x) {
    switch(op) {
        case NEG: return -x;
        case POS: return x;
        case REC: return 1.0/x;
        case SQR: return sqr(x);
        case SQRT: return sqrt(x);
        case EXP: return exp(x);
        case LOG: return log(x);
        case SIN: return sin(x);
        case COS: return cos(x);
        case TAN: return cos(x);
        case ABS: return abs(x);
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on one real argument.");
    }
}





Boolean evaluate(const Expression<Boolean>& e, const DiscreteValuation& x) {
    const ExpressionInterface<Boolean>* eptr=e._raw_pointer();
    const BinaryExpression<Boolean>* bptr=dynamic_cast<const BinaryExpression<Boolean,Operator>*>(eptr);
    if(bptr) { return _compute(bptr->_op,evaluate(bptr->_arg1,x),evaluate(bptr->_arg2,x)); }
    const UnaryExpression<Boolean>* uptr=dynamic_cast<const UnaryExpression<Boolean>*>(eptr);
    if(uptr) { return _compute(uptr->_op,evaluate(uptr->_arg,x)); }
    const ConstantExpression<Boolean>* cptr=dynamic_cast<const ConstantExpression<Boolean>*>(eptr);
    if(cptr) { return cptr->value(); }
    const BinaryExpression<Boolean,Operator,String>* bsptr=dynamic_cast<const BinaryExpression<Boolean,Operator,String>*>(eptr);
    if(bsptr) { return _compare(bsptr->_op,evaluate(bsptr->_arg1,x),evaluate(bsptr->_arg2,x)); }
    const BinaryExpression<Boolean,Operator,Integer>* bzptr=dynamic_cast<const BinaryExpression<Boolean,Operator,Integer>*>(eptr);
    if(bzptr) { return _compare(bsptr->_op,evaluate(bzptr->_arg1,x),evaluate(bzptr->_arg2,x)); }
    ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" to a Boolean using variables "<<x);
}

String evaluate(const Expression<String>& e, const DiscreteValuation& x) {
    const ExpressionInterface<String>* eptr=e._raw_pointer();
    const ConstantExpression<String>* cptr=dynamic_cast<const ConstantExpression<String>*>(eptr);
    if(cptr) { return cptr->value(); }
    const VariableExpression<String>* vptr=dynamic_cast<const VariableExpression<String>*>(eptr);
    if(vptr) { return x[vptr->variable()]; }
    ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" to a String using variables "<<x);
}

Integer evaluate(const Expression<Integer>& e, const DiscreteValuation& x) {
    const ExpressionInterface<Integer>* eptr=e._raw_pointer();
    const BinaryExpression<Integer>* bptr=dynamic_cast<const BinaryExpression<Integer>*>(eptr);
    if(bptr) { return _compute(bptr->_op,evaluate(bptr->_arg1,x),evaluate(bptr->_arg2,x)); }
    const UnaryExpression<Integer>* uptr=dynamic_cast<const UnaryExpression<Integer>*>(eptr);
    if(uptr) { return _compute(uptr->_op,evaluate(uptr->_arg,x)); }
    const ConstantExpression<Integer>* cptr=dynamic_cast<const ConstantExpression<Integer>*>(eptr);
    if(cptr) { return cptr->value(); }
    const VariableExpression<Integer>* vptr=dynamic_cast<const VariableExpression<Integer>*>(eptr);
    if(vptr) { return x[vptr->variable()]; }
    ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" to an Integer using variables "<<x);
}


template<class X> Tribool evaluate(const Expression<Tribool>& e, const ContinuousValuation<X>& x) {
    const ExpressionInterface<Tribool>* eptr=e._raw_pointer();
    const BinaryExpression<Tribool>* bptr=dynamic_cast<const BinaryExpression<Tribool>*>(eptr);
    if(bptr) { return _compute(bptr->_op,evaluate(bptr->_arg1,x),evaluate(bptr->_arg2,x)); }
    const BinaryExpression<Tribool,Operator,Real,Real>* brptr=dynamic_cast<const BinaryExpression<Tribool,Operator,Real,Real>*>(eptr);
    if(brptr) { return _compare(bptr->_op,evaluate(brptr->_arg1,x),evaluate(brptr->_arg2,x)); }
    ARIADNE_FAIL_MSG("");
}

template<class X> X evaluate(const Expression<Real>& e, const ContinuousValuation<X>& x) {
    const ExpressionInterface<Real>* eptr=e._raw_pointer();
    const BinaryExpression<Real>* bptr=dynamic_cast<const BinaryExpression<Real>*>(eptr);
    if(bptr) { return _compute(bptr->_op,evaluate(bptr->_arg1,x),evaluate(bptr->_arg2,x)); }
    const UnaryExpression<Real>* uptr=dynamic_cast<const UnaryExpression<Real>*>(eptr);
    if(uptr) { return _compute(uptr->_op,evaluate(uptr->_arg,x)); }
    const ConstantExpression<Real>* cptr=dynamic_cast<const ConstantExpression<Real>*>(eptr);
    if(cptr) { X r; _set_constant(r,cptr->value()); return r; }
    const VariableExpression<Real>* vptr=dynamic_cast<const VariableExpression<Real>*>(eptr);
    if(vptr) { return x[vptr->variable()]; }
    ARIADNE_FAIL_MSG("");
}

template<class X> X evaluate(const Expression<Real>& e, const Map<ExtendedRealVariable,X>& x) {
    const ExpressionInterface<Real>* eptr=e._raw_pointer();
    const BinaryExpression<Real>* bptr=dynamic_cast<const BinaryExpression<Real>*>(eptr);
    if(bptr) { return _compute(bptr->_op,evaluate(bptr->_arg1,x),evaluate(bptr->_arg2,x)); }
    const UnaryExpression<Real>* uptr=dynamic_cast<const UnaryExpression<Real>*>(eptr);
    if(uptr) { return _compute(uptr->_op,evaluate(uptr->_arg,x)); }
    const ConstantExpression<Real>* cptr=dynamic_cast<const ConstantExpression<Real>*>(eptr);
    if(cptr) { X r=x.begin()->second*0; _set_constant(r,cptr->value()); return r; }
    const VariableExpression<Real>* vptr=dynamic_cast<const VariableExpression<Real>*>(eptr);
    if(vptr) { ARIADNE_ASSERT_MSG(x.has_key(vptr->variable()),"Valuation "<<x<<" does not contain variable "<<vptr->variable()<<" used in expression "<<e); return x[vptr->variable()]; }
    ARIADNE_FAIL_MSG("");
}

template<class X> Tribool evaluate(const Expression<Tribool>& e, const Vector<X>& x) {
    const ExpressionInterface<Tribool>* eptr=e._raw_pointer();
    const BinaryExpression<Tribool>* bptr=dynamic_cast<const BinaryExpression<Tribool>*>(eptr);
    if(bptr) { return _compute(bptr->_op,evaluate(bptr->_arg1,x),evaluate(bptr->_arg2,x)); }
    const BinaryExpression<Tribool,Operator,Real,Real>* brptr=dynamic_cast<const BinaryExpression<Tribool,Operator,Real,Real>*>(eptr);
    if(brptr) { return _compare(bptr->_op,evaluate(brptr->_arg1,x),evaluate(brptr->_arg2,x)); }
    ARIADNE_FAIL_MSG("");
}

template<class X> X evaluate(const Expression<Real>& e, const Vector<X>& x) {
    const ExpressionInterface<Real>* eptr=e._raw_pointer();
    const BinaryExpression<Real>* bptr=dynamic_cast<const BinaryExpression<Real>*>(eptr);
    if(bptr) { return _compute(bptr->_op,evaluate(bptr->_arg1,x),evaluate(bptr->_arg2,x)); }
    const UnaryExpression<Real>* uptr=dynamic_cast<const UnaryExpression<Real>*>(eptr);
    if(uptr) { return _compute(uptr->_op,evaluate(uptr->_arg,x)); }
    const ConstantExpression<Real>* cptr=dynamic_cast<const ConstantExpression<Real>*>(eptr);
    if(cptr) { X r=x[0]*0; _set_constant(r,cptr->value()); return r; }
    const CoordinateExpression<Real>* iptr=dynamic_cast<const CoordinateExpression<Real>*>(eptr);
    if(iptr) { return x[iptr->index()]; }
    ARIADNE_FAIL_MSG("");
}

template Tribool evaluate(const Expression<Tribool>& e, const ContinuousValuation<Float>& x);
template Float evaluate(const Expression<Real>& e, const ContinuousValuation<Float>& x);
template Real evaluate(const Expression<Real>& e, const ContinuousValuation<Real>& x);

template Float evaluate(const Expression<Real>& e, const Map<ExtendedRealVariable,Float>& x);
template Interval evaluate(const Expression<Real>& e, const Map<ExtendedRealVariable,Interval>& x);
template Differential<Float> evaluate(const Expression<Real>& e, const Map< ExtendedRealVariable, Differential<Float> >& x);
template Differential<Interval> evaluate(const Expression<Real>& e, const Map< ExtendedRealVariable, Differential<Interval> >& x);
template TaylorModel evaluate(const Expression<Real>& e, const Map<ExtendedRealVariable,TaylorModel>& x);

template Float evaluate(const Expression<Real>& e, const Vector<Float>& x);
template Interval evaluate(const Expression<Real>& e, const Vector<Interval>& x);
template Differential<Float> evaluate(const Expression<Real>& e, const Vector< Differential<Float> >& x);
template Differential<Interval> evaluate(const Expression<Real>& e, const Vector< Differential<Interval> >& x);
template TaylorModel evaluate(const Expression<Real>& e, const Vector<TaylorModel>& x);

template<class X, class Y> Expression<X> substitute_variable(const VariableExpression<X>& e, const Variable<Y>& v, const Expression<Y>& c) {
    return Expression<X>(e.clone()); }
template<class X> Expression<X> substitute_variable(const VariableExpression<X>& e, const Variable<X>& v, const Expression<X>& c) {
    if(e.variable()==v) { return Expression<X>(c); } else { return Expression<X>(e.clone()); } }

template<class X> Expression<X> substitute_constant(const ConstantExpression<X>& e, const Constant<X>& con, const X& c) {
	if (e.name()==con.name()) { return Expression<X>(ConstantExpression<X>(con.name(),c)); } else	{ return Expression<X>(e.clone()); }
 }

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Expression<Y>& c) {
    const ExpressionInterface<X>* eptr=e._raw_pointer();
    const BinaryExpression<X,Operator,Y,Y>* aptr=dynamic_cast<const BinaryExpression<X,Operator,Y,Y>*>(eptr);
    if(aptr) { return make_expression<X>(aptr->_op,substitute(aptr->_arg1,v,c),substitute(aptr->_arg2,v,c)); }
    const BinaryExpression<X>* bptr=dynamic_cast<const BinaryExpression<X>*>(eptr);
    if(bptr) { return make_expression<X>(bptr->_op,substitute(bptr->_arg1,v,c),substitute(bptr->_arg2,v,c)); }
    const UnaryExpression<X>* uptr=dynamic_cast<const UnaryExpression<X>*>(eptr);
    if(uptr) { return make_expression<X>(uptr->_op,substitute(uptr->_arg,v,c)); }
    const ConstantExpression<X>* cptr=dynamic_cast<const ConstantExpression<X>*>(eptr);
    if(cptr) { return e; }
    const VariableExpression<X>* vptr=dynamic_cast<const VariableExpression<X>*>(eptr);
    if(vptr) { return substitute_variable(*vptr,v,c); }
    ARIADNE_FAIL_MSG("Cannot substitute for a named variable in an unknown expression.");
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Constant<Y>& con, const Y& c) {
    const ExpressionInterface<X>* eptr=e._raw_pointer();
    const BinaryExpression<X,Operator,Y,Y>* aptr=dynamic_cast<const BinaryExpression<X,Operator,Y,Y>*>(eptr);
    if(aptr) { return make_expression<X>(aptr->_op,substitute(aptr->_arg1,con,c),substitute(aptr->_arg2,con,c)); }
    const BinaryExpression<X>* bptr=dynamic_cast<const BinaryExpression<X>*>(eptr);
    if(bptr) { return make_expression<X>(bptr->_op,substitute(bptr->_arg1,con,c),substitute(bptr->_arg2,con,c)); }
    const UnaryExpression<X>* uptr=dynamic_cast<const UnaryExpression<X>*>(eptr);
    if(uptr) { return make_expression<X>(uptr->_op,substitute(uptr->_arg,con,c)); }
    const ConstantExpression<X>* cptr=dynamic_cast<const ConstantExpression<X>*>(eptr);
    if(cptr) { return substitute_constant(*cptr,con,c); }
    const VariableExpression<X>* vptr=dynamic_cast<const VariableExpression<X>*>(eptr);
    if(vptr) { return e; }
    ARIADNE_FAIL_MSG("Cannot substitute for a named constant in an unknown expression.");
}

template Expression<Tribool> substitute(const Expression<Tribool>& e, const Variable<Tribool>& v, const Expression<Tribool>& c);
template Expression<Tribool> substitute(const Expression<Tribool>& e, const Variable<Real>& v, const Expression<Real>& c);
template Expression<Real> substitute(const Expression<Real>& e, const Variable<Real>& v, const Expression<Real>& c);
template Expression<Real> substitute(const Expression<Real>& e, const Constant<Real>& con, const Real& c);

namespace {

template<class X> inline Expression<X> _simplify(const Expression<X>& e) {
    return e;
}

template<> inline Expression<Real> _simplify(const Expression<Real>& e) {
    // std::cout << "Simplifying expression "<< e << std::endl;
    //return e;
    const ExpressionInterface<Real>* eptr=e._raw_pointer();
    const UnaryExpression<Real>* uptr=dynamic_cast<const UnaryExpression<Real>*>(eptr);
    if(uptr) {  // Expression is unary
        Expression<Real> sarg=simplify(uptr->_arg);
        const ConstantExpression<Real>* carg=dynamic_cast<const ConstantExpression<Real>*>(sarg._raw_pointer());
        if(carg) {  // Argument is a constant
            try {
                return _compute<Real>(uptr->_op,carg->value());
            } catch (std::runtime_error) {  
                // Operator not supported by _compute, skip
            }                
        }
        return make_expression<Real>(uptr->_op,sarg);
    }
    const BinaryExpression<Real>* bptr=dynamic_cast<const BinaryExpression<Real>*>(eptr);
    if(!bptr) { return e; }     // Expression is neither unary nor binary, do not simplify
    // Simplify the two arguments
    Expression<Real> sarg1=simplify(bptr->_arg1);
    Expression<Real> sarg2=simplify(bptr->_arg2);
    // std::cout << "Expression is binary with arguments " << std::endl;
    // std::cout << "  " << sarg1 << std::endl;
    // std::cout << "  " << sarg2 << std::endl;
    // Check which arguments are constants
    const ConstantExpression<Real>* carg1=dynamic_cast<const ConstantExpression<Real>*>(sarg1._raw_pointer());
    const ConstantExpression<Real>* carg2=dynamic_cast<const ConstantExpression<Real>*>(sarg2._raw_pointer());
    if(carg1 && carg2) {    // If both arguments are constants, return the computed value of the expression
        // std::cout << "Both arguments are constants with value " << carg1->value() << " and " << carg2->value() << std::endl;
        try {
            return _compute<Real>(bptr->_op,carg1->value(),carg2->value());
        } catch(std::runtime_error) {
            // Operator not supported by _compute, skip
        }
    }
    Expression<Real> zero(0.0);
    Expression<Real> one(1.0);
    switch(eptr->type()) {
        case ADD:
            if(identical(sarg2,zero)) { return sarg1; }
            if(identical(sarg1,zero)) { return sarg2; }
            break;
        case SUB:
            if(identical(sarg2,zero)) { return sarg1; }
            if(identical(sarg1,zero)) { return -sarg2; }
            break;
        case MUL:
            if(identical(sarg1,zero)) { return sarg1; }
            if(identical(sarg2,zero)) { return sarg2; }
            if(identical(sarg1,one)) { return sarg2; }
            if(identical(sarg2,one)) { return sarg1; }
            break;
        case DIV:
            if(identical(sarg1,zero)) { return sarg1; }
            if(identical(sarg1,one)) { return rec(sarg2); }
            if(identical(sarg2,one)) { return sarg1; }
        default:
            break;
    }
    return make_expression<Real>(bptr->_op,sarg1,sarg2);

}


template<> inline Expression<Tribool> _simplify(const Expression<Tribool>& e) {
    const ExpressionInterface<Tribool>* eptr=e._ptr.operator->();
    const UnaryExpression<Tribool>* uptr=dynamic_cast<const UnaryExpression<Tribool>*>(eptr);
    if(uptr) {
        Expression<Tribool> sarg=simplify(uptr->_arg);
        if(uptr->_op==NOT) {
            const UnaryExpression<Tribool>* nuptr=dynamic_cast<const UnaryExpression<Tribool>*>(sarg._ptr.operator->());
            if(nuptr && nuptr->_op==NOT) {
                return nuptr->_arg;
            }
            const ConstantExpression<Tribool>* ncptr=dynamic_cast<const ConstantExpression<Tribool>*>(sarg._ptr.operator->());
            if(ncptr) {
                return Expression<Tribool>(!ncptr->value());
            }
        }
        return make_expression<Tribool>(uptr->_op,sarg);
    }
    const BinaryExpression<Tribool>* bptr=dynamic_cast<const BinaryExpression<Tribool>*>(eptr);
    if(bptr) {
        Expression<Tribool> sarg1=simplify(bptr->_arg1);
        Expression<Tribool> sarg2=simplify(bptr->_arg2);
        const ConstantExpression<Tribool>* carg1ptr=dynamic_cast<const ConstantExpression<Tribool>*>(sarg1._ptr.operator->());
        const ConstantExpression<Tribool>* carg2ptr=dynamic_cast<const ConstantExpression<Tribool>*>(sarg2._ptr.operator->());
        if(carg1ptr && carg2ptr) {
            if(bptr->_op==AND) { return Expression<Tribool>(carg1ptr->value() && carg2ptr->value()); }
            if(bptr->_op==OR) { return Expression<Tribool>(carg1ptr->value() || carg2ptr->value()); }
        } else if(carg1ptr) {
            if(bptr->_op==AND && carg1ptr->value()==true) { return sarg2; }
            if(bptr->_op==AND && carg1ptr->value()==false) { return sarg1; }
            if(bptr->_op==OR && carg1ptr->value()==true) { return sarg1; }
            if(bptr->_op==OR && carg1ptr->value()==false) { return sarg2; }
        } else if(carg2ptr) {
            if(bptr->_op==AND && carg2ptr->value()==true) { return sarg1; }
            if(bptr->_op==AND && carg2ptr->value()==false) { return sarg2; }
            if(bptr->_op==OR && carg2ptr->value()==true) { return sarg2; }
            if(bptr->_op==OR && carg2ptr->value()==false) { return sarg1; }
        } else {
            return make_expression<Tribool>(bptr->_op,sarg1,sarg2);
        }
    }
    return e;
}

} // namespace

template<class X> Expression<X> simplify(const Expression<X>& e) {
    return Ariadne::_simplify(e);
}

template Expression<Real> simplify(const Expression<Real>& e);
template Expression<Tribool> simplify(const Expression<Tribool>& e);


template<class X>
Polynomial<X> polynomial(const Expression<Real>& e, const Space<Real>& s)
{
    const ExpressionInterface<Real>* const eptr=e._raw_pointer();
    const Operator op=eptr->type();

    const ConstantExpression<Real>* cptr;
    const VariableExpression<Real>* vptr;
    const CoordinateExpression<Real>* iptr;
    const UnaryExpression<Real,Operator>* uptr;
    const BinaryExpression<Real,Operator>* bptr;

    switch(op) {
        case CNST:
            cptr=static_cast<const ConstantExpression<Real>*>(eptr);
            return Polynomial<X>::constant(s.size(),cptr->value()); break;
        case VAR:
            vptr=static_cast<const VariableExpression<Real>*>(eptr);
            return Polynomial<X>::variable(s.size(),s.index(vptr->variable())); break;
        case IND:
            iptr=static_cast<const CoordinateExpression<Real>*>(eptr);
            return Polynomial<X>::variable(s.size(),iptr->index()); break;
        case NEG:
            uptr=static_cast<const UnaryExpression<Real>*>(eptr);
            return -polynomial<X>(uptr->_arg,s); break;
        case ADD:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            return polynomial<X>(bptr->_arg1,s)+polynomial<X>(bptr->_arg2,s); break;
        case SUB:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            return polynomial<X>(bptr->_arg1,s)-polynomial<X>(bptr->_arg2,s); break;
        case MUL:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            return polynomial<X>(bptr->_arg1,s)*polynomial<X>(bptr->_arg2,s); break;
        case DIV:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            cptr=dynamic_cast<const ConstantExpression<Real>*>(bptr->_arg2._raw_pointer());
            ARIADNE_ASSERT_MSG(cptr,"Cannot convert expression "<<e<<" to polynomial form.");
            return polynomial<X>(bptr->_arg1,s)/cptr->value(); break;
        default:
            ARIADNE_FAIL_MSG("Cannot convert expression "<<e<<" to polynomial form.");
    }
}


template<class X>
Affine<X> affine(const Expression<Real>& e, const Space<Real>& s) {

    const ExpressionInterface<Real>* eptr=e._raw_pointer();
    Operator op=eptr->type();

    const ConstantExpression<Real>* cptr;
    const VariableExpression<Real>* vptr;
    const CoordinateExpression<Real>* iptr;
    const UnaryExpression<Real,Operator>* uptr;
    const BinaryExpression<Real,Operator>* bptr;
    const ConstantExpression<Real>* cptr1;
    const ConstantExpression<Real>* cptr2;

    switch(op) {
        case CNST:
            cptr=dynamic_cast<const ConstantExpression<Real>*>(eptr);
            return Affine<X>::constant(s.size(),cptr->value()); break;
        case VAR:
            vptr=dynamic_cast<const VariableExpression<Real>*>(eptr);
            return Affine<X>::variable(s.size(),s.index(vptr->variable())); break;
        case IND:
            iptr=dynamic_cast<const CoordinateExpression<Real>*>(eptr);
            return Affine<X>::variable(s.size(),iptr->coordinate()); break;
        case NEG:
            uptr=static_cast<const UnaryExpression<Real>*>(eptr);
            return -affine<X>(uptr->_arg,s); break;
        case ADD:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            return affine<X>(bptr->_arg1,s)+affine<X>(bptr->_arg2,s); break;
        case SUB:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            return affine<X>(bptr->_arg1,s)-affine<X>(bptr->_arg2,s); break;
        case DIV:
            bptr=dynamic_cast<const BinaryExpression<Real>*>(eptr);
            cptr=dynamic_cast<const ConstantExpression<Real>*>(bptr->_arg2._raw_pointer());
            ARIADNE_ASSERT_MSG(cptr,"Cannot convert expression "<<e<<" to affine form.");
            return affine<X>(bptr->_arg1,s)/static_cast<X>(cptr->value()); break;
        case MUL:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            cptr1=dynamic_cast<const ConstantExpression<Real>*>(bptr->_arg1._raw_pointer());
            cptr2=dynamic_cast<const ConstantExpression<Real>*>(bptr->_arg2._raw_pointer());
            ARIADNE_ASSERT_MSG(cptr1 || cptr2,"Cannot convert expression "<<e<<" to affine form.");
            if(cptr1) { return static_cast<X>(cptr1->value()) * affine<X>(bptr->_arg2,s); }
            else { return affine<X>(bptr->_arg1,s) * static_cast<X>(cptr2->value()); }
            break;
        default:
            ARIADNE_FAIL_MSG("Cannot convert expression "<<e<<" to affine form.");
    }
}

template Affine<Interval> affine(const Expression<Real>&, const Space<Real>&);
template Polynomial<Interval> polynomial(const Expression<Real>&, const Space<Real>&);
template Polynomial<Real> polynomial(const Expression<Real>&, const Space<Real>&);

Expression<Real> function(const Expression<Real>& e,  const Space<Real>& s)
{
    const ExpressionInterface<Real>* eptr=e._raw_pointer();
    const BinaryExpression<Real,Operator>* bptr=dynamic_cast<const BinaryExpression<Real,Operator>*>(eptr);
    if(bptr) { return make_expression<Real>(bptr->_op,function(bptr->_arg1,s),function(bptr->_arg2,s)); }
    const UnaryExpression<Real,Operator>* uptr=dynamic_cast<const UnaryExpression<Real,Operator>*>(eptr);
    if(uptr) { return make_expression<Real>(uptr->_op,function(uptr->_arg,s)); }
    const ConstantExpression<Real>* cptr=dynamic_cast<const ConstantExpression<Real>*>(eptr);
    if(cptr) { return Expression<Real>(*cptr); }
    const VariableExpression<Real>* vptr=dynamic_cast<const VariableExpression<Real>*>(eptr);
    if(vptr) { return Expression<Real>(new CoordinateExpression<Real>(s.index(vptr->variable()))); }
    const CoordinateExpression<Real>* iptr=dynamic_cast<const CoordinateExpression<Real>*>(eptr);
    if(iptr) { return e; }
    ARIADNE_FAIL_MSG("Cannot convert numbered variable");
}

Expression<Real> derivative(const Expression<Real>& e, const Variable<Real>& v)
{
    const ExpressionInterface<Real>* eptr=e._raw_pointer();
    const VariableExpression<Real>* vptr=static_cast<const VariableExpression<Real>*>(eptr);
    const UnaryExpression<Real>* uptr=static_cast<const UnaryExpression<Real>*>(eptr);
    const BinaryExpression<Real>* bptr=static_cast<const BinaryExpression<Real>*>(eptr);
    switch(eptr->type()) {
        case CNST:
            return Expression<Real>(0.0);
        case VAR:
            if(vptr->variable()==v) { return Expression<Real>(1.0); }
            else { return Expression<Real>(0.0); }
        case ADD:
            return simplify( derivative(bptr->_arg1,v)+derivative(bptr->_arg2,v) );
        case SUB:
            return simplify( derivative(bptr->_arg1,v)-derivative(bptr->_arg2,v) );
        case MUL:
            return simplify( bptr->_arg1*derivative(bptr->_arg2,v)+derivative(bptr->_arg1,v)*bptr->_arg2 );
        case DIV:
            return simplify( derivative(bptr->_arg1 * rec(bptr->_arg2),v) );
        case NEG:
            return simplify( - derivative(uptr->_arg,v) );
        case REC:
            return simplify( - derivative(uptr->_arg,v) * rec(sqr(uptr->_arg)) );
        case SQR:
            return simplify( 2 * derivative(uptr->_arg,v) * uptr->_arg );
        case EXP:
            return derivative(uptr->_arg,v) * uptr->_arg;
        case LOG:
            return derivative(uptr->_arg,v) * rec(uptr->_arg);
        case SIN:
            return derivative(uptr->_arg,v) * cos(uptr->_arg);
        case COS:
            return -derivative(uptr->_arg,v) * sin(uptr->_arg);
        case TAN:
            return derivative(uptr->_arg,v) * (1-sqr(uptr->_arg));
        default:
            ARIADNE_THROW(std::runtime_error,"derivative(RealExpression,RealVariable)",
                          "Cannot compute derivative of "<<e<<"\n");
    }
}


Expression<Real> indicator(Expression<tribool> e, Sign sign) {
    ExpressionInterface<tribool>* eptr=const_cast<ExpressionInterface<tribool>*>(e._ptr.operator->());
    ConstantExpression<Tribool>* cnptr;
    BinaryExpression<Tribool,Operator,Real,Real>* cptr;
    BinaryExpression<Tribool,Operator,Tribool,Tribool>* bptr;
    UnaryExpression<tribool,Operator,Tribool>* nptr;
    tribool value;
    switch(eptr->type()) {
        case CNST:
            cnptr=dynamic_cast<ConstantExpression<Tribool>*>(eptr);
            value=( sign==positive ? cnptr->value() : !cnptr->value() );
            if(value==true) { return RealExpression(+1.0); }
            else if(value==false) {  return RealExpression(-1.0); }
            else { return RealExpression(0.0); }
        case GEQ: case GT:
            cptr=dynamic_cast<BinaryExpression<Tribool,Operator,Real,Real>*>(eptr);
            assert(cptr);
            if(sign==positive) { return cptr->_arg1-cptr->_arg2; }
            else { return cptr->_arg2-cptr->_arg1; }
        case LEQ: case LT:
            cptr=dynamic_cast<BinaryExpression<Tribool,Operator,Real,Real>*>(eptr);
            assert(cptr);
            if(sign==positive) { return cptr->_arg2-cptr->_arg1; }
            else { return cptr->_arg1-cptr->_arg2; }
        case AND:
            bptr=dynamic_cast<BinaryExpression<Tribool,Operator,Tribool,Tribool>*>(eptr);
            assert(bptr);
            return min(indicator(bptr->_arg1,sign),indicator(bptr->_arg2,sign));
        case OR:
            bptr=dynamic_cast<BinaryExpression<tribool,Operator,Tribool,Tribool>*>(eptr);
            assert(bptr);
            return max(indicator(bptr->_arg1,sign),indicator(bptr->_arg2,sign));
        case NOT:
            nptr=dynamic_cast<UnaryExpression<tribool,Operator,Tribool>*>(eptr);
            assert(nptr);
            return neg(indicator(nptr->_arg,sign));
        default:
            ARIADNE_FAIL_MSG("Cannot compute indicator function of expression " << *eptr);
    }
}


tribool opposite(Expression<tribool> e1, Expression<tribool> e2) {
    // Simple test if two expressions are negations of each other.
    // Current

    typedef BinaryExpression<Tribool,Operator,Real,Real> ComparisonInterface;

    ExpressionInterface<tribool>* e1ptr=const_cast<ExpressionInterface<tribool>*>(e1._raw_pointer());
    ExpressionInterface<tribool>* e2ptr=const_cast<ExpressionInterface<tribool>*>(e2._raw_pointer());

    Operator e1op, e2op;
    switch(e1ptr->type()) {
        case GEQ: case GT:
            e1op=GEQ; break;
        case LEQ: case LT:
            e1op=LEQ; break;
        default:
            return indeterminate;
    }
    switch(e2ptr->type()) {
        case GEQ: case GT:
            e2op=GEQ; break;
        case LEQ: case LT:
            e2op=LEQ; break;
        default:
            return indeterminate;
    }
    std::cerr<<e1op<<" "<<e2op<<"\n";

    // Both expressions are <=,<,>=,> comparisons
    ComparisonInterface* c1ptr=dynamic_cast<ComparisonInterface*>(e1ptr);
    ComparisonInterface* c2ptr=dynamic_cast<ComparisonInterface*>(e2ptr);
    assert(c1ptr); assert(c2ptr);
    Expression<Real> e1arg1=c1ptr->_arg1;
    Expression<Real> e1arg2=c1ptr->_arg2;
    Expression<Real> e2arg1=c2ptr->_arg1;
    Expression<Real> e2arg2=c2ptr->_arg2;
    std::cerr<<e1arg1<<" "<<e1arg2<<" "<<e2arg1<<" "<<e2arg2<<"\n";

    // Test if the expressions are of the form a1<=a2; a1>=a2 or a1<=a2; a2<=a1
    if(e1op==e2op) {
        if(identical(e1arg1,e2arg2) && identical(e1arg2,e2arg1)) { return true; }
        else { return indeterminate; }
    } else {
        if(identical(e1arg1,e2arg1) && identical(e1arg2,e2arg2)) { return true; }
        else { return indeterminate; }
    }

}


} // namespace Ariadne
