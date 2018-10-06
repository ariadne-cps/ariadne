/***************************************************************************
 *            expression.cpp
 *
 *  Copyright 2009--17  Pieter Collins
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

#include "../utility/standard.hpp"
#include "../config.hpp"

#include "../numeric/operators.tpl.hpp"

#include "../algebra/algebra.hpp"
#include "../algebra/algebra_wrapper.hpp"

#include "../symbolic/constant.hpp"
#include "../symbolic/variables.hpp"
#include "../symbolic/expression.hpp"
#include "../symbolic/assignment.hpp"
#include "../symbolic/space.hpp"
#include "../symbolic/valuation.hpp"

#include "../function/formula.hpp"

namespace Ariadne {

template<class T>
struct ExpressionNode {
    mutable Nat count;
    Operator op;
    virtual ~ExpressionNode() = default;
    explicit ExpressionNode(const Operator& o) : count(0u), op(o) { }
    explicit ExpressionNode(OperatorCode cd, OperatorKind knd) : count(0u), op(cd,knd) { }
};



template<> class SymbolicAlgebraWrapper<Expression<Real>,Real>
    : public virtual SymbolicAlgebraInterface<Real>
    , public Expression<Real>
{
    typedef Expression<Real> A; typedef Real X;
  private:
    static A const& _cast(AlgebraInterface<X> const& a) { return static_cast<A const&>(dynamic_cast<SymbolicAlgebraWrapper<A,X>const&>(a)); }
    static A const& _cast(SymbolicAlgebraWrapper<A,X> const& a) { return static_cast<A const&>(a); }
    static SymbolicAlgebraInterface<X>* _make(A&& a) { return new SymbolicAlgebraWrapper<A,X>(std::move(a)); }
    template<class OP> static SymbolicAlgebraInterface<X>* _eval(OP op, SymbolicAlgebraWrapper<A,X> const& aw1, AlgebraInterface<X> const& ai2) {
        return _make(op(_cast(aw1),_cast(ai2))); }
    template<class OP> static SymbolicAlgebraInterface<X>* _eval(OP op, SymbolicAlgebraWrapper<A,X> const& aw1, X const& c2) {
        return _make(op(_cast(aw1),c2)); }
    template<class OP> static SymbolicAlgebraInterface<X>* _eval(OP op, SymbolicAlgebraWrapper<A,X> const& aw) {
        return _make(op(_cast(aw))); }
  public:
    SymbolicAlgebraWrapper(A const& a) : A(a) { }
    virtual SymbolicAlgebraInterface<X>* _create_zero() const { return new SymbolicAlgebraWrapper<A>(A()); }
    virtual SymbolicAlgebraInterface<X>* _create_constant(X const& c) const { return new SymbolicAlgebraWrapper<A>(A::constant(c)); }
    virtual SymbolicAlgebraInterface<X>* _create_copy() const { return new SymbolicAlgebraWrapper<A>(*this); }
    virtual SymbolicAlgebraInterface<X>* _neg() const {
        return _make(-_cast(*this)); }
    virtual SymbolicAlgebraInterface<X>* _add(AlgebraInterface<X> const& other) const {
        return _make(_cast(*this) + _cast(other)); }
    virtual SymbolicAlgebraInterface<X>* _sub(AlgebraInterface<X> const& other) const {
        return _make(_cast(*this) - _cast(other)); }
    virtual SymbolicAlgebraInterface<X>* _mul(AlgebraInterface<X> const& other) const {
        return _make(_cast(*this) * _cast(other)); }
    virtual SymbolicAlgebraInterface<X>* _add(const X& cnst) const { return _make(_cast(*this) + cnst); }
    virtual SymbolicAlgebraInterface<X>* _sub(X const& cnst) const { return _make(_cast(*this) - cnst); }
    virtual SymbolicAlgebraInterface<X>* _mul(X const& cnst) const { return _make(_cast(*this) * cnst); }
    virtual SymbolicAlgebraInterface<X>* _div(X const& cnst) const { return _make(_cast(*this) / cnst); }
    virtual SymbolicAlgebraInterface<X>* _radd(X const& cnst) const { return _make(cnst + _cast(*this)); }
    virtual SymbolicAlgebraInterface<X>* _rsub(X const& cnst) const { return _make(cnst - _cast(*this)); }
    virtual SymbolicAlgebraInterface<X>* _rmul(X const& cnst) const { return _make(cnst * _cast(*this)); }
    virtual SymbolicAlgebraInterface<X>* _pow(Nat m) const { return _make(pow(_cast(*this),static_cast<Int>(m))); }
    virtual Void _iadd(const X& c) { (*this) = (*this) + c; }
    virtual Void _imul(const X& c) { (*this) = (*this) * c; }
    virtual Void _isma(const X& c, const AlgebraInterface<X>& x) {
        (*this) = (*this) + c * _cast(x); }
    virtual Void _ifma(const AlgebraInterface<X>& x1, const AlgebraInterface<X>& x2)  {
        (*this) = (*this) + _cast(x1) * _cast(x2); }
    virtual AlgebraInterface<X>* _apply(Neg op) const { return _eval(Minus(),*this); }
    virtual AlgebraInterface<X>* _apply(Add op, AlgebraInterface<X>const& other) const { return _eval(Plus(),*this,other); }
    virtual AlgebraInterface<X>* _apply(Sub op, AlgebraInterface<X>const& other) const { return _eval(Minus(),*this,other); }
    virtual AlgebraInterface<X>* _apply(Mul op, AlgebraInterface<X>const& other) const { return _eval(Times(),*this,other); }
    virtual AlgebraInterface<X>* _apply(Add op, X const& cnst) const { return _eval(Plus(),*this,cnst); }
    virtual AlgebraInterface<X>* _apply(Mul op, X const& cnst) const { return _eval(Times(),*this,cnst); }
    virtual OutputStream& _write(OutputStream& os) const { return os << _cast(*this); }
    virtual SymbolicAlgebraInterface<X>* _apply(OperatorCode op);
  private:
    template<class OP> SymbolicAlgebraInterface<X>* _apply(OP op) { return new SymbolicAlgebraWrapper<A>(op(static_cast<A const&>(*this))); }
};

auto SymbolicAlgebraWrapper<Expression<Real>,Real>::_apply(OperatorCode op) -> SymbolicAlgebraInterface<X>* {
    switch(op) {
        case OperatorCode::SQRT: return this->_apply(Sqrt());
        case OperatorCode::EXP: return this->_apply(Exp());
        case OperatorCode::LOG: return this->_apply(Log());
        case OperatorCode::SIN: return this->_apply(Sin());
        case OperatorCode::COS: return this->_apply(Cos());
        case OperatorCode::TAN: return this->_apply(Tan());
        case OperatorCode::ATAN: return this->_apply(Atan());
        default: ARIADNE_FAIL_MSG("Unknown operator "<<op<<"\n");
    }
}

template<> template<>
Expression<Real>::operator Algebra<Real>() const {
    return Algebra<Real>(new SymbolicAlgebraWrapper<Expression<Real>,Real>(*this));
}


template<class T> struct ConstantExpressionNode : public ExpressionNode<T> {
    T val;
    ConstantExpressionNode(const T& v) : ExpressionNode<T>(OperatorCode::CNST,OperatorKind::NULLARY), val(v) { }
};
template<class T> struct VariableExpressionNode : public ExpressionNode<T> {
    Identifier var;
    VariableExpressionNode(const Identifier& v) : ExpressionNode<T>(OperatorCode::VAR,OperatorKind::VARIABLE), var(v) { }
};
template<class T, class A=T> struct UnaryExpressionNode : public ExpressionNode<T> {
    Expression<A> arg;
    UnaryExpressionNode(const Operator& oper, Expression<A> const& a)
        : ExpressionNode<T>(oper), arg(a) { }
};
template<class T, class A1=T, class A2=A1> struct BinaryExpressionNode : public ExpressionNode<T> {
    Expression<T> arg1; Expression<T> arg2;
};
template<class T> struct BinaryExpressionNode<T> : public ExpressionNode<T> {
    Expression<T> arg1; Expression<T> arg2;
    BinaryExpressionNode(const Operator& oper, Expression<T> const& a1, Expression<T> const& a2)
        : ExpressionNode<T>(oper), arg1(a1), arg2(a2) { }
};
template<class T> struct BinaryExpressionNode<typename Logic<T>::Type,T,T> : public ExpressionNode<typename Logic<T>::Type> {
    typedef typename Logic<T>::Type R; typedef T A;
    Expression<A> arg1; Expression<A> arg2;
    BinaryExpressionNode(const Operator& oper, Expression<A> const& a1, Expression<A> const& a2)
        : ExpressionNode<R>(oper), arg1(a1), arg2(a2) { }
};
template<class R, class A=R, class N=Int> struct ScalarExpressionNode : public UnaryExpressionNode<R,A> {
    N num;
    ScalarExpressionNode(const Operator& oper, Expression<R> const& a, N n)
        : UnaryExpressionNode<R,A>(oper,a), num(n) { }
};

template<class T> inline OutputStream& operator<<(OutputStream& os, const ExpressionNode<T>* e) {
    return os << (const Void*)(e);
}


template<class T> Expression<T>::Expression() : _root(new ConstantExpressionNode<T>(T())) { }
template<class T> Expression<T>::Expression(const T& c): _root(new ConstantExpressionNode<T>(c)) { }
template<class T> Expression<T>::Expression(const Constant<T>& c): _root(new ConstantExpressionNode<T>(c.value())) { }
template<class T> Expression<T>::Expression(const Variable<T>& v) : _root(new VariableExpressionNode<T>(v.name())) { }
template<class T> Expression<T> Expression<T>::constant(const T& c) {
    return Expression<T>(SharedPointer<const ExpressionNode<T>>(new ConstantExpressionNode<T>(c))); }
template<class T> Expression<T> Expression<T>::variable(const Identifier& v) {
    return Expression<T>(SharedPointer<const ExpressionNode<T>>(new VariableExpressionNode<T>(v))); }

template<class T> const Operator& Expression<T>::op() const {
    return node_ptr()->op; }
template<class T> OperatorCode Expression<T>::code() const {
    return node_ptr()->op.code(); }
template<class T> OperatorKind Expression<T>::kind() const {
    return node_ptr()->op.kind(); }
template<class T> const T& Expression<T>::val() const {
    return static_cast<const ConstantExpressionNode<T>*>(node_raw_ptr())->val; }
template<class T> const Identifier& Expression<T>::var() const {
    return static_cast<const VariableExpressionNode<T>*>(node_raw_ptr())->var; }
template<class T> const Expression<T>& Expression<T>::arg() const {
    return static_cast<const UnaryExpressionNode<T>*>(node_raw_ptr())->arg; }
template<class T> const Int& Expression<T>::num() const {
    return static_cast<const ScalarExpressionNode<T>*>(node_raw_ptr())->num; }
template<class T> const Expression<T>& Expression<T>::arg1() const {
    return static_cast<const BinaryExpressionNode<T>*>(node_raw_ptr())->arg1; }
template<class T> const Expression<T>& Expression<T>::arg2() const {
    return static_cast<const BinaryExpressionNode<T>*>(node_raw_ptr())->arg2; }
template<class R> template<class A> const Expression<A>& Expression<R>::cmp1(A*) const {
    return static_cast<const BinaryExpressionNode<R,A>*>(node_raw_ptr())->arg1; }
template<class R> template<class A> const Expression<A>& Expression<R>::cmp2(A*) const {
    return static_cast<const BinaryExpressionNode<R,A>*>(node_raw_ptr())->arg2; }

template<class T> Set<UntypedVariable> Expression<T>::arguments() const {
    const Expression<T>& e=*this;
    switch(e.kind()) {
        case OperatorKind::VARIABLE: return Set<UntypedVariable>{Variable<T>(e.var())};
        case OperatorKind::NULLARY: return Set<UntypedVariable>();
        case OperatorKind::UNARY: return e.arg().arguments();
        case OperatorKind::BINARY: return join(e.arg1().arguments(),e.arg2().arguments());
        case OperatorKind::COMPARISON: {
            const BinaryExpressionNode<T,Real>* rlp = dynamic_cast<const BinaryExpressionNode<T,Real>*>(e.node_raw_ptr());
            if(rlp) { return join(rlp->arg1.arguments(),rlp->arg2.arguments()); }
            const BinaryExpressionNode<T,String>* strp = dynamic_cast<const BinaryExpressionNode<T,String>*>(e.node_raw_ptr());
            if(strp) { return join(strp->arg1.arguments(),strp->arg2.arguments()); }
            [[fallthrough]];
        }
        default: ARIADNE_FAIL_MSG("Cannot compute arguments of expression "<<e<<" of kind "<<e.kind()<<"\n");
    }
}

template<class T> inline OutputStream& _write_comparison(OutputStream& os, const Expression<T>& f) {
    ARIADNE_FAIL_MSG("Comparison must return a logical type."); }
template<> inline OutputStream& _write_comparison(OutputStream& os, const Expression<Kleenean>& f) {
    Real* real_ptr=0; return os << "(" << f.cmp1(real_ptr) << symbol(f.op()) << f.cmp2(real_ptr) << ")"; }
template<> inline OutputStream& _write_comparison(OutputStream& os, const Expression<Boolean>& f) {
    String* string_ptr=0; return os << "(" << f.cmp1(string_ptr) << symbol(f.op()) << f.cmp2(string_ptr) << ")"; }
//FIXME: Distinguish String and Integer comparisons

template<class T> OutputStream& Expression<T>::_write(OutputStream& os) const {
    const Expression<T>& f=*this;
    switch(f.op()) {
        case OperatorCode::CNST:
            os << f.val(); return os;
        case OperatorCode::VAR:
            return os << f.var();
        case OperatorCode::ADD:
            return os << f.arg1() << '+' << f.arg2();
        case OperatorCode::SUB:
            os << f.arg1() << '-';
            switch(f.arg2().op()) { case OperatorCode::ADD: case OperatorCode::SUB: os << '(' << f.arg2() << ')'; break; default: os << f.arg2(); }
            return os;
        case OperatorCode::MUL:
            switch(f.arg1().op()) { case OperatorCode::ADD: case OperatorCode::SUB: case OperatorCode::DIV: os << '(' << f.arg1() << ')'; break; default: os << f.arg1(); }
            os << '*';
            switch(f.arg2().op()) { case OperatorCode::ADD: case OperatorCode::SUB: os << '(' << f.arg2() << ')'; break; default: os << f.arg2(); }
            return os;
        case OperatorCode::DIV:
            switch(f.arg1().op()) { case OperatorCode::ADD: case OperatorCode::SUB: case OperatorCode::DIV: os << '(' << f.arg1() << ')'; break; default: os << f.arg1(); }
            os << '/';
            switch(f.arg2().op()) { case OperatorCode::ADD: case OperatorCode::SUB: case OperatorCode::MUL: case OperatorCode::DIV: os << '(' << f.arg2() << ')'; break; default: os << f.arg2(); }
            return os;
        case OperatorCode::POW:
            return os << "pow" << '(' << f.arg() << ',' << f.num() << ')';
        default:
            switch(f.kind()) {
                case OperatorKind::UNARY: return os << f.op() << "(" << f.arg() << ")";
                case OperatorKind::BINARY: return os << f.op() << "(" << f.arg1() << "," << f.arg2() << ")";
                // FIXME: Type-cast comparison arguments correctly
                case OperatorKind::COMPARISON: return _write_comparison(os,f);
                default: ARIADNE_FAIL_MSG("Cannot output expression with operator "<<f.op()<<" of kind "<<f.kind()<<"\n");
            }
    }
}

template class Expression<Boolean>;
template class Expression<Kleenean>;
template class Expression<String>;
template class Expression<Integer>;
template class Expression<Real>;



template<class R> inline
Expression<R> make_expression(const Constant<R>& c) {
    return Expression<R>(std::make_shared<ConstantExpressionNode<R>>(c)); }
template<class R> inline
Expression<R> make_expression(const R& c) {
    return Expression<R>(std::make_shared<ConstantExpressionNode<R>>(c)); }
template<class R, class A> inline
Expression<R> make_expression(OperatorCode op, const Expression<A>& e) {
    return Expression<R>(std::make_shared<UnaryExpressionNode<R,A>>(op,e)); }
template<class R, class A, class N> inline
Expression<R> make_expression(OperatorCode op, const Expression<A>& e, N n) {
    return Expression<R>(std::make_shared<ScalarExpressionNode<R,A,N>>(op,e,n)); }
template<class R, class A1, class A2> inline
Expression<R> make_expression(OperatorCode op, const Expression<A1>& e1, Expression<A2> e2) {
    return Expression<R>(std::make_shared<BinaryExpressionNode<R,A1,A2>>(op,e1,e2)); }

template<class R, class Op, class A> inline
Expression<R> make_expression(Op op, const Expression<A> e) {
    return make_expression<R,A>(op.code(),e); }
template<class R, class Op, class A, class N> inline
Expression<R> make_expression(Op op, const Expression<A>& e, N n) {
    return make_expression<R>(op.code(),e,n); }
template<class R, class Op, class A1, class A2> inline
Expression<R> make_expression(Op op, const Expression<A1>& e1, Expression<A2> e2) {
    return make_expression<R>(op.code(),e1,e2); }

Expression<Boolean> operator&&(Expression<Boolean> const& e1, Expression<Boolean> const& e2) {
    return make_expression<Boolean>(AndOp(),e1,e2); }
Expression<Boolean> operator||(Expression<Boolean> const& e1, Expression<Boolean> const& e2) {
    return make_expression<Boolean>(OrOp(),e1,e2); }
Expression<Boolean> operator!(Expression<Boolean> const& e) {
    return make_expression<Boolean>(NotOp(),e); }


Expression<Kleenean> operator&&(Expression<Kleenean> const& e1, Expression<Kleenean> const& e2) {
    return make_expression<Kleenean>(AndOp(),e1,e2); }
Expression<Kleenean> operator||(Expression<Kleenean> const& e1, Expression<Kleenean> const& e2) {
    return make_expression<Kleenean>(OrOp(),e1,e2); }
Expression<Kleenean> operator!(Expression<Kleenean> const& e) {
    return make_expression<Kleenean>(NotOp(),e); }


Expression<Boolean> operator==(Variable<String> v1, const String& s2) {
    return make_expression<Boolean>(Equal(),Expression<String>(v1),Expression<String>::constant(s2)); }
Expression<Boolean> operator!=(Variable<String> v1, const String& s2) {
    return make_expression<Boolean>(Unequal(),Expression<String>(v1),Expression<String>::constant(s2)); }


Expression<Boolean> operator==(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Boolean>(Equal(),e1,e2); }
Expression<Boolean> operator!=(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Boolean>(Unequal(),e1,e2); }
Expression<Boolean> operator>=(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Boolean>(Geq(),e1,e2); }
Expression<Boolean> operator<=(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Boolean>(Leq(),e1,e2); }
Expression<Boolean> operator>(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Boolean>(Gtr(),e1,e2); }
Expression<Boolean> operator<(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Boolean>(Less(),e1,e2); }



Expression<Integer> operator+(Expression<Integer> const& e) {
    return make_expression<Integer>(Pos(),e); }
Expression<Integer> operator-(Expression<Integer> const& e) {
    return make_expression<Integer>(Neg(),e); }
Expression<Integer> operator+(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Integer>(Add(),e1,e2); }
Expression<Integer> operator-(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Integer>(Sub(),e1,e2); }
Expression<Integer> operator*(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Integer>(Mul(),e1,e2); }

Expression<Integer>& operator+=(Expression<Integer>& e1, Expression<Integer> const& e2) {
    return e1=e1+e2; }
Expression<Integer>& operator-=(Expression<Integer>& e1, Expression<Integer> const& e2) {
    return e1=e1-e2; }
Expression<Integer>& operator*=(Expression<Integer>& e1, Expression<Integer> const& e2) {
    return e1=e1*e2; }

Expression<Kleenean> sgn(Expression<Real> const& e) {
    return make_expression<Kleenean>(Sgn(),e); }

Expression<Kleenean> operator==(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Kleenean>(Equal(),e1,e2); }
Expression<Kleenean> operator!=(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Kleenean>(Unequal(),e1,e2); }
Expression<Kleenean> operator>=(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Kleenean>(Geq(),e1,e2); }
Expression<Kleenean> operator<=(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Kleenean>(Leq(),e1,e2); }
Expression<Kleenean> operator>(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Kleenean>(Gtr(),e1,e2); }
Expression<Kleenean> operator<(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Kleenean>(Less(),e1,e2); }


Expression<Real> operator+(Expression<Real> const& e) {
    return make_expression<Real>(Pos(),e); }
Expression<Real> operator-(Expression<Real> const& e) {
    return make_expression<Real>(Neg(),e); }
Expression<Real> operator+(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Add(),e1,e2); }
Expression<Real> operator-(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Sub(),e1,e2); }
Expression<Real> operator*(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Mul(),e1,e2); }
Expression<Real> operator/(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Div(),e1,e2); }

Expression<Real>& operator+=(Expression<Real>& e1, Expression<Real> const& e2) {
    return e1=e1+e2; }
Expression<Real>& operator-=(Expression<Real>& e1, Expression<Real> const& e2) {
    return e1=e1-e2; }
Expression<Real>& operator*=(Expression<Real>& e1, Expression<Real> const& e2) {
    return e1=e1*e2; }
Expression<Real>& operator/=(Expression<Real>& e1, Expression<Real> const& e2) {
    return e1=e1/e2; }

Expression<Real> pow(Expression<Real> const& e, Int n) {
    return make_expression<Real>(Pow(),e,n); }

Expression<Real> neg(Expression<Real> const& e) {
    return make_expression<Real>(Neg(),e); }
Expression<Real> rec(Expression<Real> const& e) {
    return make_expression<Real>(Rec(),e); }
Expression<Real> sqr(Expression<Real> const& e) {
    return make_expression<Real>(Sqr(),e); }
Expression<Real> sqrt(Expression<Real> const& e) {
    return make_expression<Real>(Sqrt(),e); }
Expression<Real> exp(Expression<Real> const& e) {
    return make_expression<Real>(Exp(),e); }
Expression<Real> log(Expression<Real> const& e) {
    return make_expression<Real>(Log(),e); }
Expression<Real> sin(Expression<Real> const& e) {
    return make_expression<Real>(Sin(),e); }
Expression<Real> cos(Expression<Real> const& e) {
    return make_expression<Real>(Cos(),e); }
Expression<Real> tan(Expression<Real> const& e) {
    return make_expression<Real>(Tan(),e); }
Expression<Real> atan(Expression<Real> const& e) {
    return make_expression<Real>(Atan(),e); }

Expression<Real> max(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Max(),e1,e2); }
Expression<Real> min(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Min(),e1,e2); }
Expression<Real> abs(Expression<Real> const& e) {
    return make_expression<Real>(Abs(),e); }


Expression<Real> operator+(Expression<Real> const& e, Real const& c) { return e + Expression<Real>::constant(c); }
Expression<Real> operator-(Expression<Real> const& e, Real const& c) { return e - Expression<Real>::constant(c); }
Expression<Real> operator*(Expression<Real> const& e, Real const& c) { return e * Expression<Real>::constant(c); }
Expression<Real> operator/(Expression<Real> const& e, Real const& c) { return e / Expression<Real>::constant(c); }
Expression<Real> operator+(Real const& c, Expression<Real> const& e) { return Expression<Real>::constant(c) + e; }
Expression<Real> operator-(Real const& c, Expression<Real> const& e) { return Expression<Real>::constant(c) - e; }
Expression<Real> operator*(Real const& c, Expression<Real> const& e) { return Expression<Real>::constant(c) * e; }
Expression<Real> operator/(Real const& c, Expression<Real> const& e) { return Expression<Real>::constant(c) / e; }

Expression<Kleenean> operator<=(Expression<Real> const& e, Real const& c) { return e <= Expression<Real>::constant(c); }
Expression<Kleenean> operator< (Expression<Real> const& e, Real const& c) { return e <  Expression<Real>::constant(c); }
Expression<Kleenean> operator>=(Expression<Real> const& e, Real const& c) { return e >= Expression<Real>::constant(c); }
Expression<Kleenean> operator> (Expression<Real> const& e, Real const& c) { return e >  Expression<Real>::constant(c); }
Expression<Kleenean> operator<=(Real const& c, Expression<Real> const& e) { return Expression<Real>::constant(c) <= e; }
Expression<Kleenean> operator< (Real const& c, Expression<Real> const& e) { return Expression<Real>::constant(c) <  e; }
Expression<Kleenean> operator>=(Real const& c, Expression<Real> const& e) { return Expression<Real>::constant(c) >= e; }
Expression<Kleenean> operator> (Real const& c, Expression<Real> const& e) { return Expression<Real>::constant(c) >  e; }




template<class A>
typename Logic<A>::Type
evaluate(const Expression<typename Logic<A>::Type>& e, const Map<Identifier,A>& x)
{
    typedef typename Logic<A>::Type R;
    A* aptr=0;
    switch(e.kind()) {
        case OperatorKind::NULLARY: return static_cast<R>(e.val());
        case OperatorKind::UNARY: return compute(e.op(),evaluate(e.arg(),x));
        case OperatorKind::BINARY: return compute(e.op(),evaluate(e.arg1(),x),evaluate(e.arg2(),x));
        case OperatorKind::COMPARISON: return compare(e.op(),evaluate(e.cmp1(aptr),x),evaluate(e.cmp2(aptr),x));
        default: ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" on "<<x<<"\n");
    }
}

template<class T>
T
evaluate(const Expression<T>& e, const Map<Identifier,T>& x)
{
    switch(e.kind()) {
        case OperatorKind::VARIABLE: return x[e.var()];
        case OperatorKind::NULLARY: return e.val();
        case OperatorKind::UNARY: return compute(e.op(),evaluate(e.arg(),x));
        case OperatorKind::BINARY: return compute(e.op(),evaluate(e.arg1(),x),evaluate(e.arg2(),x));
        case OperatorKind::GRADED: return compute(e.op(),evaluate(e.arg(),x),e.num());
        default: ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" on "<<x<<"\n");
    }
}

template Real evaluate(const Expression<Real>& e, const Map<Identifier,Real>& x);

String evaluate(const Expression<String>& e, const StringValuation& x) {
    switch(e.op()) {
        case OperatorCode::CNST: return e.val();
        case OperatorCode::VAR: return x[e.var()];
        default: ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" to a String using variables "<<x);
    }
}

Integer evaluate(const Expression<Integer>& e, const Valuation<Integer>& x) {
    return evaluate(e,x.values());
}

Boolean evaluate(const Expression<Boolean>& e, const Valuation<String>& x) {
    return evaluate(e,x.values());
}

Kleenean evaluate(const Expression<Kleenean>& e, const Valuation<Real>& x) {
    return evaluate(e,x.values());
}

Real evaluate(const Expression<Real>& e, const Valuation<Real>& x) {
    return evaluate(e,x.values());
}



template<class T> Set<Identifier> arguments(const Expression<T>& e)
{
    switch(e.kind()) {
        case OperatorKind::VARIABLE: return Set<Identifier>{e.var()};
        case OperatorKind::NULLARY: return Set<Identifier>();
        case OperatorKind::UNARY: return arguments(e.arg());
        case OperatorKind::BINARY: return join(arguments(e.arg1()),arguments(e.arg2()));
        case OperatorKind::COMPARISON: {
            const BinaryExpressionNode<T,Real>* rlp = dynamic_cast<const BinaryExpressionNode<T,Real>*>(e.node_raw_ptr());
            if(rlp) { return join(arguments(rlp->arg1),arguments(rlp->arg2)); }
            const BinaryExpressionNode<T,String>* strp = dynamic_cast<const BinaryExpressionNode<T,String>*>(e.node_raw_ptr());
            if(strp) { return join(arguments(strp->arg1),arguments(strp->arg2)); }
            [[fallthrough]];
        }
        default: ARIADNE_FAIL_MSG("Cannot compute arguments of expression "<<e<<" of kind "<<e.kind()<<"\n");
    }
}

template Set<Identifier> arguments(const Expression<Boolean>& e);
template Set<Identifier> arguments(const Expression<Kleenean>& e);
template Set<Identifier> arguments(const Expression<Real>& e);






namespace {
template<class I, class X, class Y> inline const Expression<X>& _substitute_variable(const I& ie, const I& is, const Expression<X>& e, const Expression<Y>& s) {
    ARIADNE_ASSERT_MSG(ie!=is,"Cannot substitute expression "<<s<<" for variable "<<ie<<"\n");
    return e; }
template<class I, class X> inline const Expression<X>& _substitute_variable(const I& ie, const I& is, const Expression<X>& e, const Expression<X>& s) {
    return ie==is ? s : e; }
} // namespace

template<class X, class I, class Y> Expression<X> substitute(const Expression<X>& e, const I& v, const Expression<Y>& s) {
    switch(e.kind()) {
        case OperatorKind::COMPARISON: {
            Y* yptr=0;
            const Expression<Y>& c1=e.cmp1(yptr);
            const Expression<Y>& c2=e.cmp2(yptr);
            return make_expression<X>(e.op(),substitute(c1,v,s),substitute(c2,v,s)); }
        case OperatorKind::BINARY: return make_expression<X>(e.op(),substitute(e.arg1(),v,s),substitute(e.arg2(),v,s));
        case OperatorKind::UNARY: return make_expression<X>(e.op(),substitute(e.arg(),v,s));
        case OperatorKind::NULLARY: return make_expression<X>(e.val());
        case OperatorKind::VARIABLE: return _substitute_variable(e.var(),v,e,s);
        case OperatorKind::GRADED: return make_expression<X>(e.op(),substitute(e.arg(),v,s),e.num());
        default: ARIADNE_FAIL_MSG("Cannot substitute "<<s<<" for a named variable "<<v<<" in an unknown expression "<<e<<"\n");
    }
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Expression<Y>& s) {
    return substitute(e,v.name(),s);
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Y& c) {
    return substitute(e,v.name(),Expression<Y>::constant(c));
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Assignment<Variable<Y>,Expression<Y>>& a) {
    return substitute(e,a.lhs,a.rhs);
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const List<Assignment<Variable<Y>,Expression<Y>>>& a) {
    Expression<X> r=e;
    for(SizeType i=0; i!=a.size(); ++i) {
        r=substitute(r,a[i]);
    }
    return r;
}

template<class X, class Y> SizeType substitute(const Vector<Expression<X>>& e, const List< Assignment< Variable<Y>, Expression<Y> > >& a) {
    Vector<Expression<X>> r(e.size());
    for(SizeType i=0; i!=e.size(); ++i) {
        r[i]=substitute(e[i],a);
    }
    return r;
}

template Expression<Kleenean> substitute(const Expression<Kleenean>& e, const Variable<Kleenean>& v, const Kleenean& c);
template Expression<Kleenean> substitute(const Expression<Kleenean>& e, const Variable<Real>& v, const Real& c);
template Expression<Real> substitute(const Expression<Real>& e, const Variable<Real>& v, const Real& c);
template Expression<Real> substitute(const Expression<Real>& e, const Variable<Real>& v, const Expression<Real>& c);
template Expression<Real> substitute(const Expression<Real>& e, const List< Assignment< Variable<Real>, Expression<Real> > >& c);
template Expression<Kleenean> substitute(const Expression<Kleenean>& e, const List< Assignment< Variable<Real>, Expression<Real> > >& c);


namespace {

template<class X> inline Expression<X> _simplify(const Expression<X>& e) {
    return e;
}


inline Expression<Real> _simplify(const Expression<Real>& e) {
    typedef Real R;

    if(e.kind() == OperatorKind::UNARY) {
        Expression<R> sarg=simplify(e.arg());
        if(sarg.op()==OperatorCode::CNST) {
            return Expression<R>(compute(e.op(),sarg.val()));
        } else {
            return make_expression<R>(e.op(),sarg);
        }
    }

    if(e.kind() == OperatorKind::GRADED) {
        Expression<R> sarg=simplify(e.arg());
        Expression<R> one(static_cast<R>(1));
        switch(e.op()) {
            case OperatorCode::POW:
                switch (e.num()) {
                case 0: return one;
                case 1: return sarg;
                default: return make_expression<R>(OperatorCode::POW,sarg,e.num());
                }
            default:
                return make_expression<R>(e.op(),sarg,e.num());
        }
    }

    if(e.kind() != OperatorKind::BINARY) { return e; }

    Expression<R> sarg1=simplify(e.arg1());
    Expression<R> sarg2=simplify(e.arg2());
    Expression<R> zero(static_cast<R>(0));
    Expression<R> one(static_cast<R>(1));
    switch(e.op()) {
        case OperatorCode::ADD:
            if(identical(sarg2,zero)) { return sarg1; }
            if(identical(sarg1,zero)) { return sarg2; }
            break;
        case OperatorCode::SUB:
            if(identical(sarg2,zero)) { return sarg1; }
            if(identical(sarg1,zero)) { return -sarg2; }
            break;
        case OperatorCode::MUL:
            if(identical(sarg1,zero)) { return zero; }
            if(identical(sarg2,zero)) { return zero; }
            if(identical(sarg1,one)) { return sarg2; }
            if(identical(sarg2,one)) { return sarg1; }
            break;
        case OperatorCode::DIV:
            if(identical(sarg1,zero)) { return sarg1; }
            if(identical(sarg1,one)) { return rec(sarg2); }
            if(identical(sarg2,one)) { return sarg1; }
        default:
            break;
    }
    return make_expression<R>(e.op(),sarg1,sarg2);
}

template<class I> inline Expression<Kleenean> _simplify(const Expression<Kleenean>& e) {
    typedef Kleenean T;

    if( e.kind()==OperatorKind::UNARY ) {
        Expression<T> sarg=simplify(e.arg());
        if(e.op()==OperatorCode::NOT) {
            if( sarg.op()==OperatorCode::NOT ) {
                return sarg.arg();
            }
            if( sarg.op()==OperatorCode::CNST ) {
                return Expression<T>(compute(e.op(),sarg.val()));
            }
        }
        return make_expression<T>(e.op(),sarg);
    }

    if( e.kind()==OperatorKind::BINARY ) {
        Expression<T> sarg1=simplify(e.arg1());
        Expression<T> sarg2=simplify(e.arg2());
        if( sarg1.op()==OperatorCode::CNST && sarg2.op()==OperatorCode::CNST ) {
            if(e.op()==OperatorCode::AND) { return Expression<T>(sarg1.val() && sarg2.val()); }
            if(e.op()==OperatorCode::OR) { return Expression<T>(sarg1.val() || sarg2.val()); }
            return Expression<T>(compute(e.op(),sarg1.val(),sarg2.val()));
        } else if(sarg1.op()==OperatorCode::CNST) {
            if(e.op()==OperatorCode::AND && definitely(sarg1.val())) { return sarg2; }
            if(e.op()==OperatorCode::AND && definitely(not sarg1.val())) { return sarg1; }
            if(e.op()==OperatorCode::OR && definitely(sarg1.val())) { return sarg1; }
            if(e.op()==OperatorCode::OR && definitely(not sarg1.val())) { return sarg2; }
        } else if(sarg2.op()==OperatorCode::CNST) {
            if(e.op()==OperatorCode::AND && definitely(sarg2.val())) { return sarg1; }
            if(e.op()==OperatorCode::AND && definitely(not sarg2.val())) { return sarg2; }
            if(e.op()==OperatorCode::OR && definitely(sarg2.val())) { return sarg2; }
            if(e.op()==OperatorCode::OR && definitely(not sarg2.val())) { return sarg1; }
        } else {
            return make_expression<T>(e.op(),sarg1,sarg2);
        }
    }
    return e;
}


} // namespace

template<class X> Expression<X> simplify(const Expression<X>& e) {
    return Ariadne::_simplify(e);
}

template Expression<Real> simplify(const Expression<Real>& e);
template Expression<Kleenean> simplify(const Expression<Kleenean>& e);



Expression<Real> indicator(Expression<Kleenean> e, Sign sign) {
    switch(e.op()) {
        case OperatorCode::CNST: {
            Kleenean value=( sign==Sign::POSITIVE ? e.val() : !e.val() );
            ValidatedKleenean checked_value = value.check(Effort::get_default());
            if(definitely(checked_value)) { return Expression<Real>::constant(+1); }
            else if(not possibly(checked_value)) {  return Expression<Real>::constant(-1); }
            else { return Expression<Real>::constant(0); }
        }
        case OperatorCode::VAR:
            return Expression<Real>(Variable<Real>(e.var()));
        case OperatorCode::GEQ: case OperatorCode::GT:
            if(sign==Sign::POSITIVE) { return e.cmp1<Real>()-e.cmp2<Real>(); }
            else { return e.cmp2<Real>()-e.cmp1<Real>(); }
        case OperatorCode::LEQ: case OperatorCode::LT:
            if(sign==Sign::POSITIVE) { return e.cmp2<Real>()-e.cmp1<Real>(); }
            else { return e.cmp1<Real>()-e.cmp2<Real>(); }
        case OperatorCode::AND:
            return min(indicator(e.arg1(),sign),indicator(e.arg2(),sign));
        case OperatorCode::OR:
            return max(indicator(e.arg1(),sign),indicator(e.arg2(),sign));
        case OperatorCode::NOT:
            return neg(indicator(e.arg(),sign));
        default:
            ARIADNE_FAIL_MSG("Cannot compute indicator function of expression " << e);
    }
}


template<class T> Bool is_constant(const Expression<T>& e, const typename Expression<T>::ValueType& c) {
    switch(e.op()) {
        case OperatorCode::CNST: return decide(e.val()==c);
        default: return false;
    }
}

template Bool is_constant(const Expression<Real>&, const Real&);
template Bool is_constant(const Expression<Kleenean>&, const Kleenean&);


template<class T> Bool is_variable(const Expression<T>& e, const Variable<T>& v) {
    switch(e.op()) {
        case OperatorCode::VAR: return e.var()==v.name();
        default: return false;
    }
}

template Bool is_variable(const Expression<Real>&, const Variable<Real>&);

Bool is_constant_in(const Expression<Real>& e, const Set<Variable<Real>>& vs) {
    switch(e.kind()) {
        case OperatorKind::VARIABLE: return not vs.contains(Variable<Real>(e.var()));
        case OperatorKind::NULLARY: return true;
        case OperatorKind::UNARY: case OperatorKind::SCALAR: case OperatorKind::GRADED: return is_constant_in(e.arg(),vs);
        case OperatorKind::BINARY: return is_constant_in(e.arg1(),vs) and is_constant_in(e.arg2(),vs);
        default: ARIADNE_FAIL_MSG("Cannot evaluate if expression "<<e<<" is constant on "<<vs<<"\n");
    }
}

Bool is_affine_in(const Expression<Real>& e, const Set<Variable<Real>>& vs) {
    switch(e.op()) {
        case OperatorCode::CNST: return true;
        case OperatorCode::VAR: return true;
        case OperatorCode::ADD: case OperatorCode::SUB: return is_affine_in(e.arg1(),vs) and is_affine_in(e.arg2(),vs);
        case OperatorCode::MUL: return (is_affine_in(e.arg1(),vs) and is_constant_in(e.arg2(),vs)) or (is_constant_in(e.arg1(),vs) and is_affine_in(e.arg2(),vs));
        case OperatorCode::DIV: return (is_affine_in(e.arg1(),vs) and is_constant_in(e.arg2(),vs));
        case OperatorCode::POS: case OperatorCode::NEG: return is_affine_in(e.arg(),vs);
        case OperatorCode::POW: case OperatorCode::SQR: case OperatorCode::COS: case OperatorCode::SIN: case OperatorCode::TAN: return is_constant_in(e.arg(),vs);
        default: ARIADNE_FAIL_MSG("Not currently supporting code '"<<e.op()<<"' for evaluation of affinity in given variables\n");
    }
}

Bool is_affine_in(const Vector<Expression<Real>>& e, const Set<Variable<Real>>& vs) {
    for (auto i : range(e.size()))
        if (not is_affine_in(e[i],vs)) return false;
    return true;
}

Bool is_additive_in(const Vector<Expression<Real>>& ev, const Set<Variable<Real>>& vs) {
    // We treat the vector of expressions as additive in vs if each variable in vs appears at most once in all expressions,
    // with a constant value of 1
    // (FIXME: this simplifies the case of a constant multiplier, for which would need to rescale the variable)
    // (FIXME: more generally, this simplifies the case of a diagonalisable matrix of constant multipliers)

    auto one = Expression<Real>::constant(1);
    auto zero = Expression<Real>::constant(0);

    for (auto v : vs) {
        Bool already_found = false;
        Bool already_found_one = false;
        for (auto i : range(ev.size())) {
            const Expression<Real>& e = ev[i];
            auto der = simplify(derivative(e, v));
            if (not identical(der,zero)) {
                if (already_found) {
                    return false;
                } else {
                    already_found = true;
                    if (identical(der,one)) {
                        if (already_found_one) {
                            return false;
                        } else {
                            already_found_one = true;
                        }
                    } else {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}


template<class R> Bool identical(const Expression<R>& e1, const Expression<R>& e2)
{
    if(e1.node_raw_ptr()==e2.node_raw_ptr()) { return true; }
    if(e1.op()!=e2.op()) { return false; }
    switch(e1.kind()) {
        case OperatorKind::VARIABLE:
            return e1.var()==e2.var();
        case OperatorKind::NULLARY:
            return same(e1.val(),e2.val());
        case OperatorKind::UNARY:
            return identical(e1.arg(),e2.arg());
        case OperatorKind::GRADED:
            return identical(e1.arg(),e2.arg()) && e1.num() == e2.num();
        case OperatorKind::BINARY:
            switch(e1.op()) {
            case OperatorCode::MUL: case OperatorCode::ADD:
                return (identical(e1.arg1(),e2.arg1()) && identical(e1.arg2(),e2.arg2())) ||
                       (identical(e1.arg1(),e2.arg2()) && identical(e1.arg2(),e2.arg1()));
            default:
                return identical(e1.arg1(),e2.arg1()) && identical(e1.arg2(),e2.arg2());
            }
        default:
            return false;
    }
}

template Bool identical(const Expression<Real>&, const Expression<Real>&);



Bool opposite(Expression<Kleenean> e1, Expression<Kleenean> e2) {

    OperatorCode e1op;
    OperatorCode e2op;
    switch(e1.op()) {
        case OperatorCode::GEQ: case OperatorCode::GT: e1op=OperatorCode::GEQ; break;
        case OperatorCode::LEQ: case OperatorCode::LT: e1op=OperatorCode::LEQ; break;
        default: return false;
    }
    switch(e2.op()) {
        case OperatorCode::GEQ: case OperatorCode::GT: e2op=OperatorCode::GEQ; break;
        case OperatorCode::LEQ: case OperatorCode::LT: e2op=OperatorCode::LEQ; break;
        default: return false;
    }

    // Both expressions are <=,<,>=,> comparisons
    Expression<Real> const& e1arg1=e1.cmp1<Real>();
    Expression<Real> const& e1arg2=e2.cmp2<Real>();
    Expression<Real> const& e2arg1=e1.cmp1<Real>();
    Expression<Real> const& e2arg2=e2.cmp2<Real>();

    // Test if the expressions are of the form a1<=a2; a1>=a2 or a1<=a2; a2<=a1
    if(e1op==e2op) {
        if(identical(e1arg1,e2arg2) && identical(e1arg2,e2arg1)) { return true; }
        else { return false; }
    } else {
        if(identical(e1arg1,e2arg1) && identical(e1arg2,e2arg2)) { return true; }
        else { return false; }
    }

}


Expression<Real> derivative(const Expression<Real>& e, Variable<Real> v)
{
    /*
    switch(e.op()) {
        case OperatorCode::CNST:
            return Expression<Real>::constant(0);
        case OperatorCode::VAR:
            if(e.var()==v.name()) { return Expression<Real>::constant(1); }
            else { return Expression<Real>::constant(0); }
        case OperatorCode::ADD:
            return derivative(e.arg1(),v)+derivative(e.arg2(),v);
        case OperatorCode::SUB:
            return derivative(e.arg1(),v)-derivative(e.arg2(),v);
        case OperatorCode::MUL:
            return e.arg1()*derivative(e.arg2(),v)+derivative(e.arg1(),v)*e.arg2();
        case OperatorCode::DIV:
            return derivative(e.arg1() * rec(e.arg2()),v);
        case OperatorCode::NEG:
            return  - derivative(e.arg(),v);
        case OperatorCode::REC:
            return  - derivative(e.arg(),v) * rec(sqr(e.arg()));
        case OperatorCode::SQR:
            return static_cast<Real>(2) * derivative(e.arg(),v) * e.arg();
        case OperatorCode::POW:
            return Expression<Real>::constant(e.num())*make_expression<Real>(OperatorCode::POW,e.arg(),e.num()-1)*derivative(e.arg(),v);
        case OperatorCode::EXP:
            return derivative(e.arg(),v) * e.arg();
        case OperatorCode::LOG:
            return derivative(e.arg(),v) * rec(e.arg());
        case OperatorCode::SIN:
            return derivative(e.arg(),v) * cos(e.arg());
        case OperatorCode::COS:
            return -derivative(e.arg(),v) * sin(e.arg());
        case OperatorCode::TAN:
            return derivative(e.arg(),v) * (static_cast<Real>(1)-sqr(e.arg()));*/
    switch(e.kind()) {
        case OperatorKind::NULLARY: return Expression<Real>::constant(0);
        case OperatorKind::VARIABLE: return Expression<Real>::constant(e.var()==v.name()?1:0);
        case OperatorKind::UNARY: return derivative(e.op(),e.arg(),derivative(e.arg(),v));
        case OperatorKind::BINARY: return derivative(e.op(),e.arg1(),derivative(e.arg1(),v),e.arg2(),derivative(e.arg2(),v));
        case OperatorKind::GRADED: return derivative(e.op(), e.arg(), derivative(e.arg(),v), e.num());
        default:
            ARIADNE_THROW(std::runtime_error,"derivative(Expression<Real> e, Variable<Real> v)",
                          "Cannot compute derivative of "<<e<<"\n");
    }
}



const Formula<EffectiveNumber>& cached_make_formula(const Expression<Real>& e, const Map<Identifier,Nat>& v, Map< const Void*, Formula<EffectiveNumber> >& cache)
{
    typedef EffectiveNumber Y;
    const ExpressionNode<Real>* eptr=e.node_ptr().operator->();
    if(cache.has_key(eptr)) { return cache.get(eptr); }
    switch(e.kind()) {
        case OperatorKind::VARIABLE: return insert( cache, eptr, make_formula<Y>(v[e.var()]) );
        case OperatorKind::NULLARY: return insert( cache, eptr, make_formula<Y>(e.val()) );
        case OperatorKind::UNARY: return insert( cache, eptr, make_formula<Y>(e.op(),cached_make_formula(e.arg(),v,cache)));
        case OperatorKind::BINARY: return insert( cache, eptr, make_formula<Y>(e.op(),cached_make_formula(e.arg1(),v,cache),cached_make_formula(e.arg2(),v,cache)) );
        case OperatorKind::GRADED: return insert( cache, eptr, make_formula<Y>(e.op(),cached_make_formula(e.arg(),v,cache),e.num()) );
        default: ARIADNE_FAIL_MSG("Cannot convert expression "<<e<<" to use variables "<<v<<"\n");
    }
}

Formula<EffectiveNumber> make_formula(const Expression<Real>& e, const Map<Identifier,Nat>& v)
{
    Map< const Void*, Formula<EffectiveNumber> > cache;
    return cached_make_formula(e,v,cache);
}

Formula<EffectiveNumber> make_formula(const Expression<Real>& e, const Space<Real>& spc)
{
    typedef EffectiveNumber Y;
    switch(e.kind()) {
        case OperatorKind::GRADED: return make_formula<Y>(e.op(),make_formula(e.arg(),spc),e.num());
        case OperatorKind::BINARY: return make_formula<Y>(e.op(),make_formula(e.arg1(),spc),make_formula(e.arg2(),spc));
        case OperatorKind::UNARY: return make_formula<Y>(e.op(),make_formula(e.arg(),spc));
        case OperatorKind::NULLARY: return Formula<Y>::constant(e.val());
        case OperatorKind::VARIABLE: return Formula<Y>::coordinate(spc.index(e.var()));
        default: ARIADNE_FAIL_MSG("Cannot compute formula for expression "<<e.op()<<"of kind "<<e.kind()<<" in space "<<spc);
    }
}

Formula<EffectiveNumber> make_formula(const Expression<Real>& e, const Variable<Real>& var)
{
    return make_formula(e,Space<Real>({var}));
}

Vector<Formula<EffectiveNumber>> make_formula(const Vector<Expression<Real>>& e, const Space<Real>& spc)
{
    Vector<Formula<EffectiveNumber>> res(e.size());
    for(SizeType i=0; i!=e.size(); ++i) {
        res[i]=make_formula(e[i],spc);
    }
    return res;
}

Formula<EffectiveNumber> make_formula(const Expression<Real>& e, const List<Variable<Real>>& vars)
{
    return make_formula(e,Space<Real>(vars));
}

Formula<EffectiveNumber> make_formula(const Expression<Real>& out, const List<Assignment<Variable<Real>,Expression<Real>>>& aux, const Space<Real> spc)
{
    return make_formula(substitute(out,aux),spc);
}

Vector<Formula<EffectiveNumber>> make_formula(const Vector<Expression<Real>>& out, const List<Assignment<Variable<Real>,Expression<Real>>>& aux, const Space<Real> spc)
{
    Vector<Formula<EffectiveNumber>> res(out.size());
    for(SizeType i=0; i!=out.size(); ++i) {
        res[i]=make_formula(out[i],aux,spc);
    }
    return res;
}


Expression<Real> make_expression(const Formula<Real>& f, const Space<Real>& s) {
    const List<RealVariable>& vars=s.variables();
    typedef Algebra<Real> RealAlgebra;
    RealAlgebra az(RealExpression::constant(0));
    Vector<RealAlgebra> va(vars.size(),az);
    for(SizeType i=0; i!=va.size(); ++i) { va[i]=RealAlgebra(RealExpression(vars[i])); }
    RealAlgebra fa=evaluate(f,va);
    return fa.template extract<RealExpression>();
}

Formula<Real> make_formula(const EffectiveScalarFunction& f);
Expression<Real> make_expression(const ScalarFunction<EffectiveTag>& f, const Space<Real>& s) {
    return make_expression(make_formula(f),s); }

} // namespace Ariadne
