/***************************************************************************
 *            expression.tpl.hpp
 *
 *  Copyright 2009--17  Pieter Collins
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
 *  MERCHANTABILITY or FITNESS FOperatorCode::OR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "utility/standard.hpp"

#include "numeric/operators.tpl.hpp"

#include "algebra/algebra.hpp"
#include "algebra/algebra_wrapper.hpp"

#include "symbolic/constant.hpp"
#include "symbolic/variables.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/space.hpp"
#include "symbolic/valuation.hpp"

#include "function/formula.hpp"

namespace Ariadne {

template<class T>
struct ExpressionNode {
    mutable Nat count;
    Operator op;
    virtual ~ExpressionNode();
    explicit ExpressionNode(const Operator& o) : count(0u), op(o) { }
    explicit ExpressionNode(OperatorCode cd, OperatorKind knd) : count(0u), op(cd,knd) { }
};


template<class T> struct ConstantExpressionNode : public ExpressionNode<T> {
    T val;
    ConstantExpressionNode(const T& v) : ExpressionNode<T>(OperatorCode::CNST,OperatorKind::NULLARY), val(v) { }
};
template<class T> struct NamedConstantExpressionNode : public ConstantExpressionNode<T> {
    Identifier name;
    NamedConstantExpressionNode(const Constant<T>& v) : ConstantExpressionNode<T>(v), name(v.name()) { }
};
template<class T> struct VariableExpressionNode : public ExpressionNode<T> {
    Identifier var;
    VariableExpressionNode(const Identifier& v) : ExpressionNode<T>(OperatorCode::VAR,OperatorKind::VARIABLE), var(v) { }
};
template<class T, class A=T> struct UnaryExpressionNode : public ExpressionNode<T> {
    Expression<A> arg;
    UnaryExpressionNode(const Operator& op, Expression<A> const& a) : ExpressionNode<T>(op), arg(a) { }
};
template<class T, class A1=T, class A2=A1> struct BinaryExpressionNode : public ExpressionNode<T> {
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

template<class T> inline OutputStream& operator<<(OutputStream& os, const ExpressionNode<T>* e) {
    return os << static_cast<Void const*>(e);
}


template<class T> Expression<T>::Expression() : _root(new ConstantExpressionNode<T>(T())) { };
template<class T> Expression<T>::Expression(const T& c): _root(new ConstantExpressionNode<T>(c)) { };
template<class T> Expression<T>::Expression(const Constant<T>& c): _root(new NamedConstantExpressionNode<T>(c)) { };
template<class T> Expression<T>::Expression(const Variable<T>& v) : _root(new VariableExpressionNode<T>(v.name())) { }
template<class T> Expression<T> Expression<T>::constant(const T& c) {
    return Expression<T>(SharedPointer<const ExpressionNode<T>>(new ConstantExpressionNode<T>(c))); }
template<class T> Expression<T> Expression<T>::constant(const Constant<T>& c) {
    return Expression<T>(SharedPointer<const ExpressionNode<T>>(new NamedConstantExpressionNode<T>(c))); }
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
        }
        default: ARIADNE_FAIL_MSG("Cannot compute arguments of expression "<<e<<" of kind "<<e.kind()<<"\n");
    }
}

namespace {

template<class T> inline OutputStream& _write_comparison(OutputStream& os, const Expression<T>& f) {
    ARIADNE_FAIL_MSG("Comparison must return a logical type."); }
template<> inline OutputStream& _write_comparison(OutputStream& os, const Expression<Kleenean>& f) {
    Real* real_ptr=0; return os << "(" << f.cmp1(real_ptr) << symbol(f.op()) << f.cmp2(real_ptr) << ")"; }
template<> inline OutputStream& _write_comparison(OutputStream& os, const Expression<Boolean>& f) {
    String* string_ptr=0; return os << "(" << f.cmp1(string_ptr) << symbol(f.op()) << f.cmp2(string_ptr) << ")"; }
//FIXME: Distinguish String and Integer comparisons

}

template<class T> OutputStream& Expression<T>::_write(OutputStream& os) const {
    const Expression<T>& f=*this;
    switch(f.op()) {
        //case OperatorCode::CNST: return os << std::fixed << std::setprecision(4) << fptr->val;
        case OperatorCode::CNST:
            if(auto fp=std::dynamic_pointer_cast<NamedConstantExpressionNode<T>const>(f._root)) { return os << fp->name; }
            else if(auto afp=std::dynamic_pointer_cast<ConstantExpressionNode<T>const>(f._root)) { os <<
            "c"<<afp->val; return os; }
            else { return os << "dcf"; }
            //if(f.val()==0.0) { return os << 0.0; } if(abs(f.val())<1e-4) { os << std::fixed << f.val(); } else { os << f.val(); } return os;
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

template<class R> inline
Expression<R> make_expression(const Constant<R>& c) {
    return Expression<R>(std::make_shared<NamedConstantExpressionNode<R>>(c)); }
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



template<class A> typename Logic<A>::Type evaluate(const Expression<typename Logic<A>::Type>& e, const Map<Identifier,A>& x) {
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

template<class T> T evaluate(const Expression<T>& e, const Map<Identifier,T>& x) {
    switch(e.kind()) {
        case OperatorKind::VARIABLE: return x[e.var()];
        case OperatorKind::NULLARY: return e.val();
        case OperatorKind::UNARY: return compute(e.op(),evaluate(e.arg(),x));
        case OperatorKind::BINARY: return compute(e.op(),evaluate(e.arg1(),x),evaluate(e.arg2(),x));
        case OperatorKind::GRADED: return compute(e.op(),evaluate(e.arg(),x),e.num());
        default: ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" on "<<x<<"\n");
    }
}

template<class T> T evaluate(const Expression<T>& e, const Valuation<T>& x) {
    return evaluate(e,x.values());
}

template<class A> typename Logic<A>::Type evaluate(const Expression<typename Logic<A>::Type>& e, const Valuation<A>& x) {
    return evaluate(e,x.values());
}


template<class T> Set<Identifier> arguments(const Expression<T>& e) {
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
        }
        default: ARIADNE_FAIL_MSG("Cannot compute arguments of expression "<<e<<" of kind "<<e.kind()<<"\n");
    }
}



namespace {
template<class I, class X, class Y> inline const Expression<X>& _substitute_variable(const I& ie, const I& is, const Expression<X>& e, const Expression<Y>& s) {
    ARIADNE_ASSERT_MSG(ie!=is,"Cannot substitute expression "<<s<<" for variable "<<ie<<"\n");
    return e; }
template<class I, class X> inline const Expression<X>& _substitute_variable(const I& ie, const I& is, const Expression<X>& e, const Expression<X>& s) {
    return ie==is ? s : e; }
} // namespace

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Expression<Y>& s) {
    switch(e.kind()) {
        case OperatorKind::COMPARISON: {
            Y* yptr=0;
            const Expression<Y>& c1=e.cmp1(yptr);
            const Expression<Y>& c2=e.cmp2(yptr);
            return make_expression<X>(e.op(),substitute(c1,v,s),substitute(c2,v,s)); }
        case OperatorKind::BINARY: return make_expression<X>(e.op(),substitute(e.arg1(),v,s),substitute(e.arg2(),v,s));
        case OperatorKind::UNARY: return make_expression<X>(e.op(),substitute(e.arg(),v,s));
        case OperatorKind::GRADED: return make_expression<X>(e.op(),substitute(e.arg(),v,s),e.num());
        case OperatorKind::NULLARY: return make_expression<X>(e.val());
        case OperatorKind::VARIABLE: return _substitute_variable(e.var(),v.name(),e,s);
        default: ARIADNE_FAIL_MSG("Cannot substitute "<<s<<" for a named variable "<<v<<" in an unknown expression "<<e<<"\n");
    }
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Y& c) {
    return substitute(e,v,Expression<Y>::constant(c));
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

template<class T> Bool is_constant(const Expression<T>& e, const SelfType<T>& c) {
    switch(e.op()) {
        case OperatorCode::CNST: return decide(e.val()==c);
        default: return false;
    }
}

template<class T> Bool is_variable(const Expression<T>& e, const Variable<T>& v) {
    switch(e.op()) {
        case OperatorCode::VAR: return e.var()==v.name();
        default: return false;
    }
}


template<class T> Bool identical(const Expression<T>& e1, const Expression<T>& e2)
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

} // namespace Ariadne
