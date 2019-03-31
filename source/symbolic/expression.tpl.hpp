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

using Eq = Equal;
using Neq = Unequal;

template<class SIG> struct OperatorTypedef;
template<class SIG> using OperatorType = typename OperatorTypedef<SIG>::Type;

struct CatOp {
    static constexpr OperatorCode code() { return OperatorCode::ADD; }
    static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    String operator()(String const& s1, String const& s2) const { return s1+s2; }
    friend OutputStream& operator<<(OutputStream& os, CatOp op) { return os << "cat"; }
};
inline String pos(String s) { return s; }


template<> struct OperatorTypedef<String(String)> { typedef OperatorVariant<> Type; };

template<> struct OperatorTypedef<Integer(Integer)> { typedef OperatorVariant<Nul,Pos,Neg,Sqr,Abs> Type; };
template<> struct OperatorTypedef<Integer(Integer,Integer)> { typedef OperatorVariant<Add,Sub,Mul,Max,Min> Type; };

template<> struct OperatorTypedef<Real(Real)> { typedef UnaryElementaryOperator Type; };
template<> struct OperatorTypedef<Real(Real,Real)> { typedef BinaryElementaryOperator Type; };
template<> struct OperatorTypedef<Real(Real,Integer)> { typedef GradedElementaryOperator Type; };
template<> struct OperatorTypedef<Real(Real,Int)> { typedef GradedElementaryOperator Type; };

template<> struct OperatorTypedef<Boolean(String,String)> { typedef OperatorVariant<Eq,Neq> Type; };
template<> struct OperatorTypedef<Boolean(Integer,Integer)> { typedef BinaryComparisonOperator Type; };
template<> struct OperatorTypedef<Kleenean(Real)> { typedef OperatorVariant<Sgn> Type; };
//template<> struct OperatorTypedef<Kleenean(Real,Real)> { typedef OperatorVariant<Less,Gtr,Leq,Geq> Type; };
template<> struct OperatorTypedef<Kleenean(Real,Real)> { typedef BinaryComparisonOperator Type; };

template<> struct OperatorTypedef<Boolean(Boolean)> { typedef UnaryLogicalOperator Type; };
template<> struct OperatorTypedef<Boolean(Boolean,Boolean)> { typedef BinaryLogicalOperator Type; };
template<> struct OperatorTypedef<Kleenean(Kleenean)> { typedef UnaryLogicalOperator Type; };
template<> struct OperatorTypedef<Kleenean(Kleenean,Kleenean)> { typedef BinaryLogicalOperator Type; };

template<class SIG> struct SymbolicTypedef;
template<class SIG> using SymbolicType = typename SymbolicTypedef<SIG>::Type;

template<class R, class A> struct SymbolicTypedef<R(A)> {
    typedef Symbolic<OperatorType<R(A)>,Expression<A>> Type; };
template<class R, class A1, class A2> struct SymbolicTypedef<R(A1,A2)> {
    typedef Symbolic<OperatorType<R(A1,A2)>,Expression<A1>,Expression<A2>> Type; };
template<class R, class... AS> struct SymbolicTypedef<R(AS...)> {
    typedef Symbolic<OperatorType<R(AS...)>,Expression<AS>...> Type; };


template<class T> struct ExpressionVariantTypedef;
template<class T> using ExpressionVariantType = typename ExpressionVariantTypedef<T>::Type;

template<> struct ExpressionVariantTypedef<String> {
    typedef Variant< Constant<String>, Variable<String> > Type;
};

template<> struct ExpressionVariantTypedef<Integer> {
    typedef Variant< Constant<Integer>, Variable<Integer>,
                     SymbolicType<Integer(Integer)>, SymbolicType<Integer(Integer,Integer)> > Type;
};

template<> struct ExpressionVariantTypedef<Real> {
    typedef Variant< Constant<Real>, Variable<Real>,
                     SymbolicType<Real(Real)>, SymbolicType<Real(Real,Real)>,
                     Symbolic<OperatorType<Real(Real,Integer)>,Expression<Real>,Int> > Type;
};

template<> struct ExpressionVariantTypedef<Boolean> {
    typedef Variant< Constant<Boolean>, Variable<Boolean>,
                     SymbolicType<Boolean(Boolean)>, SymbolicType<Boolean(Boolean,Boolean)>,
                     SymbolicType<Boolean(String,String)>, SymbolicType<Boolean(Integer,Integer)> > Type;
};

template<> struct ExpressionVariantTypedef<Kleenean> {
    typedef Variant< Constant<Kleenean>, Variable<Kleenean>,
                     SymbolicType<Kleenean(Kleenean)>, SymbolicType<Kleenean(Kleenean,Kleenean)>,
                     SymbolicType<Kleenean(Real)>, SymbolicType<Kleenean(Real,Real)> > Type;
};

namespace {
template<class T> inline Cnst _op_impl(Constant<T> const&) { Cnst op; return op; }
template<class T> inline Var _op_impl(Variable<T> const&) { return Var(); }
template<class OP, class... AS> inline OP _op_impl(Symbolic<OP,AS...> const& s) { return s._op; }
}

template<class T> struct ExpressionNode : public ExpressionVariantType<T> {
  public:
    template<class... AS, EnableIf<IsConstructible<ExpressionVariantType<T>,AS...>> =dummy>
        ExpressionNode(AS... as) : ExpressionVariantType<T>(as...) { }

    ExpressionVariantType<T> const& base() const { return *this; }

    template<class VIS> decltype(auto) accept(VIS&& vis) const {
        return std::visit(std::forward<VIS>(vis),static_cast<ExpressionVariantType<T>const&>(*this)); }

    Operator op() const { return this->accept([](auto s){return Operator(_op_impl(s));}); }
};

template<class T> inline OutputStream& operator<<(OutputStream& os, const ExpressionNode<T>* e) {
    return os << static_cast<Void const*>(e);
}

template<class T> using ConstantExpressionNode = Constant<T>;
template<class T> using NamedConstantExpressionNode = ConstantExpressionNode<T>;
template<class T> using VariableExpressionNode = Variable<T>;

template<class R, class A=R> using UnaryExpressionNode = Symbolic<OperatorType<R(A)>,Expression<A>>;
template<class R, class A1=R, class A2=A1> using BinaryExpressionNode = Symbolic<OperatorType<R(A1,A2)>,Expression<A1>,Expression<A2>>;
template<class R, class A=R, class N=Int> using GradedExpressionNode = Symbolic<OperatorType<R(A,N)>,Expression<A>,N>;


template<class T> Expression<T>::Expression() : Expression(T()) { }
template<class T> Expression<T>::Expression(const T& c) : Expression(Constant<T>(c)) { }
template<class T> Expression<T>::Expression(const Constant<T>& c): _root(new ExpressionNode<T>(c)) { }
template<class T> Expression<T>::Expression(const Variable<T>& v) : _root(new ExpressionNode<T>(v)) { }

template<> Expression<String>::Expression(const String& c): _root(new ExpressionNode<String>(Constant<String>(c))) { }

template<class T> Expression<T> Expression<T>::constant(const T& c) {
    return Expression<T>(c); }
template<class T> Expression<T> Expression<T>::constant(const Constant<T>& c) {
    return Expression<T>(c); }
template<class T> Expression<T> Expression<T>::variable(const Identifier& v) {
    return Expression<T>(Variable<T>(v)); }

template<class T> Operator Expression<T>::op() const {
    return this->node_ptr()->op(); }
template<class T> OperatorCode Expression<T>::code() const {
    return node_ptr()->op().code(); }
template<class T> OperatorKind Expression<T>::kind() const {
    return node_ptr()->op().kind(); }
template<class T> const T& Expression<T>::val() const {
    return std::get<Constant<T>>(node_ref()); }
template<class T> const Identifier& Expression<T>::var() const {
    return std::get<VariableExpressionNode<T>>(node_ref()).name(); }
template<class T> const Expression<T>& Expression<T>::arg() const {
    if constexpr (IsSame<T,Real>::value) { if (auto* gn = std::get_if<GradedExpressionNode<T>>(&node_ref())) { return gn->_arg; } }
    if constexpr (not IsSame<T,String>::value) { return std::get<UnaryExpressionNode<T>>(node_ref())._arg; } else { assert(false); } }
template<class T> const Int& Expression<T>::num() const {
    if constexpr (IsSame<T,Real>::value) { return std::get<GradedExpressionNode<T>>(node_ref())._num; }
    else { assert(false); } }
template<class T> const Expression<T>& Expression<T>::arg1() const {
    if constexpr (not IsSame<T,String>::value) { return std::get<BinaryExpressionNode<T>>(node_ref())._arg1; } else { assert(false); } }
template<class T> const Expression<T>& Expression<T>::arg2() const {
    if constexpr (not IsSame<T,String>::value) { return std::get<BinaryExpressionNode<T>>(node_ref())._arg2; } else { assert(false); } }
template<class R> template<class A> const Expression<A>& Expression<R>::cmp1(A*) const {
    return std::get<Symbolic<OperatorType<R(A,A)>,Expression<A>,Expression<A>>>(node_ref())._arg1; }
template<class R> template<class A> const Expression<A>& Expression<R>::cmp2(A*) const {
    return std::get<Symbolic<OperatorType<R(A,A)>,Expression<A>,Expression<A>>>(node_ref())._arg2; }


template<class T> Set<UntypedVariable> Expression<T>::arguments() const {
    return this->node_ref().accept([](auto s){return _arguments(s);});
}

template<class T> OutputStream& Expression<T>::_write(OutputStream& os) const {
    this->node_ref().accept([&os](auto e){_write_impl(os,e);}); return os;
/*
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
                case OperatorKind::UNARY: case OperatorKind::GRADED: return os << f.op() << "(" << f.arg() << ")";
                case OperatorKind::BINARY: return os << f.op() << "(" << f.arg1() << "," << f.arg2() << ")";
                // FIXME: Type-cast comparison arguments correctly
                case OperatorKind::COMPARISON: return _write_comparison(os,f);
                default: ARIADNE_FAIL_MSG("Cannot output expression with operator "<<f.op()<<" of kind "<<f.kind()<<"\n");
            }
    }
*/
    return os;
}

template<class R> inline
Expression<R> make_expression(const Constant<R>& c) {
    return Expression<R>(std::make_shared<NamedConstantExpressionNode<R>>(c)); }
template<class R> inline
Expression<R> make_expression(const R& c) {
    return Expression<R>(std::make_shared<ExpressionNode<R>>(ConstantExpressionNode<R>(Constant<R>(c)))); }
template<class R, class A> inline
Expression<R> make_expression(OperatorType<R(A)> op, const Expression<A>& e) {
    return Expression<R>(std::make_shared<ExpressionNode<R>>(UnaryExpressionNode<R,A>(op,e))); }
template<class R, class A, class N> inline
Expression<R> make_expression(OperatorType<R(A,N)> op, const Expression<A>& e, N n) {
    return Expression<R>(std::make_shared<ExpressionNode<R>>(GradedExpressionNode<R,A,N>(op,e,n))); }
template<class R, class A1, class A2> inline
Expression<R> make_expression(OperatorType<R(A1,A2)> op, const Expression<A1>& e1, Expression<A2> e2) {
    return Expression<R>(std::make_shared<ExpressionNode<R>>(BinaryExpressionNode<R,A1,A2>(op,e1,e2))); }

template<class R, class OP, class A> inline
Expression<R> make_expression(OP op, const Expression<A> e) {
    return make_expression<R>(OperatorType<R(A)>(op),e); }
template<class R, class A, class N> inline
Expression<R> make_expression(Pow op, const Expression<A>& e, N n) {
    return make_expression<R>(OperatorType<R(A,Int)>(op),e,n); }
template<class R, class OP, class A1, class A2> inline
Expression<R> make_expression(OP op, const Expression<A1>& e1, Expression<A2> e2) {
    return make_expression<R>(OperatorType<R(A1,A2)>(op),e1,e2); }



template<class A> A evaluate(const Expression<A>& e, const Map<Identifier,A>& x);

inline Integer evaluate(const Expression<Integer>& e, const Map<Identifier,String>& x) {
    return e.node_ref().accept([&x](auto en){return evaluate_as<Integer>(en,x);});
}
inline String evaluate(const Expression<String>& e, const Map<Identifier,Integer>& x) {
    return e.node_ref().accept([&x](auto en){return evaluate_as<String>(en,x);});
}


template<class A> typename Logic<A>::Type evaluate(const Expression<typename Logic<A>::Type>& e, const Map<Identifier,A>& x) {
    typedef typename Logic<A>::Type R;
    return e.node_ref().accept([&x](auto en){return evaluate_as<R>(en,x);});
}

template<class T> T evaluate(const Expression<T>& e, const Map<Identifier,T>& x) {
    return e.node_ref().accept([&x](auto en){return evaluate_as<T>(en,x);});
}

template<class T> T evaluate(const Expression<T>& e, const Valuation<T>& x) {
    return evaluate(e,x.values());
}

template<class A> typename Logic<A>::Type evaluate(const Expression<typename Logic<A>::Type>& e, const Valuation<A>& x) {
    return evaluate(e,x.values());
}


template<class T> Set<Identifier> arguments(const Expression<T>& e) {
    Set<UntypedVariable> arg_vars = e.arguments();
    Set<Identifier> arg_names;
    for (auto var : arg_vars) { arg_names.insert(var.name()); }
    return arg_names;
}


template<class T> Bool is_constant_in(const Expression<T>& e, const Set<Variable<T>>& spc) {
    return e.node_ref().accept([&spc](auto en){return is_constant_in(en,spc);});
}


namespace {
template<class X, class Y> Expression<X> _substitute(const Constant<X>& c, const Variable<Y>& v, const Expression<Y>& s) { return c; }
template<class X, class Y> Expression<X> _substitute(const Variable<X>& var, const Variable<Y>& v, const Expression<Y>& s) { return var; }
template<class X> Expression<X> _substitute(const Variable<X>& var, const Variable<X>& v, const Expression<X>& s) {
    if (var==v) { return s; } else { return var; } }
template<class X, class Y, class OP, class E> Expression<X> _substitute(const Symbolic<OP,E>& e, const Variable<Y>& v, const Expression<Y>& s) {
    return e._op(substitute(e._arg,v,s)); }
template<class X, class Y, class OP, class E1, class E2> Expression<X> _substitute(const Symbolic<OP,E1,E2>& e, const Variable<Y>& v, const Expression<Y>& s) {
    return e._op(substitute(e._arg1,v,s),substitute(e._arg2,v,s)); }
template<class X, class Y, class OP, class E> Expression<X> _substitute(const Symbolic<OP,E,Int>& e, const Variable<Y>& v, const Expression<Y>& s) {
    return e._op(substitute(e._arg,v,s),e._num); }
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Expression<Y>& s) {
    return e.node_ref().accept([&v,&s](auto en){return _substitute<X>(en,v,s);});
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

template<class X, class Y> Vector<Expression<X>> substitute(const Vector<Expression<X>>& e, const List< Assignment< Variable<Y>, Expression<Y> > >& a) {
    Vector<Expression<X>> r(e.size());
    for(SizeType i=0; i!=e.size(); ++i) {
        r[i]=substitute(e[i],a);
    }
    return r;
}



inline Bool same(Kleenean k1, Kleenean k2) { return definitely(k1==k2); }


template<class T> Bool identical(const Expression<T>& e1, const Expression<T>& e2)
{
    if(e1.node_raw_ptr()==e2.node_raw_ptr()) { return true; }
    if(e1.op()!=e2.op()) { return false; }
    return e1.node_ref().accept([&e2](auto e1n){return e2.node_ref().accept([&e1n](auto e2n){return _identical_dispatch(e1n,e2n);});});
}



/*
template<class OP1, class OP2> constexpr Bool are_inverses(OP1 op1, OP2 op2) { return false; }

constexpr Bool are_inverses(Pos,Pos) { return true; }
constexpr Bool are_inverses(Neg,Neg) { return true; }
constexpr Bool are_inverses(Rec,Rec) { return true; }
constexpr Bool are_inverses(Sqr,Sqrt) { return true; }
constexpr Bool are_inverses(Exp,Log) { return true; }
constexpr Bool are_inverses(Log,Exp) { return true; }
constexpr Bool are_inverses(Tan,Atan) { return true; }
constexpr Bool are_inverses(Atan,Tan) { return true; }
*/

constexpr Pos inverse(Pos) { return Pos(); }
constexpr Neg inverse(Neg) { return Neg(); }
constexpr Rec inverse(Rec) { return Rec(); }
constexpr Sqrt inverse(Sqr) { return Sqrt(); }
constexpr Sqr inverse(Sqrt) { return Sqr(); }
constexpr Log inverse(Exp) { return Log(); }
constexpr Exp inverse(Log) { return Exp(); }
constexpr Atan inverse(Tan) { return Atan(); }
constexpr Tan inverse(Atan) { return Tan(); }

template<class OP> using InverseType = decltype(inverse(declval<OP>()));

template<class OP> struct HasInverse {
    template<class O, class=InverseType<O>> static std::true_type test(int);
    template<class O> static std::false_type test(...);
    static const bool value = decltype(test<OP>(1))::value;
};

template<class OP, class... OPS> Bool _are_inverses(OP const& op1, OperatorVariant<OPS...> const& ops2) {
    if constexpr(HasInverse<decltype(op1)>::value) { return inverse(op1).code() == ops2.code(); }
    else { return false; }
}

template<class... OPS> Bool are_inverses(OperatorVariant<OPS...> const& ops1, OperatorVariant<OPS...> const& ops2) {
    return ops1.accept([&ops2](auto op1){ return _are_inverses(op1,ops2); });
}

namespace {

template<class X> inline Expression<X> _simplify(const Expression<X>& e) {
    return e;
}

template<class T> using C = Constant<T>;
template<class T> using V = Variable<T>;
template<class T> using E = Expression<T>;

using R = Real; using RE = E<R>; using REcr=E<R> const&; using RCcr = const C<R>&;
using K = Kleenean; using KE = E<K>; using KEcr=E<K> const&; using KCcr = const C<K>&;

template<class OP, class A1, class A2> decltype(auto) _simpl(OP op, Constant<A1> a1, Constant<A2> a2) {
    auto r=op(a1.value(),a2.value()); return Constant<decltype(r)>(r); }
template<class OP, class A> decltype(auto) _simpl(OP op, Constant<A> a) {
    auto r=op(a.value()); return Constant<decltype(r)>(r); }

template<class E> inline RE _simpl(Add op, E const& e1, E const& e2) {
    if (is_constant(e1,0)) { return e2; }
    else if (is_constant(e2,0)) { return e1; }
    else { return op(e1,e2); } }
template<class E> inline RE _simpl(Sub op, E const& e1, E const& e2) {
    if (is_constant(e1,0)) { return neg(e2); }
    else if (is_constant(e2,0)) { return e1; }
    else { return op(e1,e2); } }
template<class E> inline RE _simpl(Mul op, E const& e1, E const& e2) {
    if (is_constant(e1,0) or is_constant(e2,1)) { return e1; }
    else if (is_constant(e2,0) or is_constant(e1,1)) { return e2; }
    else { return op(e1,e2); } }
template<class E> inline RE _simpl(Div op, E const& e1, E const& e2) {
    if (is_constant(e1,0) or is_constant(e2,1)) { return e1; }
    else if (is_constant(e1,1)) { return rec(e2); } else { return op(e1,e2); } }
template<class E> inline RE _simpl(Max op, E const& e1, E const& e2) { return op(e1,e2); }
template<class E> inline RE _simpl(Min op, E const& e1, E const& e2) { return op(e1,e2); }

template<class E> inline E _simpl(Pow op, E const& e, Int n) {
    switch (n) {
        case -1: return rec(e);
        case 0: return E::constant(1);
        case 1: return e;
        case 2: return sqr(e);
        default: return pow(e,n);
    }
}

template<class E> inline E _simpl(AndOp op, E const& e1, E const& e2) {
    if (is_constant(e1,false) or is_constant(e2,true)) { return e1; }
    else if (is_constant(e1,true) or is_constant(e2,false)) { return e2; }
    else { return op(e1,e2); } }
template<class E> inline E _simpl(OrOp op, E const& e1, E const& e2) {
    if (is_constant(e1,true) or is_constant(e2,false)) { return e1; }
    else if (is_constant(e1,false) or is_constant(e2,true)) { return e2; }
    else { return op(e1,e2); } }
template<class E> inline E _simpl(NotOp op, E const& e) {
    if (e.op().code() == op.code()) { return e.arg(); }
    else { return op(e); }
}

template<class T> inline Expression<T> _simplify_node(const Constant<T>& c) { return c; }
template<class T> inline Expression<T> _simplify_node(const Variable<T>& v) { return v; }

template<class T> inline Expression<T> _simplify_node(const UnaryExpressionNode<T>& e) {
    auto op = e._op;
    auto sarg = simplify(e._arg);
    if (holds_alternative<Pos>(op)) {
        return sarg;
    } else if (auto* cp = std::get_if<Constant<T>>(&sarg.node_ref().base())) {
        return op(*cp);
    } else {
        auto* up = std::get_if<UnaryExpressionNode<T>>(&sarg.node_ref());
        if (up && are_inverses(op,up->_op)) { return up->_arg; }
        else { return op(sarg); }
    }
}
template<class T> inline Expression<T> _simplify_node(const BinaryExpressionNode<T>& e) {
    Expression<T> sarg1=simplify(e._arg1);
    Expression<T> sarg2=simplify(e._arg2);
    return e.op().accept([&](auto op){return _simpl(op,sarg1,sarg2);});
}
inline Expression<R> _simplify_node(const BinaryExpressionNode<R>& e) {
    Expression<R> sarg1=simplify(e._arg1);
    Expression<R> sarg2=simplify(e._arg2);
    return e.op().accept([&](auto op){return _simpl(op,sarg1,sarg2);});
}
template<class T> inline Expression<T> _simplify_node(const GradedExpressionNode<T>& e) {
    Expression<R> sarg=simplify(e.arg());
    return e.op().accept([&](auto op){return _simpl(op,sarg,e.num());});
}

//inline Expression<Real> _simplify(const Expression<Real>& e) {
Expression<Real> _simplify(const Expression<Real>& e) {
    return e.node_ref().accept([](auto en){return _simplify_node(en);});
}

inline Expression<Kleenean> _simplify(const UnaryExpressionNode<Kleenean,Real>& e) {
    return make_expression<Kleenean>(e.op(),simplify(e.arg()));
}
inline Expression<Kleenean> _simplify(const BinaryExpressionNode<Kleenean,Real,Real>& e) {
    return make_expression<Kleenean>(e.op(),simplify(e.arg1()),simplify(e.arg2()));
}

template<class I> inline Expression<Kleenean> _simplify(const Expression<Kleenean>& e) {
    return e.node_ref().accept([](auto en){return _simpliy(en);});
}

} // namespace

template<class X> Expression<X> simplify(const Expression<X>& e) {
    return Ariadne::_simplify(e);
}



template<class T> Bool before(Expression<T> const& e1, Expression<T> const& e2) {
    if(e1.op()==e2.op()) {
        switch(e1.kind()) {
            case OperatorKind::VARIABLE:
                return e1.var() < e2.var();
            case OperatorKind::NULLARY:
                return decide(e1.val() < e2.val());
            case OperatorKind::UNARY:
                return before(e1.arg(),e2.arg());
            case OperatorKind::GRADED:
                return (e1.num() == e2.num() ? before(e1.arg(),e2.arg()) : (e1.num() < e2.num()));
            case OperatorKind::BINARY:
                return identical(e1.arg1(),e2.arg1()) ? before(e1.arg2(),e2.arg2()) : before(e1.arg1(),e2.arg1());
            default:
                return false;
        }
    } else {
        return e1.op() < e2.op();
    }
}

struct ExpressionComparator {
    template<class A1, class A2> decltype(auto) operator() (A1&& a1, A2&& a2) const { return before(std::forward<A1>(a1),std::forward<A2>(a2)); }
};


template<class T> Expression<T> eliminate_common_subexpressions(const Expression<T>& e, std::set<Expression<T>,ExpressionComparator>& cache)
{
    auto iter=cache.find(e);
    if (iter!=cache.end()) {
        return *iter; }
    else {
        switch(e.kind()) {
            case OperatorKind::VARIABLE: case OperatorKind::NULLARY:
                cache.insert(e); return e;
            case OperatorKind::UNARY: {
                auto new_arg = eliminate_common_subexpressions(e.arg(),cache);
                if(new_arg.node_raw_ptr() == e.arg().node_raw_ptr()) {
                    cache.insert(e); return e;
                } else {
                    Expression<T> new_e=make_expression<T>(e.op(),new_arg);
                    cache.insert(new_e); return new_e;
                }
            }
            case OperatorKind::GRADED: {
                auto new_arg = eliminate_common_subexpressions(e.arg(),cache);
                if(new_arg.node_raw_ptr() == e.arg().node_raw_ptr()) {
                    cache.insert(e); return e;
                } else {
                    Expression<T> new_e=make_expression<T>(GradedElementaryOperator(e.op()),new_arg,e.num());
                    cache.insert(new_e); return new_e;
                }
            }
            case OperatorKind::BINARY: {
                auto new_arg1 = eliminate_common_subexpressions(e.arg1(),cache);
                auto new_arg2 = eliminate_common_subexpressions(e.arg2(),cache);
                if(new_arg1.node_raw_ptr() == e.arg1().node_raw_ptr() && new_arg2.node_raw_ptr() == e.arg2().node_raw_ptr()) {
                    cache.insert(e); return e;
                } else {
                    Expression<T> new_e=make_expression<T>(e.op(),new_arg1,new_arg2);
                    cache.insert(new_e); return new_e;
                }
            }
            default:
                abort();
        }
    }
}

template<class T> Void eliminate_common_subexpressions(Expression<T>& e)
{
    std::set<Expression<T>,ExpressionComparator> cache;
    e=eliminate_common_subexpressions(e,cache);
}

template<class T> Void eliminate_common_subexpressions(Vector<Expression<T>>& es)
{
    std::set<Expression<T>,ExpressionComparator> cache;
    for(SizeType i=0; i!=es.size(); ++i) { es[i]=eliminate_common_subexpressions(es[i],cache); }
}

struct ExpressionPtrComparator {
    template<class A1, class A2> decltype(auto) operator() (A1&& a1, A2&& a2) const { return a1.node_raw_ptr() < a2.node_raw_ptr(); }
};


template<class T, class C> Void _count_nodes(Expression<T> const& e, std::set<Expression<T>,C>& cache, Bool const& distinct, Nat& count) {

    auto iter=cache.find(e);
    if (!distinct || iter==cache.end()) {
        cache.insert(e);
        switch(e.kind()) {
            case OperatorKind::VARIABLE: case OperatorKind::NULLARY:
                break;
            case OperatorKind::UNARY:
                _count_nodes(e.arg(),cache,distinct,count);
                break;
            case OperatorKind::GRADED:
                _count_nodes(e.arg(),cache,distinct,count);
                break;
            case OperatorKind::BINARY:
                _count_nodes(e.arg1(),cache,distinct,count);
                _count_nodes(e.arg2(),cache,distinct,count);
                break;
            default:
                abort();
        }
        count++;
    }
}


template<class T, class C> Nat count_nodes(Expression<T> const& e, std::set<Expression<T>,C>& cache, Bool const& distinct) {
    Nat result = 0;

    _count_nodes(e, cache, distinct, result);

    return result;
}

template<class T> Nat count_nodes(Expression<T> const& e) {
    std::set<Expression<T>,ExpressionComparator> cache;
    return count_nodes(e,cache,false);
}

template<class T> Nat count_distinct_nodes(Expression<T> const& e) {
    std::set<Expression<T>,ExpressionComparator> cache;
    return count_nodes(e,cache,true);
}

template<class T> Nat count_distinct_node_ptrs(Expression<T> const& e) {
    std::set<Expression<T>,ExpressionPtrComparator> cache;
    return count_nodes(e,cache,true);
}


template<class T> Bool is_constant(const Expression<T>& e, const SelfType<T>& c) {
    auto* ce = std::get_if<Constant<T>>(&e.node_ref()); return ce && decide(ce->value()==c);
}

template<class T> Bool is_variable(const Expression<T>& e, const Variable<T>& v) {
    auto* ve = std::get_if<Variable<T>>(&e.node_ref()); return ve && *ve==v;
}


} // namespace Ariadne
