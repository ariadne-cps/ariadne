/***************************************************************************
 *            symbolic/templates.tpl.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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


namespace Ariadne {

template<class T> class Variable;

template<class OP, class... AS> using Symbolic = ExpressionTemplate<OP,AS...>;

namespace {

template<class R, class OP, class E, class V> R _evaluate_as_impl(const OP& op, const E& e, const V& v) {
    return op(evaluate(e,v)); }
template<class R, class OP, class E1, class E2, class V> R _evaluate_as_impl(const OP& op, const E1& e1, const E2& e2, const V& v) {
    return op(evaluate(e1,v),evaluate(e2,v)); }
template<class R, class OP, class E, class N, class V> R _graded_evaluate_as_impl(const OP& op, const E& e, N n, const V& v) {
    return op(evaluate(e,v),n); }

template<class R, class E, class V> R _evaluate_as_impl(const UnaryElementaryOperator& op, const E& e, const V& v) {
    return op.call_as<R>(evaluate(e,v)); }
template<class R, class E1, class E2, class V> R _evaluate_as_impl(const BinaryElementaryOperator& op, const E1& e1, const E2& e2, const V& v) {
    return op.call_as<R>(evaluate(e1,v),evaluate(e2,v)); }
template<class R, class E1, class E2, class V> R _evaluate_as_impl(const BinaryComparisonOperator& op, const E1& e1, const E2& e2, const V& v) {
    return op.call_as<R>(evaluate(e1,v),evaluate(e2,v)); }
template<class R, class E, class N, class V> R _graded_evaluate_as_impl(const GradedElementaryOperator& op, const E& e, const N& n, const V& v) {
    return op.call_as<R>(evaluate(e,v),n); }

template<class R, class A> R evaluate_as(const Constant<R>& c, const Map<Identifier,A>& x) { return c; }
template<class R, class A> R evaluate_as(const Variable<R>& v, const Map<Identifier,A>& x) { assert(false); }
template<class R> R evaluate_as(const Variable<R>& v, const Map<Identifier,R>& x) { return x[v.name()]; }
template<class R, class OP, class E, class V> R evaluate_as(const Symbolic<OP,E>& e, const V& x) {
    return _evaluate_as_impl<R>(e._op,e._arg,x); }
template<class R, class OP, class E1, class E2, class V> R evaluate_as(const Symbolic<OP,E1,E2>& e, const V& x) {
    return _evaluate_as_impl<R>(e._op,e._arg1,e._arg2,x); }
template<class R, class OP, class E, class V> R evaluate_as(const Symbolic<OP,E,Int>& e, const V& x) {
    return _graded_evaluate_as_impl<R>(e._op,e._arg,e._num,x); }




template<class Y, template<class>class E> void _write_impl(OutputStream& os, Add, E<Y> e1, E<Y> e2) {
    os << e1 << '+' << e2; }
template<class Y, template<class>class E> void _write_impl(OutputStream& os, Sub, E<Y> e1, E<Y> e2) {
    os << e1 << '+'; switch(e2.op().code()) { case Add::code(): case Sub::code(): os << '(' << e2 << ')'; break; default: os << e2; } }
template<class Y, template<class>class E> void _write_impl(OutputStream& os, Mul, E<Y> e1, E<Y> e2) {
    switch(e1.op().code()) { case Add::code(): case Sub::code(): case Div::code(): os << '(' << e1 << ')'; break; default: os << e1; } os << '*';
    switch(e2.op().code()) { case Add::code(): case Sub::code(): os << '(' << e2 << ')'; break; default: os << e2; } }
template<class Y, template<class>class E> void _write_impl(OutputStream& os, Div, E<Y> e1, E<Y> e2) {
        switch(e1.op()) { case Add::code(): case Sub::code(): case Div::code(): os << '(' << e1 << ')'; break; default: os << e1; } os << '/';
        switch(e2.op()) { case Add::code(): case Sub::code(): case Mul::code(): case Div::code(): os << '(' << e2 << ')'; break; default: os << e2; } }
template<class E1, class E2> void _write_impl(OutputStream& os, Max, E1 const& e1, E2 const& e2) { os << "max" << '(' << e1 << ',' << e2 << ')'; }
template<class E1, class E2> void _write_impl(OutputStream& os, Min, E1 const& e1, E2 const& e2) { os << "min" << '(' << e1 << ',' << e2 << ')'; }

template<class OP, class E1, class E2> void _write_impl(OutputStream& os, OP op, E1 const& e1, E2 const& e2) { os << op << '(' << e1 << ',' << e2 << ')'; }

template<class Y> void _write_impl(OutputStream& os, Symbolic<Cnst,Y> const& c) {
    os << c._val; }
template<class I> void _write_impl(OutputStream& os, Symbolic<Var,I> const& v) {
    os << v._ind; }

template<class T> inline Void _write_impl(OutputStream& os, Constant<T> const& c) {
    os << c; }
template<class T> inline Void _write_impl(OutputStream& os, Variable<T> const& v) {
    os << v; }


#warning TODO: Change to allow templated operators here
template<class A1, class A2, class... OPS> void _write_impl(OutputStream& os, Symbolic<OperatorVariant<OPS...>,A1,A2> const& s) {
    os << s._op << '(' << s._arg1 << ',' << s._arg2 << ')'; } // s._op.accept([&os,&s](auto op){_write_impl(os,op,s._arg1,s._arg2);}); }
template<class A1, class A2> void _write_impl(OutputStream& os, Symbolic<BinaryLogicalOperator,A1,A2> const& s) {
    os << s._op << '(' << s._arg1 << ',' << s._arg2 << ')'; } // s._op.accept([&os,&s](auto op){_write_impl(os,op,s._arg1,s._arg2);}); }
template<class A1, class A2> void _write_impl(OutputStream& os, Symbolic<BinaryComparisonOperator,A1,A2> const& s) {
    os << s._op << '(' << s._arg1 << ',' << s._arg2 << ')'; } // \s._op.accept([&os,&s](auto op){_write_impl(os,op,s._arg1,s._arg2);}); }
template<class A1, class A2> void _write_impl(OutputStream& os, Symbolic<BinaryElementaryOperator,A1,A2> const& s) {
    os << s._op << '(' << s._arg1 << ',' << s._arg2 << ')'; } //  { s._op.accept([&os,&s](auto op){_write_impl(os,op,s._arg1,s._arg2);}); }
template<class OP, class A> void _write_impl(OutputStream& os, Symbolic<OP,A> const& s) {
    os << s._op << '(' << s._arg << ')'; }
template<class A, template<class>class E> void _write_impl(OutputStream& os, Symbolic<BinaryElementaryOperator,A,E<A>> const& s) {
    os << '('; _write_impl(os,Symbolic<BinaryElementaryOperator,E<A>,E<A>>(s._op,E<A>(s._cnst),s._arg)); os << ')'; }
template<class A, class N> void _write_impl(OutputStream& os, Symbolic<GradedElementaryOperator,A,N> const& s) {
    os << s._op << '(' << s._arg << ',' << s._num << ')'; }

//template<class A, class... OPS> void _write_impl(OutputStream& os, Symbolic<OperatorVariant<OPS...>,A> const& s) {
//    s._op.accept([&os,&s](auto op){_write_impl(os,op,s._arg);}); }
//template<class A1, class A2, class... OPS> void _write_impl(OutputStream& os, Symbolic<OperatorVariant<OPS...>,A1,A2> const& s) {
//    s._op.accept([&os,&s](auto op){_write_impl(os,op,s._arg1,s._arg2);}); }
}


namespace {

template<class OP, class F, class J> inline decltype(auto) _derivative_impl(OP op, F const& f, J j) {
    return op.derivative(f,derivative(f,j)); }
template<class F, class J> inline auto _derivative_impl(Abs op, F const& f, J j) -> F {
    ARIADNE_THROW(std::runtime_error,"derivative(abs(f))","Cannot take derivative of non-smooth function"); }

template<class F1, class F2, class J> inline decltype(auto) _derivative_impl(Add, F1 const& f1, F2 const& f2, J j) {
    return derivative(f1,j)+derivative(f2,j); }
template<class F1, class F2, class J> inline decltype(auto) _derivative_impl(Sub, F1 const& f1, F2 const& f2, J j) {
    return derivative(f1,j)-derivative(f2,j); }
template<class F1, class F2, class J> inline decltype(auto) _derivative_impl(Mul, F1 const& f1, F2 const& f2, J j) {
    return derivative(f1,j)*f2+f1*derivative(f2,j); }
template<class F1, class F2, class J> inline decltype(auto) _derivative_impl(Div, F1 const& f1, F2 const& f2, J j) {
    return (derivative(f1,j)+derivative(f2,j)*(f1/f2))/f2; }
template<class F1, class F2, class J> inline auto _derivative_impl(Max, F1 const& f1, F2 const& f2, J j) -> decltype(max(f1,f2)){
    ARIADNE_THROW(std::runtime_error,"derivative(max(f1,f2))","Cannot take derivative of non-smooth function."); }
template<class F1, class F2, class J> inline auto _derivative_impl(Min, F1 const& f1, F2 const& f2, J j) -> decltype(min(f1,f2)) {
    ARIADNE_THROW(std::runtime_error,"derivative(min(f1,f2))","Cannot take derivative of non-smooth function."); }

template<class F, class N, class J> inline decltype(auto) _derivative_impl(Pow op, F const& f, N const& n, J j) {
    return op.derivative(f,derivative(f,j),n); }

#warning Should not need to explicitly use this
template<class F, class J, class... OPS> decltype(auto) _derivative_impl(UnaryElementaryOperator ops, F const& f, J j) {
    return ops.accept([&f,j](auto op){return _derivative_impl(op,f,j);}); }

template<class F, class J, class... OPS> decltype(auto) _derivative_impl(OperatorVariant<OPS...> ops, F const& f, J j) {
    return ops.accept([&f,j](auto op){return _derivative_impl(op,f,j);}); }
template<class F1, class F2, class J, class... OPS> decltype(auto) _derivative_impl(OperatorVariant<OPS...> ops, F1 const& f1, F2 const& f2, J j) {
    return ops.accept([&f1,&f2,j](auto op){return _derivative_impl(op,f1,f2,j);}); }
template<class F, class N, class J> decltype(auto) _derivative_impl(OperatorVariant<Pow> ops, F const& f, N n, J j) {
    return ops.accept([&f,j,n](auto op){return _derivative_impl(op,f,n,j);}); }
template<class X, template<class>class A, class J, class... OPS> decltype(auto) _derivative_impl(OperatorVariant<OPS...> ops, A<X> const& f1, X const& c2, J j) {
    return _derivative_impl(ops,f1,A<X>(c2),j); }
template<class X, template<class>class A, class J, class... OPS> decltype(auto) _derivative_impl(OperatorVariant<OPS...> ops, X const& c1, A<X> const& f2, J j) {
    return _derivative_impl(ops,A<X>(c1),f2,j); }


template<class Y, class J> decltype(auto) derivative(ExpressionTemplate<Cnst,Y> const& s, J j) {
    return ExpressionTemplate<Cnst,Y>(Y(0)); }

template<class OP, class A, class J> decltype(auto) derivative(ExpressionTemplate<OP,A> const& s, J j) {
    return _derivative_impl(s._op,s._arg,j); }
template<class OP, class A1, class A2, class J> decltype(auto) derivative(ExpressionTemplate<OP,A1,A2> const& s, J j) {
    return _derivative_impl(s._op,s._arg1,s._arg2,j); }
template<class OP, class A, template<class>class E, class J> decltype(auto) derivative(ExpressionTemplate<OP,A,E<A>> const& s, J j) {
    return _derivative_impl(s._op,E<A>(s._cnst),s._arg,j); }
template<class OP, class A, template<class>class E, class J> decltype(auto) derivative(ExpressionTemplate<OP,E<A>,A> const& s, J j) {
    return _derivative_impl(s._op,E<A>(s._arg),s._cnst,j); }
template<class OP, class A, class J> decltype(auto) derivative(ExpressionTemplate<OP,A,Int> const& s, J j) {
    return _derivative_impl(s._op,s._arg,s._num,j); }

} // namespace


namespace {

template<class T> Bool identical(const Constant<T>& c1, const Constant<T>& c2) {
    return same(c1.value(),c2.value()); }
template<class T> Bool identical(const Variable<T>& v1, const Variable<T>& v2) {
    return v1==v2; }
template<class OP, class A> Bool identical(const Symbolic<OP,A>& s1, const Symbolic<OP,A>& s2) {
    return s1._op.code()==s2._op.code() && identical(s1._arg,s2._arg); }
template<class OP, class A1, class A2> Bool identical(const Symbolic<OP,A1,A2>& s1, const Symbolic<OP,A1,A2>& s2) {
    return s1._op.code()==s2._op.code() && identical(s1._arg1,s2._arg1) && identical(s1._arg2,s2._arg2); }
template<class OP, class A> Bool identical(const Symbolic<OP,A,Int>& s1, const Symbolic<OP,A,Int>& s2) {
    return s1._op.code()==s2._op.code() && identical(s1._arg,s2._arg) && s1._num==s2._num; }

template<class E> Bool _identical_dispatch(const E& e1, const E& e2) { return identical(e1,e2); }
template<class E1, class E2> Bool _identical_dispatch(const E1& e1, const E2& e2) { return false; }

//template<class E1, class E2> Bool identical(const E1& e1, const E2& e2) { return _identical_dispatch(e1,e2); }

} // namespace



class UntypedVariable;

namespace {
template<class T> Set<UntypedVariable> _arguments(Constant<T> const& var) { return {}; }
template<class T> Set<UntypedVariable> _arguments(Variable<T> const& var) { return {var}; }
template<class OP, class E> Set<UntypedVariable> _arguments(Symbolic<OP,E> const& s) {
    return s._arg.arguments(); }
template<class OP, class E1, class E2> Set<UntypedVariable> _arguments(Symbolic<OP,E1,E2> const& s) {
    return join(s._arg1.arguments(),s._arg2.arguments()); }
template<class OP, class E1> Set<UntypedVariable> _arguments(Symbolic<OP,E1,Int> const& s) {
    return s._arg.arguments(); }
}

namespace {

template<class OP, class E, class VARS> inline Bool _is_constant_in_impl(OP op, E const& e, VARS const& vars) {
    return is_constant_in(e,vars); }
template<class OP, class E1, class E2, class VARS> inline Bool _is_constant_in_impl(OP op, E1 const& e1, E2 const& e2, VARS const& vars) {
    return is_constant_in(e1,vars) && is_constant_in(e2,vars); }
template<class E, class N, class VARS> inline Bool _is_constant_in_graded_impl(Pow op, E const& e, N n, VARS const& vars) {
    return n==0 || is_constant_in(e,vars); }

template<class E, class VARS, class... OPS> Bool _is_constant_in_impl(OperatorVariant<OPS...> ops, E const& e, VARS const& vars) {
    return ops.accept([&e,&vars](auto op){return _is_constant_in_impl(op,e, vars);}); }
template<class E1, class E2, class VARS, class... OPS> Bool _is_constant_in_impl(OperatorVariant<OPS...> ops, E1 const& e1, E2 const& e2, VARS const& vars) {
    return ops.accept([&e1,&e2,&vars](auto op){return _is_constant_in_impl(op,e1,e2, vars);}); }
template<class E, class VARS, class... OPS> Bool _is_constant_in_graded_impl(OperatorVariant<OPS...> ops, E const& e, Int n, VARS const& vars) {
    return ops.accept([&e,n,&vars](auto op){return _is_constant_in_graded_impl(op,e,n, vars);}); }
template<class OP, class E, class N, class VARS> inline Bool _is_constant_in_graded_impl(GradedElementaryOperator ops, E const& e, N n, VARS const& vars) {
    return ops.accept([&e,n,&vars](auto op){return _is_constant_in_graded_impl(op,e,n,vars);}); }

template<class T, class VARS> inline Bool is_constant_in(Constant<T> const&, VARS const&) {
    return true; }
template<class T, class VARS> inline Bool is_constant_in(Variable<T> const& v, VARS const& vars) {
    return not vars.contains(v); }
template<class OP, class A, class VARS> inline Bool is_constant_in(Symbolic<OP,A> const& s, VARS const& vars) {
    return _is_constant_in_impl(s._op,s._arg,vars); }
template<class OP, class A1, class A2, class VARS> inline Bool is_constant_in(Symbolic<OP,A1,A2> const& s, VARS const& vars) {
    return _is_constant_in_impl(s._op,s._arg1,s._arg2,vars); }
template<class OP, class A, class VARS> inline Bool is_constant_in(Symbolic<OP,A,Int> const& s, VARS const& vars) {
    return _is_constant_in_graded_impl(s._op,s._arg,s._num,vars); }



template<class E, class VARS> inline Bool _is_affine_in_impl(Neg, E const& e, VARS const& vars) {
    return is_affine_in(e,vars); }
template<class OP, class E, class VARS> inline Bool _is_affine_in_impl(OP, E const& e, VARS const& vars) {
    return is_constant_in(e,vars); }
template<class E1, class E2, class VARS> inline Bool _is_affine_in_impl(Add, E1 const& e1, E2 const& e2, VARS const& vars) {
    return is_affine_in(e1,vars) && is_affine_in(e2,vars); }
template<class E1, class E2, class VARS> inline Bool _is_affine_in_impl(Sub, E1 const& e1, E2 const& e2, VARS const& vars) {
    return is_affine_in(e1,vars) && is_affine_in(e2,vars); }
template<class E1, class E2, class VARS> inline Bool _is_affine_in_impl(Mul, E1 const& e1, E2 const& e2, VARS const& vars) {
    return (is_affine_in(e1,vars) && is_constant_in(e2,vars)) || (is_constant_in(e1,vars) && is_affine_in(e2,vars)); }
template<class E1, class E2, class VARS> inline Bool _is_affine_in_impl(Div, E1 const& e1, E2 const& e2, VARS const& vars) {
    return (is_affine_in(e1,vars) && is_constant_in(e2,vars)); }
template<class E1, class E2, class VARS> inline Bool _is_affine_in_impl(Max, E1 const& e1, E2 const& e2, VARS const& vars) {
    return is_constant_in(e1,vars) && is_constant_in(e2,vars); }
template<class E1, class E2, class VARS> inline Bool _is_affine_in_impl(Min, E1 const& e1, E2 const& e2, VARS const& vars) {
    return is_constant_in(e1,vars) && is_constant_in(e2,vars); }
template<class E, class N, class VARS> inline Bool _is_affine_in_impl(Pow, E const& e1, N n2, VARS const& vars) {
    return n2 == 0 || (n2==1 && is_affine_in(e1,vars)) || is_constant_in(e1,vars); }

template<class E, class... OPS, class VARS> Bool _is_affine_in_impl(OperatorVariant<OPS...> ops, E const& e, VARS const& vars) {
    return ops.accept([&e,&vars](auto op){return _is_affine_in_impl(op,e, vars);}); }
template<class E1, class E2, class... OPS, class VARS> Bool _is_affine_in_impl(OperatorVariant<OPS...> ops, E1 const& e1, E2 const& e2, VARS const& vars) {
    return ops.accept([&e1,&e2,&vars](auto op){return _is_affine_in_impl(op,e1,e2, vars);}); }
template<class E, class VARS, class... OPS> Bool _is_affine_in_impl(OperatorVariant<OPS...> ops, E const& e, Int n, VARS const& vars) {
    return ops.accept([&e,n,&vars](auto op){return _is_affine_in_impl(op,e,n, vars);}); }
template<class OP, class E, class N, class VARS> inline Bool _is_affine_in_impl(GradedElementaryOperator ops, E const& e, N n, VARS const& vars) {
    return ops.accept([&e,n,&vars](auto op){return _is_affine_in_impl(op,e,n,vars);}); }

template<class T, class VARS> constexpr Bool is_affine_in(Constant<T> const&, VARS const& vars) {
    return true; }
template<class T, class VARS> constexpr Bool is_affine_in(Variable<T> const&, VARS const& vars) {
    return true; }
template<class OP, class A, class VARS> Bool is_affine_in(Symbolic<OP,A> const& s, VARS const& vars) {
    return _is_affine_in_impl(s._op,s._arg,vars); }
template<class OP, class A1, class A2, class VARS> Bool is_affine_in(Symbolic<OP,A1,A2> const& s, VARS const& vars) {
    return _is_affine_in_impl(s._op,s._arg1,s._arg2,vars); }
template<class OP, class A, class VARS> Bool is_affine_in(Symbolic<OP,A,Int> const& s, VARS const& vars) {
    return _is_affine_in_impl(s._op,s._arg,s._num,vars); }


template<class E, class VARS> Bool _is_polynomial_in_impl(Neg, E const& e, VARS const& vars) { return is_polynomial_in(e, vars); }
template<class OP, class E, class VARS> Bool _is_polynomial_in_impl(OP, E const& e, VARS const& vars) { return is_constant_in(e, vars); }
template<class E1, class E2, class VARS> Bool _is_polynomial_in_impl(Add, E1 const& e1, E2 const& e2, VARS const& vars) {
    return is_polynomial_in(e1, vars) && is_polynomial_in(e2, vars); }
template<class E1, class E2, class VARS> Bool _is_polynomial_in_impl(Sub, E1 const& e1, E2 const& e2, VARS const& vars) {
    return is_polynomial_in(e1, vars) && is_polynomial_in(e2, vars); }
template<class E1, class E2, class VARS> Bool _is_polynomial_in_impl(Mul, E1 const& e1, E2 const& e2, VARS const& vars) { return is_polynomial_in(e1, vars) && is_polynomial_in(e2, vars); }
template<class E1, class E2, class VARS> Bool _is_polynomial_in_impl(Div, E1 const& e1, E2 const& e2, VARS const& vars) { return is_polynomial_in(e1, vars) && is_constant_in(e2, vars); }
template<class E1, class E2, class VARS> Bool _is_polynomial_in_impl(Max, E1 const& e1, E2 const& e2, VARS const& vars) { return is_constant_in(e1, vars) && is_constant_in(e2, vars); }
template<class E1, class E2, class VARS> Bool _is_polynomial_in_impl(Min, E1 const& e1, E2 const& e2, VARS const& vars) { return is_constant_in(e1, vars) && is_constant_in(e2, vars); }
template<class E, class N, class VARS> Bool _is_polynomial_in_impl(Pow, E const& e, N const& n, VARS const& vars) { return n==0 || (n>0 && is_polynomial_in(e, vars)) || (n<0 && is_constant_in(e,vars)); }

template<class E, class VARS, class... OPS> Bool _is_polynomial_in_impl(OperatorVariant<OPS...> ops, E const& e, VARS const& vars) {
    return ops.accept([&e,&vars](auto op){return _is_polynomial_in_impl(op,e, vars);}); }
template<class E1, class E2, class VARS, class... OPS> Bool _is_polynomial_in_impl(OperatorVariant<OPS...> ops, E1 const& e1, E2 const& e2, VARS const& vars) {
    return ops.accept([&e1,&e2,&vars](auto op){return _is_polynomial_in_impl(op,e1,e2, vars);}); }
template<class E1, class E2, class VARS, class... OPS> Bool _is_polynomial_in_impl(OperatorVariant<OPS...> ops, E1 const& e1, Int n2, VARS const& vars) {
    return ops.accept([&e1,n2,&vars](auto op){return _is_polynomial_in_impl(op,e1,n2, vars);}); }


template<class T, class VARS> constexpr Bool is_polynomial_in(Constant<T> const&, VARS const& vars) {
    return true; }
template<class T, class VARS> constexpr Bool is_polynomial_in(Variable<T> const&, VARS const& vars) {
    return true; }
template<class OP, class A, class VARS> Bool is_polynomial_in(Symbolic<OP,A> const& s, VARS const& vars) {
    return _is_polynomial_in_impl(s._op,s._arg,vars); }
template<class OP, class A1, class A2, class VARS> Bool is_polynomial_in(Symbolic<OP,A1,A2> const& s, VARS const& vars) {
    return _is_polynomial_in_impl(s._op,s._arg1,s._arg2,vars); }
template<class OP, class A, class VARS> Bool is_polynomial_in(Symbolic<OP,A,Int> const& s, VARS const& vars) {
    return _is_polynomial_in_impl(s._op,s._arg,s._num,vars); }

}

} // namespace Ariadne
