/***************************************************************************
 *            symbolic/templates.hpp
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

/*! \file symbolic/templates.hpp
 *  \brief
 */


#ifndef ARIADNE_EXPRESSION_TEMPLATES_HPP
#define ARIADNE_EXPRESSION_TEMPLATES_HPP

#include "../numeric/operators.hpp"

namespace Ariadne {

/************ Symbolic ************************************************/

template<class OP, class... AS> struct Symbolic;

template<class C> struct Symbolic<Cnst,C> {
    Cnst _op; C _val;
    explicit Symbolic(C c) : _op(), _val(c) { }
    Symbolic(Cnst o, C c) : _op(o), _val(c) { }
    template<class T> explicit operator T() const { return static_cast<T>(_val); }
    template<class... AS> auto operator() (AS... vals) const -> C { return _val; }
    friend OutputStream& operator<<(OutputStream& os, Symbolic expr) { return os << expr._val; }
};

template<class I> struct Symbolic<Var,I> {
    Var _op; I _ind;
    explicit Symbolic(I i) : _op(), _ind(i) { }
    Symbolic(Var o, I i) : _op(o), _ind(i) { }
    template<class AS> auto operator() (AS vals) const -> decltype(vals[_ind]) { return vals[_ind]; }
    friend OutputStream& operator<<(OutputStream& os, Symbolic expr) { return os << expr._ind; }
};

template<class O, class A> struct Symbolic<O,A> {
    O _op; A _arg;
    O op() const { return _op; } A arg() const { return _arg; }
    Symbolic(O o, A a) : _op(o), _arg(a) { }
//    template<class T> explicit operator T() const { return _op(static_cast<T>(_arg)); }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_arg(vals...))) {
        return _op(_arg(vals...)); }
    friend OutputStream& operator<<(OutputStream& os, Symbolic expr) {
        return os << expr._op.code() << "(" << expr._arg << ")"; }
};

template<class O, class A1, class A2> struct Symbolic<O,A1,A2> {
    O _op; A1 _arg1; A2 _arg2;
    O op() const { return _op; } A1 arg1() const { return _arg1; } A2 arg2() const { return _arg2; }
    Symbolic(O o, A1 a1, A2 a2) : _op(o), _arg1(a1), _arg2(a2) { }
//    template<class T> explicit operator T() const { return _op(static_cast<T>(_arg1),static_cast<T>(_arg2)); }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_arg1(vals...),_arg2(vals...))) {
        return _op(_arg1(vals...),_arg2(vals...)); }
    friend OutputStream& operator<<(OutputStream& os, Symbolic expr) {
        return os << expr._op.code() << "(" << expr._arg1 << "," << expr._arg2 << ")"; }
};

template<class A, class N> struct Symbolic<Pow,A,N> {
    typedef Pow O;
    O _op; A _arg; N _num;
    O op() const { return _op; } A arg() const { return _arg; } N num() const { return _num; }
    Symbolic(O o, A a, N n) : _op(o), _arg(a), _num(n) { }
//    template<class T> explicit operator T() const { return _op(static_cast<T>(_arg),_num); }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_arg(vals...),_num)) {
        return _op(_arg(vals...),_num); }
    friend OutputStream& operator<<(OutputStream& os, Symbolic expr) {
        return os << expr._op.code() << "(" << expr._arg << "," << expr._num << ")"; }
};

template<class A, class N> struct Symbolic<GradedElementaryOperator,A,N> {
    typedef GradedElementaryOperator O;
    O _op; A _arg; N _num;
    O op() const { return _op; } A arg() const { return _arg; } N num() const { return _num; }
    Symbolic(O o, A a, N n) : _op(o), _arg(a), _num(n) { }
    template<class T> explicit operator T() const { return _op(static_cast<T>(_arg),_num); }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_arg(vals...),_num)) {
        return _op(_arg(vals...),_num); }
    friend OutputStream& operator<<(OutputStream& os, Symbolic expr) {
        return os << expr._op.code() << "(" << expr._arg << "," << expr._num << ")"; }
};

template<class O, class A1, class A2, class A3> struct Symbolic<O,A1,A2,A3> {
    O _op; A1 _arg1; A2 _arg2; A3 _arg3;
    Symbolic(O o, A1 a1, A2 a2, A3 a3) : _op(o), _arg1(a1), _arg3(a3) { }
    template<class T> explicit operator T() const { return _op(static_cast<T>(_arg1),static_cast<T>(_arg2),static_cast<T>(_arg3)); }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_arg1(vals...),_arg2(vals...),_arg3(vals...))) {
        return _op(_arg1(vals...),_arg2(vals...)); }
    friend OutputStream& operator<<(OutputStream& os, Symbolic expr) {
        return os << expr._op.code() << "(" << expr._arg1 << "," << expr._arg2 << "," << expr._arg3 << ")"; }
};

template<class O, class... Args> Symbolic<O,Args...> make_expression_template(O o, Args... as) {
    return Symbolic<O,Args...>(o, as...); }

template<class O, class... Args> struct TemporaryExpression;

template<class O, class A> struct TemporaryExpression<O,A> {
    O _op; A const& _arg;
    TemporaryExpression(O o, A const& a) : _op(o), _arg(a) { }
};

template<class O, class A1, class A2> struct TemporaryExpression<O,A1,A2> {
    O _op; A1 const& _arg1; A2 const& _arg2;
    TemporaryExpression(O o, A1 const& a1, A2 const& a2) : _op(o), _arg1(a1), _arg2(a2) { }
};

template<class O, class A1, class A2, class A3> struct TemporaryExpression<O,A1,A2,A3> {
    O _op; A1 const& _arg1; A2 const& _arg2; A3 const& _arg3;
    TemporaryExpression(O o, A1 const& a1, A2 const& a2, A3 const& a3) : _op(o), _arg1(a1), _arg3(a3) { }
};


} // namespace Ariadne

#endif
