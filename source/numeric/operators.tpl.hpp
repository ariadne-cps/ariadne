/***************************************************************************
 *            operators.tcc
 *
 *  Copyright 2008-17 Pieter Collins
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

#include "../numeric/operators.hpp"

namespace Ariadne {


// Declare functions and specialisations
template<class X> LogicType<X> compare(OperatorCode cmp, const X& x1, const X& x2);
template<> Boolean compare(OperatorCode cmp, const String& s1, const String& s2);
template<> Boolean compare(OperatorCode cmp, const Integer& z1, const Integer& z2);

template<class X1, class X2> X2 compute(OperatorCode op, const X1& x1, const X2& x2);
template<> Boolean compute(OperatorCode op, const Boolean& b1, const Boolean& b2);
template<> Kleenean compute(OperatorCode op, const Kleenean& b1, const Kleenean& b2);
template<> String compute(OperatorCode op, const String& s1, const String& s2);
template<> Integer compute(OperatorCode op, const Integer& x1, const Integer& x2);

template<class X> X compute(OperatorCode op, const X& x);
template<> Boolean compute(OperatorCode op, const Boolean& b);
template<> Kleenean compute(OperatorCode op, const Kleenean& b);
template<> String compute(OperatorCode op, const String& s);
template<> Integer compute(OperatorCode op, const Integer& z);

template<class X> X compute(OperatorCode op, const X& x, Int n);
template<> String compute(OperatorCode op, const String& s, Int n);
template<> Integer compute(OperatorCode op, const Integer& z, Int n);



template<class X> LogicType<X> compare(OperatorCode cmp, const X& x1, const X& x2) {
    switch(cmp) {
        case OperatorCode::GT: case OperatorCode::GEQ: return x1>x2;
        case OperatorCode::LT: case OperatorCode::LEQ: return x1<x2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate comparison "<<cmp<<" on real arguments.");
    }
}


template<class X1, class X2> X2 compute(OperatorCode op, const X1& x1, const X2& x2) {
    switch(op) {
        case OperatorCode::ADD: case OperatorCode::SADD: return x1+x2;
        case OperatorCode::SUB: case OperatorCode::SSUB: return x1-x2;
        case OperatorCode::MUL: case OperatorCode::SMUL: return x1*x2;
        case OperatorCode::DIV: case OperatorCode::SDIV: return x1/x2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on two real arguments.");
    }
}


template<class X> X compute(OperatorCode op, const X& x) {
    switch(op) {
        case OperatorCode::POS: return pos(x);
        case OperatorCode::NEG: return neg(x);
        case OperatorCode::SQR: return sqr(x);
        case OperatorCode::REC: return rec(x);
        case OperatorCode::SQRT: return sqrt(x);
        case OperatorCode::EXP: return exp(x);
        case OperatorCode::LOG: return log(x);
        case OperatorCode::SIN: return sin(x);
        case OperatorCode::COS: return cos(x);
        case OperatorCode::TAN: return tan(x);
        case OperatorCode::ATAN: return atan(x);
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on one real argument.");
    }
}


template<class X> X compute(OperatorCode op, const X& x, Int n) {
    switch(op) {
        case OperatorCode::POW: return pow(x,n);
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on one real argument and an integer.");
    }
}

template<class X> X derivative(OperatorCode op, const X& x) {
    switch(op) {
        case OperatorCode::POS: return Pos().derivative(x);
        case OperatorCode::NEG: return Neg().derivative(x);
        case OperatorCode::SQR: return Sqr().derivative(x);
        case OperatorCode::REC: return Rec().derivative(x);
        case OperatorCode::SQRT: return Sqrt().derivative(x);
        case OperatorCode::EXP: return Exp().derivative(x);
        case OperatorCode::LOG: return Log().derivative(x);
        case OperatorCode::SIN: return Sin().derivative(x);
        case OperatorCode::COS: return Cos().derivative(x);
        case OperatorCode::TAN: return Tan().derivative(x);
        case OperatorCode::ATAN: return Atan().derivative(x);
        default: ARIADNE_FAIL_MSG("Cannot differentiate operator "<<op<<" on one real argument.");
    }
}

template<class X, class D> D derivative(OperatorCode op, const X& x, const D& dx) {
    switch(op) {
        case OperatorCode::POS: return Pos().derivative(x,dx);
        case OperatorCode::NEG: return Neg().derivative(x,dx);
        case OperatorCode::SQR: return Sqr().derivative(x,dx);
        case OperatorCode::REC: return Rec().derivative(x,dx);
        case OperatorCode::SQRT: return Sqrt().derivative(x,dx);
        case OperatorCode::EXP: return Exp().derivative(x,dx);
        case OperatorCode::LOG: return Log().derivative(x,dx);
        case OperatorCode::SIN: return Sin().derivative(x,dx);
        case OperatorCode::COS: return Cos().derivative(x,dx);
        case OperatorCode::TAN: return Tan().derivative(x,dx);
        case OperatorCode::ATAN: return Atan().derivative(x,dx);
        default: ARIADNE_FAIL_MSG("Cannot differentiate operator "<<op<<" on one real argument.");
    }
}

template<class X, class D> D derivative(OperatorCode op, const X& x1, const D& dx1, const X& x2, const D& dx2) {
    switch(op) {
        case OperatorCode::ADD: return Add().derivative(x1,dx1,x2,dx2);
        case OperatorCode::SUB: return Sub().derivative(x1,dx1,x2,dx2);
        case OperatorCode::MUL: return Mul().derivative(x1,dx1,x2,dx2);
        case OperatorCode::DIV: return Div().derivative(x1,dx1,x2,dx2);
        default: ARIADNE_FAIL_MSG("Cannot differentiate operator "<<op<<" on two real arguments.");
    }
}

template<class X, class D> D derivative(OperatorCode op, const X& x, const D& dx, Int n) {
    switch(op) {
        case OperatorCode::POW: return Pow().derivative(x,dx,n);
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on one real argument and an integer.");
    }
}


} // namespace Ariadne

