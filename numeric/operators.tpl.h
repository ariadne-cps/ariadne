/***************************************************************************
 *            operators.tcc
 *
 *  Copyright 2008-15 Pieter Collins
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

#include "utility/standard.h"

#include "numeric/operators.h"

namespace Ariadne {


// Declare functions and specialisations
template<class X> typename Logic<X>::Type compare(OperatorCode cmp, const X& x1, const X& x2);
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



template<class X> typename Logic<X>::Type compare(OperatorCode cmp, const X& x1, const X& x2) {
    switch(cmp) {
        case OperatorCode::GT: case OperatorCode::GEQ: return x1>x2;
        case OperatorCode::LT: case OperatorCode::LEQ: return x1<x2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate comparison "<<cmp<<" on real arguments.");
    }
}


template<class X1, class X2> X2 compute(OperatorCode op, const X1& x1, const X2& x2) {
    // FIXME: Should return ArithmeticType<X1,X2>
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


} // namespace Ariadne

