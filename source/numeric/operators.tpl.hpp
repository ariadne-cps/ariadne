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
template<> Real compute(OperatorCode op, const Real& x);
template<> Boolean compute(OperatorCode op, const Boolean& b);
template<> Kleenean compute(OperatorCode op, const Kleenean& b);
template<> String compute(OperatorCode op, const String& s);
template<> Integer compute(OperatorCode op, const Integer& z);

template<class X> X compute(OperatorCode op, const X& x, Int n);
template<> Real compute(OperatorCode op, const Real& x, Int n);
template<> String compute(OperatorCode op, const String& s, Int n);
template<> Integer compute(OperatorCode op, const Integer& z, Int n);


template<class L> decltype(not declval<L>()) UnaryLogicalOperator::operator()(L&& l) const {
    switch(this->_code) {
        case Code::NOT: return not std::forward<L>(l);
    }
}

template<class X> X UnaryElementaryOperator::operator()(X&& x) const {
    switch(this->_code) {
        case Code::POS: return pos(std::forward<X>(x));
        case Code::NEG: return neg(std::forward<X>(x));
        case Code::SQR: return sqr(std::forward<X>(x));
        case Code::REC: return rec(std::forward<X>(x));
        case Code::SQRT: return sqrt(std::forward<X>(x));
        case Code::EXP: return exp(std::forward<X>(x));
        case Code::LOG: return log(std::forward<X>(x));
        case Code::SIN: return sin(std::forward<X>(x));
        case Code::COS: return cos(std::forward<X>(x));
        case Code::TAN: return tan(std::forward<X>(x));
        case Code::ASIN: return asin(std::forward<X>(x));
        case Code::ACOS: return acos(std::forward<X>(x));
        case Code::ATAN: return atan(std::forward<X>(x));
        default: assert(false);
    }
}
template<class X> X compute(OperatorCode op, X const& x) {
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
//        case OperatorCode::ASIN: return asin(x);
//        case OperatorCode::ACOS: return acos(x);
        case OperatorCode::ATAN: return atan(x);
        default: assert(false);
    }
}


template<class X1, class X2> X2 compute(OperatorCode op, const X1& x1, const X2& x2) {
    switch(op) {
        case OperatorCode::ADD: case OperatorCode::SADD: return x1+x2;
        case OperatorCode::SUB: case OperatorCode::SSUB: return x1-x2;
        case OperatorCode::MUL: case OperatorCode::SMUL: return x1*x2;
        case OperatorCode::DIV: case OperatorCode::SDIV: return x1/x2;
//        case OperatorCode::MAX: return max(x1,x2);
//        case OperatorCode::MIN: return min(x1,x2);
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on two real arguments.");
    }
}

template<class X> X UnaryLatticeOperator::operator()(X&& x) const {
    switch(this->_code) {
        case Code::ABS: return abs(std::forward<X>(x));
    }
};

template<class X> LogicType<X> BinaryComparisonOperator::operator()(const X& x1, const X& x2) const {
    switch(this->_code) {
        // FIXME: Should be able to use == and != directly
        case Code::EQ:  return x1<=x2 && x2<=x2;
        case Code::NEQ: return !(x1<=x2 && x2<=x2);
        case Code::GT:  return x1> x2;
        case Code::GEQ: return x1>=x2;
        case Code::LT:  return x1< x2;
        case Code::LEQ: return x1<=x2;
        default: assert(false);
    }
}

template<class X> LogicType<X> compare(OperatorCode cmp, const X& x1, const X& x2) {
    return BinaryComparisonOperator(cmp).operator()(x1,x2);
}

template<class L> L BinaryLogicalOperator::operator()(L&& l1, L&& l2) const {
    switch(this->_code) {
        case Code::AND: return std::forward<L>(l1) && std::forward<L>(l2);
        case Code::OR: return std::forward<L>(l1) || std::forward<L>(l2);
        case Code::XOR: return std::forward<L>(l1) ^ std::forward<L>(l2);
        case Code::IMPL: return !std::forward<L>(l1) || std::forward<L>(l2);
        default: assert(false);
    }
}

template<class X1, class X2> ArithmeticType<X1,X2> BinaryArithmeticOperator::operator()(X1&& x1, X2&& x2) const {
    switch(this->_code) {
        case Code::ADD: return std::forward<X1>(x1)+std::forward<X2>(x2);
        case Code::SUB: return std::forward<X1>(x1)-std::forward<X2>(x2);
        case Code::MUL: return std::forward<X1>(x1)*std::forward<X2>(x2);
        case Code::DIV: return std::forward<X1>(x1)/std::forward<X2>(x2);
        default: assert(false);
    }
}



template<class X> X compute(OperatorCode op, const X& x, Int n) {
    switch(op) {
        case OperatorCode::POW: return pow(x,n);
        case OperatorCode::ROOT:
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on one real argument and an integer.");
    }
}

template<class X> X compute_derivative(OperatorCode op, const X& x) {
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

template<class X, class D> D compute_derivative(OperatorCode op, const X& x, const D& dx) {
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

template<class X, class D> D compute_derivative(OperatorCode op, const X& x1, const D& dx1, const X& x2, const D& dx2) {
    switch(op) {
        case OperatorCode::ADD: return Add().derivative(x1,dx1,x2,dx2);
        case OperatorCode::SUB: return Sub().derivative(x1,dx1,x2,dx2);
        case OperatorCode::MUL: return Mul().derivative(x1,dx1,x2,dx2);
        case OperatorCode::DIV: return Div().derivative(x1,dx1,x2,dx2);
        default: ARIADNE_FAIL_MSG("Cannot differentiate operator "<<op<<" on two real arguments.");
    }
}

template<class X, class D> D compute_derivative(OperatorCode op, const X& x, const D& dx, Int n) {
    switch(op) {
        case OperatorCode::POW: return Pow().derivative(x,dx,n);
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on one real argument and an integer.");
    }
}


} // namespace Ariadne

