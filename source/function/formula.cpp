/***************************************************************************
 *            formula.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "formula.hpp"

#include "../numeric/operators.tpl.hpp"

namespace Ariadne {

template<class Y> Formula<Y> Formula<Y>::_derivative(SizeType j) const
{
    const Formula<Y>& f = *this;
    switch(f.kind()) {
        case OperatorKind::NULLARY: return Formula<Y>::constant(0);
        case OperatorKind::COORDINATE: return Formula<Y>::constant(f.ind()==j?1:0);
        case OperatorKind::UNARY: return derivative(f.op(),f.arg(),derivative(f.arg(),j));
        case OperatorKind::BINARY: return derivative(f.op(),f.arg1(),derivative(f.arg1(),j),f.arg2(),derivative(f.arg2(),j));
        case OperatorKind::GRADED: return derivative(f.op(),f.arg(),derivative(f.arg(),j),f.num());
        default:
            ARIADNE_THROW(std::runtime_error,"derivative(Formula<Y>)",
                          "Cannot compute derivative of "<<f<<"\n");
    }
}

//! \brief Write to an output stream
template<class Y> OutputStream& Formula<Y>::_write(OutputStream& os) const
{
    const Formula<Y>& f = *this;
    switch(f.op()) {
        //case OperatorCode::CNST: return os << std::fixed << std::setprecision(4) << fptr->val;
        case OperatorCode::CNST:
            os << f.val(); return os;
            //if(f.val()==0.0) { return os << 0.0; } if(abs(f.val())<1e-4) { os << std::fixed << f.val(); } else { os << f.val(); } return os;
        case OperatorCode::IND:
            return os << "x" << f.ind();
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
                case OperatorKind::COMPARISON: return os << "(" << f.arg1() << symbol(f.op()) << f.arg2() << ")";
                case OperatorKind::SCALAR: return os << "(" << f.cnst() << symbol(f.op()) << f.arg() << ")";
                default: ARIADNE_FAIL_MSG("Cannot output formula with operator "<<f.op()<<" of kind "<<f.kind()<<"\n");
            }
    }
}

template class Formula<ApproximateNumber>;
template class Formula<ValidatedNumber>;
template class Formula<EffectiveNumber>;
template class Formula<ExactNumber>;

template class Formula<FloatDPApproximation>;
template class Formula<Real>;

template ApproximateNumber compute(OperatorCode op, const ApproximateNumber& x);
template ValidatedNumber compute(OperatorCode op, const ValidatedNumber& x);
template EffectiveNumber compute(OperatorCode op, const EffectiveNumber& x);
template ExactNumber compute(OperatorCode op, const ExactNumber& x);

} // namespace Ariadne
