/***************************************************************************
 *            formula.cc
 *
 *  Copyright 2008-16  Pieter Collins
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

#include "formula.h"

namespace Ariadne {

template<class Y> Formula<Y> Formula<Y>::_derivative(SizeType j) const
{
    const Formula<Y>& f = *this;
    switch(f.op()) {
        case OperatorCode::CNST:
            return Formula<Y>::constant(0);
        case OperatorCode::IND:
            if(f.ind()==j) { return Formula<Y>::constant(1); }
            else { return Formula<Y>::constant(0); }
        case OperatorCode::ADD:
            return derivative(f.arg1(),j)+derivative(f.arg2(),j);
        case OperatorCode::SUB:
            return derivative(f.arg1(),j)-derivative(f.arg2(),j);
        case OperatorCode::MUL:
            return f.arg1()*derivative(f.arg2(),j)+derivative(f.arg1(),j)*f.arg2();
        case OperatorCode::DIV:
            return derivative(f.arg1() * rec(f.arg2()),j);
        case OperatorCode::NEG:
            return  - derivative(f.arg(),j);
        case OperatorCode::REC:
            return  - derivative(f.arg(),j) * rec(sqr(f.arg()));
        case OperatorCode::SQR:
            return static_cast<Y>(2) * derivative(f.arg(),j) * f.arg();
        case OperatorCode::EXP:
            return derivative(f.arg(),j) * f.arg();
        case OperatorCode::LOG:
            return derivative(f.arg(),j) * rec(f.arg());
        case OperatorCode::SIN:
            return derivative(f.arg(),j) * cos(f.arg());
        case OperatorCode::COS:
            return -derivative(f.arg(),j) * sin(f.arg());
        case OperatorCode::TAN:
            return derivative(f.arg(),j) * (static_cast<Y>(1)-sqr(f.arg()));
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

template class Formula<Float64Approximation>;
template class Formula<Real>;

} // namespace Ariadne
