/***************************************************************************
 *            operators.cpp
 *
 *  Copyright 2008--17 Pieter Collins
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
double add(double,double);
double sub(double,double);
double mul(double,double);
double div(double,double);
}

#include "../utility/standard.hpp"

#include "../utility/string.hpp"
#include "../numeric/logical.hpp"
#include "../numeric/integer.hpp"
#include "../numeric/real.hpp"
#include "../numeric/operators.hpp"

namespace Ariadne {

template<class OP, class... AS> class Expression;

template<class Y> class Expression<Sgn,Y> : public LogicalInterface {
    Y _y;
  public:
    Expression(Sgn, Y y) : _y(y) { }
    virtual OutputStream& _write(OutputStream& os) const override;
    virtual LogicalValue _check(Effort eff) const override;
};

template<class Y> OutputStream& Expression<Sgn,Y>::_write(OutputStream& os) const {
    return os << "sgn(" << _y << ")";
}
template<class Y> LogicalValue Expression<Sgn,Y>::_check(Effort eff) const {
    return static_cast<LogicalValue>(check_sgn(_y,eff));
}
template class Expression<Sgn,Real>;


Kleenean Sgn::operator()(const Real& a) const {
    return Kleenean(LogicalHandle(std::make_shared<Expression<Sgn,Real>>(Sgn(),a)));
}

OutputStream& operator<<(OutputStream& os, const OperatorKind& knd) {
    switch(knd) {
        case OperatorKind::VARIABLE: return os << "VARIABLE";
        case OperatorKind::COORDINATE: return os << "COORDINATE";
        case OperatorKind::NULLARY: return os << "NULLARY";
        case OperatorKind::UNARY: return os << "UNARY";
        case OperatorKind::BINARY: return os << "BINARY";
        case OperatorKind::TERNARY: return os << "TERNARY";
        case OperatorKind::SCALAR: return os << "SCALAR";
        case OperatorKind::COMPARISON: return os << "COMPARISON";
        default: return os << "UNKNOWN";
    }
}

const char* name(const OperatorCode& op) {
    switch(op) {
        case OperatorCode::CNST: return "cnst"; break;
        case OperatorCode::VAR:  return "var"; break;
        case OperatorCode::IND:  return "ind"; break;
        case OperatorCode::POS:  return "pos"; break;
        case OperatorCode::NEG:  return "neg"; break;
        case OperatorCode::REC:  return "rec"; break;
        case OperatorCode::ADD:  return "add"; break;
        case OperatorCode::SUB:  return "sub"; break;
        case OperatorCode::MUL:  return "mul"; break;
        case OperatorCode::DIV:  return "div"; break;
        case OperatorCode::SADD:  return "sadd"; break;
        case OperatorCode::SSUB:  return "ssub"; break;
        case OperatorCode::SMUL:  return "smul"; break;
        case OperatorCode::SDIV:  return "sdiv"; break;
        case OperatorCode::POW:  return "pow"; break;
        case OperatorCode::NOT:  return "not"; break;
        case OperatorCode::AND:  return "and"; break;
        case OperatorCode::OR:   return "or"; break;
        case OperatorCode::XOR:  return "xor"; break;
        case OperatorCode::IMPL: return "impl"; break;
        case OperatorCode::ABS:  return "abs"; break;
        case OperatorCode::MAX:  return "max"; break;
        case OperatorCode::MIN:  return "min"; break;
        case OperatorCode::SQR:  return "sqr"; break;
        case OperatorCode::SQRT: return "sqrt"; break;
        case OperatorCode::EXP:  return "exp"; break;
        case OperatorCode::LOG:  return "log"; break;
        case OperatorCode::SIN:  return "sin"; break;
        case OperatorCode::COS:  return "cos"; break;
        case OperatorCode::TAN:  return "tan"; break;
        case OperatorCode::ASIN:  return "asin"; break;
        case OperatorCode::ACOS:  return "acos"; break;
        case OperatorCode::ATAN:  return "atan"; break;
        case OperatorCode::ITOR:  return "itor"; break;
        case OperatorCode::PULL: return "pull"; break;
        case OperatorCode::PUSH: return "push"; break;
        case OperatorCode::SGN:  return "sgn"; break;
        case OperatorCode::EQ:   return "eq"; break;
        case OperatorCode::NEQ:  return "neq"; break;
        case OperatorCode::GEQ:  return "geq"; break;
        case OperatorCode::LEQ:  return "leq"; break;
        case OperatorCode::GT:   return "gt"; break;
        case OperatorCode::LT:   return "lt"; break;
        case OperatorCode::SUBS:   return "subs"; break;
        case OperatorCode::HLF:  return "hlf"; break;
        default: return "UNKNOWN";
    }
}

const char* symbol(const OperatorCode& op) {
    switch(op) {
        case OperatorCode::POS:  return "+"; break;
        case OperatorCode::NEG:  return "-"; break;
        case OperatorCode::ADD:  return "+"; break;
        case OperatorCode::SUB:  return "-"; break;
        case OperatorCode::MUL:  return "*"; break;
        case OperatorCode::DIV:  return "/"; break;
        case OperatorCode::POW:  return "^"; break;
        case OperatorCode::SADD: return "+"; break;
        case OperatorCode::SMUL: return "*"; break;
        case OperatorCode::NOT:  return "!"; break;
        case OperatorCode::AND:  return "&"; break;
        case OperatorCode::OR:   return "|"; break;
        case OperatorCode::EQ:  return "=="; break;
        case OperatorCode::NEQ: return "!="; break;
        case OperatorCode::LEQ: return "<="; break;
        case OperatorCode::GEQ: return ">="; break;
        case OperatorCode::LT:  return "<"; break;
        case OperatorCode::GT:  return ">"; break;
        default: return "???";
    }
}

OutputStream& operator<<(OutputStream& os, const OperatorCode& op) {
    return os << name(op);
}

OperatorKind kind(OperatorCode op) {
    switch(op) {
        case OperatorCode::CNST:
            return OperatorKind::NULLARY;
        case OperatorCode::IND:
            return OperatorKind::COORDINATE;
        case OperatorCode::VAR:
            return OperatorKind::VARIABLE;
        case OperatorCode::ADD: case OperatorCode::SUB: case OperatorCode::MUL: case OperatorCode::DIV:
        case OperatorCode::MAX: case OperatorCode::MIN:
            return OperatorKind::BINARY;
        case OperatorCode::SADD: case OperatorCode::SSUB: case OperatorCode::SMUL: case OperatorCode::SDIV:
            return OperatorKind::SCALAR;
        case OperatorCode::POS: case OperatorCode::NEG: case OperatorCode::REC: case OperatorCode::SQR:
        case OperatorCode::SQRT: case OperatorCode::EXP: case OperatorCode::LOG:
        case OperatorCode::SIN: case OperatorCode::COS: case OperatorCode::TAN: case OperatorCode::ATAN:
        case OperatorCode::ABS:
            return OperatorKind::UNARY;
        case OperatorCode::POW:
            return OperatorKind::GRADED;
        case OperatorCode::AND: case OperatorCode::OR:
            return OperatorKind::BINARY;
        case OperatorCode::NOT:
            return OperatorKind::UNARY;
        case OperatorCode::EQ: case OperatorCode::NEQ: case OperatorCode::LEQ: case OperatorCode::GEQ: case OperatorCode::LT: case OperatorCode::GT:
            return OperatorKind::COMPARISON;
        default:
            ARIADNE_FAIL_MSG("Cannot deduce kind of operator "<<op<<"\n");
    }
}



template<> Boolean compare(OperatorCode cmp, const String& s1, const String& s2) {
    switch(cmp) {
        case OperatorCode::EQ:  return s1==s2;
        case OperatorCode::NEQ: return s1!=s2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate comparison "<<cmp<<" on string arguments.");
    }
}

template<> Boolean compare(OperatorCode cmp, const Integer& z1, const Integer& z2) {
    switch(cmp) {
        case OperatorCode::EQ:  return z1==z2;
        case OperatorCode::NEQ: return z1!=z2;
        case OperatorCode::LEQ: return z1<=z2;
        case OperatorCode::GEQ: return z1>=z2;
        case OperatorCode::LT:  return z1< z2;
        case OperatorCode::GT:  return z1> z2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate comparison "<<cmp<<" on integer arguments.");
    }
}


template<> Boolean compute(OperatorCode op, const Boolean& b1, const Boolean& b2) {
    switch(op) {
        case OperatorCode::AND: return b1 && b2;
        case OperatorCode::OR: return b1 || b2;
        case OperatorCode::XOR: return b1 ^ b2;
        case OperatorCode::IMPL: return !b1 || b2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on a boolean arguments.");
    }
}

template<> Kleenean compute(OperatorCode op, const Kleenean& b1, const Kleenean& b2) {
    switch(op) {
        case OperatorCode::AND: return b1 && b2;
        case OperatorCode::OR: return b1 || b2;
        case OperatorCode::XOR: return b1 ^ b2;
        case OperatorCode::IMPL: return !b1 || b2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on two boolean arguments.");
    }
}

template<> String compute(OperatorCode op, const String& s1, const String& s2) {
    switch(op) {
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on two string arguments.");
    }
}

template<> Integer compute(OperatorCode op, const Integer& x1, const Integer& x2) {
    switch(op) {
        case OperatorCode::ADD: return x1+x2;
        case OperatorCode::SUB: return x1-x2;
        case OperatorCode::MUL: return x1*x2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on two integer arguments.");
    }
}


template<> Boolean compute(OperatorCode op, const Boolean& b) {
    switch(op) {
        case OperatorCode::NOT: return !b;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on a boolean arguments.");
    }
}

template<> Kleenean compute(OperatorCode op, const Kleenean& b) {
    switch(op) {
        case OperatorCode::NOT: return !b;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on a boolean arguments.");
    }
}

template<> String compute(OperatorCode op, const String& s) {
    switch(op) {
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on one string argument.");
    }
}

template<> Integer compute(OperatorCode op, const Integer& z) {
    switch(op) {
        case OperatorCode::POS: return +z;
        case OperatorCode::NEG: return -z;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on one integer argument.");
    }
}


template<> String compute(OperatorCode op, const String& s, Int n) {
    switch(op) {
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on one string argument and an integer.");
    }
}

template<> Integer compute(OperatorCode op, const Integer& z, Int n) {
    switch(op) {
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on one integer argument and a builtin.");
    }
}

} // namespace Ariadne

