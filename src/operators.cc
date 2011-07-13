/***************************************************************************
 *            operators.cc
 *
 *  Copyright 2008-11 Pieter Collins
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


#include "operators.h"
#include <include/operators.h>

namespace Ariadne {

OutputStream& operator<<(OutputStream& os, const OperatorKind& knd) {
    switch(knd) {
        case VARIABLE: return os << "VARIABLE";
        case COORDINATE: return os << "VARIABLE";
        case NULLARY: return os << "NULLARY";
        case UNARY: return os << "UNARY";
        case BINARY: return os << "BINARY";
        case TERNARY: return os << "TERNARY";
        case SCALAR: return os << "SCALAR";
        case COMPARISON: return os << "COMPARISON";
    }
}

const char* name(const OperatorCode& op) {
    switch(op) {
        case CNST: return "cnst"; break;
        case VAR:  return "var"; break;
        case IND:  return "ind"; break;
        case POS:  return "pos"; break;
        case NEG:  return "neg"; break;
        case REC:  return "rec"; break;
        case ADD:  return "add"; break;
        case SUB:  return "sub"; break;
        case MUL:  return "mul"; break;
        case DIV:  return "div"; break;
        case POW:  return "pow"; break;
        case NOT:  return "not"; break;
        case AND:  return "and"; break;
        case OR:   return "or"; break;
        case XOR:  return "xor"; break;
        case IMPL: return "impl"; break;
        case ABS:  return "abs"; break;
        case MAX:  return "max"; break;
        case MIN:  return "min"; break;
        case SQR:  return "sqr"; break;
        case SQRT: return "sqrt"; break;
        case EXP:  return "exp"; break;
        case LOG:  return "log"; break;
        case SIN:  return "sin"; break;
        case COS:  return "cos"; break;
        case TAN:  return "tan"; break;
        case ASIN:  return "asin"; break;
        case ACOS:  return "acos"; break;
        case ATAN:  return "atan"; break;
        case ITOR:  return "itor"; break;
        case PULL: return "pull"; break;
        case PUSH: return "push"; break;
        case SGN:  return "sgn"; break;
        case EQ:   return "eq"; break;
        case NEQ:  return "neq"; break;
        case GEQ:  return "geq"; break;
        case LEQ:  return "leq"; break;
        case GT:   return "lt"; break;
        case LT:   return "gt"; break;
        case SUBS:   return "subs"; break;
        default: return "UNKNOWN";
    }
}

const char* symbol(const OperatorCode& op) {
    switch(op) {
        case POS:  return "+"; break;
        case NEG:  return "-"; break;
        case ADD:  return "+"; break;
        case SUB:  return "-"; break;
        case MUL:  return "*"; break;
        case DIV:  return "/"; break;
        case POW:  return "^"; break;
        case SADD: return "+"; break;
        case SMUL: return "*"; break;
        case NOT:  return "!"; break;
        case AND:  return "&"; break;
        case OR:   return "|"; break;
        case EQ:  return "=="; break;
        case NEQ: return "!="; break;
        case LEQ: return "<="; break;
        case GEQ: return ">="; break;
        case LT:  return "<"; break;
        case GT:  return ">"; break;
        default: return "???";
    }
}

OutputStream& operator<<(OutputStream& os, const OperatorCode& op) {
    return os << name(op);
}

OperatorKind kind(OperatorCode op) {
    switch(op) {
        case CNST:
            return NULLARY;
        case IND:
            return COORDINATE;
        case VAR:
            return VARIABLE;
        case ADD: case SUB: case MUL: case DIV:
            return BINARY;
        case POS: case NEG: case REC: case SQR:
        case SQRT: case EXP: case LOG:
        case SIN: case COS: case TAN: case ATAN:
            return UNARY;
        case POW:
            return SCALAR;
        case AND: case OR:
            return BINARY;
        case NOT:
            return UNARY;
        case EQ: case NEQ: case LEQ: case GEQ: case LT: case GT:
            return COMPARISON;
        default:
            ARIADNE_FAIL_MSG("Cannot deduce kind of operator "<<op<<"\n");;
    }
}

} // namespace Ariadne

