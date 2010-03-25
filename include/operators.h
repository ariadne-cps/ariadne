/***************************************************************************
 *            operators.h
 *
 *  Copyright 2008-9  Pieter Collins
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


/*! \file operators.h
 *  \brief Numerical operator classes
 */

#ifndef ARIADNE_OPERATORS_H
#define ARIADNE_OPERATORS_H

#include <cstdarg>
#include <cassert>
#include <iostream>
#include <string>

#include "tribool.h"
#include "numeric.h"

namespace Ariadne {




enum Operator {
    CNST,  // A constant value
    VAR,   // A named variable
    IND,   // A numbered index
    ADD,   // Addition
    SUB,   // Subtraction
    MUL,   // Multiplication
    DIV,   // Division
    POW,   // Integer power
    POS,   // Unary plus
    NEG,   // Unary negation
    REC,   // Reciprocal
    SQR,   // Square
    SQRT,  // Square root
    EXP,   // Natural exponential
    LOG,   // Natural logarithm
    SIN,   // Sine
    COS,   // Cosine
    TAN,   // Tangent
    ASIN,   // ArcSine
    ACOS,   // ArcCosine
    ATAN,   // ArcTangent
    ABS,   // Absolute value
    MAX,   // Maximum
    MIN,   // Minimum
    NOT,   // Logical not
    AND,   // Logical and
    OR,    // Logical or
    XOR,   // Logical exclusive or
    IMPL,  // Logical implication
    ITOR,   // Conversion of Int to Real
    EQ=-1,    // Equal
    NEQ=-2,   // Not equal
    GEQ=-3,   // Greater or equal
    LEQ=-4,   // Less or equal
    GT=-5,    // Greater than
    LT=-6,    // Less than
    SGN=-7,   // Compare with zero
    SUBS=-8,   // Compare as a subset
    DISJ=-9    // Compare disjointness
};

std::ostream& operator<<(std::ostream&, const Operator& op);


inline const char* symbol(const Operator& op) {
    switch(op) {
        case POS:  return "+"; break;
        case NEG:  return "-"; break;
        case ADD:  return "+"; break;
        case SUB:  return "-"; break;
        case MUL:  return "*"; break;
        case DIV:  return "/"; break;
        case POW:  return "^"; break;
        case NOT:  return "!"; break;
        case AND:  return "&"; break;
        case OR:   return "|"; break;
        case EQ:  return "=="; break;
        case NEQ: return "!="; break;
        case LEQ: return "<="; break;
        case GEQ: return ">="; break;
        case LT:  return "<"; break;
        case GT:  return ">"; break;
        default: assert(false);
    }
}

inline const char* name(const Operator& op) {
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
        case ITOR:  return "itor"; break;
        case SGN:  return "sgn"; break;
        case EQ:   return "eq"; break;
        case NEQ:  return "neq"; break;
        case GEQ:  return "geq"; break;
        case LEQ:  return "leq"; break;
        case GT:   return "lt"; break;
        case LT:   return "gt"; break;
        case SUBS:   return "subs"; break;
        default: assert(false);
    }
}

inline std::ostream& operator<<(std::ostream& os, const Operator& op) {
    return os << name(op);
}

struct GtrZero {}; struct LessZero {};

struct Gtr {
    tribool operator()(const Float& x1, const Float& x2) const {
        return (x1==x2) ? indeterminate : tribool(x1>x2); }
    tribool operator()(const Interval& x1, const Interval& x2) const {
        if(x1.lower()>x2.upper()) { return true; } else if(x1.upper()<x2.lower()) { return false; } else { return indeterminate; } }
};

struct Less {
    tribool operator()(const Float& x1, const Float& x2) const {
        return (x1==x2) ? indeterminate : tribool(x1<x2); }
    tribool operator()(const Interval& x1, const Interval& x2) const {
        if(x1.lower()>x2.upper()) { return false; } else if(x1.upper()<x2.lower()) { return true; } else { return indeterminate; } }
};

struct Equal {
    template<class T1, class T2> bool operator()(const T1& a1, const T2& a2) const { return a1 == a2; }
};

struct And {
    template<class T> T operator()(const T& a1, const T& a2) const { return a1 && a2; }
    Operator code() const { return AND; } };
struct Or {
    template<class T> T operator()(const T& a1, const T& a2) const { return a1 || a2; }
    Operator code() const { return OR; } };
struct Not {
    template<class T> T operator()(const T& a) const { return !a; }
    Operator code() const { return NOT; } };

struct Add {
    template<class T> T operator()(const T& a1, const T& a2) const { return a1+a2; }
    Operator code() const { return ADD; } };
struct Sub {
    template<class T> T operator()(const T& a1, const T& a2) const { return a1-a2; }
    Operator code() const { return SUB; } };
struct Mul {
    template<class T> T operator()(const T& a1, const T& a2) const { return a1*a2; }
    Operator code() const { return MUL; } };
struct Div {
    template<class T> T operator()(const T& a1, const T& a2) const { return a1/a2; }
    Operator code() const { return DIV; } };

struct Pow {
    template<class T, class N> T operator()(const T& a, const N& n) const { return pow(a,n); }
    Operator code() const { return DIV; } };
//struct Pow { Pow(int n) : n(n) { } template<class T> T operator()(const T& a) const { return Ariadne::pow(a,n); } int n; };

struct Neg {
    template<class T> T operator()(const T& a) const { return neg(a); }
    Operator code() const { return NEG; } };
struct Rec {
    template<class T> T operator()(const T& a) const { return rec(a); }
    Operator code() const { return REC; } };
struct Sqr {
    template<class T> T operator()(const T& a) const { return sqr(a); }
    Operator code() const { return SQR; } };
struct Sqrt {
    template<class T> T operator()(const T& a) const { return sqrt(a); }
    Operator code() const { return SQRT; } };

struct Exp {
    template<class T> T operator()(const T& a) const { return exp(a); }
    Operator code() const { return EXP; } };
struct Log {
    template<class T> T operator()(const T& a) const { return log(a); }
    Operator code() const { return LOG; } };
struct Sin {
    template<class T> T operator()(const T& a) const { return sin(a); }
    Operator code() const { return SIN; } };
struct Cos {
    template<class T> T operator()(const T& a) const { return cos(a); }
    Operator code() const { return COS; } };
struct Tan {
    template<class T> T operator()(const T& a) const { return tan(a); }
    Operator code() const { return TAN; } };


inline std::ostream& operator<<(std::ostream& os, const Less& v) { return os << "<="; }
inline std::ostream& operator<<(std::ostream& os, const Gtr& v) { return os << ">="; }
inline std::ostream& operator<<(std::ostream& os, const Equal& v) { return os << "=="; }

inline std::ostream& operator<<(std::ostream& os, const And& v) { return os << "&&"; }
inline std::ostream& operator<<(std::ostream& os, const Or& v) { return os << "||"; }
inline std::ostream& operator<<(std::ostream& os, const Not& v) { return os << "!"; }

inline std::ostream& operator<<(std::ostream& os, const Add& v) { return os << "+"; }
inline std::ostream& operator<<(std::ostream& os, const Sub& v) { return os << "-"; }
inline std::ostream& operator<<(std::ostream& os, const Mul& v) { return os << "*"; }
inline std::ostream& operator<<(std::ostream& os, const Div& v) { return os << "/"; }

inline std::ostream& operator<<(std::ostream& os, const Neg& op) { return os << "neg"; }
inline std::ostream& operator<<(std::ostream& os, const Rec& op) { return os << "rec"; }
inline std::ostream& operator<<(std::ostream& os, const Pow& op) { return os << "pow"; }
inline std::ostream& operator<<(std::ostream& os, const Sqr& op) { return os << "sqr"; }
inline std::ostream& operator<<(std::ostream& os, const Sqrt& op) { return os << "sqrt"; }

inline std::ostream& operator<<(std::ostream& os, const Exp& op) { return os << "exp"; }
inline std::ostream& operator<<(std::ostream& os, const Log& op) { return os << "log"; }
inline std::ostream& operator<<(std::ostream& os, const Sin& op) { return os << "sin"; }
inline std::ostream& operator<<(std::ostream& os, const Cos& op) { return os << "cos"; }
inline std::ostream& operator<<(std::ostream& os, const Tan& op) { return os << "tan"; }




} // namespace Ariadne

#endif // ARIADNE_OPERATORS_H
