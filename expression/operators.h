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

#include "utility/pointer.h"

#include "utility/tribool.h"
#include "utility/string.h"
#include "numeric/logical.h"
#include "numeric/numeric.h"

namespace Ariadne {

typedef void Void;
typedef bool Bool;
typedef unsigned int Nat;
typedef int Int;

typedef std::ostream OutputStream;

class String;
class Boolean;
class Tribool;

enum class OperatorKind : char {
    VARIABLE,
    COORDINATE,
    NULLARY,
    UNARY,
    BINARY,
    TERNARY,
    SCALAR,
    COMPARISON
};

OutputStream& operator<<(OutputStream& os, const OperatorKind& knd);

enum class OperatorCode : char {
    CNST,  // A constant value
    VAR,   // A named variable
    IND,   // A numbered index
    ADD,   // Addition
    SUB,   // Subtraction
    MUL,   // Multiplication
    DIV,   // Division
    FMA,   // Fused multiply-and-add
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
    SADD,  // Addition of a scalar constant
    SMUL,  // Multiplication by a scalar constant
    SFMA,   // Fused multiply-and-add with scalar constants
    ABS,   // Absolute value
    MAX,   // Maximum
    MIN,   // Minimum
    NOT,   // Logical not
    AND,   // Logical and
    OR,    // Logical or
    XOR,   // Logical exclusive or
    IMPL,  // Logical implication
    ITOR,   // Conversion of Int to Real
    PUSH,
    PULL,
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

OutputStream& operator<<(OutputStream&, const OperatorCode& op);
const char* name(const OperatorCode& op);
const char* symbol(const OperatorCode& op);

OperatorKind kind(OperatorCode op);

class Operator {
    OperatorCode _code;
    OperatorKind _kind;
  public:
    Operator(OperatorCode c, OperatorKind k) : _code(c), _kind(k) { }
    Operator(OperatorCode c) : _code(c), _kind(Ariadne::kind(_code)) { }
    template<class OP> Operator(const OP& op) : _code(op.code()), _kind(op.kind()) { }

    operator OperatorCode () const { return _code; }
    OperatorCode code() const { return _code; }
    OperatorKind kind() const { return _kind; }
};

template<class R, class A> inline R compare(OperatorCode op, const A& x1, const A& x2) {
    switch(op) {
        case OperatorCode::EQ: return static_cast<R>(x1==x2);
        case OperatorCode::NEQ: return static_cast<R>(x1!=x2);
        case OperatorCode::LEQ: return x1<=x2;
        case OperatorCode::GEQ: return x1>=x2;
        case OperatorCode::LT: return x1< x2;
        case OperatorCode::GT: return x1> x2;
        default: ARIADNE_FAIL_MSG("Cannot compute comparison "<<op<<" on "<<x1<<" and "<<x2);
    }
}
template<> inline Boolean compare(OperatorCode op, const String& x1, const String& x2) {
    switch(op) {
        case OperatorCode::EQ: return x1==x2;
        case OperatorCode::NEQ: return x1!=x2;
        default: ARIADNE_FAIL_MSG("Cannot compute comparison "<<op<<" on "<<x1<<" and "<<x2);
    }
}
template<class X> inline X compute(OperatorCode op, const X& x) {
    switch(op) {
        case OperatorCode::POS: return +x;
        case OperatorCode::NEG: return -x;
        case OperatorCode::REC: return 1/x;
        case OperatorCode::SQR: return sqr(x);
        case OperatorCode::SQRT: return sqrt(x);
        case OperatorCode::EXP: return exp(x);
        case OperatorCode::LOG: return log(x);
        case OperatorCode::SIN: return sin(x);
        case OperatorCode::COS: return cos(x);
        case OperatorCode::TAN: return tan(x);
//        case OperatorCode::ABS: return abs(x);
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x);
    }
}
template<> inline Boolean compute(OperatorCode op, const Boolean& x) {
    switch(op) {
        case OperatorCode::NOT: return !x;
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x);
    }
}
template<> inline Tribool compute(OperatorCode op, const Tribool& x) {
    switch(op) {
        case OperatorCode::NOT: return !x;
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x);
    }
}
template<> inline Integer compute(OperatorCode op, const Integer& x) {
    switch(op) {
        case OperatorCode::POS: return +x;
        case OperatorCode::NEG: return -x;
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x);
    }
}
template<> inline String compute(OperatorCode op, const String& x) {
    switch(op) {
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x);
    }
}
template<class X> inline X compute(OperatorCode op, const X& x1, const X& x2) {
    switch(op) {
        case OperatorCode::ADD: return x1+x2;
        case OperatorCode::SUB: return x1-x2;
        case OperatorCode::MUL: return x1*x2;
        case OperatorCode::DIV: return x1/x2;
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x1<<" and "<<x2);
    }
}
template<> inline Boolean compute(OperatorCode op, const Boolean& x1, const Boolean& x2) {
    switch(op) {
        case OperatorCode::AND: return x1 && x2;
        case OperatorCode::OR: return x1 || x2;
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x1<<" and "<<x2);
    }
}
template<> inline Tribool compute(OperatorCode op, const Tribool& x1, const Tribool& x2) {
    switch(op) {
        case OperatorCode::AND: return x1 && x2;
        case OperatorCode::OR: return x1 || x2;
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x1<<" and "<<x2);
    }
}
template<> inline String compute(OperatorCode op, const String& x1, const String& x2) {
    switch(op) {
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x1<<" and "<<x2);
    }
}
template<> inline Integer compute(OperatorCode op, const Integer& x1, const Integer& x2) {
    switch(op) {
        case OperatorCode::ADD: return x1+x2;
        case OperatorCode::SUB: return x1-x2;
        case OperatorCode::MUL: return x1*x2;
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x1<<" and "<<x2);
    }
}

template<class X> inline X compute(OperatorCode op, const X& x, Int n) {
    switch(op) {
        case OperatorCode::POW: return pow(x,n);
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x<<" and "<<n);
    }
}
template<> inline Boolean compute(OperatorCode op, const Boolean& x, Int n) {
    switch(op) {
        //case OperatorCode::POW: if(n>=0) return pow(x,Nat(n));
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x<<" and "<<n);
    }
}
template<> inline String compute(OperatorCode op, const String& x, Int n) {
    switch(op) {
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x<<" and "<<n);
    }
}
template<> inline Integer compute(OperatorCode op, const Integer& x, Int n) {
    switch(op) {
        //case OperatorCode::POW: if(n>=0) return pow(x,Nat(n));
        default: ARIADNE_FAIL_MSG("Cannot compute "<<op<<" on "<<x<<" and "<<n);
    }
}

struct GtrZero {}; struct LessZero {};

struct Gtr {
    template<class T1, class T2> auto operator()(const T1& a1, const T2& a2) const -> decltype(a1> a2) { return a1 >  a2; }
    OperatorCode code() const { return OperatorCode::GT; } OperatorKind kind() const { return OperatorKind::COMPARISON; }
};

struct Less {
    template<class T1, class T2> auto operator()(const T1& a1, const T2& a2) const -> decltype(a1< a2) { return a1 <  a2; }
    OperatorCode code() const { return OperatorCode::LT; } OperatorKind kind() const { return OperatorKind::COMPARISON; }
};

struct Geq {
    template<class T1, class T2> auto operator()(const T1& a1, const T2& a2) const -> decltype(a1>=a2) { return a1 >= a2; }
    OperatorCode code() const { return OperatorCode::GEQ; } OperatorKind kind() const { return OperatorKind::COMPARISON; }
};

struct Leq {
    template<class T1, class T2> auto operator()(const T1& a1, const T2& a2) const -> decltype(a1<=a2) { return a1 <= a2; }
    OperatorCode code() const { return OperatorCode::LEQ; } OperatorKind kind() const { return OperatorKind::COMPARISON; }
};

struct Equal {
    template<class T1, class T2> auto operator()(const T1& a1, const T2& a2) const -> decltype(a1==a2) { return a1 == a2; }
    OperatorCode code() const { return OperatorCode::EQ; } OperatorKind kind() const { return OperatorKind::COMPARISON; }
};

struct Unequal {
    template<class T1, class T2> auto operator()(const T1& a1, const T2& a2) const -> decltype(a1!=a2) { return a1 != a2; }
    OperatorCode code() const { return OperatorCode::NEQ; } OperatorKind kind() const { return OperatorKind::COMPARISON; }
};

struct AndOp {
    template<class T> T operator()(const T& a1, const T& a2) const { return a1 && a2; }
    OperatorCode code() const { return OperatorCode::AND; } OperatorKind kind() const { return OperatorKind::BINARY; }
};
struct OrOp {
    template<class T> T operator()(const T& a1, const T& a2) const { return a1 || a2; }
    OperatorCode code() const { return OperatorCode::OR; } OperatorKind kind() const { return OperatorKind::BINARY; }
};
struct NotOp {
    template<class T> T operator()(const T& a) const { return !a; }
    OperatorCode code() const { return OperatorCode::NOT; } OperatorKind kind() const { return OperatorKind::UNARY; }
};

struct Cnst {
    OperatorCode code() const { return OperatorCode::CNST; } OperatorKind kind() const { return OperatorKind::NULLARY; }
};

struct Ind {
    OperatorCode code() const { return OperatorCode::IND; } OperatorKind kind() const { return OperatorKind::COORDINATE; }
};

struct Var {
    OperatorCode code() const { return OperatorCode::VAR; } OperatorKind kind() const { return OperatorKind::COORDINATE; }
};

struct Add {
    OperatorCode code() const { return OperatorCode::ADD; } OperatorKind kind() const { return OperatorKind::BINARY; }
    template<class T> T operator()(const T& a1, const T& a2) const { return a1+a2; }
    template<class X,class D> D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const {
        return d1+d2; }
};
struct Sub {
    OperatorCode code() const { return OperatorCode::SUB; } OperatorKind kind() const { return OperatorKind::BINARY; }
    template<class T> T operator()(const T& a1, const T& a2) const { return a1-a2; }
    template<class X,class D> D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const {
        return d1-d2; }
};
struct Mul {
    OperatorCode code() const { return OperatorCode::MUL; } OperatorKind kind() const { return OperatorKind::BINARY; }
    template<class T> T operator()(const T& a1, const T& a2) const { return a1*a2; }
    template<class X,class D> D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const {
        return a2*d1+a1*d2; }
};
struct Div {
    OperatorCode code() const { return OperatorCode::DIV; } OperatorKind kind() const { return OperatorKind::BINARY; }
    template<class T> T operator()(const T& a1, const T& a2) const { return a1/a2; }
    template<class X,class D> D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const {
        return (d1-a1*(a1/a2))/a2; }
};

struct Pow {
    OperatorCode code() const { return OperatorCode::POW; } OperatorKind kind() const { return OperatorKind::SCALAR; }
    template<class T, class N> T operator()(const T& a, const N& n) const { return pow(a,n); }
};
struct Fma {
    OperatorCode code() const { return OperatorCode::FMA; } OperatorKind kind() const { return OperatorKind::TERNARY; }
    template<class T, class S> T operator()(const T& a1, const S& a2, const S& a3) const { return a1*a2+a3; }
};

struct Pos {
    OperatorCode code() const { return OperatorCode::POS; } OperatorKind kind() const { return OperatorKind::UNARY; }
    template<class T> T operator()(const T& a) const { return a; }
    template<class X,class D> D derivative(const X& a, const D& d) const { return d; }
};
struct Neg {
    OperatorCode code() const { return OperatorCode::NEG; } OperatorKind kind() const { return OperatorKind::UNARY; }
    template<class T> T operator()(const T& a) const { return -(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return -d; }
};
struct Rec {
    OperatorCode code() const { return OperatorCode::REC; } OperatorKind kind() const { return OperatorKind::UNARY; }
    template<class T> T operator()(const T& a) const { return rec(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return d/(-sqr(a)); }
};
struct Sqr {
    OperatorCode code() const { return OperatorCode::SQR; } OperatorKind kind() const { return OperatorKind::UNARY; }
    template<class T> T operator()(const T& a) const { return sqr(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return (2*a)*d; }
};
struct Sqrt {
    OperatorCode code() const { return OperatorCode::SQRT; } OperatorKind kind() const { return OperatorKind::UNARY; }
    template<class T> T operator()(const T& a) const { return sqrt(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return d/(2*sqrt(a)); }
};
struct Exp {
    OperatorCode code() const { return OperatorCode::EXP; } OperatorKind kind() const { return OperatorKind::UNARY; }
    template<class T> T operator()(const T& a) const { return exp(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return exp(a)*d; }
};
struct Log {
    OperatorCode code() const { return OperatorCode::LOG; } OperatorKind kind() const { return OperatorKind::UNARY; }
    template<class T> T operator()(const T& a) const { return log(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return d/a; }
};
struct Sin {
    OperatorCode code() const { return OperatorCode::SIN; } OperatorKind kind() const { return OperatorKind::UNARY; }
    template<class T> T operator()(const T& a) const { return sin(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return cos(a)*d; }
};
struct Cos {
    OperatorCode code() const { return OperatorCode::COS; } OperatorKind kind() const { return OperatorKind::UNARY; }
    template<class T> T operator()(const T& a) const { return cos(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return (-sin(a))*d; }
};
struct Tan {
    OperatorCode code() const { return OperatorCode::TAN; } OperatorKind kind() const { return OperatorKind::UNARY; }
    template<class T> T operator()(const T& a) const { return tan(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return rec(sqr(cos(a)))*d; }
};
struct Atan {
    OperatorCode code() const { return OperatorCode::ATAN; } OperatorKind kind() const { return OperatorKind::UNARY; }
    template<class T> T operator()(const T& a) const { return atan(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return rec(1+sqr(a))*d; }
};

struct Max {
    OperatorCode code() const { return OperatorCode::MAX; } OperatorKind kind() const { return OperatorKind::BINARY; }
    template<class T> T operator()(const T& a1, const T& a2) const { return max(a1,a2); }
    template<class X,class D> D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const { return a1>=a2 ? d1 : d2; }
};

struct Min {
    OperatorCode code() const { return OperatorCode::MIN; } OperatorKind kind() const { return OperatorKind::BINARY; }
    template<class T> T operator()(const T& a1, const T& a2) const { return min(a1,a2); }
    template<class X,class D> D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const { return a1<=a2 ? d1 : d2; }
};

struct Abs {
    OperatorCode code() const { return OperatorCode::MAX; } OperatorKind kind() const { return OperatorKind::BINARY; }
    template<class T> T operator()(const T& a) const { return abs(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return a>=0 ? a : -a; }
};

struct Sgn {
    OperatorCode code() const { return OperatorCode::SGN; }
    Tribool operator()(const Real& a) const { if(definitely(a>0)) { return true; } else if(definitely(a<0)) { return false; } else { return indeterminate; } }
};


inline OutputStream& operator<<(OutputStream& os, const Less& v) { return os << "<"; }
inline OutputStream& operator<<(OutputStream& os, const Gtr& v) { return os << ">"; }
inline OutputStream& operator<<(OutputStream& os, const Leq& v) { return os << "<="; }
inline OutputStream& operator<<(OutputStream& os, const Geq& v) { return os << ">="; }
inline OutputStream& operator<<(OutputStream& os, const Equal& v) { return os << "=="; }
inline OutputStream& operator<<(OutputStream& os, const Unequal& v) { return os << "!="; }

inline OutputStream& operator<<(OutputStream& os, const AndOp& v) { return os << "&&"; }
inline OutputStream& operator<<(OutputStream& os, const OrOp& v) { return os << "||"; }
inline OutputStream& operator<<(OutputStream& os, const NotOp& v) { return os << "!"; }

inline OutputStream& operator<<(OutputStream& os, const Add& v) { return os << "+"; }
inline OutputStream& operator<<(OutputStream& os, const Sub& v) { return os << "-"; }
inline OutputStream& operator<<(OutputStream& os, const Mul& v) { return os << "*"; }
inline OutputStream& operator<<(OutputStream& os, const Div& v) { return os << "/"; }

inline OutputStream& operator<<(OutputStream& os, const Pos& op) { return os << "pos"; }
inline OutputStream& operator<<(OutputStream& os, const Neg& op) { return os << "neg"; }
inline OutputStream& operator<<(OutputStream& os, const Rec& op) { return os << "rec"; }
inline OutputStream& operator<<(OutputStream& os, const Pow& op) { return os << "pow"; }
inline OutputStream& operator<<(OutputStream& os, const Sqr& op) { return os << "sqr"; }
inline OutputStream& operator<<(OutputStream& os, const Sqrt& op) { return os << "sqrt"; }

inline OutputStream& operator<<(OutputStream& os, const Exp& op) { return os << "exp"; }
inline OutputStream& operator<<(OutputStream& os, const Log& op) { return os << "log"; }
inline OutputStream& operator<<(OutputStream& os, const Sin& op) { return os << "sin"; }
inline OutputStream& operator<<(OutputStream& os, const Cos& op) { return os << "cos"; }
inline OutputStream& operator<<(OutputStream& os, const Tan& op) { return os << "tan"; }

inline OutputStream& operator<<(OutputStream& os, const Max& op) { return os << "max"; }
inline OutputStream& operator<<(OutputStream& os, const Min& op) { return os << "min"; }
inline OutputStream& operator<<(OutputStream& os, const Abs& op) { return os << "abs"; }

inline OutputStream& operator<<(OutputStream& os, const Sgn& op) { return os << "sgn"; }


} // namespace Ariadne

#endif // ARIADNE_OPERATORS_H
