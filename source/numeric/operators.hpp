/***************************************************************************
 *            operators.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file operators.hpp
 *  \brief Numerical operator classes
 */

#ifndef ARIADNE_OPERATORS_HPP
#define ARIADNE_OPERATORS_HPP

#include <iostream>

#include "logical.decl.hpp"
#include "number.decl.hpp"

namespace Ariadne {

typedef void Void;
typedef bool Bool;
typedef unsigned int Nat;
typedef int Int;

typedef std::ostream OutputStream;


enum class OperatorKind : char {
    VARIABLE,
    COORDINATE,
    NULLARY,
    UNARY,
    BINARY,
    TERNARY,
    SCALAR,
    GRADED,
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
    NUL,   // Zero
    POS,   // Unary plus
    NEG,   // Unary negation
    HLF,   // Halve
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
    SSUB,  // Subtraction from a scalar constant
    SMUL,  // Multiplication by a scalar constant
    SDIV,  // Division into a scalar constant
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

template<class X> struct Logic;
template<> struct Logic<String> { typedef Boolean Type; };
template<> struct Logic<Integer> { typedef Boolean Type; };
template<> struct Logic<Real> { typedef Kleenean Type; };
template<class X> using LogicType = typename Logic<X>::Type;

template<class X> LogicType<X> compare(OperatorCode op, const X& x1, const X& x2);
template<class X> X compute(OperatorCode op, const X& x);
template<class X1, class X2> X2 compute(OperatorCode op, const X1& x1, const X2& x2);
template<class X> X compute(OperatorCode op, const X& x, Int n);

template<class X> X derivative(OperatorCode op, const X& x);
template<class X, class D> D derivative(OperatorCode op, const X& x, const D& dx);
template<class X, class D> D derivative(OperatorCode op, const X& x1, const D& dx1, const X& x2, const D& dx2);
template<class X, class D> D derivative(OperatorCode op, const X& x, const D& dx, Int n);

template<class OBJ> struct Object { OBJ const& upcast() const { return static_cast<OBJ const&>(*this); } };
template<class OP> struct OperatorObject : Object<OP> { };
template<class CMP> struct ComparisonObject : Object<CMP> { };

template<class OP> inline OutputStream& operator<<(OutputStream& os, OperatorObject<OP> const& op) {
    return os << op.upcast().code(); }
template<class CMP> inline OutputStream& operator<<(OutputStream& os, ComparisonObject<CMP> const& cmp) {
    return os << cmp.upcast().code(); }

struct GtrZero {}; struct LessZero {};

struct Gtr : ComparisonObject<Gtr> {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(a1> a2) { return a1 >  a2; }
    static constexpr OperatorCode code() { return OperatorCode::GT; } OperatorKind kind() { return OperatorKind::COMPARISON; }
};

struct Less : ComparisonObject<Less> {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(a1< a2) { return a1 <  a2; }
    static constexpr OperatorCode code() { return OperatorCode::LT; } static constexpr OperatorKind kind() { return OperatorKind::COMPARISON; }
};

struct Geq : ComparisonObject<Geq> {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(a1>=a2) { return a1 >= a2; }
    static constexpr OperatorCode code() { return OperatorCode::GEQ; } static constexpr OperatorKind kind() { return OperatorKind::COMPARISON; }
};

struct Leq : ComparisonObject<Leq> {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(a1<=a2) { return a1 <= a2; }
    static constexpr OperatorCode code() { return OperatorCode::LEQ; } static constexpr OperatorKind kind() { return OperatorKind::COMPARISON; }
};

struct Equal : ComparisonObject<Equal> {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(a1==a2) { return a1 == a2; }
    static constexpr OperatorCode code() { return OperatorCode::EQ; } static constexpr OperatorKind kind() { return OperatorKind::COMPARISON; }
};

struct Unequal : ComparisonObject<Unequal> {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(a1!=a2) { return a1 != a2; }
    static constexpr OperatorCode code() { return OperatorCode::NEQ; } static constexpr OperatorKind kind() { return OperatorKind::COMPARISON; }
};

struct XOrOp : OperatorObject<XOrOp> {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 xor a2) { return a1 xor a2; }
    static constexpr OperatorCode code() { return OperatorCode::XOR; } static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
};
struct AndOp : OperatorObject<AndOp> {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 and a2) { return a1 and a2; }
    static constexpr OperatorCode code() { return OperatorCode::AND; } static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
};
struct OrOp : OperatorObject<OrOp> {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 or a2) { return a1 or a2; }
    static constexpr OperatorCode code() { return OperatorCode::OR; } static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
};
struct NotOp : OperatorObject<NotOp> {
    template<class A> auto operator()(A&& a) const -> decltype(not a) { return not a; }
    static constexpr OperatorCode code() { return OperatorCode::NOT; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
};

struct Cnst : OperatorObject<Cnst> {
    static constexpr OperatorCode code() { return OperatorCode::CNST; } static constexpr OperatorKind kind() { return OperatorKind::NULLARY; }
};

struct Ind : OperatorObject<Ind> {
    static constexpr OperatorCode code() { return OperatorCode::IND; } static constexpr OperatorKind kind() { return OperatorKind::COORDINATE; }
};

struct Var : OperatorObject<Var> {
    static constexpr OperatorCode code() { return OperatorCode::VAR; } static constexpr OperatorKind kind() { return OperatorKind::COORDINATE; }
};

struct Plus {
    template<class A> auto operator()(A&& a) const -> decltype(+a) { return +a; }
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(a1+a2) { return a1+a2; }
};
struct Minus {
    template<class A> auto operator()(A&& a) const -> decltype(-a) { return -a; }
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(a1-a2) { return a1-a2; }
};
struct Times {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(a1*a2) { return a1*a2; }
};
struct Divides {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(a1/a2) { return a1/a2; }
};

struct RAdd : OperatorObject<RAdd> {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(add(a2,a1)) { return add(a2,a1); }
};
struct RSub : OperatorObject<RSub> {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(sub(a2,a1)) { return sub(a2,a1); }
};
struct RMul : OperatorObject<RMul> {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(mul(a2,a1)) { return mul(a2,a1); }
};
struct RDiv : OperatorObject<RDiv> {
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(div(a2,a1)) { return div(a2,a1); }
};

struct FirstArgumentTag { };
struct SecondArgumentTag { };
struct ThirdArgumentTag { };

struct Add : OperatorObject<Add> {
    static constexpr OperatorCode code() { return OperatorCode::ADD; } static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(add(a1,a2)) { return add(a1,a2); }
    template<class X> X derivative(const X& a1, const X& a2, FirstArgumentTag) const { return a1.create(1); }
    template<class X> X derivative(const X& a1, const X& a2, SecondArgumentTag) const { return a2.create(1); }
    template<class X,class D> D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const { return d1+d2; }
};
struct Sub : OperatorObject<Sub> {
    static constexpr OperatorCode code() { return OperatorCode::SUB; } static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(sub(a1,a2)) { return sub(a1,a2); }
    template<class X> X derivative(const X& a1, const X& a2, FirstArgumentTag) const { return a1.create(+1); }
    template<class X> X derivative(const X& a1, const X& a2, SecondArgumentTag) const { return a2.create(-1); }
    template<class X,class D> D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const { return d1-d2; }
};
struct Mul : OperatorObject<Mul> {
    static constexpr OperatorCode code() { return OperatorCode::MUL; } static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(mul(a1,a2)) { return mul(a1,a2); }
    template<class X> X derivative(const X& a1, const X& a2, FirstArgumentTag) const { return a2; }
    template<class X> X derivative(const X& a1, const X& a2, SecondArgumentTag) const { return a1; }
    template<class X,class D> D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const { return a2*d1+a1*d2; }
};
struct Div : OperatorObject<Div> {
    static constexpr OperatorCode code() { return OperatorCode::DIV; } static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(div(a1,a2)) { return div(a1,a2); }
    template<class X> X derivative(const X& a1, const X& a2, FirstArgumentTag) const { return rec(a2); }
    template<class X> X derivative(const X& a1, const X& a2, SecondArgumentTag) const { return div(neg(a1),sqr(a2)); }
    template<class X,class D> D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const { return (d1-d2*(a1/a2))/a2; }
};

struct Pow : OperatorObject<Pow> {
    static constexpr OperatorCode code() { return OperatorCode::POW; } static constexpr OperatorKind kind() { return OperatorKind::GRADED; }
    template<class A, class N> auto operator()(A&& a, N&& n) const -> decltype(pow(a,n)) { return pow(a,n); }
    template<class X, class N> X derivative(const X& a, N n) const { return n*pow(a,n-1); }
    template<class X,class D, class N> D derivative(const X& a, const D& d, N n) const { return derivative(a,n)*d; }
};
struct Fma : OperatorObject<Fma> {
    static constexpr OperatorCode code() { return OperatorCode::FMA; } static constexpr OperatorKind kind() { return OperatorKind::TERNARY; }
    template<class A1, class A2, class A3> auto operator()(A1&& a1, A2&& a2, A3&& a3) const -> decltype(a1*a2+a3) { return a1*a2+a3; }
    template<class X> X derivative(const X& a1, const X& a2, const X& a3, FirstArgumentTag) const { return a2; }
    template<class X> X derivative(const X& a1, const X& a2, const X& a3, SecondArgumentTag) const { return a1; }
    template<class X> X derivative(const X& a1, const X& a2, const X& a3, ThirdArgumentTag) const { return nul(a3)+1; }
    template<class X,class D> D derivative(const X& a1, const D& d1, const X& a2, const D& d2, const X& a3, const D& d3) const {
        return a2*d1+a1*d2+d3; }
};

struct Nul : OperatorObject<Nul> {
    static constexpr OperatorCode code() { return OperatorCode::NUL; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(nul(a)) { return nul(a); }
    template<class X> X derivative(const X& a) const { return nul(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return nul(d); }
};
struct Pos : OperatorObject<Pos> {
    static constexpr OperatorCode code() { return OperatorCode::POS; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(pos(a)) { return pos(a); }
    template<class X> X derivative(const X& a) const { return nul(a)+1; }
    template<class X,class D> D derivative(const X& a, const D& d) const { return d; }
};
struct Neg : OperatorObject<Neg> {
    static constexpr OperatorCode code() { return OperatorCode::NEG; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(neg(a)) { return neg(a); }
    template<class X> X derivative(const X& a) const { return nul(a)-1; }
    template<class X,class D> D derivative(const X& a, const D& d) const { return -d; }
};
struct Hlf : OperatorObject<Neg> {
    static constexpr OperatorCode code() { return OperatorCode::HLF; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(hlf(a)) { return hlf(a); }
    template<class X> X derivative(const X& a) const { return hlf(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return hlf(d); }
};
struct Rec : OperatorObject<Rec> {
    static constexpr OperatorCode code() { return OperatorCode::REC; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(rec(a)) { return rec(a); }
    template<class X> X derivative(const X& a) const { return rec(neg(sqr(a))); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return d/(neg(sqr(a))); }
};
struct Sqr : OperatorObject<Sqr> {
    static constexpr OperatorCode code() { return OperatorCode::SQR; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(sqr(a)) { return sqr(a); }
    template<class X> X derivative(const X& a) const { return 2*a; }
    template<class X,class D> D derivative(const X& a, const D& d) const { return (2*a)*d; }
};
struct Sqrt : OperatorObject<Sqrt> {
    static constexpr OperatorCode code() { return OperatorCode::SQRT; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(sqrt(a)) { return sqrt(a); }
    template<class X> X derivative(const X& a) const { return rec(2*sqrt(a)); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return d/(2*sqrt(a)); }
};
struct Exp : OperatorObject<Exp> {
    static constexpr OperatorCode code() { return OperatorCode::EXP; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(exp(a)) { return exp(a); }
    template<class X> X derivative(const X& a) const { return exp(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return exp(a)*d; }
};
struct Log : OperatorObject<Log> {
    static constexpr OperatorCode code() { return OperatorCode::LOG; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(log(a)) { return log(a); }
    template<class X> X derivative(const X& a) const { return rec(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return d/a; }
};
struct Sin : OperatorObject<Sin> {
    static constexpr OperatorCode code() { return OperatorCode::SIN; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(sin(a)) { return sin(a); }
    template<class X> X derivative(const X& a) const { return cos(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return derivative(a)*d; }
};
struct Cos : OperatorObject<Cos> {
    static constexpr OperatorCode code() { return OperatorCode::COS; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(cos(a)) { return cos(a); }
    template<class X> X derivative(const X& a) const { return neg(sin(a)); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return derivative(a)*d; }
};
struct Tan : OperatorObject<Tan> {
    static constexpr OperatorCode code() { return OperatorCode::TAN; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(tan(a)) { return tan(a); }
    template<class X> X derivative(const X& a) const { return rec(sqr(cos(a))); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return derivative(a)*d; }
};
struct Asin : OperatorObject<Asin> {
    static constexpr OperatorCode code() { return OperatorCode::ASIN; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(asin(a)) { return asin(a); }
    template<class X> X derivative(const X& a) const { return rec(1-sqrt(a)); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return derivative(a)*d; }
};
struct Acos : OperatorObject<Acos> {
    static constexpr OperatorCode code() { return OperatorCode::ACOS; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(acos(a)) { return acos(a); }
    template<class X> X derivative(const X& a) const { return rec(1-sqrt(a)); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return derivative(a)*d; }
};
struct Atan : OperatorObject<Atan> {
    static constexpr OperatorCode code() { return OperatorCode::ATAN; } static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    template<class A> auto operator()(A&& a) const -> decltype(atan(a)) { return atan(a); }
    template<class X> X derivative(const X& a) const { return rec(1+sqr(a)); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return derivative(a)*d; }
};

struct Max : OperatorObject<Max> {
    static constexpr OperatorCode code() { return OperatorCode::MAX; } static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(max(a1,a2)) { return max(a1,a2); }
    template<class X,class D> D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const { return a1>=a2 ? d1 : d2; }
};

struct Min : OperatorObject<Min> {
    static constexpr OperatorCode code() { return OperatorCode::MIN; } static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    template<class A1, class A2> auto operator()(A1&& a1, A2&& a2) const -> decltype(min(a1,a2)) { return min(a1,a2); }
    template<class X,class D> D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const { return a1<=a2 ? d1 : d2; }
};

struct Abs : OperatorObject<Abs> {
    static constexpr OperatorCode code() { return OperatorCode::ABS; } static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    template<class A> auto operator()(A&& a) const -> decltype(abs(a)) { return abs(a); }
    template<class X,class D> D derivative(const X& a, const D& d) const { return a>=0 ? a : -a; }
};

struct Sgn : ComparisonObject<Sgn> {
    static constexpr OperatorCode code() { return OperatorCode::SGN; }
    Kleenean operator()(const Real& a) const;
};


class UnaryOperator {
    OperatorCode _op;
  public:
    template<class OP> UnaryOperator(OP op) : _op(op) { }
    template<class X> X operator() (X const& x) { return compute(_op,x); }
    template<class X> X derivative(X const& x) { return derivative(_op,x); }
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

inline OutputStream& operator<<(OutputStream& os, const Plus& v) { return os << "+"; }
inline OutputStream& operator<<(OutputStream& os, const Minus& v) { return os << "-"; }
inline OutputStream& operator<<(OutputStream& os, const Times& v) { return os << "*"; }
inline OutputStream& operator<<(OutputStream& os, const Divides& v) { return os << "/"; }


} // namespace Ariadne

#endif // ARIADNE_OPERATORS_HPP
