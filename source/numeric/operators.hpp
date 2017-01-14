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
#include <cassert>

#include "logical.decl.hpp"
#include "number.decl.hpp"

namespace Ariadne {

typedef void Void;
typedef bool Bool;
typedef unsigned int Nat;
typedef int Int;

typedef std::ostream OutputStream;

template<class X> struct Logic;
template<> struct Logic<String> { typedef Boolean Type; };
template<> struct Logic<Integer> { typedef Boolean Type; };
template<> struct Logic<Real> { typedef Kleenean Type; };
template<class X> using LogicType = typename Logic<X>::Type;


class Operator {
  public:
    enum class Code : char;
    enum class Kind : char;
    friend Kind get_kind(Code);
  private:
    Code _code;
    Kind _kind;
  public:
    Operator(Code c, Kind k) : _code(c), _kind(k) { }
    Operator(Code c) : _code(c), _kind(get_kind(_code)) { }
    template<class OP> Operator(const OP& op) : _code(op.code()), _kind(op.kind()) { }

    operator Code () const { return _code; }
    Code code() const { return _code; }
    Kind kind() const { return _kind; }
};

using OperatorCode=Operator::Code;
using OperatorKind=Operator::Kind;



enum class Operator::Kind : char {
    VARIABLE,
    COORDINATE,
    NULLARY,
    UNARY,
    BINARY,
    TERNARY,
    GRADED,
    SCALAR,
    COMPARISON
};

enum class Operator::Code : char {
    CNST,  // A constant value
    VAR,   // A named variable
    IND,   // A numbered index
    ADD,   // Addition
    SUB,   // Subtraction
    MUL,   // Multiplication
    DIV,   // Division
    FMA,   // Fused multiply-and-add
    POW,   // Integer power
    ROOT,   // Integer root
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
    SSUB,   // Scalar Subtraction
    SMUL,  // Multiplication by a scalar constant
    SDIV,   // Scalar Division
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
    EQ=-6,    // Equal
    NEQ=-5,   // Not equal
    GEQ=-4,   // Greater or equal
    LEQ=-3,   // Less or equal
    GT=-2,    // Greater than
    LT=-1,    // Less than
    SGN=-7,   // Compare with zero
    SUBS=-8,   // Compare as a subset
    DISJ=-9    // Compare disjointness
};

OutputStream& operator<<(OutputStream& os, const OperatorKind& knd);
OutputStream& operator<<(OutputStream&, const OperatorCode& op);
const char* name(const OperatorCode& op);
const char* symbol(const OperatorCode& op);


template<class X> LogicType<X> compare(OperatorCode op, const X& x1, const X& x2);
template<class X> X compute(OperatorCode op, const X& x);
template<class X1, class X2> X2 compute(OperatorCode op, const X1& x1, const X2& x2);
template<class X> X compute(OperatorCode op, const X& x, Int n);
template<class X1, class X2> X2 compute(OperatorCode op, const X1& x1, const X2& x2);

template<class X> X compute_derivative(OperatorCode op, const X& x);
template<class X, class D> D compute_derivative(OperatorCode op, const X& x, const D& dx);
template<class X, class D> D compute_derivative(OperatorCode op, const X& x1, const D& dx1, const X& x2, const D& dx2);
template<class X, class D> D compute_derivative(OperatorCode op, const X& x, const D& dx, Int n);

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
struct Root : OperatorObject<Root> {
    OperatorCode code() const { return OperatorCode::ROOT; } OperatorKind kind() const { return OperatorKind::GRADED; }
    template<class A, class N> friend auto root(A&& a, N&& n) -> decltype(exp(log(a)/n)) { return exp(log(a)/n); }
    template<class A, class N> auto operator()(A&& a, N&& n) const -> decltype(root(a,n)) { return root(a,n); }
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


/*
template<class OP> struct CodedOperator {
    typedef OP CodeType;
    operator Operator() const { return Operator(static_cast<Operator::Code>(_code)); }
    CodeType code() const { return _code; }
  protected:
    CodeType _code;
};
struct SpecialOperator : CodedOperator<SpecialOperator::Code> {
    enum class Code : char;
};
*/

struct SpecialOperator {
    enum class Code : char;
    operator Operator() const { return Operator(static_cast<Operator::Code>(_code)); }
    Code code() const { return _code; }
  private:
    Code _code;
};

struct UnaryLogicalOperator {
    enum class Code : char;
    template<class OP, EnableIf<IsOneOf<OP,NotOp>> =dummy> UnaryLogicalOperator(OP op) : _code(op) { }
    template<class L> decltype(not declval<L>()) operator()(L&& l) const;
    Code _code;
};
struct UnaryArithmeticOperator {
    enum class Code : char;
    template<class OP, EnableIf<IsOneOf<OP,Nul,Pos,Neg,Sqr,Rec>> =dummy> UnaryArithmeticOperator(OP op) : _code(op) { }
    template<class X> X operator()(X&& x) const;
    Code _code;
};
struct UnaryElementaryOperator {
    enum class Code : char;
    explicit UnaryElementaryOperator(OperatorCode cd);
    UnaryElementaryOperator(Code cd) : _code(cd) { }
    template<class OP, EnableIf<IsOneOf<OP,Pos,Neg,Sqr,Hlf,Rec,Sqrt,Exp,Log,Sin,Cos,Tan,Atan,Abs>> =dummy> UnaryElementaryOperator(OP op) : _code(op) { }
    template<class X> X operator()(X&& x) const;
    Code _code;
};
struct UnaryLatticeOperator {
    enum class Code : char;
    explicit UnaryLatticeOperator(OperatorCode cd);
    template<class OP, EnableIf<IsOneOf<OP,Abs>> =dummy> UnaryLatticeOperator(OP op) : _code(op) { }
    template<class X> X operator()(X&& x) const;
    Code _code;
};

struct BinaryLogicalOperator {
    enum class Code : char;
    explicit BinaryLogicalOperator(OperatorCode cd);
    template<class OP, EnableIf<IsOneOf<OP,AndOp,OrOp>> =dummy> BinaryLogicalOperator(OP op) : _code(op) { }
    template<class L> L operator()(L&& l1, L&& l2) const;
    Code _code;
};
struct BinaryComparisonOperator {
    enum struct Code : char;
    explicit BinaryComparisonOperator(OperatorCode cd);
    template<class OP, EnableIf<IsOneOf<OP,Equal,Unequal,Less,Gtr,Leq,Geq>> =dummy> BinaryComparisonOperator(OP op) : _code(op) { }
    template<class X> typename Logic<X>::Type operator()(const X& x1, const X& x2) const;
    Code _code;
};
struct BinaryLatticeOperator {
    enum struct Code : char;
    explicit BinaryLatticeOperator(OperatorCode cd);
    template<class OP, EnableIf<IsOneOf<OP,Max,Min>> =dummy> BinaryLatticeOperator(OP op) : _code(op) { }
    template<class X> X operator()(X&& x1, X&& x2) const;
    Code _code;
};
struct BinaryArithmeticOperator {
    enum struct Code : char;
    explicit BinaryArithmeticOperator(OperatorCode cd);
    template<class OP, EnableIf<IsOneOf<OP,Add,Sub,Mul,Div>> =dummy> BinaryArithmeticOperator(OP op) : _code(op) { }
    template<class X1, class X2> ArithmeticType<X1,X2> operator()(X1&& x1, X2&& x2) const;
    Code _code;
};
struct BinaryElementaryOperator {
    enum struct Code : char;
    BinaryElementaryOperator(Code cd) : _code(cd) { }
    explicit BinaryElementaryOperator(OperatorCode cd);
    template<class OP, EnableIf<IsOneOf<OP,Add,Sub,Mul,Div,Max,Min>> =dummy> BinaryElementaryOperator(OP op) : _code(op) { }
    template<class X1, class X2> ArithmeticType<X1,X2> operator()(X1&& x1, X2&& x2) const;
    Code _code;
};

struct GradedElementaryOperator {
    enum struct Code : char;
    explicit GradedElementaryOperator(OperatorCode cd);
    template<class OP, EnableIf<IsOneOf<OP,Pow>> =dummy> GradedElementaryOperator(OP op) : _code(op) { }
    template<class X> X operator()(X&& x, Integer const& n) const;
    Code _code;
};

struct TernaryArithmeticOperator {
    enum struct Code : char;
    template<class OP, EnableIf<IsOneOf<OP,Fma>> =dummy> TernaryArithmeticOperator(OP op) : _code(op) { }
    template<class X> X operator()(X&& x1, X&& x2, X&& x3) const;
    Code _code;
};

enum class SpecialOperator::Code : char { CNST=(char)OperatorCode::CNST, VAR, IND };
enum class UnaryLogicalOperator::Code : char { NOT=(char)OperatorCode::NOT };
enum class UnaryLatticeOperator::Code : char { ABS=(char)OperatorCode::ABS };
enum class UnaryArithmeticOperator::Code : char { POS=(char)OperatorCode::POS, NEG, REC, SQR };
enum class BinaryLogicalOperator::Code : char { AND=(char)OperatorCode::AND, OR, XOR, IMPL };
enum class BinaryComparisonOperator::Code : char { EQ=(char)OperatorCode::EQ, NEQ, GEQ, LEQ, GT, LT };
enum class BinaryLatticeOperator::Code : char { MAX=(char)OperatorCode::MAX, MIN };
enum class BinaryArithmeticOperator::Code : char { ADD=(char)OperatorCode::ADD, SUB, MUL, DIV };

enum class UnaryElementaryOperator::Code : char { POS=(char)OperatorCode::POS, NEG, HLF, REC, SQR, SQRT, EXP, LOG, SIN, COS, TAN, ASIN, ACOS, ATAN, ABS=(char)OperatorCode::ABS };
enum class BinaryElementaryOperator::Code : char { ADD=(char)OperatorCode::ADD, SUB, MUL, DIV, MAX=(char)OperatorCode::MAX, MIN };
enum class GradedElementaryOperator::Code : char { POW=(char)OperatorCode::POW, ROOT };

inline BinaryComparisonOperator::BinaryComparisonOperator(OperatorCode cd) : _code(static_cast<Code>(cd)) { assert(Code::EQ <= _code && _code <= Code::LT); }

static_assert((char)UnaryElementaryOperator::Code::ATAN==(char)Operator::Code::ATAN,"");
static_assert((char)BinaryLogicalOperator::Code::IMPL==(char)Operator::Code::IMPL,"");
static_assert((char)BinaryLatticeOperator::Code::MIN==(char)Operator::Code::MIN,"");
static_assert((char)BinaryComparisonOperator::Code::LT==(char)Operator::Code::LT,"");
static_assert((char)BinaryArithmeticOperator::Code::DIV==(char)Operator::Code::DIV,"");
static_assert((char)BinaryElementaryOperator::Code::DIV==(char)Operator::Code::DIV,"");
static_assert((char)BinaryElementaryOperator::Code::MIN==(char)Operator::Code::MIN,"");


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
