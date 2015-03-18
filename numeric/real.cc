/***************************************************************************
 *            real.cc
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file real.cc
 *  \brief
 */

#include "utility/module.h"
#include "expression/operators.h"
#include "expression/templates.h"

#include "real.h"
#include "logical.h"
#include "integer.h"
#include "rational.h"

#include "dyadic.h"
#include "decimal.h"

#include "float.h"

#include "number_wrapper.h"

namespace Ariadne {

TwoExp Accuracy::error() const {
    return TwoExp(-(Int)this->bits());
}

typedef Real::Interface RealInterface;

struct Real::Interface {
  public:
    virtual ~Interface() = default;
    virtual BoundedFloat64 _value() const = 0;
    virtual BoundedFloat64 _evaluate(Precision64) const = 0;
    virtual BoundedFloatMP _evaluate(PrecisionMP) const = 0;
  public:
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

template<class O, class... AS> struct RealWrapper;

template<class O, class A> struct RealWrapper<O,A> : virtual RealInterface, ExpressionTemplate<O,A>, BoundedFloat64 {
    RealWrapper(O o, A a) : ExpressionTemplate<O,A>(o,a)
        , BoundedFloat64(static_cast<ExpressionTemplate<O,A>const&>(*this).operator BoundedFloat64()) { }
    virtual BoundedFloat64 _value() const { return static_cast<BoundedFloat64 const&>(*this); }
    virtual BoundedFloat64 _evaluate(Precision64 pr) const {  return static_cast<BoundedFloat64>(*this); }
    virtual BoundedFloatMP _evaluate(PrecisionMP pr) const {  return this->_op(this->_arg.get(pr)); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<ExpressionTemplate<O,A> const&>(*this); }
};

template<class O, class A1, class A2> struct RealWrapper<O,A1,A2> : virtual RealInterface, ExpressionTemplate<O,A1,A2>, BoundedFloat64 {
    RealWrapper(O o, A1 a1, A2 a2) : ExpressionTemplate<O,A1,A2>(o,a1,a2)
        , BoundedFloat64(static_cast<ExpressionTemplate<O,A1,A2>const&>(*this).operator BoundedFloat64()) { }
    virtual BoundedFloat64 _value() const { return static_cast<BoundedFloat64 const&>(*this); }
    virtual BoundedFloat64 _evaluate(Precision64 pr) const {  return static_cast<BoundedFloat64>(*this); }
    virtual BoundedFloatMP _evaluate(PrecisionMP pr) const {  return this->_op(this->_arg1.get(pr),this->_arg2.get(pr)); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<ExpressionTemplate<O,A1,A2> const&>(*this); }
};

template<class A, class N> struct RealWrapper<Pow,A,N> : virtual RealInterface, ExpressionTemplate<Pow,A,N>, BoundedFloat64 {
    RealWrapper(Pow o, A a, N n) : ExpressionTemplate<Pow,A,N>(o,a,n)
        , BoundedFloat64(static_cast<ExpressionTemplate<Pow,A,N>const&>(*this).operator BoundedFloat64()) { }
    virtual BoundedFloat64 _value() const { return static_cast<BoundedFloat64 const&>(*this); }
    virtual BoundedFloat64 _evaluate(Precision64 pr) const {  return static_cast<BoundedFloat64>(*this); }
    virtual BoundedFloatMP _evaluate(PrecisionMP pr) const {  return this->_op(this->_arg.get(pr),this->_n); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<ExpressionTemplate<Pow,A,N> const&>(*this); }
};

template<class X> struct RealConstant : RealInterface, BoundedFloat64 {
    X _c;
  public:
    RealConstant(X const& x) : BoundedFloat64(x), _c(x) { }
    virtual BoundedFloat64 _value() const { return static_cast<BoundedFloat64 const&>(*this); }
    virtual BoundedFloat64 _evaluate(Precision64 pr) const { return static_cast<BoundedFloat64>(*this); }
    virtual BoundedFloatMP _evaluate(PrecisionMP pr) const { return BoundedFloatMP(this->_c,pr); }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_c; }
};

template<> struct RealConstant<Integer> : RealInterface, BoundedFloat64 {
    typedef Integer X;
    X _c;
  public:
    RealConstant(X const& x) : BoundedFloat64(x), _c(x) { }
    virtual BoundedFloat64 _value() const { return static_cast<BoundedFloat64 const&>(*this); }
    virtual BoundedFloat64 _evaluate(Precision64 pr) const { return static_cast<BoundedFloat64>(*this); }
    virtual BoundedFloatMP _evaluate(PrecisionMP pr) const { return BoundedFloatMP(this->_c,pr); }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_c; }
};

template<> struct RealConstant<BoundedFloat64> : RealInterface, BoundedFloat64 {
    typedef BoundedFloat64 X;
  public:
    RealConstant(X const& x) : BoundedFloat64(x) { }
    virtual BoundedFloat64 _value() const { return static_cast<BoundedFloat64 const&>(*this); }
    virtual BoundedFloat64 _evaluate(Precision64 pr) const { return static_cast<BoundedFloat64>(*this); }
    virtual BoundedFloatMP _evaluate(PrecisionMP pr) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<BoundedFloat64 const&>(*this); }
};

template<class O, class... A> inline Real make_real(O o, A... a) {
    return Real(std::make_shared<RealWrapper<O,A...>>(o,a...));
}

inline Real::Real(SharedPointer<RealInterface> p) : _ptr(p) { }


// FIXME: Is this necessary?
Real::Real(double l, double a, double u)
    : Real(std::make_shared<RealConstant<BoundedFloat64>>(BoundedFloat64(l,u)))
{
}

// FIXME: Is this necessary?
Real::Real(double x)
    : Real(std::make_shared<RealConstant<BoundedFloat64>>(BoundedFloat64(x)))
{
}

Real::Real(Dyadic const& d) : Real(Rational(d)) { }
Real::Real(Decimal const& d) : Real(Rational(d)) { }

UpperFloat64 Real::upper() const { return this->_ptr->_value(); }
LowerFloat64 Real::lower() const { return this->_ptr->_value(); }
ApproximateFloat64 Real::approx() const { return this->_ptr->_value(); }

double Real::get_d() const { return this->approx().get_d(); }

/*
template<class PR> Float<Metric,PR>::Float(Real const& x) : Float<Metric,PR>(x.lower(),x.upper()) { }
template<class PR> Float<Bounded,PR>::Float(Real const& x) : Float<Bounded,PR>(x.lower(),x.upper()) { }
template<class PR> Float<Upper,PR>::Float(Real const& x) : Float<Upper,PR>(x.upper()) { }
template<class PR> Float<Lower,PR>::Float(Real const& x) : Float<Lower,PR>(x.lower()) { }
template<class PR> Float<Approximate,PR>::Float(Real const& x) : Float<Approximate,PR>(x.approx()) { }
*/

template<> Float<Metric,Precision64>::Float(Real const& x) : Float<Metric,Precision64>(x.approx().raw(),max(sub_up(x.upper().raw(),x.approx().raw()),sub_up(x.approx().raw(),x.lower().raw()))) { }
template<> Float<Bounded,Precision64>::Float(Real const& x) : Float<Bounded,Precision64>(x.lower(),x.upper()) { }
template<> Float<Upper,Precision64>::Float(Real const& x) : Float<Upper,Precision64>(x.upper()) { }
template<> Float<Lower,Precision64>::Float(Real const& x) : Float<Lower,Precision64>(x.lower()) { }
template<> Float<Approximate,Precision64>::Float(Real const& x) : Float<Approximate,Precision64>(x.approx()) { }

Real::Real(std::uint64_t m, Void*) : Real(std::make_shared<RealConstant<Integer>>(m)) { }
Real::Real(std::int64_t n, Void*) : Real(std::make_shared<RealConstant<Integer>>(n)) { }

Real::Real() : Real(std::make_shared<RealConstant<Integer>>(0)) { }
Real::Real(Integer const& x) : Real(std::make_shared<RealConstant<Integer>>(x)) { }
Real::Real(Rational const& x) : Real(std::make_shared<RealConstant<Rational>>(x)) { }
Real::Real(ExactFloat64 x) : Real(std::make_shared<RealConstant<ExactFloat64>>(x)) { }

Real add(Real x1, Real x2) { return make_real(Add(),x1,x2); }
Real sub(Real x1, Real x2) { return make_real(Sub(),x1,x2); }
Real mul(Real x1, Real x2) { return make_real(Mul(),x1,x2); }
Real div(Real x1, Real x2) { return make_real(Div(),x1,x2); }
Real pow(Real x1, Nat m2) { return make_real(Pow(),x1,Int(m2)); }
Real pow(Real x1, Int n2) { return make_real(Pow(),x1,n2); }
Real pos(Real x) { return make_real(Pos(),x); }
Real neg(Real x) { return make_real(Neg(),x); }
Real sqr(Real x) { return make_real(Sqr(),x); }
Real rec(Real x) { return make_real(Rec(),x); }
Real sqrt(Real x) { return make_real(Sqrt(),x); }
Real exp(Real x) { return make_real(Exp(),x); }
Real log(Real x) { return make_real(Log(),x); }
Real sin(Real x) { return make_real(Sin(),x); }
Real cos(Real x) { return make_real(Cos(),x); }
Real tan(Real x) { return make_real(Tan(),x); }
Real atan(Real x) { return make_real(Atan(),x); }

PositiveReal abs(Real x) { return PositiveReal(make_real(Abs(),x)); }
Real max(Real x1, Real x2) { return make_real(Max(),x1,x2); }
Real min(Real x1, Real x2) { return make_real(Min(),x1,x2); }

ErrorFloat64 mag(Real x) { return mag(static_cast<BoundedFloat64>(x)); }

Real operator+(Real x) { return make_real(Pos(),x); }
Real operator-(Real x) { return make_real(Neg(),x); }
Real operator+(Real x1, Real x2) { return make_real(Add(),x1,x2); }
Real operator-(Real x1, Real x2) { return make_real(Sub(),x1,x2); }
Real operator*(Real x1, Real x2) { return make_real(Mul(),x1,x2); }
Real operator/(Real x1, Real x2) { return make_real(Div(),x1,x2); }
Real& operator+=(Real& x1, Real x2) { return x1=make_real(Add(),x1,x2); }
Real& operator-=(Real& x1, Real x2) { return x1=make_real(Sub(),x1,x2); }
Real& operator*=(Real& x1, Real x2) { return x1=make_real(Mul(),x1,x2); }
Real& operator/=(Real& x1, Real x2) { return x1=make_real(Div(),x1,x2); }

OutputStream& operator<<(OutputStream& os, Real const& x) { return x._ptr->_write(os); }

Bool same(Real x1, Real x2) { ARIADNE_NOT_IMPLEMENTED; }

NegSierpinski eq(Real x1, Real x2) { return BoundedFloat64(x1)==BoundedFloat64(x2); }
Kleenean lt(Real x1, Real x2) { return BoundedFloat64(x1)< BoundedFloat64(x2); }

PositiveReal dist(Real x1, Real x2) { return abs(sub(x1,x2)); }

template<class O, class... ARGS> struct LogicalWrapper;

template<class O> struct LogicalWrapper<O,Real,Real> : virtual LogicalInterface, ExpressionTemplate<O,Real,Real> {
    LogicalWrapper(O o, Real a1, Real a2)
        : ExpressionTemplate<O,Real,Real>(o,a1,a2) { }
    virtual LogicalValue _check(Effort e) const;
    virtual OutputStream& _write(OutputStream& os) const {
        return os << static_cast<ExpressionTemplate<O,Real,Real> const&>(*this); }
};

template<class O> LogicalValue LogicalWrapper<O,Real,Real>::_check(Effort e) const {
    if(e==0u) { Precision64 p; return static_cast<LogicalValue>(this->_op(this->_arg1(p),this->_arg2(p))); }
    else { PrecisionMP p(e*64); return static_cast<LogicalValue>(this->_op(this->_arg1(p),this->_arg2(p))); }
}

template<class P, class O, class... ARGS> Logical<P> make_logical(O op, ARGS ...args) {
    return Logical<P>(std::make_shared<LogicalWrapper<O,ARGS...>>(op,args...));
}


Falsifyable operator==(Real x1, Real x2) { return make_logical<EffectiveLower>(Equal(),x1,x2); }
Verifyable operator!=(Real x1, Real x2) { return make_logical<EffectiveUpper>(Unequal(),x1,x2); }
Quasidecidable operator< (Real x1, Real x2) { return make_logical<Effective>(Less(),x1,x2); }
Quasidecidable operator> (Real x1, Real x2) { return make_logical<Effective>(Gtr(),x1,x2); }
Quasidecidable operator<=(Real x1, Real x2) { return make_logical<Effective>(Leq(),x1,x2); }
Quasidecidable operator>=(Real x1, Real x2) { return make_logical<Effective>(Geq(),x1,x2); }

NegSierpinski operator==(Real x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Sierpinski operator!=(Real x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Kleenean operator< (Real x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Kleenean operator> (Real x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Kleenean operator<=(Real x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Kleenean operator>=(Real x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }

template<> String class_name<Real>() { return "Real"; }
template<> String class_name<PositiveReal>() { return "PositiveReal"; }

const Real pi = Real(3.1415926535897930, 3.141592653589793238, 3.1415926535897936);
const Real infinity = Real(std::numeric_limits<double>::infinity());

BoundedFloat64 Real::operator() (Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

BoundedFloatMP Real::operator() (PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}

BoundedFloat64 Real::get(Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

BoundedFloatMP Real::get(PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}

BoundedFloatMP Real::evaluate(Accuracy accuracy) const {
    Nat effort=1;
    Nat acc=accuracy.bits();
    PrecisionMP precision(effort*64);
    ErrorFloatMP error_bound(FloatMP(Rational(two_exp(-acc).get_d()),FloatMP::upward,precision));
    ErrorFloatMP error=2u*error_bound;
    BoundedFloatMP res;
    while (!(error.raw()<error_bound.raw())) {
        res=(*this)(precision);
        error=res.error();
        effort+=1;
        precision=PrecisionMP(effort*64);
    }
    return res;
}




LowerReal::LowerReal(Real r) : _ptr(r._ptr) {
}

LowerFloat64 LowerReal::operator() (Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

LowerFloatMP LowerReal::operator() (PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}

LowerFloat64 LowerReal::get(Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

LowerFloatMP LowerReal::get(PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}

UpperReal::UpperReal(Real r) : _ptr(r._ptr) {
}

UpperFloat64 UpperReal::operator() (Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

UpperFloatMP UpperReal::operator() (PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}

UpperFloat64 UpperReal::get(Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

UpperFloatMP UpperReal::get(PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}

inline Real const& cast_real(LowerReal const& lr) { return reinterpret_cast<Real const&>(lr); }
inline Real const& cast_real(UpperReal const& ur) { return reinterpret_cast<Real const&>(ur); }
inline Real const& make_signed(PositiveReal const& pr) { return pr; }
inline PositiveReal const& cast_positive(Real const& pr) { return static_cast<PositiveReal const&>(pr); }
inline LowerReal const& make_lower(Real const& r) { return reinterpret_cast<LowerReal const&>(r); }
inline UpperReal const& make_upper(Real const& r) { return reinterpret_cast<UpperReal const&>(r); }

LowerReal max(LowerReal lr1, LowerReal lr2) { return make_lower(max(cast_real(lr1),cast_real(lr2))); }
LowerReal min(LowerReal lr1, LowerReal lr2) { return make_lower(min(cast_real(lr1),cast_real(lr2))); }
Real min(LowerReal lr1, Real r2) { return min(cast_real(lr1),r2); }
Real min(Real r1, LowerReal lr2) { return min(r1,cast_real(lr2)); }

UpperReal max(UpperReal ur1, UpperReal ur2) { return make_upper(max(cast_real(ur1),cast_real(ur2))); }
Real max(UpperReal ur1, Real r2) { return max(cast_real(ur1),r2); }
Real max(Real r1, UpperReal ur2) { return max(r1,cast_real(ur2)); }
UpperReal min(UpperReal ur1, UpperReal ur2) { return make_upper(min(cast_real(ur1),cast_real(ur2))); }

LowerReal neg(UpperReal ur) { return make_lower(neg(cast_real(ur))); }
UpperReal neg(LowerReal lr) { return make_upper(neg(cast_real(lr))); }
LowerReal add(LowerReal lr1, LowerReal lr2) { return make_lower(add(cast_real(lr1),cast_real(lr2))); }
UpperReal add(UpperReal ur1, UpperReal ur2) { return make_upper(add(cast_real(ur1),cast_real(ur2))); }
LowerReal add(LowerReal lr1, UpperReal ur2) { return make_lower(add(cast_real(lr1),cast_real(ur2))); }
UpperReal add(UpperReal ur1, LowerReal lr2) { return make_upper(add(cast_real(ur1),cast_real(lr2))); }

PositiveBoundedFloat64 PositiveReal::get(Precision64 pr) const {
    return PositiveBoundedFloat64(this->_ptr->_evaluate(pr));
}

PositiveBoundedFloatMP PositiveReal::get(PrecisionMP pr) const {
    return PositiveBoundedFloatMP(this->_ptr->_evaluate(pr));
}

PositiveReal max(PositiveReal pr1, PositiveReal pr2) { return cast_positive(max(make_signed(pr1),make_signed(pr2))); }
PositiveReal min(PositiveReal pr1, PositiveReal pr2) { return cast_positive(min(make_signed(pr1),make_signed(pr2))); }
PositiveReal rec(PositiveReal pr) { return cast_positive(rec(make_signed(pr))); }
PositiveReal add(PositiveReal pr1, PositiveReal pr2) { return cast_positive(add(make_signed(pr1),make_signed(pr2))); }
PositiveReal mul(PositiveReal pr1, PositiveReal pr2) { return cast_positive(mul(make_signed(pr1),make_signed(pr2))); }
PositiveReal div(PositiveReal pr1, PositiveReal pr2) { return cast_positive(div(make_signed(pr1),make_signed(pr2))); }


PositiveLowerFloat64 PositiveLowerReal::get(Precision64 pr) const {
    return PositiveLowerFloat64(this->_ptr->_evaluate(pr));
}

PositiveLowerFloatMP PositiveLowerReal::get(PrecisionMP pr) const {
    return PositiveLowerFloatMP(this->_ptr->_evaluate(pr));
}

PositiveUpperFloat64 PositiveUpperReal::get(Precision64 pr) const {
    return PositiveUpperFloat64(this->_ptr->_evaluate(pr));
}

PositiveUpperFloatMP PositiveUpperReal::get(PrecisionMP pr) const {
    return PositiveUpperFloatMP(this->_ptr->_evaluate(pr));
}

PositiveUpperReal rec(PositiveLowerReal plr) { return cast_positive(rec(cast_real(plr))); }
PositiveLowerReal rec(PositiveUpperReal pur) { return cast_positive(rec(cast_real(pur))); }
PositiveLowerReal add(PositiveLowerReal plr1, PositiveLowerReal plr2) { return cast_positive(add(cast_real(plr1),cast_real(plr2))); }
PositiveUpperReal add(PositiveUpperReal pur1, PositiveUpperReal pur2) { return cast_positive(add(cast_real(pur1),cast_real(pur2))); }
PositiveLowerReal mul(PositiveLowerReal plr1, PositiveLowerReal plr2) { return cast_positive(mul(cast_real(plr1),cast_real(plr2))); }
PositiveUpperReal mul(PositiveUpperReal pur1, PositiveUpperReal pur2) { return cast_positive(mul(cast_real(pur1),cast_real(pur2))); }
PositiveLowerReal div(PositiveLowerReal plr1, PositiveUpperReal pur2) { return cast_positive(div(cast_real(plr1),cast_real(pur2))); }
PositiveUpperReal div(PositiveUpperReal pur1, PositiveLowerReal plr2) { return cast_positive(div(cast_real(pur1),cast_real(plr2))); }

} // namespace Ariadne

