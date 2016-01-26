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
 *  You should have received a copy of the GNU G3c767e04cec413f9afb4c30b521ca71ceb5b0409eneral Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file real.cc
 *  \brief
 */

#include "utility/module.h"
#include "numeric/operators.h"
#include "expression/templates.h"

#include "logical.h"
#include "real.h"
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
    virtual Float64Bounds _value() const = 0;
    virtual Float64Bounds _evaluate(Precision64) const = 0;
    virtual FloatMPBounds _evaluate(PrecisionMP) const = 0;
  public:
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

template<class O, class... AS> struct RealWrapper;

template<class O, class A> struct RealWrapper<O,A> : virtual RealInterface, ExpressionTemplate<O,A>, Float64Bounds {
    RealWrapper(O o, A a) : ExpressionTemplate<O,A>(o,a)
        , Float64Bounds(static_cast<ExpressionTemplate<O,A>const&>(*this).operator Float64Bounds()) { }
    virtual Float64Bounds _value() const { return static_cast<Float64Bounds const&>(*this); }
    virtual Float64Bounds _evaluate(Precision64 pr) const {  return static_cast<Float64Bounds>(*this); }
    virtual FloatMPBounds _evaluate(PrecisionMP pr) const {  return this->_op(this->_arg.get(pr)); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<ExpressionTemplate<O,A> const&>(*this); }
};

template<class O, class A1, class A2> struct RealWrapper<O,A1,A2> : virtual RealInterface, ExpressionTemplate<O,A1,A2>, Float64Bounds {
    RealWrapper(O o, A1 a1, A2 a2) : ExpressionTemplate<O,A1,A2>(o,a1,a2)
        , Float64Bounds(static_cast<ExpressionTemplate<O,A1,A2>const&>(*this).operator Float64Bounds()) { }
    virtual Float64Bounds _value() const { return static_cast<Float64Bounds const&>(*this); }
    virtual Float64Bounds _evaluate(Precision64 pr) const {  return static_cast<Float64Bounds>(*this); }
    virtual FloatMPBounds _evaluate(PrecisionMP pr) const {  return this->_op(this->_arg1.get(pr),this->_arg2.get(pr)); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<ExpressionTemplate<O,A1,A2> const&>(*this); }
};

template<class A, class N> struct RealWrapper<Pow,A,N> : virtual RealInterface, ExpressionTemplate<Pow,A,N>, Float64Bounds {
    RealWrapper(Pow o, A a, N n) : ExpressionTemplate<Pow,A,N>(o,a,n)
        , Float64Bounds(static_cast<ExpressionTemplate<Pow,A,N>const&>(*this).operator Float64Bounds()) { }
    virtual Float64Bounds _value() const { return static_cast<Float64Bounds const&>(*this); }
    virtual Float64Bounds _evaluate(Precision64 pr) const {  return static_cast<Float64Bounds>(*this); }
    virtual FloatMPBounds _evaluate(PrecisionMP pr) const {  return this->_op(this->_arg.get(pr),this->_n); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<ExpressionTemplate<Pow,A,N> const&>(*this); }
};

template<class X> struct RealConstant : RealInterface, Float64Bounds {
    X _c;
  public:
    RealConstant(X const& x) : Float64Bounds(x), _c(x) { }
    virtual Float64Bounds _value() const { return static_cast<Float64Bounds const&>(*this); }
    virtual Float64Bounds _evaluate(Precision64 pr) const { return static_cast<Float64Bounds>(*this); }
    virtual FloatMPBounds _evaluate(PrecisionMP pr) const { return FloatMPBounds(this->_c,pr); }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_c; }
};

template<> struct RealConstant<Integer> : RealInterface, Float64Bounds {
    typedef Integer X;
    X _c;
  public:
    RealConstant(X const& x) : Float64Bounds(x), _c(x) { }
    virtual Float64Bounds _value() const { return static_cast<Float64Bounds const&>(*this); }
    virtual Float64Bounds _evaluate(Precision64 pr) const { return static_cast<Float64Bounds>(*this); }
    virtual FloatMPBounds _evaluate(PrecisionMP pr) const { return FloatMPBounds(this->_c,pr); }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_c; }
};

template<> struct RealConstant<Float64Bounds> : RealInterface, Float64Bounds {
    typedef Float64Bounds X;
  public:
    RealConstant(X const& x) : Float64Bounds(x) { }
    virtual Float64Bounds _value() const { return static_cast<Float64Bounds const&>(*this); }
    virtual Float64Bounds _evaluate(Precision64 pr) const { return static_cast<Float64Bounds>(*this); }
    virtual FloatMPBounds _evaluate(PrecisionMP pr) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<Float64Bounds const&>(*this); }
};

template<class O, class... A> inline Real make_real(O o, A... a) {
    return Real(std::make_shared<RealWrapper<O,A...>>(o,a...));
}

inline Real::Real(SharedPointer<RealInterface> p) : _ptr(p) { }


// FIXME: Is this necessary?
Real::Real(double l, double a, double u)
    : Real(std::make_shared<RealConstant<Float64Bounds>>(Float64Bounds(l,u)))
{
}

// FIXME: Is this necessary?
Real::Real(double x)
    : Real(std::make_shared<RealConstant<Float64Bounds>>(Float64Bounds(x)))
{
}

Real::Real(Dyadic const& d) : Real(Rational(d)) { }
Real::Real(Decimal const& d) : Real(Rational(d)) { }

Float64UpperBound Real::upper() const { return this->_ptr->_value(); }
Float64LowerBound Real::lower() const { return this->_ptr->_value(); }
Float64Approximation Real::approx() const { return this->_ptr->_value(); }

double Real::get_d() const { return this->approx().get_d(); }

/*
template<class PR> Float<MetricTag,PR>::Float(Real const& x) : Float<MetricTag,PR>(x.lower(),x.upper()) { }
template<class PR> Float<BoundedTag,PR>::Float(Real const& x) : Float<BoundedTag,PR>(x.lower(),x.upper()) { }
template<class PR> Float<UpperTag,PR>::Float(Real const& x) : Float<UpperTag,PR>(x.upper()) { }
template<class PR> Float<LowerTag,PR>::Float(Real const& x) : Float<LowerTag,PR>(x.lower()) { }
template<class PR> Float<ApproximateTag,PR>::Float(Real const& x) : Float<ApproximateTag,PR>(x.approx()) { }
*/

template<> Float<MetricTag,Precision64>::Float(Real const& x) : Float<MetricTag,Precision64>(x.approx().raw(),max(sub_up(x.upper().raw(),x.approx().raw()),sub_up(x.approx().raw(),x.lower().raw()))) { }
template<> Float<BoundedTag,Precision64>::Float(Real const& x) : Float<BoundedTag,Precision64>(x.lower(),x.upper()) { }
template<> Float<UpperTag,Precision64>::Float(Real const& x) : Float<UpperTag,Precision64>(x.upper()) { }
template<> Float<LowerTag,Precision64>::Float(Real const& x) : Float<LowerTag,Precision64>(x.lower()) { }
template<> Float<ApproximateTag,Precision64>::Float(Real const& x) : Float<ApproximateTag,Precision64>(x.approx()) { }

Real::Real(std::uint64_t m, Void*) : Real(std::make_shared<RealConstant<Integer>>(m)) { }
Real::Real(std::int64_t n, Void*) : Real(std::make_shared<RealConstant<Integer>>(n)) { }

Real::Real() : Real(std::make_shared<RealConstant<Integer>>(0)) { }
Real::Real(Integer const& x) : Real(std::make_shared<RealConstant<Integer>>(x)) { }
Real::Real(Rational const& x) : Real(std::make_shared<RealConstant<Rational>>(x)) { }
Real::Real(Float64Value x) : Real(std::make_shared<RealConstant<Float64Value>>(x)) { }

Real add(Real const& x1, Real const& x2) { return make_real(Add(),x1,x2); }
Real sub(Real const& x1, Real const& x2) { return make_real(Sub(),x1,x2); }
Real mul(Real const& x1, Real const& x2) { return make_real(Mul(),x1,x2); }
Real div(Real const& x1, Real const& x2) { return make_real(Div(),x1,x2); }
Real pow(Real const& x1, Nat m2) { return make_real(Pow(),x1,Int(m2)); }
Real pow(Real const& x1, Int n2) { return make_real(Pow(),x1,n2); }
Real pos(Real const& x) { return make_real(Pos(),x); }
Real neg(Real const& x) { return make_real(Neg(),x); }
Real sqr(Real const& x) { return make_real(Sqr(),x); }
Real rec(Real const& x) { return make_real(Rec(),x); }
Real sqrt(Real const& x) { return make_real(Sqrt(),x); }
Real exp(Real const& x) { return make_real(Exp(),x); }
Real log(Real const& x) { return make_real(Log(),x); }
Real sin(Real const& x) { return make_real(Sin(),x); }
Real cos(Real const& x) { return make_real(Cos(),x); }
Real tan(Real const& x) { return make_real(Tan(),x); }
Real atan(Real const& x) { return make_real(Atan(),x); }

PositiveReal abs(Real const& x) { return PositiveReal(make_real(Abs(),x)); }
Real max(Real const& x1, Real const& x2) { return make_real(Max(),x1,x2); }
Real min(Real const& x1, Real const& x2) { return make_real(Min(),x1,x2); }

Float64Error mag(Real const& x) { return mag(static_cast<Float64Bounds>(x)); }

OutputStream& operator<<(OutputStream& os, Real const& x) { return x._ptr->_write(os); }

Bool same(Real const& x1, Real const& x2) {
    return same(Float64Bounds(x1),Float64Bounds(x2));
}

PositiveReal dist(Real const& x1, Real const& x2) { return abs(sub(x1,x2)); }

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

NegatedSierpinskian eq(Real const& x1, Real const& x2) { return make_logical<EffectiveLowerTag>(Equal(),x1,x2); }
Kleenean lt(Real const& x1, Real const& x2) { return make_logical<EffectiveTag>(Less(),x1,x2); }

Falsifyable operator==(Real const& x1, Real const& x2) { return make_logical<EffectiveLowerTag>(Equal(),x1,x2); }
Verifyable operator!=(Real const& x1, Real const& x2) { return make_logical<EffectiveUpperTag>(Unequal(),x1,x2); }
Quasidecidable operator< (Real const& x1, Real const& x2) { return make_logical<EffectiveTag>(Less(),x1,x2); }
Quasidecidable operator> (Real const& x1, Real const& x2) { return make_logical<EffectiveTag>(Gtr(),x1,x2); }
Quasidecidable operator<=(Real const& x1, Real const& x2) { return make_logical<EffectiveTag>(Leq(),x1,x2); }
Quasidecidable operator>=(Real const& x1, Real const& x2) { return make_logical<EffectiveTag>(Geq(),x1,x2); }

ValidatedNegatedSierpinskian operator==(Real const& x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
ValidatedSierpinskian operator!=(Real const& x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Kleenean operator< (Real const& x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Kleenean operator> (Real const& x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Kleenean operator<=(Real const& x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Kleenean operator>=(Real const& x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }

template<> String class_name<Real>() { return "Real"; }
template<> String class_name<PositiveReal>() { return "PositiveReal"; }

const Real pi = Real(3.1415926535897930, 3.141592653589793238, 3.1415926535897936);
const Real infinity = Real(std::numeric_limits<double>::infinity());

Float64Bounds Real::operator() (Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

FloatMPBounds Real::operator() (PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}

Float64Bounds Real::get(Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

FloatMPBounds Real::get(PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}

FloatMPBounds Real::evaluate(Accuracy accuracy) const {
    Nat effort=1;
    Nat acc=accuracy.bits();
    PrecisionMP precision(effort*64);
    FloatMPError error_bound(FloatMP(Rational(two_exp(-acc).get_d()),FloatMP::upward,precision));
    FloatMPError error=2u*error_bound;
    FloatMPBounds res;
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

Float64LowerBound LowerReal::operator() (Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

FloatMPLowerBound LowerReal::operator() (PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}

Float64LowerBound LowerReal::get(Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

FloatMPLowerBound LowerReal::get(PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}

UpperReal::UpperReal(Real r) : _ptr(r._ptr) {
}

Float64UpperBound UpperReal::operator() (Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

FloatMPUpperBound UpperReal::operator() (PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}

Float64UpperBound UpperReal::get(Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

FloatMPUpperBound UpperReal::get(PrecisionMP pr) const {
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

PositiveFloat64Bounds PositiveReal::get(Precision64 pr) const {
    return PositiveFloat64Bounds(this->_ptr->_evaluate(pr));
}

PositiveFloatMPBounds PositiveReal::get(PrecisionMP pr) const {
    return PositiveFloatMPBounds(this->_ptr->_evaluate(pr));
}

PositiveReal max(PositiveReal pr1, PositiveReal pr2) { return cast_positive(max(make_signed(pr1),make_signed(pr2))); }
PositiveReal min(PositiveReal pr1, PositiveReal pr2) { return cast_positive(min(make_signed(pr1),make_signed(pr2))); }
PositiveReal rec(PositiveReal pr) { return cast_positive(rec(make_signed(pr))); }
PositiveReal add(PositiveReal pr1, PositiveReal pr2) { return cast_positive(add(make_signed(pr1),make_signed(pr2))); }
PositiveReal mul(PositiveReal pr1, PositiveReal pr2) { return cast_positive(mul(make_signed(pr1),make_signed(pr2))); }
PositiveReal div(PositiveReal pr1, PositiveReal pr2) { return cast_positive(div(make_signed(pr1),make_signed(pr2))); }


PositiveFloat64LowerBound PositiveLowerReal::get(Precision64 pr) const {
    return PositiveFloat64LowerBound(this->_ptr->_evaluate(pr));
}

PositiveFloatMPLowerBound PositiveLowerReal::get(PrecisionMP pr) const {
    return PositiveFloatMPLowerBound(this->_ptr->_evaluate(pr));
}

PositiveFloat64UpperBound PositiveUpperReal::get(Precision64 pr) const {
    return PositiveFloat64UpperBound(this->_ptr->_evaluate(pr));
}

PositiveFloatMPUpperBound PositiveUpperReal::get(PrecisionMP pr) const {
    return PositiveFloatMPUpperBound(this->_ptr->_evaluate(pr));
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
