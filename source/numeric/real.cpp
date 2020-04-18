/***************************************************************************
 *            numeric/real.cpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file numeric/real.cpp
 *  \brief
 */

#include "../utility/module.hpp"
#include "../numeric/operators.hpp"
#include "../symbolic/templates.hpp"

#include "logical.hpp"
#include "real.hpp"
#include "integer.hpp"
#include "rational.hpp"

#include "twoexp.hpp"
#include "dyadic.hpp"
#include "decimal.hpp"

#include "float_bounds.hpp"

#include "real_interface.hpp"
#include "sequence.hpp"
#include "number_wrapper.hpp"

namespace Ariadne {

TwoExp Accuracy::error() const {
    return TwoExp(-(Int)this->bits());
}

template<class X> struct ValidatedRealWrapper;

template<> struct ValidatedRealWrapper<DyadicBounds> : public ValidatedRealInterface, public DyadicBounds {
    typedef DyadicBounds X;
    ValidatedRealWrapper<X>(X const& x) : X(x) { }
    virtual DyadicBounds _get() const override final { return *this; }
    virtual FloatDPBounds _get(DoublePrecision pr) const override final { return FloatDPBounds(static_cast<DyadicBounds const&>(*this),pr); }
    virtual FloatMPBounds _get(MultiplePrecision pr) const override final { return FloatMPBounds(static_cast<DyadicBounds const&>(*this),pr); }
    virtual OutputStream& _write(OutputStream& os) const override final { return os << static_cast<DyadicBounds const&>(*this); }
};

template<class O, class... AS> struct RealWrapper;

template<class O, class A> struct RealWrapper<O,A> : virtual RealInterface, Symbolic<O,A>, FloatDPBounds {
    RealWrapper(O o, A a) : Symbolic<O,A>(o,a)
        , FloatDPBounds(this->_op(this->_arg.get(dp))) { }
    virtual ValidatedReal _compute(Effort eff) const { return ValidatedReal(this->_compute_get(MP(eff.work()+2))); }
    virtual FloatDPBounds _compute_get(DoublePrecision pr) const {  return static_cast<FloatDPBounds>(*this); }
    virtual FloatMPBounds _compute_get(MultiplePrecision pr) const {  return this->_op(this->_arg.get(pr)); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<Symbolic<O,A> const&>(*this); }
};

template<class O, class A1, class A2> struct RealWrapper<O,A1,A2> : virtual RealInterface, Symbolic<O,A1,A2>, FloatDPBounds {
    RealWrapper(O o, A1 a1, A2 a2) : Symbolic<O,A1,A2>(o,a1,a2)
        , FloatDPBounds(this->_op(this->_arg1.get(dp),this->_arg2.get(dp))) { }
    virtual ValidatedReal _compute(Effort eff) const { return ValidatedReal(this->_compute_get(MP(eff.work()+2))); }
    virtual FloatDPBounds _compute_get(DoublePrecision pr) const {  return static_cast<FloatDPBounds>(*this); }
    virtual FloatMPBounds _compute_get(MultiplePrecision pr) const {  return this->_op(this->_arg1.get(pr),this->_arg2.get(pr)); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<Symbolic<O,A1,A2> const&>(*this); }
};

template<class A, class N> struct RealWrapper<Pow,A,N> : virtual RealInterface, Symbolic<Pow,A,N>, FloatDPBounds {
    RealWrapper(Pow o, A a, N n) : Symbolic<Pow,A,N>(o,a,n)
        , FloatDPBounds(this->_op(this->_arg.get(dp),n)) { }
    virtual ValidatedReal _compute(Effort eff) const { return ValidatedReal(this->_compute_get(MP(eff.work()+2))); }
    virtual FloatDPBounds _compute_get(DoublePrecision pr) const {  return static_cast<FloatDPBounds>(*this); }
    virtual FloatMPBounds _compute_get(MultiplePrecision pr) const {  return this->_op(this->_arg.get(pr),this->_num); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<Symbolic<Pow,A,N> const&>(*this); }
};

template<class X> struct RealWrapper<Cnst,X> : RealInterface, FloatDPBounds {
    X _c;
  public:
    RealWrapper(X const& x) : FloatDPBounds(x,dp), _c(x) { }
    virtual ValidatedReal _compute(Effort eff) const { return ValidatedReal(this->_compute_get(MP(eff.work()+2))); }
    virtual FloatDPBounds _compute_get(DoublePrecision pr) const { return static_cast<FloatDPBounds const&>(*this); }
    virtual FloatMPBounds _compute_get(MultiplePrecision pr) const { return FloatMPBounds(this->_c,pr); }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_c; }
};

template<> struct RealWrapper<Cnst,FloatDPBounds> : RealInterface, FloatDPBounds {
    typedef FloatDPBounds X;
  public:
    RealWrapper(X const& x) : FloatDPBounds(x,dp) { }
    virtual ValidatedReal _compute(Effort eff) const { return ValidatedReal(static_cast<FloatDPBounds const&>(*this)); }
    virtual FloatDPBounds _compute_get(DoublePrecision pr) const { return static_cast<FloatDPBounds const&>(*this); }
    virtual FloatMPBounds _compute_get(MultiplePrecision pr) const { return FloatMPBounds(FloatMP(this->lower_raw(),down,pr),FloatMP(this->upper_raw(),up,pr)); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<FloatDPBounds const&>(*this); }
};

template<> struct RealWrapper<Cnst,EffectiveNumber> : RealInterface, FloatDPBounds {
    typedef EffectiveNumber X;
    X _c;
  public:
    RealWrapper(X const& x) : FloatDPBounds(x,dp) { }
    virtual ValidatedReal _compute(Effort eff) const { return ValidatedReal(this->_compute_get(MP(eff.work()+2))); }
    virtual FloatDPBounds _compute_get(DoublePrecision pr) const { return static_cast<FloatDPBounds const&>(*this); }
    virtual FloatMPBounds _compute_get(MultiplePrecision pr) const { return FloatMPBounds(this->_c,pr); }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_c; }
};

template<class Y> struct RealLimit;

template<> struct RealLimit<Real> : RealInterface {
    Sequence<Real> _seq;
  public:
    RealLimit(Sequence<Real> const& seq) : _seq(seq) { }
    virtual FloatDPBounds _compute_get(DoublePrecision pr) const {
        return FloatDPBounds(_seq[53u],pr).pm(FloatDPError(two^-53,DoublePrecision())); }
    virtual FloatMPBounds _compute_get(MultiplePrecision pr) const {
        Effort eff(pr.bits()+1u); return this->_compute(eff).get(pr); }
    virtual ValidatedReal _compute(Effort eff) const {
        Nat n = eff.work()+1u; Accuracy acc(n+1u); return ValidatedReal(_seq[n].compute(Accuracy(n)).get().pm(two^(-n))); }
    virtual OutputStream& _write(OutputStream& os) const {
        return os << "{" << _seq[0u] << ", " << _seq[1u] << ", " <<_seq[2u] << ", ... }"; }
};

template<> struct RealLimit<Dyadic> : RealInterface {
    Sequence<Dyadic> _seq;
  public:
    RealLimit(Sequence<Dyadic> const& seq) : _seq(seq) { }
    virtual FloatDPBounds _compute_get(DoublePrecision pr) const {
        Dyadic w=_seq[53u]; return FloatDPBounds(w,pr).pm(FloatDPError(two^-53,DoublePrecision())); }
    virtual FloatMPBounds _compute_get(MultiplePrecision pr) const {
        Nat n = pr.bits()+1u; Accuracy acc(n+1u); return FloatMPBounds(_seq[n],pr).pm(FloatMPError(two^-n,pr)); }
    virtual ValidatedReal _compute(Effort eff) const {
        Nat n = eff.work()+1u; Accuracy acc(n+1u); return DyadicBounds(_seq[n]).pm(two^-n); }
    virtual OutputStream& _write(OutputStream& os) const {
        return os << "{" << _seq[0u] << ", " << _seq[1u] << ", " <<_seq[2u] << ", ... }"; }
};

template<> struct RealLimit<DyadicBounds> : RealInterface {
    Sequence<DyadicBounds> _seq;
  public:
    RealLimit(ConvergentSequence<DyadicBounds> const& seq) : _seq(seq) { }
    virtual ValidatedReal _compute(Effort eff) const {
        Nat n = eff.work()+1u; return _seq[n]; }
    virtual FloatDPBounds _compute_get(DoublePrecision pr) const {
        return FloatDPBounds(_seq[0u],pr); }
    virtual FloatMPBounds _compute_get(MultiplePrecision pr) const {
        return FloatMPBounds(_seq[pr.bits()],pr); }
    virtual OutputStream& _write(OutputStream& os) const {
        return os << "{" << _seq[0u] << ", " << _seq[1u] << ", " <<_seq[2u] << ", ... }"; }
};

template<class O, class... A> inline Real make_real(O o, A... a) {
    return Real(std::make_shared<RealWrapper<O,A...>>(o,a...));
}

Real::Real(SharedPointer<RealInterface> p) : _ptr(p) { }

Real::Real(ConvergentSequence<DyadicBounds> const& seq) : Real(std::make_shared<RealLimit<DyadicBounds>>(seq)) { }
Real::Real(FastCauchySequence<Dyadic> const& seq) : Real(std::make_shared<RealLimit<Dyadic>>(seq)) { }

// FIXME: Is this necessary?
Real::Real(double l, double a, double u)
    : Real(std::make_shared<RealWrapper<Cnst,FloatDPBounds>>(FloatDPBounds(l,u)))
{
}

// FIXME: Is this necessary?
Real::Real(double x)
    : Real(std::make_shared<RealWrapper<Cnst,FloatDPBounds>>(FloatDPBounds(x)))
{
}

ValidatedNegatedSierpinskian operator==(Real const& x1, Int64 n2);
ValidatedSierpinskian operator!=(Real const& x1, Int64 n2);
Kleenean operator< (Real const& x1, Int64 n2);
Kleenean operator> (Real const& x1, Int64 n2);
Kleenean operator<=(Real const& x1, Int64 n2);
Kleenean operator>=(Real const& x1, Int64 n2);


UpperReal Real::upper() const { return UpperReal(this->_ptr); }
LowerReal Real::lower() const { return LowerReal(this->_ptr); }

double Real::get_d() const { return this->get(dp).get_d(); }

/*
template<class PR> FloatBall<PR>::Float(Real const& x) : FloatBall<PR>(x.lower(),x.upper()) { }
template<class PR> FloatBounds<PR>::Float(Real const& x) : FloatBounds<PR>(x.lower(),x.upper()) { }
template<class PR> FloatUpperBound<PR>::Float(Real const& x) : FloatUpperBound<PR>(x.upper()) { }
template<class PR> FloatLowerBound<PR>::Float(Real const& x) : FloatLowerBound<PR>(x.lower()) { }
template<class PR> FloatApproximation<PR>::Float(Real const& x) : FloatApproximation<PR>(x.approx()) { }
*/

Real::Real(std::uint64_t m, Void*) : Real(Integer(m)) { }
Real::Real(std::int64_t n, Void*) : Real(Integer(n)) { }

Real::Real() : Real(Integer(0)) { }
Real::Real(ExactDouble d) : Real(std::make_shared<RealWrapper<Cnst,ExactDouble>>(d)) { }
//Real::Real(ExactDouble d) : Real(Dyadic(d)) { }
Real::Real(Integer const& z) : Real(std::make_shared<RealWrapper<Cnst,Integer>>(z)) { }
Real::Real(Dyadic const& w) : Real(std::make_shared<RealWrapper<Cnst,Dyadic>>(w)) { }
Real::Real(Decimal const& d) : Real(std::make_shared<RealWrapper<Cnst,Decimal>>(d)) { }
Real::Real(Rational const& q) : Real(std::make_shared<RealWrapper<Cnst,Rational>>(q)) { }
Real::Real(EffectiveNumber q) : Real(std::make_shared<RealWrapper<Cnst,EffectiveNumber>>(q)) { }
Real::Real(FloatDPValue x) : Real(Dyadic(x.get_d())) { ARIADNE_DEPRECATED("Real::Real(FloatDPValue)","Use Real([Exact]Double) or Real(Dyadic) instead."); }


Real add(Real const& x1, Real const& x2) { return make_real(Add(),x1,x2); }
Real sub(Real const& x1, Real const& x2) { return make_real(Sub(),x1,x2); }
Real mul(Real const& x1, Real const& x2) { return make_real(Mul(),x1,x2); }
Real div(Real const& x1, Real const& x2) { return make_real(Div(),x1,x2); }
Real pow(Real const& x1, Nat m2) { return make_real(Pow(),x1,Int(m2)); }
Real pow(Real const& x1, Int n2) { return make_real(Pow(),x1,n2); }
Real nul(Real const& x) { return Real(0); }
Real pos(Real const& x) { return make_real(Pos(),x); }
Real neg(Real const& x) { return make_real(Neg(),x); }
Real hlf(Real const& x) { return make_real(Hlf(),x); }
Real sqr(Real const& x) { return make_real(Sqr(),x); }
Real rec(Real const& x) { return make_real(Rec(),x); }
Real sqrt(Real const& x) { return make_real(Sqrt(),x); }
Real exp(Real const& x) { return make_real(Exp(),x); }
Real log(Real const& x) { return make_real(Log(),x); }
Real sin(Real const& x) { return make_real(Sin(),x); }
Real cos(Real const& x) { return make_real(Cos(),x); }
Real tan(Real const& x) { return make_real(Tan(),x); }
Real asin(Real const& x) { return make_real(Asin(),x); }
Real acos(Real const& x) { return make_real(Acos(),x); }
Real atan(Real const& x) { return make_real(Atan(),x); }

PositiveReal abs(Real const& x) { return PositiveReal(make_real(Abs(),x)); }
Real max(Real const& x1, Real const& x2) { return make_real(Max(),x1,x2); }
Real min(Real const& x1, Real const& x2) { return make_real(Min(),x1,x2); }

Real limit(FastCauchySequence<Real> const& seq) { return Real(std::make_shared<RealLimit<Real>>(seq)); }
Real limit(FastCauchySequence<Dyadic> const& seq) { return Real(std::make_shared<RealLimit<Dyadic>>(seq)); }

Real choose(Case<LowerKleenean,Real> const& c1, Case<LowerKleenean,Real> const& c2) {
    if (choose(c1.condition(),c2.condition())) { return c1.term(); } else { return c2.term(); } }

class WhenRealExpression : public RealInterface {
    UpperKleenean _p1, _p2; Real _r1, _r2;
  public:
    WhenRealExpression(UpperKleenean const& p1, Real const& r1, UpperKleenean const& p2, Real const& r2) : _p1(p1), _p2(p2), _r1(r1), _r2(r2) { }
    WhenRealExpression(Case<UpperKleenean,Real> const& c1, Case<UpperKleenean,Real> const& c2)
        : WhenRealExpression(c1.condition(),c1.term(),c2.condition(),c2.term()) { }
    virtual ValidatedReal _compute(Effort eff) const;
    virtual FloatDPBounds _compute_get(DoublePrecision pr) const { return this->_compute(Effort(0u)).get(pr); }
    virtual FloatMPBounds _compute_get(MultiplePrecision pr) const { return this->_compute(Effort(pr.bits())).get(pr); }
  public:
    virtual OutputStream& _write(OutputStream& os) const;
};
Real when(Case<UpperKleenean,Real> const& c1, Case<UpperKleenean,Real> const& c2) {
    return Real(std::make_shared<WhenRealExpression>(c1,c2)); }

ValidatedReal WhenRealExpression::_compute(Effort eff) const {
    if(not possibly(_p1.check(eff))) { return _r2.compute(eff); }
    if(not possibly(_p2.check(eff))) { return _r1.compute(eff); }

    ValidatedReal vr1=_r1.compute(eff);
    ValidatedReal vr2=_r2.compute(eff);

    DyadicBounds w1=vr1.get();
    DyadicBounds w2=vr2.get();

    ARIADNE_ASSERT(w1.lower_raw()<=w2.upper_raw() && w1.upper_raw()>=w2.lower_raw());
    return DyadicBounds(min(w1.lower_raw(),w2.lower_raw()),max(w1.upper_raw(),w2.upper_raw()));
}

OutputStream& WhenRealExpression::_write(OutputStream& os) const {
    return os << "when(" << _p1 << " => " << _r1 << " & " << _p2 << " => " << _r2 << ")";
}


// May be true if r>a; may be false if r<b
// Equivalent to choose(r>a,r<b)
Boolean nondeterministic_greater(Real const& r, Rational const& a, Rational const& b) {
    ARIADNE_PRECONDITION(a<b);
    FloatDPBounds x0=r.get(dp);
    if(x0.lower_raw()>a) { return true; } else if(x0.upper_raw()<b) { return false; }
    Nat bits=64;
    while(true) {
        MultiplePrecision pr(bits);
        FloatMPBounds x=r.get(pr);
        if(x.lower_raw()>a) { return true; } else if(x.upper_raw()<b) { return false; }
        bits+=64;
    }
}


PositiveUpperReal mag(Real const& r) { return abs(r); }
FloatDPError mag(Real const& r, DoublePrecision pr) { return mag(r.get(pr)); }

OutputStream& operator<<(OutputStream& os, Real const& x) { return x._ptr->_write(os); }

Bool same(Real const& r1, Real const& r2) {
    // FIXME: Use symbolic approach
    DoublePrecision pr;
    FloatDPBounds x1(r1,pr);
    FloatDPBounds x2(r2,pr);
    return x1.lower_raw()==x2.upper_raw() && x1.upper_raw() == x2.lower_raw();
}

PositiveReal dist(Real const& r1, Real const& r2) { return abs(sub(r1,r2)); }

template<class O, class... ARGS> struct LogicalWrapper;

template<class O> struct LogicalWrapper<O,Real> : virtual LogicalInterface, Symbolic<O,Real> {
    LogicalWrapper(O o, Real a)
        : Symbolic<O,Real>(o,a) { }
    virtual LogicalValue _check(Effort e) const;
    virtual OutputStream& _write(OutputStream& os) const {
        return os << static_cast<Symbolic<O,Real> const&>(*this); }
};

template<class O> LogicalValue LogicalWrapper<O,Real>::_check(Effort e) const {
    if(e==0u) { DoublePrecision p; return static_cast<LogicalValue>(this->_op(this->_arg.get(p))); }
    else { MultiplePrecision p(e*64); return static_cast<LogicalValue>(this->_op(this->_arg.get(p))); }
}

template<class O> struct LogicalWrapper<O,Real,Real> : virtual LogicalInterface, Symbolic<O,Real,Real> {
    LogicalWrapper(O o, Real a1, Real a2)
        : Symbolic<O,Real,Real>(o,a1,a2) { }
    virtual LogicalValue _check(Effort e) const;
    virtual OutputStream& _write(OutputStream& os) const {
        return os << static_cast<Symbolic<O,Real,Real> const&>(*this); }
};

template<class O> LogicalValue LogicalWrapper<O,Real,Real>::_check(Effort e) const {
    if(e==0u) { DoublePrecision p; return static_cast<LogicalValue>(this->_op(this->_arg1.get(p),this->_arg2.get(p))); }
    else { MultiplePrecision p(e*64); return static_cast<LogicalValue>(this->_op(this->_arg1.get(p),this->_arg2.get(p))); }
}

template<class R, class O, class... ARGS> R make_logical(O op, ARGS ...args) {
    return R(std::make_shared<LogicalWrapper<O,ARGS...>>(op,args...));
}
template<class P, class O, class... ARGS> LogicalType<P> make_logical_type(O op, ARGS ...args) {
    return LogicalType<P>(std::make_shared<LogicalWrapper<O,ARGS...>>(op,args...));
}

NegatedSierpinskian eq(Real const& x1, Real const& x2) { return make_logical<NegatedSierpinskian>(Equal(),x1,x2); }
Kleenean lt(Real const& x1, Real const& x2) { return make_logical<Kleenean>(Less(),x1,x2); }

Falsifyable operator==(Real const& x1, Real const& x2) { return make_logical<NegatedSierpinskian>(Equal(),x1,x2); }
Verifyable operator!=(Real const& x1, Real const& x2) { return make_logical<Sierpinskian>(Unequal(),x1,x2); }
Quasidecidable operator< (Real const& x1, Real const& x2) { return make_logical<Kleenean>(Less(),x1,x2); }
Quasidecidable operator> (Real const& x1, Real const& x2) { return make_logical<Kleenean>(Gtr(),x1,x2); }
Quasidecidable operator<=(Real const& x1, Real const& x2) { return make_logical<Kleenean>(Leq(),x1,x2); }
Quasidecidable operator>=(Real const& x1, Real const& x2) { return make_logical<Kleenean>(Geq(),x1,x2); }

Kleenean sgn(Real const& x) { return make_logical<Kleenean>(Sgn(),x); }

ValidatedNegatedSierpinskian operator==(Real const& x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
ValidatedSierpinskian operator!=(Real const& x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Kleenean operator< (Real const& x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Kleenean operator> (Real const& x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Kleenean operator<=(Real const& x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Kleenean operator>=(Real const& x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }

Integer round(Real const& r) {
    DyadicBounds wb=r.compute(Effort(0)).get();
    Dyadic wc=hlf(wb.lower()+wb.upper());
    return round(wc);
}

template<> String class_name<Real>() { return "Real"; }
template<> String class_name<PositiveReal>() { return "PositiveReal"; }

const Real pi = 4*atan(1); //Real(3.1415926535897930, 3.141592653589793238, 3.1415926535897936);
const Real infinity = Real(std::numeric_limits<double>::infinity());

FloatDPBounds Real::get(DoublePrecision pr) const {
    return this->_ptr->_compute_get(pr);
}

FloatMPBounds Real::get(MultiplePrecision pr) const {
    return this->_ptr->_compute_get(pr);
}

ValidatedReal Real::compute(Effort eff) const {
    return this->_ptr->_compute(eff);
    return ValidatedReal(DyadicBounds(this->_ptr->_compute_get(MultiplePrecision(eff.work()))));
}

ValidatedReal Real::compute(Accuracy accuracy) const {
    Nat effort=1;
    Nat acc=accuracy.bits();
    MultiplePrecision precision(effort*64);
    FloatMPError error_bound(FloatMP(two^-acc,upward,precision));
    FloatMPError error=2u*error_bound;
    FloatMPBounds res;
    while (!(error.raw()<error_bound.raw())) {
        res=this->get(precision);
        error=res.error();
        effort+=1;
        precision=MultiplePrecision(effort*64);
    }
    return ValidatedReal(DyadicBounds(res));
}

ValidatedKleenean check_sgn(Real r, Effort eff) {
    auto x = r.get(MultiplePrecision(eff));
    if(definitely(x>0)) { return true; }
    else if(definitely(x<0)) { return false; }
    else { return indeterminate; }
}



LowerReal::LowerReal(SharedPointer<RealInterface> p) : _ptr(p) {
}

LowerReal::LowerReal(Real r) : _ptr(r._ptr) {
}

FloatDPLowerBound LowerReal::operator() (DoublePrecision pr) const {
    return this->_ptr->_compute_get(pr);
}

FloatMPLowerBound LowerReal::operator() (MultiplePrecision pr) const {
    return this->_ptr->_compute_get(pr);
}

FloatDPLowerBound LowerReal::get(DoublePrecision pr) const {
    return this->_ptr->_compute_get(pr);
}

FloatMPLowerBound LowerReal::get(MultiplePrecision pr) const {
    return this->_ptr->_compute_get(pr);
}

UpperReal::UpperReal(SharedPointer<RealInterface> p) : _ptr(p) {
}

UpperReal::UpperReal(Real r) : _ptr(r._ptr) {
}

FloatDPUpperBound UpperReal::operator() (DoublePrecision pr) const {
    return this->_ptr->_compute_get(pr);
}

FloatMPUpperBound UpperReal::operator() (MultiplePrecision pr) const {
    return this->_ptr->_compute_get(pr);
}

FloatDPUpperBound UpperReal::get(DoublePrecision pr) const {
    return this->_ptr->_compute_get(pr);
}

FloatMPUpperBound UpperReal::get(MultiplePrecision pr) const {
    return this->_ptr->_compute_get(pr);
}

inline Real const& cast_real(LowerReal const& lr) { return reinterpret_cast<Real const&>(lr); }
inline Real const& cast_real(UpperReal const& ur) { return reinterpret_cast<Real const&>(ur); }
inline Real const& make_signed(PositiveReal const& pr) { return pr; }
inline LowerReal const& make_lower(Real const& r) { return reinterpret_cast<LowerReal const&>(r); }
inline UpperReal const& make_upper(Real const& r) { return reinterpret_cast<UpperReal const&>(r); }

PositiveReal cast_positive(Real const& pr) { return static_cast<PositiveReal const&>(pr); }

LowerReal max(LowerReal const& lr1, LowerReal const& lr2) { return make_lower(max(cast_real(lr1),cast_real(lr2))); }
LowerReal min(LowerReal const& lr1, LowerReal const& lr2) { return make_lower(min(cast_real(lr1),cast_real(lr2))); }
Real min(LowerReal const& lr1, Real const& r2) { return min(cast_real(lr1),r2); }
Real min(Real const& r1, LowerReal const& lr2) { return min(r1,cast_real(lr2)); }

UpperReal max(UpperReal const& ur1, UpperReal const& ur2) { return make_upper(max(cast_real(ur1),cast_real(ur2))); }
Real max(UpperReal const& ur1, Real const& r2) { return max(cast_real(ur1),r2); }
Real max(Real r1, UpperReal const& ur2) { return max(r1,cast_real(ur2)); }
UpperReal min(UpperReal const& ur1, UpperReal const& ur2) { return make_upper(min(cast_real(ur1),cast_real(ur2))); }

LowerReal neg(UpperReal const& ur) { return make_lower(neg(cast_real(ur))); }
UpperReal neg(LowerReal const& lr) { return make_upper(neg(cast_real(lr))); }
LowerReal add(LowerReal const& lr1, LowerReal const& lr2) { return make_lower(add(cast_real(lr1),cast_real(lr2))); }
UpperReal add(UpperReal const& ur1, UpperReal const& ur2) { return make_upper(add(cast_real(ur1),cast_real(ur2))); }

PositiveFloatDPBounds PositiveReal::get(DoublePrecision pr) const {
    return PositiveFloatDPBounds(this->_ptr->_compute_get(pr));
}

PositiveFloatMPBounds PositiveReal::get(MultiplePrecision pr) const {
    return PositiveFloatMPBounds(this->_ptr->_compute_get(pr));
}

PositiveReal max(PositiveReal const& pr1, PositiveReal const& pr2) { return cast_positive(max(make_signed(pr1),make_signed(pr2))); }
PositiveReal min(PositiveReal const& pr1, PositiveReal const& pr2) { return cast_positive(min(make_signed(pr1),make_signed(pr2))); }
PositiveReal add(PositiveReal const& pr1, PositiveReal const& pr2) { return cast_positive(add(make_signed(pr1),make_signed(pr2))); }
PositiveReal mul(PositiveReal const& pr1, PositiveReal const& pr2) { return cast_positive(mul(make_signed(pr1),make_signed(pr2))); }
PositiveReal div(PositiveReal const& pr1, PositiveReal const& pr2) { return cast_positive(div(make_signed(pr1),make_signed(pr2))); }
PositiveReal rec(PositiveReal const& pr) { return cast_positive(rec(make_signed(pr))); }
PositiveReal sqrt(PositiveReal const& pr) { return cast_positive(sqrt(make_signed(pr))); }
PositiveReal atan(PositiveReal const& pr) { return cast_positive(atan(make_signed(pr))); }


PositiveFloatDPLowerBound PositiveLowerReal::get(DoublePrecision pr) const {
    return PositiveFloatDPLowerBound(this->_ptr->_compute_get(pr));
}

PositiveFloatMPLowerBound PositiveLowerReal::get(MultiplePrecision pr) const {
    return PositiveFloatMPLowerBound(this->_ptr->_compute_get(pr));
}

PositiveFloatDPUpperBound PositiveUpperReal::get(DoublePrecision pr) const {
    return PositiveFloatDPUpperBound(this->_ptr->_compute_get(pr));
}

PositiveFloatMPUpperBound PositiveUpperReal::get(MultiplePrecision pr) const {
    return PositiveFloatMPUpperBound(this->_ptr->_compute_get(pr));
}

PositiveUpperReal rec(PositiveLowerReal plr) { return cast_positive(rec(cast_real(plr))); }
PositiveLowerReal rec(PositiveUpperReal pur) { return cast_positive(rec(cast_real(pur))); }
PositiveLowerReal add(PositiveLowerReal plr1, PositiveLowerReal plr2) { return cast_positive(add(cast_real(plr1),cast_real(plr2))); }
PositiveUpperReal add(PositiveUpperReal pur1, PositiveUpperReal pur2) { return cast_positive(add(cast_real(pur1),cast_real(pur2))); }
PositiveLowerReal mul(PositiveLowerReal plr1, PositiveLowerReal plr2) { return cast_positive(mul(cast_real(plr1),cast_real(plr2))); }
PositiveUpperReal mul(PositiveUpperReal pur1, PositiveUpperReal pur2) { return cast_positive(mul(cast_real(pur1),cast_real(pur2))); }
PositiveLowerReal div(PositiveLowerReal plr1, PositiveUpperReal pur2) { return cast_positive(div(cast_real(plr1),cast_real(pur2))); }
PositiveUpperReal div(PositiveUpperReal pur1, PositiveLowerReal plr2) { return cast_positive(div(cast_real(pur1),cast_real(plr2))); }

LowerReal mul(LowerReal lr1, PositiveReal pr2) { return mul(cast_real(lr1),make_signed(pr2)); }
UpperReal mul(UpperReal ur1, PositiveReal pr2) { return mul(cast_real(ur1),make_signed(pr2)); }
LowerReal mul(PositiveReal pr1, LowerReal lr2) { return mul(make_signed(pr1),cast_real(lr2)); }
UpperReal mul(PositiveReal pr1, UpperReal ur2) { return mul(make_signed(pr1),cast_real(ur2)); }
LowerReal div(LowerReal lr1, PositiveReal pr2) { return div(cast_real(lr1),make_signed(pr2)); }
UpperReal div(UpperReal ur1, PositiveReal pr2) { return div(cast_real(ur1),make_signed(pr2)); }
LowerReal div(PositiveReal pr1, UpperReal ur2) { return div(make_signed(pr1),cast_real(ur2)); }
UpperReal div(PositiveReal pr1, LowerReal lr2) { return div(make_signed(pr1),cast_real(lr2)); }


static_assert(IsConstructible<FloatDP,Dyadic,FloatDP::RoundingModeType,FloatDP::PrecisionType>::value,"");
static_assert(IsConstructible<FloatMP,Dyadic,FloatMP::RoundingModeType,FloatMP::PrecisionType>::value,"");

ValidatedReal::ValidatedReal(DyadicBounds const& y) : _ptr(std::make_shared<ValidatedRealWrapper<DyadicBounds>>(y)) { }
DyadicBounds ValidatedReal::get() const { return this->_ptr->_get(); }
FloatDPBounds ValidatedReal::get(DoublePrecision pr) const { return this->_ptr->_get(pr); }
FloatMPBounds ValidatedReal::get(MultiplePrecision pr) const { return this->_ptr->_get(pr); }
OutputStream& operator<<(OutputStream& os, ValidatedReal const& vr) { return vr._ptr->_write(os); }

} // namespace Ariadne
