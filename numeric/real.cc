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

namespace Ariadne {

typedef Real::Interface RealInterface;

struct Real::Interface {
  public:
    virtual ~Interface() = default;
    virtual BoundedFloat64 _value() const = 0;
    virtual BoundFloat64 _evaluate(Precision64) const = 0;
    virtual BoundFloatMP _evaluate(PrecisionMP) const = 0;
  public:
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

template<class O, class... AS> struct RealWrapper : virtual RealInterface, ExpressionTemplate<O,AS...>, BoundedFloat64 {
    RealWrapper(O o, AS... as) : ExpressionTemplate<O,AS...>(o,as...)
        , BoundedFloat64(static_cast<ExpressionTemplate<O,AS...>const&>(*this).operator BoundedFloat64()) { }
    virtual BoundedFloat64 _value() const { return static_cast<BoundedFloat64 const&>(*this); }
    //virtual BoundFloat64 _evaluate(Precision64 pr) const { ExpressionTemplate<O,AS...> const& self=*this; return self(pr); }
    //virtual BoundFloatMP _evaluate(PrecisionMP pr) const { ExpressionTemplate<O,AS...> const& self=*this; return self(pr); }
    virtual BoundFloat64 _evaluate(Precision64 pr) const {  return static_cast<BoundedFloat64>(*this); }
    virtual BoundFloatMP _evaluate(PrecisionMP pr) const {  ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<ExpressionTemplate<O,AS...> const&>(*this); }
};

template<class X> struct RealConstant : RealInterface, BoundedFloat64 {
    X _c;
  public:
    RealConstant(X const& x) : BoundedFloat64(x), _c(x) { }
    virtual BoundedFloat64 _value() const { return static_cast<BoundedFloat64 const&>(*this); }
    virtual BoundFloat64 _evaluate(Precision64 pr) const { return static_cast<BoundedFloat64>(*this); }
    virtual BoundFloatMP _evaluate(PrecisionMP pr) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_c; }
};

template<> struct RealConstant<Integer> : RealInterface, BoundedFloat64 {
    typedef Integer X;
     X _c;
  public:
    RealConstant(X const& x) : BoundedFloat64(x), _c(x) { }
    virtual BoundedFloat64 _value() const { return static_cast<BoundedFloat64 const&>(*this); }
    virtual BoundFloat64 _evaluate(Precision64 pr) const { return static_cast<BoundedFloat64>(*this); }
    virtual BoundFloatMP _evaluate(PrecisionMP pr) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_c; }
};

template<> struct RealConstant<BoundedFloat64> : RealInterface, BoundedFloat64 {
    typedef BoundedFloat64 X;
  public:
    RealConstant(X const& x) : BoundedFloat64(x) { }
    virtual BoundedFloat64 _value() const { return static_cast<BoundedFloat64 const&>(*this); }
    virtual BoundFloat64 _evaluate(Precision64 pr) const { return static_cast<BoundedFloat64>(*this); }
    virtual BoundFloatMP _evaluate(PrecisionMP pr) const { ARIADNE_NOT_IMPLEMENTED; }
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
template<class PR> Float<Validated,PR>::Float(Real const& x) : Float<Validated,PR>(x.lower(),x.upper()) { }
template<class PR> Float<Upper,PR>::Float(Real const& x) : Float<Upper,PR>(x.upper()) { }
template<class PR> Float<Lower,PR>::Float(Real const& x) : Float<Lower,PR>(x.lower()) { }
template<class PR> Float<Approximate,PR>::Float(Real const& x) : Float<Approximate,PR>(x.approx()) { }
*/

template<> Float<Validated,Precision64>::Float(Real const& x) : Float<Validated,Precision64>(x.lower(),x.upper()) { }
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

Real abs(Real x) { return make_real(Abs(),x); }
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
Tribool lt(Real x1, Real x2) { return BoundedFloat64(x1)< BoundedFloat64(x2); }

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
    else { ARIADNE_NOT_IMPLEMENTED; }
    // else { PrecisionMP p=e*64; return static_cast<LogicalValue>(this->_op(this->_arg1(p),this->_arg2(p))); }
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
Tribool operator< (Real x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Tribool operator> (Real x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Tribool operator<=(Real x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }
Tribool operator>=(Real x1, Int64 n2) { ARIADNE_NOT_IMPLEMENTED; }

template<> String class_name<Real>() { return "Real"; }

const Real pi = Real(3.1415926535897930, 3.141592653589793238, 3.1415926535897936);
const Real infinity = Real(std::numeric_limits<double>::infinity());

BoundFloat64 Real::operator() (Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

BoundFloatMP Real::operator() (PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}

BoundFloatMP Real::evaluate(Accuracy accuracy) const {
    ARIADNE_NOT_IMPLEMENTED;
/*
    Nat effort=1;
    Nat acc=accuracy.bits();
    PrecisionMP precision=effort*64;
    ErrorFloatMP error_bound(two_exp(-acc).get_d(),precision);
        std::cerr << "  acc="<<acc<<" max_err=="<<error_bound<<"\n";
    ErrorFloatMP error=2*error_bound;
    BoundFloatMP res;
    while (!(error.raw()<error_bound.raw())) {
        res=(*this)(precision);
        error=res.error();
        std::cerr << "  eff="<<effort<<" prec="<<precision<<" err="<<error<<" res="<<res<<"\n";
        effort+=1;
        precision=effort*64;
    }
    return res;
*/
}




LowerFloat64 LowerReal::operator() (Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

LowerFloatMP LowerReal::operator() (PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}


UpperFloat64 UpperReal::operator() (Precision64 pr) const {
    return this->_ptr->_evaluate(pr);
}

UpperFloatMP UpperReal::operator() (PrecisionMP pr) const {
    return this->_ptr->_evaluate(pr);
}


} // namespace Ariadne

