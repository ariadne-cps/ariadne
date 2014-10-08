/***************************************************************************
 *            real.cc
 *
 *  Copyright 2008-10  Pieter Collins
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

#include "standard.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <limits>



#include "config.h"
#include "real.h"
#include "operators.h"

#include "float.h"
#include "float-approximate.h"
#include "float-validated.h"
#include "float-exact.h"

namespace Ariadne {

namespace{
static const double _pi_up=3.1415926535897936;
static const double _pi_approx=3.1415926535897931;
static const double _pi_down=3.1415926535897931;
static const double _infinity=std::numeric_limits<double>::infinity();
}

class RealInterface {
  public:
    virtual ~RealInterface() { }
    virtual operator ApproximateFloat () const = 0;
    virtual operator ValidatedFloat () const = 0;
    virtual Void write(OutputStream& os) const = 0;
};

class RealBody : public RealInterface {
  private:
    ValidatedFloat _ivl;
    ApproximateFloat _flt;
  protected:
    RealBody(const ValidatedFloat& ivl, const ApproximateFloat& flt) : _ivl(ivl), _flt(flt) {
        ARIADNE_ASSERT((ivl.lower_value()<=flt.value()) && (flt.value()<=ivl.upper_value())); }
    RealBody(const ValidatedFloat& ivl) : _ivl(ivl), _flt(midpoint(ivl)) { }
  public:
    virtual operator ApproximateFloat () const final { return _flt; }
    virtual operator ValidatedFloat () const final { return _ivl; }
    virtual Void write(OutputStream& os) const = 0;
};

class RealConstant
    : public RealBody
{
    String _name;
  public:
    RealConstant(const String& name, double lower, double nearest, double upper)
        : RealBody(ValidatedFloat(lower,upper),ApproximateFloat(nearest)), _name(name) { }
    RealConstant(const String& name, ValidatedFloat bounds, ApproximateFloat approx) : RealBody(bounds,approx), _name(name) { }
    RealConstant(const String& name, ValidatedFloat bounds) : RealBody(bounds), _name(name) { }
    virtual Void write(OutputStream& os) const final { os << _name; }
};

class IntervalReal
    : public RealBody
{
  public:
    IntervalReal(double l, double u) : RealBody(ValidatedFloat(l,u),ApproximateFloat((l+u)/2)) { }
    IntervalReal(double l, double x, double u) : RealBody(ValidatedFloat(l,u),ApproximateFloat(x)) { }
    virtual Void write(OutputStream& os) const final { os << this->operator ValidatedFloat(); }
};

class IntegerReal
    : public RealBody
{
    Integer _value;
  public:
    IntegerReal(const Integer& z) : RealBody(ValidatedFloat(z),ApproximateFloat(midpoint(ValidatedFloat(z)))), _value(z) { }
    virtual Void write(OutputStream& os) const final { os << _value; }
};

#ifdef HAVE_GMPXX_H
class RationalReal
    : public RealBody
{
    Rational _value;
  public:
    RationalReal(const Rational& q) : RealBody(ValidatedFloat(q),ApproximateFloat(q)), _value(q) { }
    virtual Void write(OutputStream& os) const final { os << _value; }
};
#endif // HAVE_GMPXX_H

class DecimalReal
    : public RealBody
{
    Decimal _value;
  public:
    DecimalReal(const Decimal& d) : RealBody(ValidatedFloat(d),ApproximateFloat(d)), _value(d) { }
    virtual Void write(OutputStream& os) const final { os << _value; }
};

class DyadicReal
    : public RealBody
{
    Dyadic _value;
  public:
    DyadicReal(const Dyadic& x) : RealBody(ValidatedFloat(x),ApproximateFloat(x)), _value(x) { }
    virtual Void write(OutputStream& os) const final { os << _value; }
};

class ExactFloatReal
    : public RealBody
{
    ExactFloat _value;
  public:
    ExactFloatReal(double x) : ExactFloatReal(ExactFloat(x)) { }
    ExactFloatReal(const ExactFloat& x) : RealBody(ValidatedFloat(x),ApproximateFloat(x)), _value(x) { }
    virtual Void write(OutputStream& os) const final { os << _value; }
};

// Needed for rec
static inline ApproximateFloat operator/(int n, ApproximateFloat x) { return ApproximateFloat(n)/x; }
static inline ValidatedFloat operator/(int n, ValidatedFloat x) { return ExactFloat(n)/x; }

class UnaryReal
    : public RealBody
{
    Operator _op; Real _arg;
  public:
    UnaryReal(Operator op, const Real& arg)
        : RealBody(compute(op,ValidatedFloat(arg)),compute(op,ApproximateFloat(arg))), _op(op), _arg(arg) { }
    virtual Void write(OutputStream& os) const final { os << _op << "(" << _arg << ")"; }
};

class BinaryReal
    : public RealBody
{
    Operator _op; Real _arg1; Real _arg2;
  public:
    BinaryReal(Operator op, const Real& arg1, const Real& arg2)
        : RealBody(compute(op,ValidatedFloat(arg1),ValidatedFloat(arg2)),compute(op,ApproximateFloat(arg1),ApproximateFloat(arg2)))
        , _op(op), _arg1(arg1), _arg2(arg2) { }
    virtual Void write(OutputStream& os) const final { os << _op << "(" << _arg1 << "," << _arg2 << ")"; }
};

Real::Real(RealInterface* raw_ptr) : _ptr(raw_ptr) { }

Real::~Real() { }
Real::Real() : _ptr(new IntegerReal(0)) { }
Real::Real(double l, double u) : _ptr(new IntervalReal(l,u)) { }
Real::Real(double l, double x, double u) : _ptr(new IntervalReal(l,x,u)) { }

Real::Real(unsigned int m) : _ptr(new IntegerReal(m)) { }
Real::Real(int n) : _ptr(new IntegerReal(n)) { }
Real::Real(double x) : _ptr(new ExactFloatReal(x)) { }
Real::Real(const Dyadic& d) : _ptr(new DyadicReal(d)) { }
Real::Real(const Decimal& d) : _ptr(new DecimalReal(d)) { }
#ifdef HAVE_GMPXX_H
Real::Real(const Integer& z) : _ptr(new IntegerReal(z)) { }
Real::Real(const Rational& q) : _ptr(new RationalReal(q)) { }
#endif // HAVE_GMPXX_H
Real::Real(const ExactFloat& x) : _ptr(new ExactFloatReal(x)) { }

Real::Real(const Real& x) : _ptr(x._ptr) { }
Real& Real::operator=(const Real& x) { this->_ptr=x._ptr; return *this; }
double Real::get_d() const { return this->_ptr->operator ApproximateFloat().get_d(); }

Real::operator UpperFloat() const { return (this->_ptr->operator ValidatedFloat()).upper(); }
Real::operator ValidatedFloat() const { return this->_ptr->operator ValidatedFloat(); }
Real::operator ApproximateFloat() const { return this->_ptr->operator ApproximateFloat(); }

ApproximateFloat::ApproximateFloat(const Real& x) : ApproximateFloat(x._ptr->operator ApproximateFloat()) { }
ValidatedFloat::ValidatedFloat(const Real& x) : ValidatedFloat(x._ptr->operator ValidatedFloat()) { }

ExactInterval::ExactInterval(const Real& x) : ExactInterval(ValidatedFloat(x)) { }

Real _make_real(const ExactInterval& ivl) { return Real(ivl.lower().get_d(),ivl.upper().get_d()); }

Real operator+(const Real& x) { return Real(new UnaryReal(POS,x)); }
Real operator-(const Real& x) { return Real(new UnaryReal(NEG,x)); }
Real operator+(const Real& x, const Real& y) { return Real(new BinaryReal(ADD,x,y)); }
Real operator-(const Real& x, const Real& y) { return Real(new BinaryReal(SUB,x,y)); }
Real operator*(const Real& x, const Real& y) { return Real(new BinaryReal(MUL,x,y)); }
Real operator/(const Real& x, const Real& y) { return Real(new BinaryReal(DIV,x,y)); }

Float mag(const Real& x) { ARIADNE_NOT_IMPLEMENTED; }

const Real pi=Real(new RealConstant("pi",_pi_down,_pi_approx,_pi_up));
const Real infinity=Real(new RealConstant("inf",inf,inf,inf));

Real abs(const Real& x) { return Real(new UnaryReal(ABS,x)); }
Real pos(const Real& x) { return Real(new UnaryReal(POS,x)); }
Real neg(const Real& x) { return Real(new UnaryReal(NEG,x)); }
Real sqr(const Real& x) { return Real(new UnaryReal(SQR,x)); }
Real rec(const Real& x) { return Real(new UnaryReal(REC,x)); }
Real add(const Real& x, const Real& y) { return Real(new BinaryReal(ADD,x,y)); }
Real sub(const Real& x, const Real& y) { return Real(new BinaryReal(SUB,x,y)); }
Real mul(const Real& x, const Real& y) { return Real(new BinaryReal(MUL,x,y)); }
Real div(const Real& x, const Real& y) { return Real(new BinaryReal(DIV,x,y)); }
Real pow(const Real& x, uint m) { return Real(new BinaryReal(POW,x,m)); }
Real pow(const Real& x, int n) { return Real(new BinaryReal(POW,x,n)); }
Real sqrt(const Real& x) { return Real(new UnaryReal(SQRT,x)); }
Real exp(const Real& x) { return Real(new UnaryReal(EXP,x)); }
Real log(const Real& x) { return Real(new UnaryReal(LOG,x)); }
Real sin(const Real& x) { return Real(new UnaryReal(SIN,x)); }
Real cos(const Real& x) { return Real(new UnaryReal(COS,x)); }
Real tan(const Real& x) { return Real(new UnaryReal(TAN,x)); }
Real asin(const Real& x) { return Real(new UnaryReal(ASIN,x)); }
Real acos(const Real& x) { return Real(new UnaryReal(ACOS,x)); }
Real atan(const Real& x) { return Real(new UnaryReal(ATAN,x)); }



std::ostream& operator<<(std::ostream& os, const Real& x)
{
    x._ptr->write(os);
    return os;

}




} // namespace Ariadne

