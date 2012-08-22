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

#include <iostream>
#include <iomanip>
#include <cassert>
#include <limits>

#include "config.h"
#include "real.h"
#include "operators.h"


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
    virtual operator Float () const = 0;
    virtual operator Interval () const = 0;
    virtual Void write(OutputStream& os) const = 0;
};

class RealBody : public RealInterface {
  private:
    Interval _ivl;
    Float _flt;
  protected:
    RealBody(const Interval& ivl, const Float& flt) : _ivl(ivl), _flt(flt) { ARIADNE_ASSERT(contains(ivl,flt)); }
    RealBody(const Interval& ivl) : _ivl(ivl), _flt(midpoint(ivl)) { }
  public:
    virtual operator Float () const final { return _flt; }
    virtual operator Interval () const final { return _ivl; }
    virtual Void write(OutputStream& os) const = 0;
};

class RealConstant
    : public RealBody
{
    String _name;
  public:
    RealConstant(const String& name, double lower, double nearest, double upper)
        : RealBody(Interval(lower,upper),Float(nearest)), _name(name) { }
    RealConstant(const String& name, Interval bounds, Float approx) : RealBody(bounds,approx), _name(name) { }
    RealConstant(const String& name, Interval bounds) : RealBody(bounds), _name(name) { }
    virtual Void write(OutputStream& os) const final { os << _name; }
};

class IntervalReal
    : public RealBody
{
  public:
    IntervalReal(double l, double u) : RealBody(Interval(l,u),Float((l+u)/2)) { }
    IntervalReal(double l, double x, double u) : RealBody(Interval(l,u),Float(x)) { }
    virtual Void write(OutputStream& os) const final { os << this->operator Interval(); }
};

class IntegerReal
    : public RealBody
{
    Integer _value;
  public:
    IntegerReal(const Integer& z) : RealBody(Interval(z),Float(midpoint(Interval(z)))), _value(z) { }
    virtual Void write(OutputStream& os) const final { os << _value; }
};

class RationalReal
    : public RealBody
{
    Rational _value;
  public:
    RationalReal(const Rational& q) : RealBody(Interval(q),Float(q)), _value(q) { }
    virtual Void write(OutputStream& os) const final { os << _value; }
};

class FloatReal
    : public RealBody
{
    ExactFloat _value;
  public:
    FloatReal(double x) : RealBody(Interval(x),Float(x)), _value(x) { }
    FloatReal(const ExactFloat& x) : RealBody(Interval(x),Float(x)), _value(x) { }
    virtual Void write(OutputStream& os) const final { os << _value; }
};

class UnaryReal
    : public RealBody
{
    Operator _op; Real _arg;
  public:
    UnaryReal(Operator op, const Real& arg)
        : RealBody(compute(op,Interval(arg)),compute(op,Float(arg))), _op(op), _arg(arg) { }
    virtual Void write(OutputStream& os) const final { os << _op << "(" << _arg << ")"; }
};

class BinaryReal
    : public RealBody
{
    Operator _op; Real _arg1; Real _arg2;
  public:
    BinaryReal(Operator op, const Real& arg1, const Real& arg2)
        : RealBody(compute(op,Interval(arg1),Interval(arg2)),compute(op,Float(arg1),Float(arg2)))
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
Real::Real(double x) : _ptr(new FloatReal(x)) { }
Real::Real(const ExactFloat& x) : _ptr(new FloatReal(x)) { }
Real::Real(const Real& x) : _ptr(x._ptr) { }
Real& Real::operator=(const double& x) { *this=Real(x); return *this; }
Real& Real::operator=(const ExactFloat& x) { *this=Real(x); return *this; }
Real& Real::operator=(const Real& x) { this->_ptr=x._ptr; return *this; }
double Real::get_d() const { return this->_ptr->operator Float().get_d(); }

Float::Float(const Real& x) : v(x._ptr->operator Float().get_d()) { }
Interval::Interval(const Real& x) : l(x._ptr->operator Interval().lower()), u(x._ptr->operator Interval().lower()) { }
Float& Float::operator=(const Real& x) { *this=Float(x); return *this; }
Interval& Interval::operator=(const Real& x) { *this=Interval(x); return *this; }

Real _make_real(const Interval& ivl) { return Real(ivl.lower().get_d(),ivl.upper().get_d()); }

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



#ifdef HAVE_GMPXX_H

Real::Real(const std::string& str)
{
    Rational q;
    bool decimal_point=false;
    uint decimal_places=0;
    int sign=1;
    const char* c_ptr=str.c_str();
    while(*c_ptr != 0) {
        const char& c=*c_ptr;
        if(c=='.') {
            if(decimal_point) {
                ARIADNE_THROW(std::runtime_error,"Real(String)","real literal \""<<str<<"\" has more than one decimal point.");
            }
            else {
                decimal_point=true;
            }
        } else if(c>='0' && c<='9') {
            q=q*10+(c-'0');
            if(decimal_point) {
                ++decimal_places;
            }
        } else if(c=='-' && c_ptr==str.c_str()) {
            sign=-1;
        } else {
            ARIADNE_THROW(std::runtime_error,"Real(String)","invalid symbol '"<<c<<"' in string literal \""<<str<<"\"");
        }
        ++c_ptr;
    }
    for(uint i=0; i!=decimal_places; ++i) {
        q=q/10;
    }
    *this=Real(sign*q);
}

Real::Real(const Rational& q)
{
    rounding_mode_t rnd=get_rounding_mode();
    double x=q.get_d();
    volatile double ml=-x;
    volatile double u=x;
    set_rounding_upward();
    while(-ml>static_cast<const mpq_class&>(q)) {
        ml+=std::numeric_limits<double>::min();
    }
    while(u<static_cast<const mpq_class&>(q)) {
        u+=std::numeric_limits<double>::min();
    }
    *this=Real(-ml,x,u);
    set_rounding_mode(rnd);
}

#else

Real::Real(const std::string& str)
{
    ARIADNE_THROW(std::runtime_error,"Need GMP library to convert string literal to Real.");
}

#endif



} // namespace Ariadne

