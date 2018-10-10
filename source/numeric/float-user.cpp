/***************************************************************************
 *            float-user.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "../config.hpp"
#include "../utility/macros.hpp"
#include "../utility/exceptions.hpp"

#include "../numeric/float-user.hpp"

#include "../numeric/integer.hpp"
#include "../numeric/dyadic.hpp"
#include "../numeric/decimal.hpp"
#include "../numeric/rational.hpp"
#include "../numeric/real.hpp"

#include "../numeric/number_wrapper.hpp"

#include "float-user.tpl.hpp"

#include "float_error.tpl.hpp"
#include "float_value.tpl.hpp"
#include "float_ball.tpl.hpp"
#include "float_bounds.tpl.hpp"
#include "float_upper_bound.tpl.hpp"
#include "float_lower_bound.tpl.hpp"
#include "float_approximation.tpl.hpp"

namespace Ariadne {

FloatDP set(FloatMP const& x, RoundUpward rnd, DoublePrecision pr) { return FloatDP(Dyadic(x),rnd,pr); }
FloatDP set(FloatMP const& x, RoundDownward rnd, DoublePrecision pr) { return FloatDP(Dyadic(x),rnd,pr); }


const FloatDPValue infty = FloatDPValue(FloatDP::inf(dp));

OutputStream& operator<<(OutputStream& os, Rounding const& rnd) {
    return os << ( rnd._rbp == ROUND_TO_NEAREST ? "near" : (rnd._rbp == ROUND_DOWNWARD ? "down" : "up") ); }


FloatError<DoublePrecision> operator"" _error(long double lx) {
    double x=lx;
    assert(x==lx);
    return FloatError<DoublePrecision>(FloatDP(x));
}

FloatValue<DoublePrecision> operator"" _exact(long double lx) {
    double x=lx;
    assert(x==lx);
    return FloatValue<DoublePrecision>(x);
}

FloatBall<DoublePrecision> operator"" _near(long double lx) {
    volatile double x=lx;
    volatile long double le=std::abs((long double)x-lx);
    volatile double e=le;
    while(e<le) { e*=(1+std::numeric_limits<double>::epsilon()); }
    return FloatBall<DoublePrecision>(x,e);
}

FloatUpperBound<DoublePrecision> operator"" _upper(long double lx) {
    static const double eps = std::numeric_limits<double>::epsilon();
    static const double min = std::numeric_limits<double>::min();
    double x=lx;
    if(x<lx) { x+=min; }
    while (x<lx) { x+=std::abs(x)*eps; }
    return FloatUpperBound<DoublePrecision>(x);
}

FloatLowerBound<DoublePrecision> operator"" _lower(long double lx) {
    static const double eps = std::numeric_limits<double>::epsilon();
    static const double min = std::numeric_limits<double>::min();
    double x=lx;
    if(x>lx) { x-=min; }

    while (x>lx) { x-=std::abs(x)*eps; }
    return FloatLowerBound<DoublePrecision>(x);
}

FloatApproximation<DoublePrecision> operator"" _approx(long double lx) {
    double x=lx;
    return FloatApproximation<DoublePrecision>(x);
}


inline int log10floor(double const& x) { return std::max(std::floor(std::log10(x)),-65280.); }
inline int log10floor(FloatMP const& x) { return log10floor(x.get_d()); }
inline int abslog10floor(double const& x) { return log10floor(std::abs(x)); }

template<> OutputStream& Operations<FloatBounds<MultiplePrecision>>::_write(OutputStream& os, const FloatBounds<MultiplePrecision>& x)
{
    static const double log2ten = 3.3219280948873621817;
    using std::max; using std::min;
    FloatMP const& l=x.lower_raw();
    FloatMP const& u=x.upper_raw();
    double ldbl=l.get_d();
    double udbl=u.get_d();
    if(ldbl==0.0 && udbl==0.0) { return os << "0.0[:]"; }
    int errplc=static_cast<int>(FloatError<MultiplePrecision>::output_places);
    //int bndplc=FloatBounds<MultiplePrecision>::output_places;
    int precplc=x.precision()/log2ten;
    int log10wdth=log10floor(sub(to_nearest,u,l));
    int log10mag=log10floor(max(-ldbl,udbl));
    int dgtswdth=errplc-(log10wdth+1); // Digits appropriate given width of interval
    //int dgtsbnd=bndplc-(log10mag+1); // Digits appropriate given asked-for precision of bounded objects
    int dgtsprec=precplc-(log10mag+1); // Digits appropriate given precision of objects
    Nat dgts=static_cast<Nat>(max(min(dgtswdth,dgtsprec),1));
    DecimalPlaces plcs{dgts};
    String lstr=print(l,plcs,MPFR_RNDD);
    String ustr=print(u,plcs,MPFR_RNDU);
    auto lcstr=lstr.c_str();
    auto ucstr=ustr.c_str();
    size_t cpl=0;
    if(ldbl*udbl>=0 && abslog10floor(ldbl)==abslog10floor(udbl)) {
        while(lcstr[cpl]!='\0' && lcstr[cpl]==ustr[cpl]) { ++cpl; }
    }
    char ocstr[1024];
    ocstr[0]='\0';
    strncat(ocstr,lcstr,cpl);
    strcat(ocstr,"[");
    strcat(ocstr,lcstr+cpl);
    strcat(ocstr,":");
    strcat(ocstr,ucstr+cpl);
    strcat(ocstr,"]");
    return os << ocstr;
}

template<> OutputStream& Operations<FloatBounds<DoublePrecision>>::_write(OutputStream& os, const FloatBounds<DoublePrecision>& x)
{
    MultiplePrecision prec(64);
    return os << FloatBounds<MultiplePrecision>(FloatMP(x.lower_raw(),prec),FloatMP(x.upper_raw(),prec));
}



template<> OutputStream& Operations<FloatBall<MultiplePrecision>>::_write(OutputStream& os, FloatBall<MultiplePrecision> const& x) {
    // Write based on number of correct digits
    static const double log2ten = 3.3219280948873621817;
    static const char pmstr[] = "\u00b1";
    static const char hlfstr[] = "\u00bd";
    FloatMP const& v=x.value_raw();
    FloatMP const& e=x.error_raw();
    double edbl=e.get_d();
    // Compute the number of decimal places to be displayed
    Nat errplc = static_cast<Nat>(FloatError<MultiplePrecision>::output_places);
    Nat log10err = static_cast<Nat>(log10floor(edbl));
    Nat dgtserr = errplc-(log10err+1);
    Nat dgtsval = (x.value().raw()==0) ? dgtserr : std::floor((x.value().precision()+1-x.value().raw().exponent())/log2ten);
    Nat dgts = std::max(std::min(dgtsval,dgtserr),errplc);
    if(edbl==0.0) { dgts = dgtsval; }
    DecimalPlaces plcs{dgts};
    // Get string version of mpfr values
    String vstr=print(v,plcs,MPFR_RNDN);
    String estr=print(e,plcs,MPFR_RNDU);

    // Find position of first significan digit of error
    auto vcstr=vstr.c_str(); auto ecstr=estr.c_str();
    size_t cpl=0;
    if(edbl==0.0) {
        cpl=std::strlen(vcstr);
    } else if(edbl<1.0) {
        const char* vptr = std::strchr(vcstr,'.');
        const char* eptr = std::strchr(ecstr,'.');
        ++vptr; ++eptr;
        while((*eptr)=='0') { ++eptr; ++vptr; }
        cpl = static_cast<size_t>(vptr-vcstr);
    }

    // Chop and catenate strings
    static const size_t buf_sz = 1024;
    char ocstr[buf_sz];
    ocstr[0]='\0';
    std::strncat(ocstr,vcstr,cpl);
    std::strcat(ocstr,"[");
    std::strcat(ocstr,vcstr+cpl);
    std::strcat(ocstr,pmstr);
    std::strcat(ocstr,ecstr+cpl);
    std::strcat(ocstr,hlfstr);
    std::strcat(ocstr,"]");
    return os << ocstr;

    return os << x.value() << "\u00b1" << x.error();
}

template<> OutputStream& Operations<FloatBall<MultiplePrecision,DoublePrecision>>::_write(OutputStream& os, FloatBall<MultiplePrecision,DoublePrecision> const& x) {
    MultiplePrecision prec(64);
    return os << FloatBall<MultiplePrecision,MultiplePrecision>(x.value_raw(),FloatMP(x.error_raw(),prec));
}

template<> OutputStream& Operations<FloatBall<DoublePrecision>>::_write(OutputStream& os, FloatBall<DoublePrecision> const& x) {
    MultiplePrecision prec(64);
    return os << FloatBall<MultiplePrecision>(FloatMP(x.value_raw(),prec),FloatMP(x.error_raw(),prec));
}


template<> Nat integer_cast<Nat,FloatDPApproximation>(FloatDPApproximation const& x) {
    return std::round(x.get_d()); }

template<> Int integer_cast<Int,FloatDPApproximation>(FloatDPApproximation const& x) {
    return std::round(x.get_d()); }

template<> Nat integer_cast<Nat,FloatDPLowerBound>(FloatDPLowerBound const& x) {
    return std::round(x.get_d()); }

template<> Int integer_cast<Int,FloatDPLowerBound>(FloatDPLowerBound const& x) {
    return std::round(x.get_d()); }

template<> Int integer_cast<Int,FloatDPBounds>(FloatDPBounds const& x) {
    return std::round((x.lower().get_d()+x.upper().get_d())/2); }

template<> Int integer_cast<Int,FloatMPBounds>(FloatMPBounds const& x) {
    return std::round((x.lower().get_d()+x.upper().get_d())/2); }

template<> Nat integer_cast<Nat,FloatMPApproximation>(FloatMPApproximation const& x) {
    return std::round(x.get_d()); }
template<> Int integer_cast<Int,FloatMPApproximation>(FloatMPApproximation const& x) {
    return std::round(x.get_d()); }



template class Approximation<FloatDP>;
template class LowerBound<FloatDP>;
template class UpperBound<FloatDP>;
template class Bounds<FloatDP>;
template class Ball<FloatDP>;
template class Value<FloatDP>;

template class Approximation<FloatMP>;
template class LowerBound<FloatMP>;
template class UpperBound<FloatMP>;
template class Bounds<FloatMP>;
template class Ball<FloatMP>;
template class Value<FloatMP>;

template class Ball<FloatMP, FloatDP>;

template class Operations<FloatDPApproximation>;
template class Operations<FloatDPLowerBound>;
template class Operations<FloatDPUpperBound>;
template class Operations<FloatDPBounds>;
template class Operations<FloatDPBall>;
template class Operations<FloatDPValue>;
template class Operations<PositiveFloatDPApproximation>;
template class Operations<PositiveFloatDPLowerBound>;
template class Operations<PositiveFloatDPUpperBound>;
template class Operations<PositiveFloatDPBounds>;
template class Operations<FloatDPError>;

template class Operations<FloatMPApproximation>;
template class Operations<FloatMPLowerBound>;
template class Operations<FloatMPUpperBound>;
template class Operations<FloatMPBounds>;
template class Operations<FloatMPBall>;
template class Operations<FloatMPValue>;
template class Operations<PositiveFloatMPApproximation>;
template class Operations<PositiveFloatMPLowerBound>;
template class Operations<PositiveFloatMPUpperBound>;
template class Operations<PositiveFloatMPBounds>;
template class Operations<FloatMPError>;

template class Operations<FloatMDPBall>;


template<> String class_name<FloatDPApproximation>() { return "FloatDPApproximation"; }
template<> String class_name<FloatDPLowerBound>() { return "FloatDPLowerBound"; }
template<> String class_name<FloatDPUpperBound>() { return "FloatDPUpperBound"; }
template<> String class_name<FloatDPBounds>() { return "FloatDPBounds"; }
template<> String class_name<FloatDPBall>() { return "FloatDPBall"; }
template<> String class_name<FloatDPValue>() { return "FloatDPValue"; }
template<> String class_name<FloatMPApproximation>() { return "FloatMPApproximation"; }
template<> String class_name<FloatMPLowerBound>() { return "FloatMPLowerBound"; }
template<> String class_name<FloatMPUpperBound>() { return "FloatMPUpperBound"; }
template<> String class_name<FloatMPBounds>() { return "FloatMPBounds"; }
template<> String class_name<FloatMPBall>() { return "FloatMPBall"; }
template<> String class_name<FloatMPValue>() { return "FloatMPValue"; }

template<> String class_name<FloatMDPBall>() { return "FloatMDPBall"; }

} // namespace Ariadne
