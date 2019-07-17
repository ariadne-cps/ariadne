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

namespace Ariadne {

namespace {
inline int log10floor(double const& x) { return std::max(std::floor(std::log10(x)),-65280.); }
inline int log10floor(FloatMP const& x) { return log10floor(x.get_d()); }
inline int abslog10floor(double const& x) { return log10floor(std::abs(x)); }

template<class FE, class FLT, DisableIf<IsSame<FE,FLT>> =dummy> inline FE _make_error(FLT const& x) {
    typename FE::PrecisionType pre; return FE(Dyadic(x),upward,pre); }
template<class FE, class FLT, EnableIf<IsSame<FE,FLT>> =dummy> inline FE _make_error(FLT const& x) {
    return x; }
template<class FE, class FLT, class PRE> inline FE _make_error(FLT const& x, PRE pre) {
    return FE(Dyadic(x),upward,pre); }
}




template<class F> Nat Error<F>::output_places = 3;
template<class F> Nat Approximation<F>::output_places = 4;
template<class F> Nat Bounds<F>::output_places=8;
template<class F> Nat Value<F>::output_places = 16;

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

template<class F> LowerBound<F>::LowerBound(LowerBound<F> const& x, PR pr) : _l(x._l,downward,pr) {}
template<class F> UpperBound<F>::UpperBound(UpperBound<F> const& x, PR pr) : _u(x._u,upward,pr) {}
template<class F, class FE> Ball<F,FE>::Ball(Ball<F,FE> const& x, PR pr) : _v(x._v,near,pr), _e(x._e,up,_error_precision<PRE>(pr)) {
    F d = (this->_v>=x._v) ? sub(up,this->_v,x._v) : sub(up,x._v,this->_v); _e=add(up,_e,_make_error<FE>(d));}
template<class F, class FE> Ball<F,FE>::Ball(Bounds<F> const& x, PRE pre) : _v(x.value_raw()) , _e(_make_error<FE>(x.error_raw(),pre)) {}


template<class F> Approximation<F>::Approximation(Real const& r, PR pr) : Approximation<F>(r.get(pr)) {}
template<class F> Approximation<F>::Approximation(ApproximateNumber const& y, PR pr) : Approximation<F>(y.get(ApproximateTag(),pr)) {}
template<class F> Approximation<F>::operator ApproximateNumber() const { return ApproximateNumber(new NumberWrapper<Approximation<F>>(*this));}

template<class F> LowerBound<F>::LowerBound(Real const& r, PR pr) : LowerBound(r.get(pr)) {}
template<class F> LowerBound<F>::LowerBound(ValidatedLowerNumber const& y, PR pr) : LowerBound(y.get(LowerTag(),pr)) {}
template<class F> LowerBound<F>::operator ValidatedLowerNumber() const {
    ARIADNE_NOT_IMPLEMENTED; //return ValidatedLowerNumber(new NumberWrapper<LowerBound<F>>(*this));
}

template<class F> UpperBound<F>::UpperBound(Real const& r, PR pr) : UpperBound(r.get(pr)) {}
template<class F> UpperBound<F>::UpperBound(ValidatedUpperNumber const& y, PR pr) : UpperBound(y.get(UpperTag(),pr)) {}
template<class F> UpperBound<F>::operator ValidatedUpperNumber() const {
    ARIADNE_NOT_IMPLEMENTED; // return ValidatedUpperNumber(new NumberWrapper<UpperBound<F>>(*this));}
}

template<class F> Bounds<F>::Bounds(Real const& x, PR pr) : Bounds(x.get(pr)) {}
template<class F> Bounds<F>::Bounds(LowerBound<F> const& lower, ValidatedUpperNumber const& upper) : Bounds<F>(lower,lower.create(upper)) { }
template<class F> Bounds<F>::Bounds(ValidatedLowerNumber const& lower, UpperBound<F> const& upper) : Bounds<F>(upper.create(lower),upper) { }
template<class F> Bounds<F>::Bounds(ValidatedLowerNumber const& lower, ValidatedUpperNumber const& upper, PR pr) : Bounds<F>(lower.get(LowerTag(),pr),upper.get(UpperTag(),pr)) { }
template<class F> Bounds<F>::Bounds(ValidatedNumber const& y, PR pr) : Bounds(y.get(BoundedTag(),pr)) {}
template<class F> Bounds<F>::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<Bounds<F>>(*this));}

template<class F, class FE> Ball<F,FE>::Ball(ExactDouble const& d, PR pr) : _v(d,pr), _e(0,_error_precision<PRE>(pr)) {}
template<class F, class FE> Ball<F,FE>::Ball(TwoExp const& t, PR pr) : _v(t,pr), _e(0u,_error_precision<PRE>(pr)) {}
template<class F, class FE> Ball<F,FE>::Ball(Integer const& z, PR pr) : _v(z,near,pr), _e(abs(Dyadic(_v)-z),up,_error_precision<PRE>(pr)) {}
template<class F, class FE> Ball<F,FE>::Ball(Dyadic const& w, PR pr) : _v(w,near,pr), _e(abs(Dyadic(_v)-w),up,_error_precision<PRE>(pr)) {}
template<class F, class FE> Ball<F,FE>::Ball(Decimal const& d, PR pr) : Ball(Rational(d),pr) {}
template<class F, class FE> Ball<F,FE>::Ball(Rational const& q, PR pr) : _v(q,near,pr), _e(abs(Rational(_v)-q),up,_error_precision<PRE>(pr)) {}
template<class F, class FE> Ball<F,FE>::Ball(Dyadic const& w, PR pr, PRE pre) : _v(F(w,near,pr)), _e(abs(Dyadic(_v)-w),up,pre) {}
template<class F, class FE> Ball<F,FE>::Ball(Rational const& q, PR pr, PRE pre) : _v(F(q,near,pr)), _e(abs(Rational(_v)-q),up,pre) {}

template<class F, class FE> Ball<F,FE>::Ball(Real const& r, PR pr) : Ball(r.get(pr)) {}
template<class F, class FE> Ball<F,FE>::Ball(ValidatedNumber const& y, PR pr) : Ball(y.get(MetricTag(),pr)) {}
template<class F, class FE> Ball<F,FE>::Ball(Real const& r, PR pr, PRE pre) : Ball(r.get(pr),pre) {}
template<class F, class FE> Ball<F,FE>::Ball(ValidatedNumber const& y, PR pr, PRE pre) : Ball(y.get(MetricTag(),pr,pre)) {}
template<class F, class FE> Ball<F,FE>::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<Ball<F,FE>>(*this));}


template<class F> Value<F>::Value(ExactDouble const& d, PR pr)
    : _v(d,pr)
{
    ARIADNE_ASSERT_MSG(Dyadic(this->_v)==d,"Exact double"<<d<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

template<class F> Value<F>::Value(TwoExp const& t, PR pr)
    : _v(t,pr)
{
    ARIADNE_ASSERT_MSG(Dyadic(this->_v)==Dyadic(t),"Power-of-two "<<t<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

template<class F> Value<F>::Value(Integer const& z, PR pr)
    : _v(z,to_nearest,pr)
{
    Rational q(_v);
    ARIADNE_ASSERT_MSG(Dyadic(this->_v)==z,"Integer "<<z<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

template<class F> Value<F>::Value(Dyadic const& w, PR pr)
    : _v(w,pr)
{
    ARIADNE_ASSERT_MSG(Dyadic(this->_v)==w || is_nan(w),"Dyadic number "<<w<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

template<class F> Value<F>::Value(Value<F> const& x, PR pr)
    : _v(x._v,to_nearest,pr)
{
    ARIADNE_ASSERT_MSG(*this==x,"Exact FloatValue "<<x<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

template<class F> Value<F>::operator Dyadic() const {
    return this->_v.operator Dyadic();
}

template<class F> Value<F>::operator Rational() const {
    return Rational(this->operator Dyadic());
}

template<class F> Value<F>& Value<F>::operator=(TwoExp const& t) {
    _v=F(t,this->precision());
    return *this;
}

template<class F> Value<F>& Value<F>::operator=(Integer const& z) {
    _v=F(z,this->precision());
    ARIADNE_ASSERT_MSG(Dyadic(_v)==z,"Integer "<<z<<" cannot be assigned exactly to a floating-point number with precision "<<this->precision()<<"; nearest is "<<(*this));
    return *this;
}

template<class F> Value<F>& Value<F>::operator=(Dyadic const& w) {
    _v=F(w,this->precision());
    ARIADNE_ASSERT_MSG(Dyadic(_v)==w,"Dyadic number "<<w<<" cannot be assigned exactly to a floating-point number with precision "<<this->precision()<<"; nearest is "<<(*this));
    return *this;
}

template<class F> auto Value<F>::create(ValidatedNumber const& y) const -> Ball<F> {
    return Ball<F>(y,this->precision());
}

/*
template<class F> Value<F>::operator ExactNumber() const {
    return ExactNumber(new NumberWrapper<Value<F>>(*this));
}
*/







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
    DecimalPlaces plcs{dgts}
;
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


template<class F, class FE> InputStream& Operations<Ball<F,FE>>::_read(InputStream& is, Ball<F,FE>& x) {
    static const char pmstr[] = "\u00b1";
    char cpm[3];
    F _v; FE _e;
    auto rnd=F::get_rounding_mode();
    F::set_rounding_to_nearest();
    is >> _v;
    is >> cpm[0] >> cpm[1];
    F::set_rounding_mode(rnd);
    auto rnde=FE::get_rounding_mode();
    FE::set_rounding_upward();
    is >> _e;
    FE::set_rounding_mode(rnde);
    ARIADNE_ASSERT(not is.fail());
    ARIADNE_ASSERT(std::strcmp(cpm,pmstr));
    x._v=_v; x._e=_e;
    return is;
}


template<class F, class FE> OutputStream& Operations<Ball<F,FE>>::_write(OutputStream& os, Ball<F,FE> const& x) {
    return os << x.value() << "\u00b1" << x.error();
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

    DecimalPlaces plcs{dgts}
;
    // Get string version of mpfr values
    String vstr=print(v,plcs,MPFR_RNDN);
    String estr=print(e,plcs,MPFR_RNDU);

    // Find position of first significan digit of error
    auto vcstr=vstr.c_str(); auto ecstr=estr.c_str();
    size_t cpl=0;
    if(edbl==0.0) {
        cpl=std::strlen(vcstr);
    }
 else if(edbl<1.0) {
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






/*
Rational cast_exact(Real const& x) {
    return Rational(cast_exact(FloatApproximation<DoublePrecision>(x,dp)));
}
*/


template<class F> struct Operations<Error<F>> {
    static OutputStream& _write(OutputStream& os, Error<F> const& x) {
        return write(os,x.raw(),DecimalPrecision{Error<F>::output_places},upward);
    }
    static InputStream& _read(InputStream& is, Error<F>& x) {
        UpperBound<F> xu; is >> xu; x=Error<F>(xu); return is;
    }
};




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





/*
template<class F> Approximation<F> _make_float(Number<ApproximateTag> x) { return Approximation<F>(x); }
template<class F> LowerBound<F> _make_float(Number<ValidatedLowerTag> x) { return LowerBound<F>(x); }
template<class F> UpperBound<F> _make_float(Number<ValidatedUpperTag> x) { return UpperBound<F>(x); }
template<class F> Bounds<F> _make_float(Number<ValidatedTag> x) { return Bounds<F>(x); }
template<class F> Bounds<F> _make_float(Number<EffectiveTag> x) { return Bounds<F>(x); }
template<class F> Bounds<F> _make_float(Number<ExactTag> x) { return Bounds<F>(x); }
template<class F> Bounds<F> _make_float(Real r) { return Bounds<F>(r); }
template<class F> Bounds<F> _make_float(Rational q) { return Bounds<F>(q); }
template<class F> Value<F> _make_float(Integer z) { return Value<F>(z); }
*/


template class Approximation<FloatDP>;
template class LowerBound<FloatDP>;
template class UpperBound<FloatDP>;
template class Bounds<FloatDP>;
template class Ball<FloatDP>;
template class Value<FloatDP>;
template class Error<FloatDP>;

template class Approximation<FloatMP>;
template class LowerBound<FloatMP>;
template class UpperBound<FloatMP>;
template class Bounds<FloatMP>;
template class Ball<FloatMP>;
template class Value<FloatMP>;
template class Error<FloatMP>;

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
template<> String class_name<FloatDPError>() { return "FloatDPError"; }

template<> String class_name<FloatMPApproximation>() { return "FloatMPApproximation"; }
template<> String class_name<FloatMPLowerBound>() { return "FloatMPLowerBound"; }
template<> String class_name<FloatMPUpperBound>() { return "FloatMPUpperBound"; }
template<> String class_name<FloatMPBounds>() { return "FloatMPBounds"; }
template<> String class_name<FloatMPBall>() { return "FloatMPBall"; }
template<> String class_name<FloatMPValue>() { return "FloatMPValue"; }
template<> String class_name<FloatMPError>() { return "FloatMPError"; }

template<> String class_name<FloatMDPBall>() { return "FloatMDPBall"; }


}
 // namespace Ariadne
