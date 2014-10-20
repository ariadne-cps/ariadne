/***************************************************************************
 *            float64.cc
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

/*! \file float64.cc
 *  \brief
 */

#include "utility/module.h"

#include <limits>
#include <cmath>
#include <cfenv> // For correct rounding when printing

#include "utility/stdlib.h"

#include "float64.h"

#include "logical.h"
#include "integer.h"
#include "rational.h"

#include "rounding.h"

/******** Print rounding mode ************************************************/
#if defined ARIADNE_C99_ROUNDING
    #warning "Using standard fenv.h C header file for setting the rounding mode."
#elif defined ARIADNE_BOOST_ROUNDING
    #if defined BOOST_NUMERIC_INTERVAL_DETAIL_C99_ROUNDING_CONTROL_HPP
        #warning "Using Boost interval library standard fenv.h C header for setting the rounding mode."
    #else
        #warning "Using Boost interval library hardware rounding for setting the rounding mode."
    #endif
#elif defined ARIADNE_GCC_ROUNDING
    #warning "Using ordinary GCC inline assembler for setting the rounding mode."
#elif defined ARIADNE_EGCC_ROUNDING
    #warning "Using extended GCC inline assembler for setting the rounding mode."
#elif defined ARIADNE_SSE_ROUNDING
    #warning "Using SSE <xmmintrin.h> header file for setting the rounding mode."
#elif defined ARIADNE_MSVC_ROUNDING
    #warning "Using Microsoft Visual Studio inline assembler for setting the rounding mode."
#else
    #warning "No rounding mode defined."
#endif

namespace Ariadne {

/************ Initialisation *************************************************/

bool init() {
    set_default_rounding();
    std::cerr<<std::boolalpha<<std::fixed<<std::setprecision(6);
    std::cout<<std::boolalpha<<std::fixed<<std::setprecision(6);
    std::clog<<std::boolalpha<<std::fixed<<std::setprecision(6);
    return true;
}

static const bool initialized=init();

/************  Publicly-accessible rounding-mode changing *******************/

typedef unsigned short rounding_mode_t;

//! \ingroup NumericModule \brief Set the active rounding mode.
void set_rounding_mode(RoundingModeType rnd) { _set_rounding_mode(rnd); }
//! \ingroup NumericModule \brief Get the active rounding mode.
RoundingModeType get_rounding_mode() { return _get_rounding_mode(); }

void set_rounding_to_nearest() { _set_rounding_to_nearest(); }
void set_rounding_downward() { _set_rounding_downward(); }
void set_rounding_upward() { _set_rounding_upward(); }
void set_rounding_toward_zero() { _set_rounding_toward_zero(); }

void set_default_rounding() { _set_rounding_upward(); }

/************ Double-precision constants **********************************************************/

typedef long long int int64_t;

static const double _quarter_pi_up=0.78539816339744839;
static const double _half_pi_up=1.5707963267948968;
static const double _pi_up=3.1415926535897936;
static const double _pi_down=3.1415926535897931;

static const double _quarter_pi_approx=0.78539816339744828;
static const double _half_pi_approx=1.5707963267948966;
static const double _pi_approx=3.1415926535897931;
static const double _two_pi_approx=6.2831853071795862;
static const double _sqrt2_approx=0.70710678118654757;
static const double _log2_approx=0.6931471805599453094;

static const double pi_approx=3.1415926535897931;
static const double pi_up    =3.1415926535897936;
static const double pi_down  =3.1415926535897931;

static const BoundedFloat64 pi_ivl(pi_down,pi_up);

/************ Operations on raw builtin double type **********************************************************/

static inline char rounding_mode_char()
{
    if((get_rounding_mode() & 3072) == 0000) { return 'n'; }
    if((get_rounding_mode() & 3072) == 1024) { return 'd'; }
    if((get_rounding_mode() & 3072) == 2048) { return 'u'; }
    if((get_rounding_mode() & 3072) == 3072) { return 'z'; }
    return '?';
}


static inline double horner_rnd(int n, double x, const long long int* c)
{
    volatile double y=1./c[n];
    for(int i=n-1; i>=0; --i) {
        y=1.0/c[i]+x*y;
    }
    return y;
}

static inline double horner_opp(int n, double x, const long long int* c)
{
    volatile double y=-1./c[n];
    for(int i=n-1; i>=0; --i) {
        y=-1.0/c[i]+x*y;
    }
    return -y;
}

bool is_integer(double x) {
    double intpart;
    return std::modf(x, &intpart) == 0.0;
}

template<class Z, EnableIf<IsIntegral<Z>>> Z integer_cast(Flt64 x) {
    Z r=x.get_d();
    assert(r==x.get_d());
    return r;
}

template uint integer_cast(Flt64);
template int integer_cast(Flt64);

inline volatile double add_near(volatile double x, volatile double y) {
    _set_rounding_to_nearest(); volatile double r=x+y; _set_rounding_upward(); return r; }
inline volatile double sub_near(volatile double x, volatile double y) {
    _set_rounding_to_nearest(); volatile double r=x-y; _set_rounding_upward(); return r; }
inline volatile double mul_near(volatile double x, volatile double y) {
    _set_rounding_to_nearest(); volatile double r=x*y; _set_rounding_upward(); return r; }
inline volatile double div_near(volatile double x, volatile double y) {
    _set_rounding_to_nearest(); volatile double r=x/y; _set_rounding_upward(); return r; }

inline volatile double next_rnd(volatile double x) { x=x+1e-300; x=x-1e-300; return x; }
inline volatile double next_opp(volatile double x) { volatile double y=-x; y=y+1e-300; y=y-1e-300; return -y; }

inline volatile double add_rnd(volatile double x, volatile double y) { return x+y; }
inline volatile double sub_rnd(volatile double x, volatile double y) { return x-y; }
inline volatile double mul_rnd(volatile double x, volatile double y) { return x*y; }
inline volatile double div_rnd(volatile double x, volatile double y) { return x/y; }
inline volatile double add_opp(volatile double x, volatile double y) { volatile double mx=-x; volatile double mr=mx-y; return -mr; }
inline volatile double sub_opp(volatile double x, volatile double y) { volatile double mx=-x; volatile double mr=mx+y; return -mr; }
inline volatile double mul_opp(volatile double x, volatile double y) { volatile double mx=-x; volatile double mr=mx*y; return -mr; }
inline volatile double div_opp(volatile double x, volatile double y) { volatile double mx=-x; volatile double mr=mx/y; return -mr; }

inline volatile double med_rnd(volatile double x, volatile double y) { return 0.5*(x+y); }
inline volatile double rad_rnd(volatile double x, volatile double y) { return 0.5*(y-x); }


inline volatile double pos(volatile double x) { return +x; }
inline volatile double neg(volatile double x) { return -x; }
inline volatile double min(volatile double x1, volatile double x2) { return (x1<x2)?x1:x2; }
inline volatile double max(volatile double x1, volatile double x2) { return (x1>x2)?x1:x2; }
inline double abs(double x) { return std::fabs(x); }

// Rounded power
volatile double pow_rnd(volatile double x, uint m)
{
    if(m==0) { return 1.0; }
    if(x==0) { return 0.0; }
    //if(x>=0.0 || (m%2==0)) { double r,p; r=1.0; p=(x>=0)?x:-x; while(m) { if(m%2) { r=mul_rnd(r,p); } p=mul_rnd(p,p); m/=2; } return r; }
    //else { double r,p r=-1.0; p=x; while(m) { if(m%2) { r=mul_rnd(-r,p); } p=mul_rnd(-p,p); m/=2; } return r; }
    if(x>0.0 || (m%2==0)) { volatile double r,p; r=1.0; p=std::abs(x); while(m) { if(m%2) { r=r*p; } p=p*p; m/=2; } return r; }
    else { volatile double r,p; r=-1.0; p=x; while(m) { if(m%2) { r=(-r)*p; } p=(-p)*p; m/=2; } return r; }
}

// Rounded power
volatile double pow_rnd(volatile double x, int n)
{
    if(n>=0) { return pow_rnd(x,uint(n)); }
    ARIADNE_ASSERT(x!=0.0);
    if(x>0.0 || (n%2==-1)) { volatile double r=1.0/x; return pow_rnd(r,uint(-n)); }
    else { volatile double r=-1.0/x; return pow_rnd(r,uint(-n)); }
}


double sqrt_rnd(double x)
{
    // long int c[]={ 0, 6, -360, 15120, -604800, 23950080, -946218790, 37362124800 };
    ARIADNE_ASSERT_MSG(x>=0, " x = "<<x<<"\n");

    if(x==0.0) { return 0.0; }
    int n; volatile double y,a,b;
    y=frexp(x,&n);
    if(n%2) { y*=2; n-=1; }
    assert(y>=0.5 && y<=2.0);

    a=0.0; b=y;
    // Test for (a-b) small, rather than exactly zero
    while( std::abs(a-b)/2 > std::numeric_limits<double>::epsilon()) {
        a=b;
        b=(a+y/a)/2;
    }

    return ldexp(b,n/2);
}


double pos_sin_rnd_series(double x);
double neg_sin_rnd_series(double x);
double pos_cos_rnd_series(double x);
double neg_cos_rnd_series(double x);
double tan_rnd_series(double x);


double texp(double x) {
    double r=1.0; double t=1.0;
    for(uint i=1; i!=20; ++i) {
        t*=x; t/=i; r+=t;
    }
    return r;
}

// Correctly rounded exponential function
double exp_rnd(double x)
{
    static const long long int c[7]={ 1LL, 6LL, -360LL, 15120LL, -604800LL, 23950080LL, -946218790LL };

    // Set w=r(exp(r)+1)/(exp(r)-1).
    // Then w=2 + s/c1 + s^2/c2 + ...
    // where s=r^2, and
    // and exp(r) = y = (w+r)/(w-r)
    // The first six terms are sufficient to compute w
    // Note that we always have w>0 and y>0

    // Note that dy/dr = 2w/(w-e)^2 > 0
    // and dy/dw = -2r/(w-r)^2, so if r<0 we need to under-approximate w

    if(x==0.0) { return 1.0; }
    double n=std::floor(x/_log2_approx+0.5);
    volatile double r,s,t,w,y;
    volatile double log2;
    log2=(n>0.0) ? next_opp(_log2_approx) : next_rnd(_log2_approx);
    r=x+(-n)*log2;

    ARIADNE_ASSERT_MSG(r>=-0.4, " r = "<<r<<", x = "<<x);
    ARIADNE_ASSERT_MSG(r<=+0.4, " r = "<<r<<", x = "<<x);

    if(r<0) {
        // Compute w by standard Horner's rule gives correct rounding since w is monotone increasing in s
        s = r*r;
        w = 1+horner_rnd(6,s,c);
        //c[0] + s * (1./c[1] + s * (1./c[2] + s * (1./c[3] + s * (1./c[4] + s * (1./c[5] + s * (1./c[6]) ) ) ) ) );
    } else {
        // Compute -w by standard Horner's rule gives opposite rounding since w is monotone increasing in s
        s = (-r)*r; s=-s;
        w = horner_opp(6,s,c);
        w = -1-w; w=-w;
        //w = -c[0] + s * (-1./c[1] + s * (-1./c[2] + s * (-1./c[3] + s * (-1./c[4] + s * (-1./c[5] + s * (-1./c[6]) ) ) ) ) );
        //w=-w;
    }

    t=r-w;
    y=(w+r)/(-t);


    long int m=static_cast<long int>(n);
    int e;
    volatile double z=frexp(y,&e);
    z=ldexp(z,m+e);
    return z;
}




double log_rnd(double x) {
    static const long long int c[12]={ 1LL, 3LL, 5LL, 7LL, 9LL, 11LL, 13LL, 15LL, 17LL, 19LL, 21LL, 23LL };

    ARIADNE_ASSERT(x>0.0);
    // Write x=2^ny with 1/sqrt(2) <= y <= sqrt(2)
    // Write log(y)=log(1+z)-log(1-z) where z=(y-1)/(y+1) and y=(1+z)/(1-z),
    // Note that y is monotone increasing in z (and vice-versa)
    // The constraints on y give |z| < 0.172.
    // Use a Taylor expansion to compute log(1+z)-log(1-z), yielding
    // We have log(y) = 2(z + 1/3 z^3 + 1/5 z^5+...)
    // We obtain reasonable sufficient accuracy by taking terms up to z^19
    // Let s=z^2 and w=(1+s/3+s^2/5+...+s^9/19), so log(y)=2*z*w
    // Note that if z<0 (corresponding to y<1) then we need to use
    // opposite rounding to compute s and w.

    if(x==1.0) { return 0.0; }

    int n;
    volatile double y,z,s,t,w,ly;

    y=frexp(x,&n);
    if(y<_sqrt2_approx) { y*=2; --n; }

    if(y>=1.0) {
        t=-1-y;
        z=(y-1)/(-t);
        s=z*z;
        w=horner_rnd(10,s,c);
        ly=2*z*w;
    } else {
        t=1-y;
        z=(-t)/(1+y);
        s=(-z)*z;
        s=-s;
        w=horner_opp(10,s,c);
        ly=2*z*w;
    }

    volatile double log2rnd=next_rnd(_log2_approx);


    return log2rnd*n+ly;
}

double pi_rnd() {
    switch(get_rounding_mode()) {
        case ROUND_NEAR: return _pi_approx;
        case ROUND_DOWN: return _pi_down;
        case ROUND_UP: return _pi_up;
        default: return _pi_approx;
    }
}

double pi_opp() {
    switch(get_rounding_mode()) {
        case ROUND_NEAR: return _pi_approx;
        case ROUND_DOWN: return _pi_up;
        case ROUND_UP: return _pi_down;
        default: return _pi_approx;
    }
}

double sin_rnd(double x) {
    volatile double two_pi_rnd=2*pi_rnd();
    volatile double two_pi_opp=2*pi_opp();

    volatile double half_pi_rnd=two_pi_rnd/4;
    volatile double half_pi_opp=two_pi_opp/4;

    int q = (long int)(std::floor(x/_quarter_pi_approx)) % 8;
    if(q<-4) { q+=8; } if(q>=4) { q-=8; }
    volatile double n=-std::floor(x/_two_pi_approx+0.5);

    volatile double y,w,s;

    // Set to true if sin is decreasing so we want opposite rounding
    bool want_opposite=(q<-2||q>=2);
    // if n is negative then we need to switch rounding of two_pi
    volatile double two_pi_corr=((n>=0.0) ^ want_opposite) ? two_pi_rnd : two_pi_opp;

    // Scale onto interval from -pi to pi
    if(want_opposite) { y=-x+(-n)*(two_pi_corr); y=-y; } else { y=x+n*two_pi_corr; }

    assert(-_two_pi_approx<=y && y<=_two_pi_approx);


    switch(q) {
    case -4: { w = -y - 2*half_pi_opp; w=-w; s=neg_sin_rnd_series(w); break; }
    case -3: { w = -y - 1*half_pi_rnd; w=+w; s=neg_cos_rnd_series(w); break; }
    case -2: { w = +y + 1*half_pi_rnd; w=+w; s=neg_cos_rnd_series(w); break; }
    case -1: { w = +y + 0*half_pi_opp; w=-w; s=neg_sin_rnd_series(w); break; }
    case +0: { w = +y + 0*half_pi_opp; w=+w; s=pos_sin_rnd_series(w); break; }
    case +1: { w = +y - 1*half_pi_opp; w=-w; s=pos_cos_rnd_series(w); break; }
    case +2: { w = -y + 1*half_pi_rnd; w=-w; s=pos_cos_rnd_series(w); break; }
    case +3: { w = -y + 2*half_pi_rnd; w=+w; s=pos_sin_rnd_series(w); break; }
    default: { assert(false); }
    }

    return s;
}

double cos_rnd(double x) {
    const double pi_rnd=Ariadne::pi_rnd();
    const double pi_opp=Ariadne::pi_opp();
    const double half_pi_rnd=pi_rnd/2;

    if(x==0.0) { return 1.0; }

    // Compute a constant n such that n*pi<x<(n+1)*pi
    long int n_rnd=(long int)(std::floor(div_rnd(x,pi_opp)));
    long int n_opp=(long int)(std::floor(div_opp(x,pi_rnd)));

    if(n_rnd!=n_opp) {
        if(n_rnd>n_opp) {
            // Rounding upwards
            if(n_rnd%2==0) { return 1.0; }
            volatile double y1=sub_rnd(pi_rnd*n_rnd,x);
            volatile double y2=sub_rnd(x,pi_opp*n_rnd);
            assert(y1>=0 && y2>=0);
            volatile double w=max(y1,y2);
            return neg_cos_rnd_series(w);
        } else {
            // Rounding downwards
            if(n_rnd%2==0) { return -1.0; }
            volatile double y1=sub_opp(pi_opp*n_opp,x);
            volatile double y2=sub_opp(x,pi_rnd*n_opp);
            assert(y1>=0 && y2>=0);
            volatile double w=max(y1,y2);
            return pos_cos_rnd_series(w);
        }
    }


    // Set y=x-n*pi (with opposite rounding) if n is even
    // and to (n+1)*pi-x (with opposite rounding) if n is odd
    volatile double y;
    if(n_rnd%2==0) {
        y=sub_opp(x,mul_rnd(n_rnd,pi_rnd));
    } else {
        y=sub_opp(mul_opp(n_rnd+1,pi_opp),x);
    }

    int q = (long int)(std::floor(y/_quarter_pi_approx)) % 8;
    assert(q<=4);

    volatile double w,c;
    if(q==0) {
        w=y;
        c=pos_cos_rnd_series(w);
    } else if(q==1 || q==2) {
        w=sub_rnd(pi_rnd/2,y);
        if(w>=0.0) { c=pos_sin_rnd_series(w); }
        else { c=neg_sin_rnd_series(-w); }
    } else if(q==3 || q==4) {
        w=sub_opp(pi_opp,y);
        c=neg_cos_rnd_series(w);
    } else {
        assert(false);
    }

    return c;



    volatile double z=0.0;

    ARIADNE_ASSERT(-_two_pi_approx<=y && y<=_two_pi_approx);
    switch(q) {
    case -4: { w = +y + 2*half_pi_rnd; w=+w; w=max(w,z); c=neg_cos_rnd_series(w); break; }
    case -3: { w = +y + 1*half_pi_rnd; w=-w; w=max(w,z); c=neg_sin_rnd_series(w); break; }
    case -2: { w = +y + 1*half_pi_rnd; w=+w; w=max(w,z); c=pos_sin_rnd_series(w); break; }
    case -1: { w = +y + 0*half_pi_rnd; w=-w; w=max(w,z); c=pos_cos_rnd_series(w); break; }
    case +0: { w = -y + 0*half_pi_rnd; w=-w; w=max(w,z); c=pos_cos_rnd_series(w); break; }
    case +1: { w = -y + 1*half_pi_rnd; w=+w; w=max(w,z); c=pos_sin_rnd_series(w); break; }
    case +2: { w = -y + 1*half_pi_rnd; w=-w; w=max(w,z); c=neg_sin_rnd_series(w); break; }
    case +3: { w = -y + 2*half_pi_rnd; w=+w; w=max(w,z); c=neg_cos_rnd_series(w); break; }
    default: { assert(false); }
    }

    return c;
}


double pos_sin_rnd_series(double x) {
    ARIADNE_ASSERT(x>=0.0);
    ARIADNE_ASSERT(x<=0.7853981634);
    static const long long int c[9]={ 1LL, -6LL, 120LL, -5040LL, 362880LL, -39916800LL, 6227020800LL, -1307674368000LL, 355687428096000LL };

    // Compute sin(x) by Taylor series
    volatile double s,w,y;
    // TODO: Use Horner's scheme. Need a different algorithm for the case that x is negative
    s=x*x; w=horner_rnd(8,s,c); y=x*w;


    return y;
}

double neg_sin_rnd_series(double x) {
    ARIADNE_ASSERT(x>=0.0);
    ARIADNE_ASSERT(x<=0.7853981634);
    static const long long int c[9]={ 1LL, -6LL, 120LL, -5040LL, 362880LL, -39916800LL, 6227020800LL, -1307674368000LL, 355687428096000LL };

    // Compute sin(x) by Taylor series
    volatile double s,w,y;
    // TODO: Use Horner's scheme. Need a different algorithm for the case that x is negative
    s=(-x)*x; s=-s; w=horner_opp(8,s,c); y=x*(-w);


    return y;
}


double pos_cos_rnd_series(double x) {
    ARIADNE_ASSERT(x>=0.0);
    ARIADNE_ASSERT(x<=0.7853981634);

    static const long long int c[9]={ 1LL, -2LL, 24LL, -720LL, 40320LL, -3628800LL, 479001600LL, -87178291200LL, 20922789888000LL };

    // Compute cos(x) by Taylor series. Since cos(x) is decreasing in x^2,
    // we need to use opposite rounding for the computation of x^2

    volatile double s,y;
    s=(-x)*x;
    s=-s;
    y=0.0;
    y=horner_rnd(8,s,c);
    return y;
}

double neg_cos_rnd_series(double x) {
    ARIADNE_ASSERT(x>=0.0);
    ARIADNE_ASSERT(x<=0.7853981634);

    static const long long int c[9]={ 1LL, -2LL, 24LL, -720LL, 40320LL, -3628800LL, 479001600LL, -87178291200LL, 20922789888000LL };
    volatile double s,y;
    s=x*x;
    y=0.0;
    y=horner_opp(8,s,c);
    return -y;
}

double tan_rnd(double x) {

    volatile double y,q,r,s,t,u,v;

    double n=std::floor(x/_pi_approx+0.5);

    volatile double pi_corr=(n>=0) ? next_opp(_pi_approx) : next_rnd(_pi_approx);
    y=x-n*pi_corr;


    ARIADNE_ASSERT(y>=-_pi_up/2);
    ARIADNE_ASSERT(y<=+_pi_up/2);

    // Use the double-angle formula tan(2x) = tan(x)/(1-tan^2(x))
    // Note that the function y/(1-y^2) is monotone increasing for |y|<1
    // To get enough accuracy, we use the double angle formula twice,
    // ensuring that |x|<=pi/8.

    q=y/4;
    r=tan_rnd_series(q);
    if(y>=0) {
        u=r*r-1;
        s=(2*r)/(-u);
        v=s*s-1;
        t=(2*s)/(-v);
    } else {
        u=(-r)*r;
        u=1+u;
        s=(2*r)/u;
        v=(-s)*s;
        v=1+v;
        t=(2*s)/v;
    }
    return t;
}

double tan_rnd_series(double x) {
    // Need |x|<=pi/8
    ARIADNE_ASSERT(x>=-_pi_up/8);
    ARIADNE_ASSERT(x<=+_pi_up/8);

    // Numerators of Taylor coefficients
    static const int64_t cn[13]={
        1LL, 1LL, 2LL, 17LL, 62LL, 1382LL, 21844LL, 929569LL, 6404582LL, 443861162LL,
        18888466084LL, 113927491862LL, 58870668456604LL };

    // Denominators of Taylor coefficients
    static const int64_t cd[13]={
        1LL, 3LL, 15LL, 315LL, 2835LL, 155925LL, 6081075LL, 638512875LL, 10854718875LL, 1856156927625LL,
        194896477400625LL, 2900518163668125LL, 3698160658676859375LL };


    // To get enough accuracy, we need |x|<pi/8
    // since we can't store the rational coefficients exactly.
    // Note that we could in principle use approximations to the
    // coesfficients, but this goes against the "spirit" of the code.

    volatile double c,s,w,r;
    if(x>=0) {
        s=x*x;
        w=double(cn[12])/cd[12];
        for(int i=11; i>=0; --i) {
            c=double(cn[i])/cd[i];
            w=c+s*w;
        }
        r=x*w;
    } else {
        s=(-x)*x; s=-s;
        w=double(-cn[12])/cd[12];
        for(int i=12; i>=0; --i) {
            c=double(-cn[i])/cd[i];
            w=c+s*w;
        }
        r=x*(-w);
    }
    return r;
}

/************ Generic **********************************************************/

template<class X> inline X generic_pow(X p, uint m) {
    X r=static_cast<X>(1);
    while(m!=0) {
        if(m%2==1) { r=r*p; }
        p=sqr(p); m=m/2;
    }
    return r;
}

template<class X> inline X generic_pow(X p, int n) {
    return n>=0 ? generic_pow(p,uint(n)) : rec(generic_pow(p,uint(-n)));
}

/************ Float64 **********************************************************/

Boolean eq(Rational const& q1, Flt64 x2) { return eq(q1,Rational(x2)); }
Comparison cmp(Rational const& q1, Flt64 x2) { return cmp(q1,Rational(x2)); }


Rational::Rational(Flt64 x) : Rational(x.get_d()) { }

inline Boolean operator==(Flt64 x, Rational const& q) { return eq(q,x); }
inline Boolean operator!=(Flt64 x, Rational const& q) { return not eq(q,x); }
inline Boolean operator< (Flt64 x, Rational const& q) { return cmp(q,x)==Comparison::GREATER; }
inline Boolean operator> (Flt64 x, Rational const& q) { return cmp(q,x)==Comparison::LESS; }
inline Boolean operator<=(Flt64 x, Rational const& q) { return cmp(q,x)!=Comparison::LESS; }
inline Boolean operator>=(Flt64 x, Rational const& q) { return cmp(q,x)!=Comparison::GREATER; }

Flt64 Flt64::eps() { return Flt64(std::numeric_limits<double>::epsilon()); }
Flt64 next_up(Flt64 x) { return Flt64(next_rnd(x.get_d())); }
Flt64 next_down(Flt64 x) { return Flt64(next_opp(x.get_d())); }

Flt64 add_near(Flt64 x, Flt64 y) { return Flt64(add_near(x.get_d(),y.get_d())); }
Flt64 sub_near(Flt64 x, Flt64 y) { return Flt64(add_near(x.get_d(),y.get_d())); }
Flt64 mul_near(Flt64 x, Flt64 y) { return Flt64(add_near(x.get_d(),y.get_d())); }
Flt64 div_near(Flt64 x, Flt64 y) { return Flt64(add_near(x.get_d(),y.get_d())); }

Flt64 add_up(Flt64 x, Flt64 y) { return Flt64(add_rnd(x.get_d(),y.get_d())); }
Flt64 sub_up(Flt64 x, Flt64 y) { return Flt64(sub_rnd(x.get_d(),y.get_d())); }
Flt64 mul_up(Flt64 x, Flt64 y) { return Flt64(mul_rnd(x.get_d(),y.get_d())); }
Flt64 div_up(Flt64 x, Flt64 y) { return Flt64(div_rnd(x.get_d(),y.get_d())); }

Flt64 add_down(Flt64 x, Flt64 y) { return Flt64(add_opp(x.get_d(),y.get_d())); }
Flt64 sub_down(Flt64 x, Flt64 y) { return Flt64(sub_opp(x.get_d(),y.get_d())); }
Flt64 mul_down(Flt64 x, Flt64 y) { return Flt64(mul_opp(x.get_d(),y.get_d())); }
Flt64 div_down(Flt64 x, Flt64 y) { return Flt64(div_opp(x.get_d(),y.get_d())); }

Flt64 rad_up(Flt64 x, Flt64 y) { return Flt64(0.5*(y.get_d()-x.get_d())); }
Flt64 med_near(Flt64 x, Flt64 y) { set_rounding_to_nearest(); Flt64 r(0.5*(x.get_d()+y.get_d())); set_rounding_upward(); return r; }

Flt64 pow_up(Flt64 x, Nat m) { return Flt64(pow_rnd(x.get_d(),m)); }
Flt64 pow_down(Flt64 x, Nat m) { set_rounding_downward(); Flt64 r=pow_rnd(x.get_d(),m); set_rounding_upward(); return r; }
Flt64 pow_approx(Flt64 x, Nat m) { set_rounding_to_nearest(); Flt64 r=pow_rnd(x.get_d(),m); set_rounding_upward(); return r; }
Flt64 pow_up(Flt64 x, Int n) { return Flt64(pow_rnd(x.get_d(),n)); }
Flt64 pow_down(Flt64 x, Int n) { set_rounding_downward(); Flt64 r=pow_rnd(x.get_d(),n); set_rounding_upward(); return r; }
Flt64 pow_approx(Flt64 x, Int n) { set_rounding_to_nearest(); Flt64 r=pow_rnd(x.get_d(),n); set_rounding_upward(); return r; }

Flt64 sqrt_up(Flt64 x) { return Flt64(sqrt_rnd(x.get_d())); }
Flt64 sqrt_down(Flt64 x) { set_rounding_downward(); Flt64 r=sqrt_rnd(x.get_d()); set_rounding_upward(); return r; }
Flt64 sqrt_approx(Flt64 x) { set_rounding_to_nearest(); Flt64 r=std::sqrt(x.get_d()); set_rounding_upward(); return r; }

Flt64 exp_up(Flt64 x) { return Flt64(exp_rnd(x.get_d())); }
Flt64 exp_down(Flt64 x) { set_rounding_downward(); Flt64 r=exp_rnd(x.get_d()); set_rounding_upward(); return r; }
Flt64 exp_approx(Flt64 x) { set_rounding_to_nearest(); Flt64 r=std::exp(x.get_d()); set_rounding_upward(); return r; }
Flt64 log_up(Flt64 x) { return Flt64(log_rnd(x.get_d())); }
Flt64 log_down(Flt64 x) { set_rounding_downward(); Flt64 r=log_rnd(x.get_d()); set_rounding_upward(); return r; }
Flt64 log_approx(Flt64 x) { set_rounding_to_nearest(); Flt64 r=std::log(x.get_d()); set_rounding_upward(); return r; }

Flt64 sin_up(Flt64 x) { return Flt64(sin_rnd(x.get_d())); }
Flt64 sin_down(Flt64 x) { set_rounding_downward(); Flt64 r=sin_rnd(x.get_d()); set_rounding_upward(); return r; }
Flt64 sin_approx(Flt64 x) { set_rounding_to_nearest(); Flt64 r=std::sin(x.get_d()); set_rounding_upward(); return r; }
Flt64 cos_up(Flt64 x) { return Flt64(cos_rnd(x.get_d())); }
Flt64 cos_down(Flt64 x) { set_rounding_downward(); Flt64 r=cos_rnd(x.get_d()); set_rounding_upward(); return r; }
Flt64 cos_approx(Flt64 x) { set_rounding_to_nearest(); Flt64 r=std::cos(x.get_d()); set_rounding_upward(); return r; }

Flt64 add_rnd(Flt64 x, Flt64 y) { return Flt64(add_rnd(x.get_d(),y.get_d())); }
Flt64 sub_rnd(Flt64 x, Flt64 y) { return Flt64(sub_rnd(x.get_d(),y.get_d())); }
Flt64 mul_rnd(Flt64 x, Flt64 y) { return Flt64(mul_rnd(x.get_d(),y.get_d())); }
Flt64 div_rnd(Flt64 x, Flt64 y) { return Flt64(div_rnd(x.get_d(),y.get_d())); }
Flt64 add_opp(Flt64 x, Flt64 y) { return Flt64(add_opp(x.get_d(),y.get_d())); }
Flt64 sub_opp(Flt64 x, Flt64 y) { return Flt64(sub_opp(x.get_d(),y.get_d())); }
Flt64 mul_opp(Flt64 x, Flt64 y) { return Flt64(mul_opp(x.get_d(),y.get_d())); }
Flt64 div_opp(Flt64 x, Flt64 y) { return Flt64(div_opp(x.get_d(),y.get_d())); }

// The following does not work!
//Flt64 div_opp(Flt64 x, Flt64 y) { Flt64 mx=-x; Flt64 mr=mx/y; Flt64 r=-mr; return r; }

// The following is equivalent to the current version
//Flt64 div_opp(Flt64 x, Flt64 y) { volatile double mx=-x.d; volatile double mr=mx/y.d; return -mr; }


/************ ExactFloat64 *****************************************************/

const ExactFloat64 inf = ExactFloat64(1.0/0.0);

Float64<Exact>::Float64() : _v(0.0) { }

Float64<Exact>::Float64(Integer value) : _v(value.get_si()) { assert(_v==this->get_flt()); }

//Float64<Exact>::Float64(double value) : _v(value) { }

Float64<Exact>::Float64(Flt64 value) : _v(value) { }

Float64<Exact>::Float64(Int64 n) : _v(n.get_si()) { assert(_v==n.get_si()); }

Flt64 ExactFloat64::get_flt() const { return _v; }

double ExactFloat64::get_d() const { return _v; }

ExactFloat64 ExactFloat64::create_zero() const {
    return ExactFloat64(0.0);
}

ExactFloat64 operator+(ExactFloat64 x) {
    return ExactFloat64(x._v);
}

ExactFloat64 operator-(ExactFloat64 x) {
    return ExactFloat64(-x._v);
}

ExactFloat64 operator*(ExactFloat64 x, TwoExp y) {
    x._v *= y.get_d(); return x;
}

ExactFloat64 operator/(ExactFloat64 x, TwoExp y) {
    x._v /= y.get_d(); return x;
}

ExactFloat64& operator*=(ExactFloat64& x, TwoExp y) {
    x._v *= y.get_d(); return x;
}

ExactFloat64& operator/=(ExactFloat64& x, TwoExp y) {
    x._v /= y.get_d(); return x;
}

// The following are in the header since the return type can be chosen by a macro
// BoundFloat64 operator+(ExactFloat64 x, ExactFloat64 y);
// BoundFloat64 operator-(ExactFloat64 x, ExactFloat64 y);
// BoundFloat64 operator*(ExactFloat64 x, ExactFloat64 y);
// BoundFloat64 operator/(ExactFloat64 x, ExactFloat64 y);


MetrcFloat64 add(ExactFloat64 x, ExactFloat64 y, ValueErrorTag) {
    return MetrcFloat64(add_near(x._v,y._v),rad_rnd(add_rnd(x._v,y._v),add_opp(x._v,y._v))); }
MetrcFloat64 sub(ExactFloat64 x, ExactFloat64 y, ValueErrorTag) {
    return MetrcFloat64(sub_near(x._v,y._v),rad_rnd(sub_rnd(x._v,y._v),sub_opp(x._v,y._v))); }
MetrcFloat64 mul(ExactFloat64 x, ExactFloat64 y, ValueErrorTag) {
    return MetrcFloat64(mul_near(x._v,y._v),rad_rnd(mul_rnd(x._v,y._v),mul_opp(x._v,y._v))); }
MetrcFloat64 div(ExactFloat64 x, ExactFloat64 y, ValueErrorTag) {
    return MetrcFloat64(div_near(x._v,y._v),rad_rnd(div_rnd(x._v,y._v),div_opp(x._v,y._v))); }

BoundFloat64 add(ExactFloat64 x, ExactFloat64 y, LowerUpperTag) { return BoundFloat64(add_opp(x._v,y._v),add_rnd(x._v,y._v)); }
BoundFloat64 sub(ExactFloat64 x, ExactFloat64 y, LowerUpperTag) { return BoundFloat64(sub_opp(x._v,y._v),sub_rnd(x._v,y._v)); }
BoundFloat64 mul(ExactFloat64 x, ExactFloat64 y, LowerUpperTag) { return BoundFloat64(mul_opp(x._v,y._v),mul_rnd(x._v,y._v)); }
BoundFloat64 div(ExactFloat64 x, ExactFloat64 y, LowerUpperTag) { return BoundFloat64(div_opp(x._v,y._v),div_rnd(x._v,y._v)); }

ExactFloat64 pos(ExactFloat64 x) { return ExactFloat64(x._v); }
ExactFloat64 neg(ExactFloat64 x) { return ExactFloat64(-x._v); }
BoundFloat64 sqr(ExactFloat64 x) { return BoundFloat64(mul_opp(x._v,x._v),mul_rnd(x._v,x._v)); }
BoundFloat64 rec(ExactFloat64 x) { return BoundFloat64(div_opp(1.0,x._v),div_rnd(1.0,x._v)); }
BoundFloat64 pow(ExactFloat64 x, Nat m) { return pow(BoundFloat64(x),m); }
BoundFloat64 pow(ExactFloat64 x, Int n) { return pow(BoundFloat64(x),n); }
BoundFloat64 sqrt(ExactFloat64 x) { return sqrt(BoundFloat64(x)); }
BoundFloat64 exp(ExactFloat64 x) { return exp(BoundFloat64(x)); }
BoundFloat64 log(ExactFloat64 x) { return log(BoundFloat64(x)); }
BoundFloat64 sin(ExactFloat64 x) { return sin(BoundFloat64(x)); }
BoundFloat64 cos(ExactFloat64 x) { return cos(BoundFloat64(x)); }
BoundFloat64 tan(ExactFloat64 x) { return tan(BoundFloat64(x)); }
BoundFloat64 atan(ExactFloat64 x) { return atan(BoundFloat64(x)); }

ExactFloat64 max(ExactFloat64 x, ExactFloat64 y) { return ExactFloat64(std::max(x._v,y._v)); }
ExactFloat64 min(ExactFloat64 x, ExactFloat64 y) { return ExactFloat64(std::min(x._v,y._v)); }
ExactFloat64 abs(ExactFloat64 x) { return ExactFloat64(std::fabs(x._v)); }
ErrorFloat64 mag(ExactFloat64 x) { return ErrorFloat64(std::fabs(x._v)); }

Boolean operator==(ExactFloat64 x, ExactFloat64 y) { return x.get_d()==y.get_d(); }
Boolean operator!=(ExactFloat64 x, ExactFloat64 y) { return x.get_d()!=y.get_d(); }
Boolean operator<=(ExactFloat64 x, ExactFloat64 y) { return x.get_d()<=y.get_d(); }
Boolean operator>=(ExactFloat64 x, ExactFloat64 y) { return x.get_d()>=y.get_d(); }
Boolean operator< (ExactFloat64 x, ExactFloat64 y) { return x.get_d()< y.get_d(); }
Boolean operator> (ExactFloat64 x, ExactFloat64 y) { return x.get_d()> y.get_d(); }

inline Boolean eq(Rational const& q, ExactFloat64 x) { return eq(q,x.get_flt()); }
inline Comparison cmp(Rational const& q, ExactFloat64 x) { return cmp(q,x.get_flt()); }
Boolean operator==(Rational const& q1, ExactFloat64 x2) { return eq(q1,x2); }
Boolean operator!=(Rational const& q1, ExactFloat64 x2) { return !(eq(q1,x2)); }
Boolean operator< (Rational const& q1, ExactFloat64 x2) { return cmp(q1,x2)==Comparison::LESS; }
Boolean operator> (Rational const& q1, ExactFloat64 x2) { return cmp(q1,x2)==Comparison::GREATER; }
Boolean operator<=(Rational const& q1, ExactFloat64 x2) { return cmp(q1,x2)!=Comparison::GREATER; }
Boolean operator>=(Rational const& q1, ExactFloat64 x2) { return cmp(q1,x2)!=Comparison::LESS; }

Boolean operator==(ExactFloat64 x1, Rational const& q2) { return q2==x1; }
Boolean operator!=(ExactFloat64 x1, Rational const& q2) { return q2!=x1; }
Boolean operator< (ExactFloat64 x1, Rational const& q2) { return q2< x1; }
Boolean operator> (ExactFloat64 x1, Rational const& q2) { return q2> x1; }
Boolean operator<=(ExactFloat64 x1, Rational const& q2) { return q2<=x1; }
Boolean operator>=(ExactFloat64 x1, Rational const& q2) { return q2>=x1; }

Boolean operator==(ExactFloat64 x1, Int64 n2) { return x1.get_d() == n2.get_si(); }
Boolean operator!=(ExactFloat64 x1, Int64 n2) { return x1.get_d() != n2.get_si(); }
Boolean operator< (ExactFloat64 x1, Int64 n2) { return x1.get_d() <  n2.get_si(); }
Boolean operator> (ExactFloat64 x1, Int64 n2) { return x1.get_d() >  n2.get_si(); }
Boolean operator<=(ExactFloat64 x1, Int64 n2) { return x1.get_d() <= n2.get_si(); }
Boolean operator>=(ExactFloat64 x1, Int64 n2) { return x1.get_d() >= n2.get_si(); }

Bool same(ExactFloat64 x, ExactFloat64 y) { return x._v==y._v; }

OutputStream& operator<<(OutputStream& os, ExactFloat64 const& x) {
    rounding_mode_t rnd = get_rounding_mode();
    fesetround(FE_TONEAREST);
    os << std::showpoint << std::setprecision(ExactFloat64::output_precision) << x._v;
    set_rounding_mode(rnd);
    return os;
}

Void ExactFloat64::set_output_precision(int p) {
    output_precision=p;
}

int ExactFloat64::output_precision = 18;

Rational::Rational(const ExactFloat64& x) : Rational(x.get_d()) {
}


/************ BoundFloat64 *****************************************************/

Float64<Bound>::Float64() : _l(0.0), _u(0.0) {
}

Float64<Bound>::Float64(double x) : _l(x), _u(x) {
}

Float64<Bound>::Float64(double x, Exact) : _l(x), _u(x) {
}

Float64<Bound>::Float64(Flt64 value) : _l(value), _u(value) {
}

Float64<Bound>::Float64(const Integer& z) : BoundFloat64(Rational(z)) {
}

Float64<Bound>::Float64(const Rational& q) : _l(q.get_d()), _u(q.get_d()) {
    while(Rational(Flt64(_l))>q) { _l=sub_opp(_l,1e-300); }
    while(Rational(Flt64(_u))<q) { _u=add_rnd(_u,1e-300); }
}

Float64<Bound>::Float64(LowerFloat64 lower, UpperFloat64 upper) : _l(lower.get_d()), _u(upper.get_d()) {
}

Float64<Bound>::Float64(double lower, double upper) : _l(lower), _u(upper) {
}

BoundFloat64 BoundFloat64::create_zero() const {
    return BoundFloat64(0.0);
}

BoundFloat64 operator+(BoundFloat64 x) {
    return BoundFloat64(x._l,x._u);
}

BoundFloat64 operator-(BoundFloat64 x) {
    return BoundFloat64(-x._u,-x._l);
}

BoundFloat64 operator+(BoundFloat64 x,BoundFloat64 y) {
    return BoundFloat64(add_opp(x._l,y._l),add_rnd(x._u,y._u));
}

BoundFloat64 operator-(BoundFloat64 x,BoundFloat64 y) {
    return BoundFloat64(sub_opp(x._l,y._u),sub_rnd(x._u,y._l));
}

BoundFloat64 operator*(BoundFloat64 x,BoundFloat64 y) {
    if(x._l>=0) {
        if(y._l>=0) {
            return BoundFloat64{mul_opp(x._l,y._l),mul_rnd(x._u,y._u)};
        } else if(y._u<=0) {
            return BoundFloat64{mul_opp(x._u,y._l),mul_rnd(x._l,y._u)};
        } else {
            return BoundFloat64{mul_opp(x._u,y._l),mul_rnd(x._u,y._u)};
        }
    }
    else if(x._u<=0) {
        if(y._l>=0) {
            return BoundFloat64{mul_opp(x._l,y._u),mul_rnd(x._u,y._l)};
        } else if(y._u<=0) {
            return BoundFloat64{mul_opp(x._u,y._u),mul_rnd(x._l,y._l)};
        } else {
            return BoundFloat64{mul_opp(x._l,y._u),mul_rnd(x._l,y._l)};
        }
    } else {
        if(y._l>=0) {
            return BoundFloat64{mul_opp(x._l,y._u),mul_rnd(x._u,y._u)};
        } else if(y._u<=0) {
            return BoundFloat64{mul_opp(x._u,y._l),mul_rnd(x._l,y._l)};
        } else {
            using std::min; using std::max;
            return BoundFloat64{min(mul_opp(x._u,y._l),mul_rnd(x._l,y._u)),
                              max(mul_rnd(x._l,y._l),mul_rnd(x._u,y._u))};
        }
    }
    return BoundFloat64();
}

#ifdef ARIADNE_UNDEF
BoundFloat64 operator*(BoundFloat64 x, BoundFloat64 y) {
    BoundFloat64 r;
    if(x._l>=0) {
        if(y._l>=0) {
            r._l=mul_opp(x._l,y._l); r._u=mul_rnd(x._u,y._u);
        } else if(y._u<=0) {
            r._l=mul_opp(x._u,y._l); r._u=mul_rnd(x._l,y._u);
        } else {
            r._l=mul_opp(x._u,y._l); r._u=mul_rnd(x._u,y._u);
        }
    }
    else if(x._u<=0) {
        if(y._l>=0) {
            r._l=mul_opp(x._l,y._u); r._u=mul_rnd(x._u,y._l);
        } else if(y._u<=0) {
            r._l=mul_opp(x._u,y._u); r._u=mul_rnd(x._l,y._l);
        } else {
            r._l=mul_opp(x._l,y._u); r._u=mul_rnd(x._l,y._l);
        }
    } else {
        if(y._l>=0) {
            r._l=mul_opp(x._l,y._u); r._u=mul_rnd(x._u,y._u);
        } else if(y._u<=0) {
            r._l=mul_opp(x._u,y._l); r._u=mul_rnd(x._l,y._l);
        } else {
            r._l=min(rnd_opp(x._u,y._l),mul_rnd(x._l,y._u)),
            r._u=max(mul_rnd(x._l,y._l),mul_rnd(x._u,y._u));
        }
    }
    return BoundFloat64(r._l,r._u);
}
#endif

BoundFloat64 operator/(BoundFloat64 x,BoundFloat64 y) {
    if(y._l>0) {
        if(x._l>0) {
            return BoundFloat64{div_opp(x._l,y._u),div_rnd(x._u,y._l)};
        } else if(x._u<0) {
            return BoundFloat64{div_opp(x._l,y._l),div_rnd(x._u,y._u)};
        } else {
            return BoundFloat64{div_opp(x._l,y._l),div_rnd(x._u,y._l)};
        }
    } else if(y._u<0) {
        if(x._l>0) {
            return BoundFloat64{div_opp(x._u,y._u),div_rnd(x._l,y._l)};
        } else if(x._u<0) {
            return BoundFloat64{div_opp(x._u,y._l),div_rnd(x._l,y._u)};
        } else {
            return BoundFloat64{div_opp(x._u,y._u),div_rnd(x._l,y._u)};
        }
    } else {
        ARIADNE_THROW(DivideByZeroError,"operator/(BoundFloat64 x, BoundedFloat64 y)",x <<"/"<<y);
        return BoundFloat64();
    }
}

BoundFloat64& operator+=(BoundFloat64& x, BoundFloat64 y) {
    x._u=add_rnd(x._u,y._u);
    x._l=add_opp(x._l,y._l);
    return x;
}

BoundFloat64& operator-=(BoundFloat64& x, BoundFloat64 y) {
    x._u=sub_rnd(x._u,y._l);
    x._l=sub_opp(x._l,y._u);
    return x;
}

BoundFloat64& operator*=(BoundFloat64& x, BoundFloat64 y) {
    return x=x*y;
}

BoundFloat64& operator/=(BoundFloat64& x, BoundFloat64 y) {
    return x=x/y;
}

BoundFloat64 add(BoundFloat64 x, BoundFloat64 y) {
    return x+y;
}

BoundFloat64 sub(BoundFloat64 x, BoundFloat64 y) {
    return x-y;
}

BoundFloat64 mul(BoundFloat64 x, BoundFloat64 y) {
    return x*y;
}

BoundFloat64 div(BoundFloat64 x, BoundFloat64 y) {
    return x/y;
}

BoundFloat64 pos(BoundFloat64 x) {
    return BoundFloat64(x._l,x._u);
}

BoundFloat64 neg(BoundFloat64 x) {
    return BoundFloat64(-x._u,-x._l);
}

BoundFloat64 rec(BoundFloat64 x) {
    volatile double one=1.0;
    if(x._l>0 or x._u<0) {
        return BoundFloat64{div_opp(one,x._u),div_rnd(one,x._l)};
    } else {
        return BoundFloat64{-inf,+inf};
    }
}

BoundFloat64 sqr(BoundFloat64 x) {
    if(x._l>0) {
        return BoundFloat64{mul_opp(x._l,x._l),mul_rnd(x._u,x._u)};
    } else if(x._u<0) {
        return BoundFloat64{mul_opp(x._u,x._u),mul_rnd(x._l,x._l)};
    } else {
        return BoundFloat64{0.0,std::max(mul_rnd(x._l,x._l),mul_rnd(x._u,x._u))};
    }
}

BoundFloat64 min(BoundFloat64 x, BoundFloat64 y) {
    using std::min;
    return BoundFloat64(min(x._l,y._l),min(x._u,y._u));
}

BoundFloat64 max(BoundFloat64 x, BoundFloat64 y) {
    return BoundFloat64(max(x._l,y._l),max(x._u,y._u));
}

BoundFloat64 abs(BoundFloat64 x) {
    return BoundFloat64(max(max(x._l,-x._u),0.0),max(-x._l,x._u));
}

LowerFloat64 mig(BoundFloat64 x) {
    return LowerFloat64(max(max(x._l,-x._u),0.0));
}

ErrorFloat64 mag(BoundFloat64 x) {
    return ErrorFloat64(max(-x._l,x._u));
}

BoundFloat64 pow(BoundFloat64 x, Nat m) {
    if(m%2==0) { x=abs(x); }
    set_rounding_downward();
    volatile double r_l=pow_rnd(x._l,m);
    set_rounding_upward();
    volatile double r_u=pow_rnd(x._u,m);
    return BoundedFloat64(r_l,r_u);
}

BoundFloat64 pow(BoundFloat64 x, Int n) {
//    if(n<0) { return pow(rec(x),Nat(-n)); }
    if(n<0) { return rec(pow(x,Nat(-n))); }
    else return pow(x,Nat(n));
}

BoundFloat64 sqrt(BoundFloat64 x) {
    set_rounding_downward();
    volatile double r_l=std::sqrt(x._l);
    set_rounding_upward();
    volatile double r_u=std::sqrt(x._u);
    //volatile double tl=-x._l; tl=exp(x._l); tl=-tl; tl=1/tl; tl=-tl;
    return BoundFloat64(r_l,r_u);
}

BoundFloat64 exp(BoundFloat64 x) {
    set_rounding_downward();
    volatile double r_l=exp_rnd(x._l);
    set_rounding_upward();
    volatile double r_u=exp_rnd(x._u);
    return BoundedFloat64(r_l,r_u);
}

BoundFloat64 log(BoundFloat64 x) {
    set_rounding_downward();
    volatile double r_l=log_rnd(x._l);
    set_rounding_upward();
    volatile double r_u=log_rnd(x._u);
    return BoundedFloat64(r_l,r_u);
}

BoundFloat64 sin(BoundFloat64 x) {
    return cos(x-pi_ivl/ExactFloat64(2.0));
}

BoundFloat64 cos(BoundFloat64 x) {

    ARIADNE_ASSERT(x._l<=x._u);
    if(x.error().get_flt()>2*pi_down) {
        return BoundedFloat64(-1.0,+1.0);
    }

    Flt64 n=std::floor(x._l/(2*pi_approx)+0.5);
    x=x-ExactFloat64(2*n)*pi_ivl;
    ARIADNE_ASSERT(x._l<=x._u);

    ARIADNE_ASSERT(x._l<=pi_up);
    ARIADNE_ASSERT(x._u>=-pi_up);

    volatile double r_l,r_u;
    if(x._l<=-pi_down) {
        if(x._u<=0.0) { r_l=-1.0; r_u=cos_rnd(x._u); }
        else { r_l=-1.0; r_u=+1.0; }
    } else if(x._l<=0.0) {
        if(x._u<=0.0) { set_rounding_downward(); r_l=cos_rnd(x._l); set_rounding_upward(); r_u=cos_rnd(x._u); }
        else if(x._u<=pi_down) { set_rounding_downward(); r_l=cos_rnd(max(-x._l,x._u)); set_rounding_upward(); r_u=+1.0; }
        else { r_l=-1.0; r_u=+1.0; }
    } else if(x._l<=pi_up) {
        if(x._u<=pi_down) { set_rounding_downward(); r_l=cos_rnd(x._u); set_rounding_upward(); r_u=cos_rnd(x._l); }
        else if(x._u<=2*pi_down) { r_l=-1.0; r_u=cos_rnd(min(x._l,sub_opp(2*pi_down,x._u))); }
        else { r_l=-1.0; r_u=+1.0; }
    } else {
        assert(false);
    }
    return BoundedFloat64(r_l,r_u);
}

BoundFloat64 tan(BoundFloat64 x) {
    return mul(sin(x),rec(cos(x)));
}

BoundFloat64 atan(BoundFloat64 x) {
    ARIADNE_NOT_IMPLEMENTED;
}

Tribool operator==(BoundFloat64 x, BoundFloat64 y) {
    if(x._l>y._u || x._u<y._l) { return false; }
    else if(x._l==y._u && x._u==y._l) { return true; }
    else if(x._l+x._u && y._l+y._u) { return Tribool(LogicalValue::LIKELY); }
    else { return indeterminate; }
}

Tribool operator!=(BoundFloat64 x, BoundFloat64 y) {
    return !(x==y);
}

Tribool operator< (BoundFloat64 x, BoundFloat64 y) {
    if(x._u<y._l) { return true; }
    else if (x._l>=y._u) { return false; }
    else { return indeterminate; }
}

Tribool operator> (BoundFloat64 x, BoundFloat64 y) {
    return y<x;
}

Tribool operator<=(BoundFloat64 x, BoundFloat64 y) {
    return !(y<x);
}

Tribool operator>=(BoundFloat64 x, BoundFloat64 y) {
    return !(x<y);
}

Tribool operator==(BoundFloat64 x, Rational const& q) {
    if (x._l==x._u) { return x._l==q; }
    else if (x._l> q || x._l<q) { return false; }
    else { return indeterminate;}
}

Tribool operator!=(BoundFloat64 x, Rational const& q) {
    if (x._l==x._u) { return x._l!=q; }
    else if (x._l>q || x._u<q) { return true; }
    else { return indeterminate; }
}

Tribool operator< (BoundFloat64 x, Rational const& q) {
    if (x._u<q) { return true; }
    else if(x._l>=q) { return false; }
    else { return indeterminate; }
}

Tribool operator> (BoundFloat64 x, Rational const& q) {
    if (x._l>q) { return true; }
    else if(x._u<=q) { return false; }
    else { return indeterminate; }
}

Tribool operator<=(BoundFloat64 x, Rational const& q) {
    if (x._u<=q) { return true; }
    else if(x._l>q) { return false; }
    else { return indeterminate; }
}

Tribool operator>=(BoundFloat64 x, Rational const& q) {
    if (x._l>=q) { return true; }
    else if(x._u<q) { return false; }
    else { return indeterminate; }
}

Bool same(BoundFloat64 x, BoundFloat64 y) {
    return (x._l==y._l) && (x._u==y._u);
}

OutputStream& operator<<(OutputStream& os, BoundFloat64 const& x) {
    rounding_mode_t rnd=get_rounding_mode();
    os << '{';
    fesetround(FE_DOWNWARD);
    os << std::showpoint << std::setprecision(BoundFloat64::output_precision) << x.lower().get_d();
    os << ':';
    fesetround(FE_UPWARD);
    os << std::showpoint << std::setprecision(BoundFloat64::output_precision) << x.upper().get_d();
    set_rounding_mode(rnd);
    os << '}';
    return os;
}

LowerFloat64 BoundFloat64::lower() const {
    return LowerFloat64{_l};
}

UpperFloat64 BoundFloat64::upper() const {
    return UpperFloat64{_u};
}

ExactFloat64 BoundFloat64::value() const {
    set_rounding_to_nearest(); volatile double c=(_l+_u)/2; set_rounding_upward(); return ExactFloat64(c);
}

ErrorFloat64 BoundFloat64::error() const {
    set_rounding_to_nearest(); volatile double c=(_l+_u)/2; set_rounding_upward(); return ErrorFloat64(std::max(_u-c,c-_l));
}

ErrorFloat64 BoundFloat64::width() const {
    return ErrorFloat64(_u-_l);
}

ErrorFloat64 BoundFloat64::radius() const {
    return ErrorFloat64{sub_rnd(_u,_l)/2};
}

BoundFloat64::operator MetrcFloat64() const {
    BoundFloat64 const& x=*this;
    double v=add_near(x._l,x._u)/2;
    double e=(std::max(sub_rnd(x._u,v),sub_rnd(v,x._l)));
    assert(e>=0.0);
    return MetrcFloat64(v,e);
}
/*
Rational::operator BoundFloat64() const {
    Rational const& q=*this;
    volatile double u=q.get_d();
    volatile double nl=-u;
    while(Rational(-nl)>q) { nl+=std::numeric_limits<double>::minimum(); }
    while(Rational(u)<q) { u+=std::numeric_limits<double>::minimum(); }
    return BoundFloat64(-nl,u);
}
*/
Void BoundFloat64::set_output_precision(int p) {
    output_precision=p;
}

int BoundFloat64::get_output_precision() {
    return output_precision;
}

int BoundFloat64::output_precision = 6;


/************ Midpoint/RadiusFloat64 *******************************************/

Float64<Error>::Float64() : _e(0.0) { }
Float64<Error>::Float64(double eb) : _e(eb) {
    if(eb<0) { std::cerr<<"ERROR: eb="<<std::scientific<<std::setprecision(1)<<eb<<"\n"; } assert(eb>=-0.0); }
Float64<Error>::Float64(Flt64 eb) : _e(eb) { }
Float64<Error>::Float64(UpperFloat64 ub) : _e(ub.get_flt()) {
    if(ub.get_flt()<0.0) { std::cerr<<"ERROR: ub="<<ub<<"\n"; } assert(ub.get_flt()>=-0.0); }
ErrorFloat64::operator UpperFloat64() const { return UpperFloat64(_e); }
Flt64 ErrorFloat64::get_flt() const { return _e; }
double ErrorFloat64::get_d() const { return _e; }

UpperFloat64 operator+(ErrorFloat64 x) { return UpperFloat64(+x._e); }
LowerFloat64 operator-(ErrorFloat64 x) { return LowerFloat64(-x._e); }
ErrorFloat64 operator+(ErrorFloat64 x, ErrorFloat64 y) { return ErrorFloat64(add_rnd(x._e,y._e)); }
ErrorFloat64 operator*(ErrorFloat64 x, ErrorFloat64 y) {  return ErrorFloat64(mul_rnd(x._e,y._e)); }
ErrorFloat64 operator*(ErrorFloat64 x, UpperFloat64 y) {  assert(y>0.0); return ErrorFloat64(mul_rnd(x._e,y.get_d())); }
ErrorFloat64 operator/(ErrorFloat64 x, LowerFloat64 y) {  assert(y>0.0); return ErrorFloat64(div_rnd(x._e,y.get_d())); }
ErrorFloat64& operator+=(ErrorFloat64& x, ExactFloat64 y) { ARIADNE_ASSERT(y>=0.0_x); x._e=x._e+y.get_d(); return x; }
ErrorFloat64& operator+=(ErrorFloat64& x, ErrorFloat64 y) { return x=x+y; }
ErrorFloat64& operator*=(ErrorFloat64& x, ErrorFloat64 y) {  return x=x*y; }
ErrorFloat64 add(ErrorFloat64 x, ErrorFloat64 y) {  return ErrorFloat64(add_rnd(x._e,y._e)); }
ErrorFloat64 mul(ErrorFloat64 x, ErrorFloat64 y) {  return ErrorFloat64(mul_rnd(x._e,y._e)); }
ErrorFloat64 sqr(ErrorFloat64 x) {  return x*x; }
ErrorFloat64 pow(ErrorFloat64 x, Nat m) {  return generic_pow(x,m); }
ErrorFloat64 max(ErrorFloat64 x, ErrorFloat64 y) {  return ErrorFloat64(std::max(x._e,y._e)); }

Bool same(ErrorFloat64 x, ErrorFloat64 y) {  return x._e==y._e; }
Fuzzy operator==(ErrorFloat64 x, Flt64 y) {  return x._e==y; }

int ErrorFloat64::output_precision = 3;
Void ErrorFloat64::set_output_precision(int p) { output_precision=p; }

OutputStream& operator<<(OutputStream& os, ErrorFloat64 const& x) {
    // FIXME: Should use os << std::defaultfloat in C++11, not currently in gcc
    rounding_mode_t rnd = get_rounding_mode();
    fesetround(FE_UPWARD);
    os.unsetf(std::ios_base::floatfield);
    os << "\u00b1" << std::setprecision(ErrorFloat64::output_precision) << x._e << std::fixed;
    set_rounding_mode(rnd);
    return os;
}

/************ MetrcFloat64 ***********************************************/

Float64<Metrc>::Float64() : _v(0.0), _e(0.0) { }
Float64<Metrc>::Float64(double x) : _v(x), _e(0.0) { }
Float64<Metrc>::Float64(double x, ExactTag) : _v(x), _e(0.0) { }
Float64<Metrc>::Float64(double x, double r) : _v(x), _e(r) {
    if(!(r>=0.0)) { std::cerr<<"Error: MetrcFloat64(x,r) with x="<<x<<", r="<<r<<"\n"; } assert(r>=0.0); }
Float64<Metrc>::Float64(ExactFloat64 value, ErrorFloat64 error)
    : _v(value.get_d()), _e(error.get_d()) { assert(_e>=0.0); }
Float64<Metrc>::Float64(Flt64 value) : _v(value), _e(0.0) { }
Float64<Metrc>::Float64(Integer const& z) : MetrcFloat64(BoundFloat64(z)) { }
Float64<Metrc>::Float64(Rational const& q) : MetrcFloat64(BoundFloat64(q)) { }
MetrcFloat64 ExactFloat64::pm(ErrorFloat64 _e) const { return MetrcFloat64(*this,_e); }
ExactFloat64 MetrcFloat64::value() const { return ExactFloat64(_v); }
ErrorFloat64 MetrcFloat64::error() const { return ErrorFloat64(_e); }
UpperFloat64 MetrcFloat64::upper() const { return UpperFloat64(_v+_e); }
LowerFloat64 MetrcFloat64::lower() const { volatile double nl=_e-_v; return LowerFloat64(-nl); }
MetrcFloat64 MetrcFloat64::create_zero() const { return ValidatedFloat64(0.0); }

OutputStream& operator<<(OutputStream& os, MetrcFloat64 const& x) {
    //return os << x.value() << x.error();
    rounding_mode_t rnd=get_rounding_mode();
    fesetround(FE_TONEAREST);
    os << std::showpoint << std::setprecision(MetrcFloat64::output_precision) << x.value().get_d();
    fesetround(FE_UPWARD);
    os.unsetf(std::ios_base::floatfield);
    os << "\u00b1" << std::setprecision(3) << x.error().get_d() << std::fixed;
    set_rounding_mode(rnd);
    return os;
}

Void MetrcFloat64::set_output_precision(int p) {
    output_precision=p;
}

int MetrcFloat64::output_precision = 8;

MetrcFloat64 operator+(MetrcFloat64 x) {
    return MetrcFloat64(+x._v,x._e);
}

MetrcFloat64 operator-(MetrcFloat64 x) {
    return MetrcFloat64(-x._v,x._e);
}

MetrcFloat64 operator+(MetrcFloat64 x, MetrcFloat64 y) {
    set_rounding_to_nearest();
    volatile double rv=x._v+y._v;
    set_rounding_upward();
    volatile double ru=x._v+y._v;
    volatile double nx1v=-x._v;
    volatile double nrl=nx1v-y._v;
    volatile double re=(ru+nrl)/2+(x._e+y._e);
    return MetrcFloat64(rv,re);
}

MetrcFloat64 operator-(MetrcFloat64 x, MetrcFloat64 y) {
    set_rounding_to_nearest();
    volatile double rv=x._v-y._v;
    set_rounding_upward();
    volatile double ru=x._v-y._v;
    volatile double nx1v=-x._v;
    volatile double nrl=nx1v+y._v;
    volatile double re=(ru+nrl)/2+(x._e+y._e);
    return MetrcFloat64(rv,re);
}

MetrcFloat64 operator*(MetrcFloat64 x, MetrcFloat64 y) {
    volatile double x1v=x._v;
    volatile double x2v=y._v;
    volatile double x1e=x._e;
    volatile double x2e=y._e;
    set_rounding_to_nearest();
    volatile double rv=x1v*x2v;
    set_rounding_upward();
    volatile double rvu=x1v*x2v;
    volatile double mrvl=(-x1v)*x2v;
    volatile double re=(rvu+mrvl)/2+x1e*x2e+abs(x1v)*x2e+x1e*abs(x2v);
    return MetrcFloat64(rv,re);
}

MetrcFloat64 operator/(MetrcFloat64 x, MetrcFloat64 y) {
    return x*rec(y);
}

MetrcFloat64& operator+=(MetrcFloat64& x, MetrcFloat64 y) {
    return x=x+y;
}

MetrcFloat64& operator-=(MetrcFloat64& x, MetrcFloat64 y) {
    return x=x-y;
}

MetrcFloat64& operator*=(MetrcFloat64& x, MetrcFloat64 y) {
    return x=x*y;
}

MetrcFloat64& operator/=(MetrcFloat64& x, MetrcFloat64 y) {
    return x=x/y;
}

MetrcFloat64 add(MetrcFloat64 x, MetrcFloat64 y) {
    return x+y;
}

MetrcFloat64 sub(MetrcFloat64 x, MetrcFloat64 y) {
    return x-y;
}

MetrcFloat64 mul(MetrcFloat64 x, MetrcFloat64 y) {
    return x*y;
}

MetrcFloat64 div(MetrcFloat64 x, MetrcFloat64 y) {
    return x/y;
}

MetrcFloat64 pow(MetrcFloat64 x, Nat m) {
    return generic_pow(x,m);
}

MetrcFloat64 pow(MetrcFloat64 x, Int n) {
    return generic_pow(x,n);
}

MetrcFloat64 pos(MetrcFloat64 x) {
    return MetrcFloat64(+x._v,x._e);
}

MetrcFloat64 neg(MetrcFloat64 x) {
    return MetrcFloat64(-x._v,x._e);
}

MetrcFloat64 rec(MetrcFloat64 x) {
    volatile double xv=x._v;
    volatile double xe=x._e;
    volatile double mrxu=(-1)/(xv+xe);
    volatile double rxl=1/(xv-xe);
    volatile double re=(rxl+mrxu)/2;
    set_rounding_to_nearest();
    volatile double rv=(rxl-mrxu)/2;
    set_rounding_upward();
    return MetrcFloat64(rv,re);
}

MetrcFloat64 sqr(MetrcFloat64 x) {
    MetrcFloat64 r=x*x;
    if(r._e>r._v) {
        r._e=(r._e+r._v)/2;
        r._v=r._e;
    }
    return r;
}

MetrcFloat64 sqrt(MetrcFloat64 x) {
    return MetrcFloat64(sqrt(BoundFloat64(x)));
}

MetrcFloat64 exp(MetrcFloat64 x) {
    return MetrcFloat64(exp(BoundFloat64(x)));
}

MetrcFloat64 log(MetrcFloat64 x) {
    return MetrcFloat64(log(BoundFloat64(x)));
}

MetrcFloat64 sin(MetrcFloat64 x) {
    return MetrcFloat64(sin(BoundFloat64(x)));
}

MetrcFloat64 cos(MetrcFloat64 x) {
    return MetrcFloat64(cos(BoundFloat64(x)));
}

MetrcFloat64 tan(MetrcFloat64 x) {
    return MetrcFloat64(tan(BoundFloat64(x)));
}

MetrcFloat64 atan(MetrcFloat64 x) {
    return MetrcFloat64(atan(BoundFloat64(x)));
}


MetrcFloat64 abs(MetrcFloat64 x) {
    if(x._e<abs(x._v)) { return x; }
    else { double rv=(abs(x._v)+x._e)/2; return MetrcFloat64(rv,rv); }
}

MetrcFloat64 max(MetrcFloat64 x1, MetrcFloat64 x2) {
    return ((x1+x2)+abs(x1-x2))/2;
}

MetrcFloat64 min(MetrcFloat64 x1, MetrcFloat64 x2) {
    return ((x1+x2)-abs(x1-x2))/2;
}


ErrorFloat64 mag(MetrcFloat64 x) {
    return ErrorFloat64(abs(x._v)+x._e);
}

Bool same(MetrcFloat64 x, MetrcFloat64 y) {
    return (x._v==y._v) && (x._e==y._e);
}

Tribool operator==(MetrcFloat64 x, MetrcFloat64 y) {
    if(x._e==0.0 && y._e==0.0) { return x._v==y._v; }
    if(x._v==y._v) { return Tribool(LogicalValue::LIKELY); }

    volatile double d=x._v>=y._v ? x._v-y._v : y._v-x._v;
    volatile double e=x._e+y._e;
    if(d>e) { return false; }
    else { return indeterminate; }
}

Tribool operator!=(MetrcFloat64 x, MetrcFloat64 y) {
    return !(x==y);
}

Tribool operator< (MetrcFloat64 x, MetrcFloat64 y) {
    volatile double e=x._e+x._e;
    volatile double d=x._v-y._v;
    if(d+e<0) { return true; }
    d=y._v-x._v;
    if(d+e<=0) { return false; }
    return indeterminate;
}

Tribool operator> (MetrcFloat64 x, MetrcFloat64 y) {
    return y<x;
}

Tribool operator<=(MetrcFloat64 x, MetrcFloat64 y) {
    return !(y<x);
}

Tribool operator>=(MetrcFloat64 x, MetrcFloat64 y) {
    return !(x<y);
}



/************ Lower/UpperFloat64 ***********************************************/

Float64<Lower>::Float64() : _l(0.0) { }
Float64<Lower>::Float64(double lb, Exact) : _l(lb) { }
Float64<Lower>::Float64(Flt64 value) : _l(value) { }
Float64<Lower>::Float64(Integer const& z) : LowerFloat64(BoundFloat64(z)) { }
Float64<Lower>::Float64(Rational const& q) : LowerFloat64(BoundFloat64(q)) { }
Flt64 LowerFloat64::get_flt() const { return _l; }
double LowerFloat64::get_d() const { return _l; }
LowerFloat64 operator+(LowerFloat64 x) { return LowerFloat64(x._l); }
LowerFloat64 operator-(UpperFloat64 x) { return LowerFloat64(-x._u); }
LowerFloat64 operator+(LowerFloat64 x, LowerFloat64 y) { return LowerFloat64(add_opp(x._l,y._l)); }
LowerFloat64 operator-(LowerFloat64 x, UpperFloat64 y) {  return LowerFloat64(sub_opp(x._l,y._u)); }
LowerFloat64 operator*(LowerFloat64 x, LowerFloat64 y) { assert(x._l>=0.0 && y._l>=0.0); return LowerFloat64(mul_opp(x._l,y._l)); }
LowerFloat64 operator/(LowerFloat64 x, UpperFloat64 y) { assert(x._l>=0.0 && y._u>=0.0); return LowerFloat64(div_opp(x._l,y._u)); }
LowerFloat64 operator-(LowerFloat64 x, ErrorFloat64 y) { return LowerFloat64(sub_opp(x._l,y._e)); }

LowerFloat64 add(LowerFloat64 x, LowerFloat64 y) { return LowerFloat64(add_opp(x._l,y._l)); }
LowerFloat64 sub(LowerFloat64 x, UpperFloat64 y) { return LowerFloat64(sub_opp(x._l,y._u)); }
LowerFloat64 mul(LowerFloat64 x, LowerFloat64 y) { assert(x._l>=0.0 && y._l>=0.0); return LowerFloat64(mul_opp(x._l,y._l)); }
LowerFloat64 div(LowerFloat64 x, UpperFloat64 y) { assert(x._l>=0.0 && y._u>=0.0); return LowerFloat64(div_opp(x._l,y._u)); }
LowerFloat64 pos(LowerFloat64 x) { return LowerFloat64(+x._l); }
LowerFloat64 sqr(LowerFloat64 x) { assert(x._l>=0.0); return LowerFloat64(mul_opp(x._l,x._l)); }
LowerFloat64 neg(UpperFloat64 x) { return LowerFloat64(-x._u); }
LowerFloat64 rec(UpperFloat64 x) { assert(x._u>=0.0); return LowerFloat64(div_opp(1.0,x._u)); }
LowerFloat64 sqrt(LowerFloat64 x) { set_rounding_downward(); auto r_u=sqrt_rnd(x._l); set_default_rounding(); return LowerFloat64(r_u); }
LowerFloat64 exp(LowerFloat64 x) { set_rounding_downward(); auto r_u=exp_rnd(x._l); set_default_rounding(); return LowerFloat64(r_u); }
LowerFloat64 log(LowerFloat64 x) { set_rounding_downward(); auto r_u=log_rnd(x._l); set_default_rounding(); return LowerFloat64(r_u); }
LowerFloat64 atan(LowerFloat64 x) { ARIADNE_NOT_IMPLEMENTED; }

LowerFloat64 max(LowerFloat64 x, LowerFloat64 y) { return LowerFloat64(std::max(x._l,y._l)); }
LowerFloat64 min(LowerFloat64 x, LowerFloat64 y) { return LowerFloat64(std::min(x._l,y._l)); }

Float64<Upper>::Float64() : _u(0.0) { }
Float64<Upper>::Float64(double ub, Exact) : _u(ub) { }
Float64<Upper>::Float64(Flt64 value) : _u(value) { }
Float64<Upper>::Float64(Integer const& z) : UpperFloat64(BoundFloat64(z)) { }
Float64<Upper>::Float64(Rational const& q) : UpperFloat64(BoundFloat64(q)) { }
Flt64 UpperFloat64::get_flt() const { return _u; }
double UpperFloat64::get_d() const { return _u; }
UpperFloat64 operator+(UpperFloat64 x) { return UpperFloat64(x._u); }
UpperFloat64 operator-(LowerFloat64 x) { return UpperFloat64(-x._l); }
UpperFloat64 operator+(UpperFloat64 x, UpperFloat64 y) { return UpperFloat64(add_rnd(x._u,y._u)); }
UpperFloat64 operator-(UpperFloat64 x, LowerFloat64 y) {  return UpperFloat64(sub_rnd(x._u,y._l)); }
UpperFloat64 operator*(UpperFloat64 x, UpperFloat64 y) { assert(x._u>=0.0 && y._u>=0.0); return UpperFloat64(mul_rnd(x._u,y._u)); }
UpperFloat64 operator/(UpperFloat64 x, LowerFloat64 y) { assert(x._u>=0.0 && y._l>=0.0); return UpperFloat64(div_rnd(x._u,y._l)); }
UpperFloat64 operator+(UpperFloat64 x, ErrorFloat64 y) { return UpperFloat64(add_rnd(x._u,y._e)); }

UpperFloat64 add(UpperFloat64 x, UpperFloat64 y) { return UpperFloat64(add_rnd(x._u,y._u)); }
UpperFloat64 sub(UpperFloat64 x, LowerFloat64 y) { return UpperFloat64(sub_rnd(x._u,y._l)); }
UpperFloat64 mul(UpperFloat64 x, UpperFloat64 y) { assert(x._u>=0.0 && y._u>=0.0); return UpperFloat64(mul_rnd(x._u,y._u)); }
UpperFloat64 div(UpperFloat64 x, LowerFloat64 y) { assert(x._u>=0.0 && y._l>=0.0); return UpperFloat64(div_rnd(x._u,y._l)); }
UpperFloat64 pos(UpperFloat64 x) { return UpperFloat64(+x._u); }
UpperFloat64 sqr(UpperFloat64 x) { assert(x._u>=0.0); return UpperFloat64(mul_rnd(x._u,x._u)); }
UpperFloat64 neg(LowerFloat64 x) { return UpperFloat64(-x._l); }
UpperFloat64 rec(LowerFloat64 x) { assert(x._l>=0.0); return UpperFloat64(div_rnd(1.0,x._l)); }
UpperFloat64 sqrt(UpperFloat64 x) { return UpperFloat64(sqrt_rnd(x._u)); }
UpperFloat64 exp(UpperFloat64 x) { return UpperFloat64(exp_rnd(x._u)); }
UpperFloat64 log(UpperFloat64 x) { return UpperFloat64(log_rnd(x._u)); }
UpperFloat64 atan(UpperFloat64 x) { ARIADNE_NOT_IMPLEMENTED; }

UpperFloat64 max(UpperFloat64 x, UpperFloat64 y) { return UpperFloat64(std::max(x._u,y._u)); }
UpperFloat64 min(UpperFloat64 x, UpperFloat64 y) { return UpperFloat64(std::min(x._u,y._u)); }

NegSierpinski operator==(UpperFloat64 x, LowerFloat64 y) { return x._u>=y._l; }
NegSierpinski operator==(LowerFloat64 x, UpperFloat64 y) { return x._l<=y._u; }
Sierpinski operator!=(UpperFloat64 x, LowerFloat64 y) { return x._u< y._l; }
Sierpinski operator!=(LowerFloat64 x, UpperFloat64 y) { return x._l> y._u; }

Sierpinski operator< (UpperFloat64 x, LowerFloat64 y) { return x.get_d()< y.get_d(); }
Sierpinski operator> (LowerFloat64 x, UpperFloat64 y) { return x.get_d()> y.get_d(); }
Sierpinski operator<=(UpperFloat64 x, LowerFloat64 y) { return x.get_d()<=y.get_d(); }
Sierpinski operator>=(LowerFloat64 x, UpperFloat64 y) { return x.get_d()>=y.get_d(); }

NegSierpinski operator< (LowerFloat64 x, UpperFloat64 y) { return x.get_d()< y.get_d(); }
NegSierpinski operator> (UpperFloat64 x, LowerFloat64 y) { return x.get_d()> y.get_d(); }
NegSierpinski operator<=(LowerFloat64 x, UpperFloat64 y) { return x.get_d()<=y.get_d(); }
NegSierpinski operator>=(UpperFloat64 x, LowerFloat64 y) { return x.get_d()>=y.get_d(); }

ApprxFloat64 operator+(LowerFloat64 x, UpperFloat64 y) { return ApprxFloat64(x)+ApprxFloat64(y); }
ApprxFloat64 operator+(UpperFloat64 x, LowerFloat64 y) { return ApprxFloat64(x)+ApprxFloat64(y); }
ApprxFloat64 operator-(LowerFloat64 x, LowerFloat64 y) { return ApprxFloat64(x)-ApprxFloat64(y); }
ApprxFloat64 operator-(UpperFloat64 x, UpperFloat64 y) { return ApprxFloat64(x)-ApprxFloat64(y); }

NegSierpinski operator==(LowerFloat64 x, Rational const& y) { return x._l<=y; }
Sierpinski operator!=(LowerFloat64 x, Rational const& y) { return x._l> y; }
NegSierpinski operator< (LowerFloat64 x, Rational const& y) { return x._l< y; }
Sierpinski operator> (LowerFloat64 x, Rational const& y) { return x._l> y; }
NegSierpinski operator<=(LowerFloat64 x, Rational const& y) { return x._l<=y; }
Sierpinski operator>=(LowerFloat64 x, Rational const& y) { return x._l>=y; }

NegSierpinski operator==(UpperFloat64 x, Rational const& y) { return x._u>=y; }
Sierpinski operator!=(UpperFloat64 x, Rational const& y) { return x._u< y; }
Sierpinski operator< (UpperFloat64 x, Rational const& y) { return x._u< y; }
NegSierpinski operator> (UpperFloat64 x, Rational const& y) { return x._u> y; }
Sierpinski operator<=(UpperFloat64 x, Rational const& y) { return x._u<=y; }
NegSierpinski operator>=(UpperFloat64 x, Rational const& y) { return x._u>=y; }

OutputStream& operator<<(OutputStream& os, LowerFloat64 const& x) {
    rounding_mode_t rnd = get_rounding_mode();
    fesetround(FE_DOWNWARD);
    os << std::showpoint << std::setprecision(BoundedFloat64::get_output_precision()) << x._l;
    set_rounding_mode(rnd);
    return os;
}

OutputStream& operator<<(OutputStream& os, UpperFloat64 const& x) {
    rounding_mode_t rnd = get_rounding_mode();
    fesetround(FE_UPWARD);
    os << std::showpoint << std::setprecision(BoundFloat64::get_output_precision()) << x._u;
    set_rounding_mode(rnd);
    return os;
}



/************ ApprxFloat64 *****************************************************/

Float64<Apprx>::Float64()
    : _a(0.0) {
}

Float64<Apprx>::Float64(double value)
    : _a(value) {
}

Float64<Apprx>::Float64(Flt64 value)
    : _a(value) {
}

Float64<Apprx>::Float64(const Integer& z)
    : _a(z.get_si()) {
}

Float64<Apprx>::Float64(const Rational& q)
    : _a(q.get_d()) {
}

Flt64 ApprxFloat64::get_flt() const {
    return _a;
}

double ApprxFloat64::get_d() const {
    return _a;
}

ApprxFloat64 ApprxFloat64::create_zero() const { return ApprxFloat64(0.0); }

ApprxFloat64 operator+(ApprxFloat64 x) {
    return ApprxFloat64(x._a);
}

ApprxFloat64 operator-(ApprxFloat64 x) {
    return ApprxFloat64(-x._a);
}

ApprxFloat64 operator+(ApprxFloat64 x, ApprxFloat64 y) {
    return ApprxFloat64(add_near(x._a,y._a));
}

ApprxFloat64 operator-(ApprxFloat64 x, ApprxFloat64 y) {
    return ApprxFloat64(sub_near(x._a,y._a));
}

ApprxFloat64 operator*(ApprxFloat64 x, ApprxFloat64 y) {
    return ApprxFloat64(mul_near(x._a,y._a));
}

ApprxFloat64 operator/(ApprxFloat64 x, ApprxFloat64 y) {
    return ApprxFloat64(div_near(x._a,y._a));
}

ApprxFloat64& operator+=(ApprxFloat64& x, ApprxFloat64 y) {
    x._a=add_near(x._a,y._a); return x;
}

ApprxFloat64& operator-=(ApprxFloat64& x, ApprxFloat64 y) {
    x._a=sub_near(x._a,y._a); return x;
}

ApprxFloat64& operator*=(ApprxFloat64& x, ApprxFloat64 y) {
    x._a=mul_near(x._a,y._a); return x;
}

ApprxFloat64& operator/=(ApprxFloat64& x, ApprxFloat64 y) {
    x._a=div_near(x._a,y._a); return x;
}


ApprxFloat64 add(ApprxFloat64 x, ApprxFloat64 y) {
    return x+y;
}

ApprxFloat64 sub(ApprxFloat64 x, ApprxFloat64 y) {
    return x-y;
}

ApprxFloat64 mul(ApprxFloat64 x, ApprxFloat64 y) {
    return x*y;
}

ApprxFloat64 div(ApprxFloat64 x, ApprxFloat64 y) {
    return x/y;
}

ApprxFloat64 pos(ApprxFloat64 x) {
    return ApprxFloat64(+x._a);
}

ApprxFloat64 neg(ApprxFloat64 x) {
    return ApprxFloat64(-x._a);
}
ApprxFloat64 sqr(ApprxFloat64 x) {
    return ApprxFloat64(mul_near(x._a,x._a));
}

ApprxFloat64 rec(ApprxFloat64 x) {
    return ApprxFloat64(div_near(1.0,x._a));
}

ApprxFloat64 pow(ApprxFloat64 x, Nat m) {
    return std::pow(x._a,m);
}

ApprxFloat64 pow(ApprxFloat64 x, Int n) {
    return std::pow(x._a,n);
}

ApprxFloat64 sqrt(ApprxFloat64 x) {
    return std::sqrt(x._a);
}

ApprxFloat64 exp(ApprxFloat64 x) {
    return std::exp(x._a);
}

ApprxFloat64 log(ApprxFloat64 x) {
    return std::log(x._a);
}

ApprxFloat64 sin(ApprxFloat64 x) {
    return std::sin(x._a);
}

ApprxFloat64 cos(ApprxFloat64 x) {
    return std::cos(x._a);
}

ApprxFloat64 tan(ApprxFloat64 x) {
    return std::tan(x._a);
}

ApprxFloat64 atan(ApprxFloat64 x) {
    return std::atan(x._a);
}


ApprxFloat64 abs(ApprxFloat64 x) {
    return ApprxFloat64(fabs(x.get_d()));
}

ApprxFloat64 max(ApprxFloat64 x1, ApprxFloat64 x2) {
    return ApprxFloat64(std::max(x1.get_d(),x2.get_d()));
}

ApprxFloat64 min(ApprxFloat64 x1, ApprxFloat64 x2) {
    return ApprxFloat64(std::min(x1.get_d(),x2.get_d()));
}


ApprxFloat64 fma(ApprxFloat64 x, ApprxFloat64 y, ApprxFloat64 z) {
//    return std::fma(x,y,z);
    return x+y*z;
}

OutputStream& operator<<(OutputStream& os, ApprxFloat64 const& x) {
    os << std::setprecision(ApprxFloat64::output_precision) << x._a;
    return os;
}

Void ApprxFloat64::set_output_precision(int p) {
    output_precision=p;
}

int ApprxFloat64::output_precision = 4;


Fuzzy operator==(ApprxFloat64 x, ApprxFloat64 y) { return x.get_d()==y.get_d(); }
Fuzzy operator!=(ApprxFloat64 x, ApprxFloat64 y) { return x.get_d()!=y.get_d(); }
Fuzzy operator<=(ApprxFloat64 x, ApprxFloat64 y) { return x.get_d()<=y.get_d(); }
Fuzzy operator>=(ApprxFloat64 x, ApprxFloat64 y) { return x.get_d()>=y.get_d(); }
Fuzzy operator< (ApprxFloat64 x, ApprxFloat64 y) { return x.get_d()< y.get_d(); }
Fuzzy operator> (ApprxFloat64 x, ApprxFloat64 y) { return x.get_d()> y.get_d(); }

Bool same(ApprxFloat64 x, ApprxFloat64 y) { return x._a == y._a; }

/************ Flt64 conversions **********************************************/

ExactFloat64::operator MetrcFloat64 () const {
    return MetrcFloat64(_v,0.0);
}

ExactFloat64::operator BoundFloat64 () const {
    return BoundFloat64(_v,_v);
}

ExactFloat64::operator LowerFloat64 () const {
    return LowerFloat64(_v);
}

ExactFloat64::operator UpperFloat64 () const {
    return UpperFloat64(_v);
}

ExactFloat64::operator ApprxFloat64 () const {
    return ApprxFloat64(_v);
}

MetrcFloat64::operator BoundFloat64 () const {
    return BoundFloat64(this->lower(),this->upper());
}

MetrcFloat64::operator UpperFloat64 () const {
    return UpperFloat64(add_rnd(_v,_e));
}

MetrcFloat64::operator LowerFloat64 () const {
    return LowerFloat64(sub_opp(_v,_e));
}

MetrcFloat64::operator ApprxFloat64 () const {
    return ApprxFloat64(_v);
}

BoundFloat64::operator UpperFloat64 () const {
    return UpperFloat64(_u);
}

BoundFloat64::operator LowerFloat64 () const {
    return LowerFloat64(_l);
}

BoundFloat64::operator ApprxFloat64 () const {
    return ApprxFloat64(add_near(_l,_u)/2);
}

UpperFloat64::operator ApprxFloat64 () const {
    return ApprxFloat64(_u);
}

LowerFloat64::operator ApprxFloat64 () const {
    return ApprxFloat64(_l);
}

ExactFloat64 make_exact(const ApprxFloat64& x) {
    return ExactFloat64(x.get_d());
}

ErrorFloat64 make_error(const ApprxFloat64& x) {
    return ErrorFloat64(x.get_d());
}


/************ Flt64 creation **********************************************/

ApprxFloat64 Float64Creator::create(ApproximateNumber y, Approximate p) {
    return y.get(Approximate());
}

LowerFloat64 Float64Creator::create(LowerNumber y, Lower p) {
    return y.get(Lower());
}

UpperFloat64 Float64Creator::create(UpperNumber y, Upper p) {
    return y.get(Upper());
}

MetricFloat64 Float64Creator::create(ValidatedNumber y, Metric p) {
    return y.get(Metric());
}

BoundedFloat64 Float64Creator::create(ValidatedNumber y, Bounded p) {
    return y.get(Bounded());
}

BoundedFloat64 Float64Creator::create(ValidatedNumber y, Validated p) {
    return y.get(Bounded());
}

BoundedFloat64 Float64Creator::create(ValidatedNumber y, Effective p) {
    return y.get(Bounded());
}

BoundedFloat64 Float64Creator::create(ValidatedNumber y, Exact p) {
    return y.get(Bounded());
}



/************ Flt64 validated paradigm **********************************************/

Bool refines(ErrorFloat64 x, ErrorFloat64 y) {
    return (x._e <= y._e);
}


ExactFloat64 estimate(BoundFloat64 x) {
    return x.value();
}

ErrorFloat64 error(BoundFloat64 x) {
    return x.error();
}

BoundFloat64 refinement(BoundFloat64 x1, BoundFloat64 x2) {
    return BoundFloat64(max(x1._l,x2._l),min(x1._u,x2._u));
}

Bool refines(BoundFloat64 x, BoundFloat64 y) {
    return (x._l >= y._l) && (x._u <= y._u);
}

Bool inconsistent(BoundFloat64 x1, BoundFloat64 x2) {
    return (x1._l > x2._u) || (x1._u < x2._l);
}

Bool represents(BoundFloat64 x, ExactFloat64 y) {
    return (x._l <= y.get_d()) && (x._u >= y.get_d());
}


ExactFloat64 estimate(MetrcFloat64 x) {
    return x.value();
}

ErrorFloat64 error(MetrcFloat64 x) {
    return x.error();
}

MetrcFloat64 refinement(MetrcFloat64 x1, MetrcFloat64 x2) {
    return MetrcFloat64(refinement(BoundFloat64(x1),BoundFloat64(x2)));
}

Bool refines(MetrcFloat64 x, MetrcFloat64 y) {
    // Assumes rounding upwards
    auto d=(x._v>=y._v) ? sub_rnd(x._v,y._v) : sub_rnd(y._v,x._v);
    return add_rnd(d,x._e) <= y._e;
}

Bool inconsistent(MetrcFloat64 x1, MetrcFloat64 x2) {
    auto d=(x1._v>=x2._v) ? sub_rnd(x1._v,x2._v) : sub_rnd(x2._v,x1._v);
    return add_rnd(x1._e,x2._e)<=d;
}

Bool represents(MetrcFloat64 x, ExactFloat64 y) {
    if(x.value().get_flt()<=y.get_flt()) { return x.upper().get_flt()>=y.get_flt(); }
    else { return x.lower().get_flt()<=y.get_flt(); }
}


/************ Flt64 extended literals **********************************************/

ExactFloat64 operator"" _x (long double x) {
    return ExactFloat64(x);
}

ErrorFloat64 operator"" _e (long double x) {
    double ex=x; while(ex<x) { ex=next_rnd(ex); }
    return ErrorFloat64(ex);
}

UpperFloat64 operator"" _u (long double x) {
    double ux=x; while(ux<x) { ux=next_rnd(ux); }
    return UpperFloat64(ux);
}

LowerFloat64 operator"" _l (long double x) {
    double lx=x; while(lx>x) { lx=next_opp(lx); }
    return LowerFloat64(lx);
}

ApprxFloat64 operator"" _a (long double x) {
    return ApprxFloat64(x);
}

/************ Class name **********************************************/

template<> String class_name<Flt64>() { return "Float"; }
template<> String class_name<ExactFloat64>() { return "ExactFloat"; }
template<> String class_name<ErrorFloat64>() { return "ErrorFloat"; }
template<> String class_name<MetrcFloat64>() { return "ValidatedFloat"; }
template<> String class_name<BoundFloat64>() { return "BoundedFloat"; }
template<> String class_name<LowerFloat64>() { return "LowerFloat"; }
template<> String class_name<UpperFloat64>() { return "UpperFloat"; }
template<> String class_name<ApprxFloat64>() { return "ApproximateFloat"; }

} // namespace Ariadne
