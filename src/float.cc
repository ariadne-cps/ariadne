/***************************************************************************
 *            float.cc
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

#include "macros.h"
#include "dyadic.h"
#include "decimal.h"
#include "rational.h"
#include "rounding.h"

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

#include "float.h"
#include "float-exact.h"


namespace Ariadne {


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



static inline double next_rnd(double x) {
    volatile double y=+x; y=y+1e-300; y=y-1e-300; return +y;
}

static inline double next_opp(double x) {
    volatile double y=-x; y=y+1e-300; y=y-1e-300; return -y;
}


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

// Rounded power
double pow_rnd(double x, uint m)
{
    if(m==0) { return 1.0; }
    if(x==0) { return 0.0; }
    //if(x>=0.0 || (m%2==0)) { double r,p; r=1.0; p=(x>=0)?x:-x; while(m) { if(m%2) { r=mul_rnd(r,p); } p=mul_rnd(p,p); m/=2; } return r; }
    //else { double r,p r=-1.0; p=x; while(m) { if(m%2) { r=mul_rnd(-r,p); } p=mul_rnd(-p,p); m/=2; } return r; }
    if(x>0.0 || (m%2==0)) { volatile double r,p; r=1.0; p=std::abs(x); while(m) { if(m%2) { r=r*p; } p=p*p; m/=2; } return r; }
    else { volatile double r,p; r=-1.0; p=x; while(m) { if(m%2) { r=(-r)*p; } p=(-p)*p; m/=2; } return r; }
}

// Rounded power
double pow_rnd(double x, int n)
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
    while(a!=b) {
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
        case to_nearest: return _pi_approx;
        case downward: return _pi_down;
        case upward: return _pi_up;
        default: return _pi_approx;
    }
}

double pi_opp() {
    switch(get_rounding_mode()) {
        case to_nearest: return _pi_approx;
        case downward: return _pi_up;
        case upward: return _pi_down;
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

inline double max(double x1, double x2) { return std::max(x1,x2); }

inline double add_rnd(double x1, double x2) { return (volatile double&)x1+(volatile double&)x2; }
inline double sub_rnd(double x1, double x2) { return (volatile double&)x1-(volatile double&)x2; }
inline double mul_rnd(double x1, double x2) { return (volatile double&)x1*(volatile double&)x2; }
inline double div_rnd(double x1, double x2) { return (volatile double&)x1/(volatile double&)x2; }

inline double add_opp(double x, double y) { volatile double t=(-x)-y; return -t; }
inline double sub_opp(double x, double y) { volatile double t=(-x)+y; return -t; }
inline double mul_opp(double x, double y) { volatile double t=(-x)*y; return -t; }
inline double div_opp(double x, double y) { volatile double t=(-x)/y; return -t; }


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
            volatile double w=std::max(y1,y2);
            return neg_cos_rnd_series(w);
        } else {
            // Rounding downwards
            if(n_rnd%2==0) { return -1.0; }
            volatile double y1=sub_opp(pi_opp*n_opp,x);
            volatile double y2=sub_opp(x,pi_rnd*n_opp);
            assert(y1>=0 && y2>=0);
            volatile double w=std::max(y1,y2);
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

Float pow_rnd(Float x, int n)
{
    return pow_rnd(x.dbl,n);
}

Float sqrt_rnd(Float x)
{
    return sqrt_rnd(x.dbl);
}

Float exp_rnd(Float x)
{
    return exp_rnd(x.dbl);
}

Float log_rnd(Float x)
{
    return log_rnd(x.dbl);
}

Float sin_rnd(Float x)
{
    return sin_rnd(x.dbl);
}

Float cos_rnd(Float x)
{
    return cos_rnd(x.dbl);
}

Float tan_rnd(Float x)
{
    return tan_rnd(x.dbl);
}

} // namespace Ariadne


#ifdef ARIADNE_ENABLE_SERIALIZATION

#include "serialization.h"

namespace Ariadne {

void serialize(boost::archive::text_oarchive& a, Float& flt, const unsigned int v) {
    const double x=flt.get_d();
    a << x;
};

void serialize(boost::archive::text_iarchive& a, Float& flt, const unsigned int v) {
    flt=std::numeric_limits<double>::quiet_NaN();
    double x;
    a >> x;
    flt=x;
}

} // namespace Ariadne

#endif /* ARIADNE_ENABLE_SERIALIZATION */



