/***************************************************************************
 *            numeric/float64.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include "../utility/standard.hpp"

#include <iostream>
#include <iomanip>
#include <cassert>
#include <limits>



#include "../config.hpp"

#include "../utility/macros.hpp"
#include "../numeric/builtin.hpp"
#include "../numeric/twoexp.hpp"
#include "../numeric/dyadic.hpp"
#include "../numeric/decimal.hpp"
#include "../numeric/rational.hpp"
#include "../numeric/rounding.hpp"
#include "../numeric/floatdp.hpp"
#include "../numeric/floatmp.hpp"


namespace Ariadne {

/************  Publicly-accessible rounding-mode changing *******************/

typedef unsigned short rounding_mode_t;

Void set_rounding_mode(BuiltinRoundingModeType rnd) { _set_rounding_mode(rnd); }
BuiltinRoundingModeType get_rounding_mode() { return _get_rounding_mode(); }

Void set_rounding_to_nearest() { _set_rounding_to_nearest(); }
Void set_rounding_downward() { _set_rounding_downward(); }
Void set_rounding_upward() { _set_rounding_upward(); }
Void set_rounding_toward_zero() { _set_rounding_toward_zero(); }

Void set_default_rounding() { _set_rounding_upward(); }

static const double _pi_up=3.1415926535897936;
static const double _pi_down=3.1415926535897931;

static const double _quarter_pi_approx=0.78539816339744828;
static const double _pi_approx=3.1415926535897931;
static const double _two_pi_approx=6.2831853071795862;
static const double _sqrt2_approx=0.70710678118654757;
static const double _log2_approx=0.6931471805599453094;

static const double _pi_near=3.1415926535897931;

// The floor of log_10(|x|). abslog10floor(x)+1 is the number of digits of x before the decimal point.
// Returns the smallest representable  int  if x is zero
int abslog10floor(double x) {
    return static_cast<int>(std::floor(std::log10(std::abs(x))));
    //static const double MIN_INT=-2147483648.;
    //return std::max(std::floor(std::log10(std::abs(x))),MIN_INT);
}

static inline double next_rnd(double x) {
    volatile double y=+x; y=y+1e-300; y=y-1e-300; return +y;
}

static inline double next_opp(double x) {
    volatile double y=-x; y=y+1e-300; y=y-1e-300; return -y;
}


static inline double horner_rnd(Int n, double x, const long long int* c)
{
    volatile double y=1./c[n];
    for(Int i=n-1; i>=0; --i) {
        y=1.0/c[i]+x*y;
    }
    return y;
}

static inline double horner_opp(Int n, double x, const long long int* c)
{
    volatile double y=-1./c[n];
    for(Int i=n-1; i>=0; --i) {
        y=-1.0/c[i]+x*y;
    }
    return -y;
}

// Rounded power
double pow_rnd(double x, Nat m)
{
    if(m==0) { return 1.0; }
    if(x==0) { return 0.0; }
    //if(x>=0.0 || (m%2==0)) { double r,p; r=1.0; p=(x>=0)?x:-x; while(m) { if(m%2) { r=mul_rnd(r,p); } p=mul_rnd(p,p); m/=2; } return r; }
    //else { double r,p r=-1.0; p=x; while(m) { if(m%2) { r=mul_rnd(-r,p); } p=mul_rnd(-p,p); m/=2; } return r; }
    if(x>0.0 || (m%2==0)) { volatile double r,p; r=1.0; p=std::abs(x); while(m) { if(m%2) { r=r*p; } p=p*p; m/=2; } return r; }
    else { volatile double r,p; r=-1.0; p=x; while(m) { if(m%2) { r=(-r)*p; } p=(-p)*p; m/=2; } return r; }
}

// Rounded power
double pow_rnd(double x, Int n)
{
    if(n>=0) { return pow_rnd(x,Nat(n)); }
    ARIADNE_ASSERT(x!=0.0);
    if(x>0.0 || (n%2==-1)) { volatile double r=1.0/x; return pow_rnd(r,Nat(-n)); }
    else { volatile double r=-1.0/x; return pow_rnd(r,Nat(-n)); }
}

double sqrt_rnd(double x)
{
    // long int c[]={ 0, 6, -360, 15120, -604800, 23950080, -946218790, 37362124800 };
    ARIADNE_ASSERT_MSG(x>=0, " x = "<<x<<"\n");

    if(x==0.0) { return 0.0; }
    Int n; volatile double y,a,b;
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
    for(Nat i=1; i!=20; ++i) {
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
    Int e;
    volatile double z=frexp(y,&e);
    z=ldexp(z,m+e);
    return z;
}




double log_rnd(double x) {
    static const long long int c[12]={ 1LL, 3LL, 5LL, 7LL, 9LL, 11LL, 13LL, 15LL, 17LL, 19LL, 21LL, 23LL };

    ARIADNE_ASSERT_MSG(x>=0.0,"log(x): x="<<x);
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

    if(x==0.0) { return std::numeric_limits<double>::infinity(); }
    if(x==1.0) { return 0.0; }

    Int n;
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
        case ROUND_TO_NEAREST: return _pi_approx;
        case ROUND_DOWNWARD: return _pi_down;
        case ROUND_UPWARD: return _pi_up;
        default: return _pi_approx;
    }
}

double pi_opp() {
    switch(get_rounding_mode()) {
        case ROUND_TO_NEAREST: return _pi_approx;
        case ROUND_DOWNWARD: return _pi_up;
        case ROUND_UPWARD: return _pi_down;
        default: return _pi_approx;
    }
}

double sin_rnd(double x) {
    volatile double two_pi_rnd=2*pi_rnd();
    volatile double two_pi_opp=2*pi_opp();

    volatile double half_pi_rnd=two_pi_rnd/4;
    volatile double half_pi_opp=two_pi_opp/4;

    Int q = (long int)(std::floor(x/_quarter_pi_approx)) % 8;
    if(q<-4) { q+=8; } if(q>=4) { q-=8; }
    volatile double n=-std::floor(x/_two_pi_approx+0.5);

    volatile double y,w,s;

    // Set to true if sin is decreasing so we want opposite rounding
    Bool want_opposite=(q<-2||q>=2);
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
    default: { s=0; assert(false); }
    }

    return s;
}

inline double max(double x1, double x2) { return std::max(x1,x2); }

double sqr_rnd(double x) { return (volatile double&)x*(volatile double&)x; }
double rec_rnd(double x) { return 1.0/(volatile double&)x; }
double add_rnd(double x1, double x2) { return (volatile double&)x1+(volatile double&)x2; }
double sub_rnd(double x1, double x2) { return (volatile double&)x1-(volatile double&)x2; }
double mul_rnd(double x1, double x2) { return (volatile double&)x1*(volatile double&)x2; }
double div_rnd(double x1, double x2) { return (volatile double&)x1/(volatile double&)x2; }

double add_opp(double x, double y) { volatile double t=(-x)-y; return -t; }
double sub_opp(double x, double y) { volatile double t=(-x)+y; return -t; }
double mul_opp(double x, double y) { volatile double t=(-x)*y; return -t; }
double div_opp(double x, double y) { volatile double t=(-x)/y; return -t; }

double neg_rec_rnd(double x) { return (-1.0)/(volatile double&)x; }
double neg_rec_opp(double x) { volatile double t=1.0/x; return -t; }


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

    Int q = (long int)(std::floor(y/_quarter_pi_approx)) % 8;
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
        w=0;c=0;
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
        for(Int i=11; i>=0; --i) {
            c=double(cn[i])/cd[i];
            w=c+s*w;
        }
        r=x*w;
    } else {
        s=(-x)*x; s=-s;
        w=double(-cn[12])/cd[12];
        for(Int i=12; i>=0; --i) {
            c=double(-cn[i])/cd[i];
            w=c+s*w;
        }
        r=x*(-w);
    }
    return r;
}

double atan_rnd_series(double x) {
    // atan(x) = \sum_{n=0}^{\infty} (-1)^n/(2n+1) x^{2n+1}
    assert(std::abs(x)<0.5);
    static const long long int c[19]={ 1LL, -3LL, 5LL, -7LL, 9LL, -11LL, 13LL, -15LL, 17LL, -19LL, 21LL, -23LL, 25LL, -27LL, 29LL, -31LL, 33LL, -35LL, 37LL };

    volatile double s,w,r;
    if(x>=0) {
        s=mul_rnd(x,x);
        w=horner_rnd(18,s,c);
        r=x*w;
    } else {
        s=mul_opp(x,x);
        w=horner_opp(17,s,c);
        r=x*w;
    }
    return r;
}

double atan_rnd(double x) {
    // use range reduction
    // atan(-x) = -atan(x)
    // atan(1/x) = pi/2 - arctan(x) for x>0
    // atan(1/x) = -pi/2 - arctan(x) for x<0
    //
    // atan(x) = sgn(x)*pi/2 + arctan(-1/x)
    // atan(x) = atan(c) + atan((x-c)/(1+x*c))
    // atan(x) = 2*atan(c) where c = x/(1+sqrt(1+x^2))

    // For x in [0,1], perform range reduction twice
    // First, take c=sqrt(2)-1~=0.4142, reducing x to [0,c]
    // Then, take d=c/(1+sqrt(1+c^2))~=0.1989, reducing x to [0,d]

    if(x>1.0) { return pi_rnd()/2+atan_rnd(-1.0/x); }
    if(x<-1.0) { return -pi_opp()/2+atan_rnd(-1.0/x); }

    // The values below are closest values for c=sqrt(2)-1 and d=c/(1+sqrt(1+c^2))
    //  static const long double c=0.414213562373095048801689l;
    //  static const long double d=0.198912367379658006911598l;
    //  static const long double atan_c=0.392699081698724154807830l;
    //  static const long double atan_d=0.196349540849362077403915l;

    // TODO: Check if method below is sound i.e. error in c does not give error in result
    // The values below are can be taked as exact values for approximately given c and d
    static const double atan_c=0.392699081698724139;
    static const double atan_d=0.196349540849362077;

    static const long double c=0.41421356237309503087445222702100977585359942167997360229492187500l;
    static const long double d=0.198912367379657998957341417944899575331874075345695018768310546875l;

    volatile long double neg_c=-c;
    volatile long double neg_d=-d;

    if(x>c) {
        volatile double t=-1.0+neg_c*x;
        volatile double r=-1.0/t;
        return atan_c + atan_rnd((x-c)*r);
    } else if(x<-c) {
        volatile double t=-1.0+c*x;
        volatile double r=-1.0/t;
        return (-atan_c) + atan_rnd((x+c)*r);
    } else if(x>d) {
        volatile double t=-1.0+neg_d*x;
        volatile double r=-1.0/t;
        return atan_d + atan_rnd_series((x-d)*r);
    } else if(x<-d) {
        volatile double t=-1.0+d*x;
        volatile double r=-1.0/t;
        return (-atan_d) + atan_rnd((x+d)*r);
    }
    return atan_rnd_series(x);

}



FloatDP::FloatDP(ExactDouble const& d, PrecisionType)
    : FloatDP(d.get_d())
{
}

FloatDP::FloatDP(TwoExp const& t, PrecisionType)
    : FloatDP(std::ldexp(1.0,t.exponent()))
{
    ARIADNE_ASSERT(Dyadic(*this)==Dyadic(t));
}

FloatDP::FloatDP(Dyadic const& w, PrecisionType)
    : FloatDP(w.get_d())
{
    ARIADNE_ASSERT(Dyadic(*this)==w || is_nan(w));
}

FloatDP::FloatDP(double d, RoundingModeType rnd, PrecisionType)
    : FloatDP(d)
{
}

FloatDP::FloatDP(Integer const& x, RoundingModeType rnd, PrecisionType pr)
    : FloatDP(Dyadic(x),rnd,pr)
{
}

FloatDP::FloatDP(Dyadic const& w, RoundingModeType rnd, PrecisionType pr)
    : FloatDP(w.get_d())
{
    if (is_finite(w)) {
         RoundingModeType old_rnd=get_rounding_mode();
         if(rnd==ROUND_UPWARD) {
             set_rounding_upward();
             while (Dyadic(dbl)<w) { dbl+=std::numeric_limits<double>::min(); }
             set_rounding_mode(old_rnd);
         }
         if(rnd==ROUND_DOWNWARD) {
             set_rounding_downward();
             while (Dyadic(dbl)>w) { dbl-=std::numeric_limits<double>::min(); }
             set_rounding_mode(old_rnd);
         }
     }
}

FloatDP::FloatDP(Decimal const& dec, RoundingModeType rnd, PrecisionType pr)
    : FloatDP(Rational(dec),rnd,pr)
{
}

FloatDP::FloatDP(Rational const& q, RoundingModeType rnd, PrecisionType pr)
    : FloatDP(q.get_d())
{
    if (is_finite(q)) {
        RoundingModeType old_rnd=get_rounding_mode();
        if(rnd==ROUND_UPWARD) {
            set_rounding_upward();
            while (Rational(dbl)<q) { dbl+=std::numeric_limits<double>::min(); }
            set_rounding_mode(old_rnd);
        }
        if(rnd==ROUND_DOWNWARD) {
            set_rounding_downward();
            while (Rational(dbl)>q) { dbl-=std::numeric_limits<double>::min(); }
            set_rounding_mode(old_rnd);
        }
    }
}

FloatDP::FloatDP(FloatDP const& x, RoundingModeType rnd, PrecisionType)
    : FloatDP(x)
{
}

FloatDP::operator Dyadic () const {
    return Dyadic(this->dbl);
}

FloatDP::operator Rational () const {
    return Rational(this->dbl);
}

FloatDP pow_rnd(FloatDP x, Int n)
{
    return pow_rnd(x.dbl,n);
}

FloatDP sqrt_rnd(FloatDP x)
{
    return sqrt_rnd(x.dbl);
}

FloatDP exp_rnd(FloatDP x)
{
    return exp_rnd(x.dbl);
}

FloatDP log_rnd(FloatDP x)
{
    return log_rnd(x.dbl);
}

FloatDP sin_rnd(FloatDP x)
{
    return sin_rnd(x.dbl);
}

FloatDP cos_rnd(FloatDP x)
{
    return cos_rnd(x.dbl);
}

FloatDP tan_rnd(FloatDP x)
{
    return tan_rnd(x.dbl);
}

FloatDP atan_rnd(FloatDP x)
{
    return atan_rnd(x.dbl);
}

FloatDP FloatDP::pi(BuiltinRoundingModeType rnd, DoublePrecision pr) {
    switch(rnd) {
        case FloatDP::ROUND_UPWARD: return _pi_up;
        case FloatDP::ROUND_DOWNWARD: return _pi_down;
        case FloatDP::ROUND_TO_NEAREST: return _pi_near;
        default: assert(false); return _pi_near;
    }
}

FloatDP::RoundingModeType FloatDP::get_rounding_mode() { return Ariadne::get_rounding_mode(); }
Void FloatDP::set_rounding_mode(RoundingModeType rnd) { Ariadne::set_rounding_mode(rnd); }
Void FloatDP::set_rounding_downward() { Ariadne::set_rounding_downward(); }
Void FloatDP::set_rounding_upward() { Ariadne::set_rounding_upward(); }
Void FloatDP::set_rounding_to_nearest() { Ariadne::set_rounding_to_nearest(); }
Void FloatDP::set_rounding_toward_zero() { Ariadne::set_rounding_toward_zero(); }

FloatDP::PrecisionType FloatDP::get_default_precision() { return FloatDP::PrecisionType(); }
FloatDP::PrecisionType FloatDP::precision() const { return FloatDP::PrecisionType(); }
Void FloatDP::set_precision(FloatDP::PrecisionType) { }

FloatDP FloatDP::min(PrecisionType) { return std::numeric_limits<double>::min(); }
FloatDP FloatDP::max(PrecisionType) { return std::numeric_limits<double>::max(); }
FloatDP FloatDP::eps(PrecisionType) { return std::numeric_limits<double>::epsilon(); }

FloatDP FloatDP::inf(Sign sgn, PrecisionType pr) {
    switch (sgn) {
    case Sign::POSITIVE: return std::numeric_limits<double>::infinity();
    case Sign::NEGATIVE: return -std::numeric_limits<double>::infinity();
    default: return std::numeric_limits<double>::quiet_NaN();
    }
}
FloatDP FloatDP::inf(PrecisionType pr) { return FloatDP::inf(Sign::POSITIVE,pr); }
FloatDP FloatDP::nan(PrecisionType pr) { return FloatDP::inf(Sign::ZERO,pr); }

template<class R, class A> R integer_cast(A const&);
template<> Nat integer_cast<Nat,FloatDP>(FloatDP const& x) { return static_cast<Nat>(x.dbl); }
template<> Int integer_cast<Int,FloatDP>(FloatDP const& x) { return static_cast<Int>(x.dbl); }


Comparison cmp(FloatDP x1, Rational const& q2) {
    if(std::isfinite(x1.get_d())) { return cmp(Rational(x1),q2); }
    else { return x1.get_d()>0.0 ? Comparison::GREATER : Comparison::LESS; }
}
Comparison cmp(Rational const& q1, FloatDP x2) {
    if(std::isfinite(x2.get_d())) { return cmp(q1,Rational(x2)); }
    else { return x2.get_d()>0.0 ? Comparison::LESS : Comparison::GREATER; }
}

/*
OutputStream& write(OutputStream& os, FloatDP const& x, Nat bits, BuiltinRoundingModeType rnd) {
    Nat dgts = std::ceil(bits*std::log(2))/std::log(10);
    FloatDP::RoundingModeType old_rnd=FloatDP::get_rounding_mode();
    FloatDP::set_rounding_mode(rnd);
    os << std::setprecision(dgts) << x.get_d();
    FloatDP::set_rounding_mode(old_rnd);
    return os;
}
*/

OutputStream& write(OutputStream& os, FloatMP const& x, DecimalPlaces plcs, MPFRRoundingModeType rnd);
OutputStream& write(OutputStream& os, FloatMP const& x, DecimalPrecision figs, MPFRRoundingModeType rnd);

MPFRRoundingModeType to_mpfr_rounding_mode(BuiltinRoundingModeType rnd) {
    assert(rnd==ROUND_TO_NEAREST || rnd==ROUND_UPWARD || rnd==ROUND_DOWNWARD);
    return (rnd==FloatDP::ROUND_TO_NEAREST) ? FloatMP::ROUND_TO_NEAREST
               : (rnd==FloatDP::ROUND_DOWNWARD) ? FloatMP::ROUND_DOWNWARD : FloatMP::ROUND_UPWARD;
}

OutputStream& write(OutputStream& os, FloatDP const& x, DecimalPlaces plcs, BuiltinRoundingModeType rnd) {
    MultiplePrecision pr_mp(53);
    MPFRRoundingModeType rnd_mp=to_mpfr_rounding_mode(rnd);
    return write(os,FloatMP(x,pr_mp),plcs,rnd_mp);
}

OutputStream& write(OutputStream& os, FloatDP const& x, DecimalPrecision figs, BuiltinRoundingModeType rnd) {
    MultiplePrecision pr_mp(53);
    MPFRRoundingModeType rnd_mp=to_mpfr_rounding_mode(rnd);
    return write(os,FloatMP(x,pr_mp),figs,rnd_mp);
}


OutputStream& operator<<(OutputStream& os, FloatDP const& x) {
    return os << x.dbl;
}

InputStream& operator>>(InputStream& is, FloatDP& x) {
    double r; is >> r; x.dbl=r; return is;
}

template<> String class_name<double>() { return "double"; }

template<> String class_name<FloatDP>() { return "FloatDP"; }

} // namespace Ariadne
