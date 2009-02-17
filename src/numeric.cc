/***************************************************************************
 *            numeric.cc
 *
 *  Copyright 2008  Pieter Collins
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
#include <cassert>
#include <stdint.h>

#include "config.h"
#include "macros.h"

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

#include "numeric.h"

namespace Ariadne {


static const double quarter_pi_up=0.78539816339744839;
static const double half_pi_up=1.5707963267948968;
static const double pi_up=3.1415926535897936;
static const double pi_down=3.1415926535897931;

static const double quarter_pi_approx=0.78539816339744828;
static const double half_pi_approx=1.5707963267948966;
static const double pi_approx=3.1415926535897931;
static const double two_pi_approx=6.2831853071795862;
static const double sqrt2_approx=0.70710678118654757;
static const double log2_approx=0.6931471805599453094;
  

uint16_t 
fac(uint16_t n) 
{ 
    ARIADNE_ASSERT(n<9); // Maximum factorial in 32 bits
    uint16_t  r=1; 
    for(uint16_t i=1; i<=n; ++i) { 
        r*=i; 
    } 
    return r; 
}


uint16_t 
bin(uint16_t n, uint16_t k) 
{ 
    ARIADNE_ASSERT(n<16);  // Maximum computable bin(n,n/2) using 16 bits
                           // Note that this is shorter than the maximum representable factorial
    if(k>n+1) { std::cerr << "bin("<<n<<","<<k<<")\n"; }
    if(k==n+1) { return 0; }
    ARIADNE_ASSERT(k<=n);
    uint16_t r=1; 
    for(uint16_t i=1; i<=k; ++i) { 
        r*=(n+1-i); 
        r/=i; 
    } 
    return r; 
}

uint32_t 
fac(uint32_t n) 
{ 
    ARIADNE_ASSERT(n<13); // Maximum factorial in 32 bits
    uint32_t  r=1; 
    for(uint32_t i=1; i<=n; ++i) { 
        r*=i; 
    } 
    return r; 
}


uint32_t 
bin(uint32_t n, uint32_t k) 
{ 
    ARIADNE_ASSERT(n<31);  // Maximum computable bin(n,n/2) using 32 bits
                           // Note that this is shorter than the maximum representable factorial
    if(k>n+1) { std::cerr << "bin("<<n<<","<<k<<")\n"; }
    if(k==n+1) { return 0; }
    ARIADNE_ASSERT(k<=n);
    uint32_t r=1; 
    for(uint32_t i=1; i<=k; ++i) { 
        r*=(n+1-i); 
        r/=i; 
    } 
    return r; 
}



uint64_t
fac(uint64_t n) 
{ 
    ARIADNE_ASSERT(n<21); // Maximum factorial in 64 bits
    uint64_t  r=1; 
    for(uint64_t i=1; i<=n; ++i) { 
        r*=i; 
    } 
    return r; 
}


uint64_t
bin(uint64_t n, uint64_t k) 
{ 
    ARIADNE_ASSERT(n<63);  // Maximum computable bin(n,n/2) using 64 bits
                           // Note that this is shorter than the maximum representable factorial
    if(k>n+1) { std::cerr << "bin("<<n<<","<<k<<")\n"; }
    if(k==n+1) { return 0; }
    ARIADNE_ASSERT(k<=n);
    uint64_t r=1; 
    for(uint64_t i=1; i<=k; ++i) { 
        r*=(n+1-i); 
        r/=i; 
    } 
    return r; 
}


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
    if(x>0.0 || (m%2==0)) { volatile double r,p; r=1.0; p=abs(x); while(m) { if(m%2) { r=r*p; } p=p*p; m/=2; } return r; }
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
    ARIADNE_ASSERT(x>=0);
    
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
        //std::cerr<<i<<" "<<t<<" "<<r<<"\n";
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
    double n=floor(x/log2_approx+0.5);
    volatile double r,s,t,w,y;
    volatile double log2;
    log2=(n>0.0) ? next_opp(log2_approx) : next_rnd(log2_approx);
    r=x+(-n)*log2;

    ARIADNE_ASSERT(r>=-0.4);
    ARIADNE_ASSERT(r<=+0.4);

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
    if(y<sqrt2_approx) { y*=2; --n; }

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
    
    volatile double log2rnd=next_rnd(log2_approx);
    
  
    return log2rnd*n+ly;
}

double pi_rnd() {
    switch(get_rounding_mode()) {
        case to_nearest: return pi_approx;
        case downward: return pi_down;
        case upward: return pi_up;
        default: return pi_approx;
    } 
}

double pi_opp() {
    switch(get_rounding_mode()) {
        case to_nearest: return pi_approx;
        case downward: return pi_up;
        case upward: return pi_down;
        default: return pi_approx;
    } 
}

double sin_rnd(double x) { 
    //std::cerr<<"sin_rnd("<<x<<")\n";
    volatile double two_pi_rnd=2*pi_rnd();
    volatile double two_pi_opp=2*pi_opp();

    volatile double half_pi_rnd=two_pi_rnd/4;
    volatile double half_pi_opp=two_pi_opp/4;

    int q = (long int)(std::floor(x/quarter_pi_approx)) % 8;
    if(q<-4) { q+=8; } if(q>=4) { q-=8; }     
    volatile double n=-std::floor(x/two_pi_approx+0.5);

    volatile double y,w,s;

    // Set to true if sin is decreasing so we want opposite rounding
    bool want_opposite=(q<-2||q>=2); 
    // if n is negative then we need to switch rounding of two_pi
    volatile double two_pi_corr=((n>=0.0) ^ want_opposite) ? two_pi_rnd : two_pi_opp;

    // Scale onto interval from -pi to pi
    if(want_opposite) { y=-x+(-n)*(two_pi_corr); y=-y; } else { y=x+n*two_pi_corr; }
   
    assert(-two_pi_approx<=y && y<=two_pi_approx);


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
    //std::cerr<<"cos_rnd("<<x<<")\n";
    const double pi_rnd=Ariadne::pi_rnd();
    const double pi_opp=Ariadne::pi_opp();
    const double half_pi_rnd=pi_rnd/2;

    if(x==0.0) { return 1.0; }

    // Compute a constant n such that n*pi<x<(n+1)*pi 
    long int n_rnd=(long int)(std::floor(div_rnd(x,pi_opp)));
    long int n_opp=(long int)(std::floor(div_opp(x,pi_rnd)));
    
    if(n_rnd!=n_opp) { 
        //std::cerr<<"  n_rnd="<<n_rnd<<" n_opp="<<n_opp<<"\n";
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

    int q = (long int)(std::floor(y/quarter_pi_approx)) % 8;
    assert(q<=4);
    //std::cerr<<"  n="<<n_rnd<<" y="<<y<<" q="<<q;
   
    volatile double w,c;
    if(q==0) {
        w=y;
        //std::cerr<<" w="<<w<<"\n";
        c=pos_cos_rnd_series(w);
    } else if(q==1 || q==2) {
        w=sub_rnd(pi_rnd/2,y);
        //std::cerr<<" w="<<w<<"\n";
        if(w>=0.0) { c=pos_sin_rnd_series(w); }
        else { c=neg_sin_rnd_series(-w); }
    } else if(q==3 || q==4) {
        w=sub_opp(pi_opp,y);
        //std::cerr<<" w="<<w<<"\n";
        c=neg_cos_rnd_series(w);
    } else {
        assert(false);
    }

    return c;

        

    volatile double z=0.0;

    ARIADNE_ASSERT(-two_pi_approx<=y && y<=two_pi_approx);
    using std::cerr; using std::endl;
    switch(q) { 
    case -4: { w = +y + 2*half_pi_rnd; w=+w; w=max(w,z); cerr<<w<<endl; c=neg_cos_rnd_series(w); break; }
    case -3: { w = +y + 1*half_pi_rnd; w=-w; w=max(w,z); cerr<<w<<endl; c=neg_sin_rnd_series(w); break; }
    case -2: { w = +y + 1*half_pi_rnd; w=+w; w=max(w,z); cerr<<w<<endl; c=pos_sin_rnd_series(w); break; }
    case -1: { w = +y + 0*half_pi_rnd; w=-w; w=max(w,z); cerr<<w<<endl; c=pos_cos_rnd_series(w); break; }
    case +0: { w = -y + 0*half_pi_rnd; w=-w; w=max(w,z); cerr<<w<<endl; c=pos_cos_rnd_series(w); break; }
    case +1: { w = -y + 1*half_pi_rnd; w=+w; w=max(w,z); cerr<<w<<endl; c=pos_sin_rnd_series(w); break; }
    case +2: { w = -y + 1*half_pi_rnd; w=-w; w=max(w,z); cerr<<w<<endl; c=neg_sin_rnd_series(w); break; }
    case +3: { w = -y + 2*half_pi_rnd; w=+w; w=max(w,z); cerr<<w<<endl; c=neg_cos_rnd_series(w); break; }
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

    double n=std::floor(x/pi_approx+0.5);

    volatile double pi_corr=(n>=0) ? next_opp(pi_approx) : next_rnd(pi_approx);
    y=x-n*pi_corr;


    ARIADNE_ASSERT(y>=-pi_up/2);
    ARIADNE_ASSERT(y<=+pi_up/2);

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
    ARIADNE_ASSERT(x>=-pi_up/8);
    ARIADNE_ASSERT(x<=+pi_up/8);

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


static inline Float cos_down(Float x) { 
    set_rounding_downward(); Float y=cos_rnd(x); return y; 
}

static inline Float cos_up(Float x) { 
    set_rounding_upward(); Float y=cos_rnd(x); return y; 
}


inline double _add_down(volatile double x, volatile double y) { set_rounding_downward(); return x+y; }
inline double _add_up(volatile double x, volatile double y) { set_rounding_upward(); return x+y; }
inline double _sub_down(volatile double x, volatile double y) { set_rounding_downward(); return x-y; }
inline double _sub_up(volatile double x, volatile double y) { set_rounding_upward(); return x-y; }
inline double _mul_down(volatile double x, volatile double y) { set_rounding_downward(); return x*y; }
inline double _mul_up(volatile double x, volatile double y) { set_rounding_upward(); return x*y; }
inline double _div_down(volatile double x, volatile double y) { set_rounding_downward(); return x/y; }
inline double _div_up(volatile double x, volatile double y) { set_rounding_upward(); return x/y; }

Interval trunc(Interval x, uint n) 
{
    Interval e=Interval(pow(2.0,52-n));
    Interval y=x+e;
    return y-e;
}

Interval rec(Interval i) 
{
    volatile double rl,ru;
    if(i.l>0 || i.u<0) {
        rounding_mode_t rnd=get_rounding_mode();
        rl=_div_down(1.0,i.u); 
        ru=_div_up(1.0,i.l);
        set_rounding_mode(rnd);
    } else {
        rl=-inf();
        ru=+inf();
        ARIADNE_THROW(DivideByZeroException,"Interval rec(Interval ivl)","ivl="<<i);
    }
    return Interval(rl,ru);
}


Interval mul(Interval i1, Interval i2) 
{
    volatile double rl,ru;
    rounding_mode_t rnd=get_rounding_mode();
    if(i1.l>=0) {
        if(i2.l>=0) {
            rl=_mul_down(i1.l,i2.l); ru=_mul_up(i1.u,i2.u);
        } else if(i2.u<=0) {
            rl=_mul_down(i1.u,i2.l); ru=_mul_up(i1.l,i2.u);
        } else {
            rl=_mul_down(i1.u,i2.l); ru=_mul_up(i1.u,i2.u);
        }
    }
    else if(i1.u<=0) {
        if(i2.l>=0) {
            rl=_mul_down(i1.l,i2.u); ru=_mul_up(i1.u,i2.l);
        } else if(i2.u<=0) {
            rl=_mul_down(i1.u,i2.u); ru=_mul_up(i1.l,i2.l);
        } else {
            rl=_mul_down(i1.l,i2.u); ru=_mul_up(i1.l,i2.l);
        }
    } else {
        if(i2.l>=0) {
            rl=_mul_down(i1.l,i2.u); ru=_mul_up(i1.u,i2.u);
        } else if(i2.u<=0) {
            rl=_mul_down(i1.u,i2.l); ru=_mul_up(i1.l,i2.l);
        } else {
            set_rounding_mode(downward);
            rl=min(i1.u*i2.l,i1.l*i2.u);
            set_rounding_mode(upward);
            ru=max(i1.l*i2.l,i1.u*i2.u);
        }
    }
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}


Interval mul(Interval i1, Float x2) 
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double rl,ru;
    if(x2>=0) {
        rl=_mul_down(i1.l,x2); ru=_mul_up(i1.u,x2);
    } else {
        rl=_mul_down(i1.u,x2); ru=_mul_up(i1.l,x2);
    }
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}


Interval div(Interval i1, Interval i2) 
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double rl,ru;
    if(i2.l>=0) {
        if(i1.l>=0) {
            rl=_div_down(i1.l,i2.u); ru=_div_up(i1.u,i2.l);
        } else if(i1.u<=0) {
            rl=_div_down(i1.l,i2.l); ru=_div_up(i1.u,i2.u);
        } else {
            rl=_div_down(i1.l,i2.l); ru=_div_up(i1.u,i2.l);
        }
    }
    else if(i2.u<=0) {
        if(i1.l>=0) {
            rl=_div_down(i1.u,i2.u); ru=_div_up(i1.l,i2.l);
        } else if(i1.u<=0) {
            rl=_div_down(i1.u,i2.l); ru=_div_up(i1.l,i2.u);
        } else {
            rl=_div_down(i1.u,i2.u); ru=_div_up(i1.l,i2.u);
        } 
    }
    else {
        ARIADNE_THROW(DivideByZeroException,"Interval div(Interval ivl1, Interval ivl2)","ivl1="<<i1<<", ivl2="<<i2);
        rl=-inf(); ru=+inf();
    }
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}



Interval div(Interval i1, Float x2) 
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double rl,ru;
    if(x2>0) {
        rl=_div_down(i1.l,x2); ru=_div_up(i1.u,x2); 
    } else if(x2<0) {
        rl=_div_down(i1.u,x2); ru=_div_up(i1.l,x2);
    } else {
        rl=-inf(); ru=+inf(); 
    }
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}


Interval div(Float x1, Interval i2) 
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double rl,ru;
    if(i2.l<=0 && i2.u>=0) {
        ARIADNE_THROW(DivideByZeroException,"Interval div(Float x1, Interval ivl2)","x1="<<x1<<", ivl2="<<i2);
        rl=-inf(); ru=+inf();
    } else if(x1>=0) {
        rl=_div_down(x1,i2.u); ru=_div_up(x1,i2.l); 
    } else {
        rl=_div_down(x1,i2.l); ru=_div_up(x1,i2.u); 
    } 
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

Interval sqr(Interval i) 
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double rl,ru;
    if(i.l >=0) {
        rl=_mul_down(i.l,i.l); ru=_mul_up(i.u,i.u);
    } else if(i.u<=0) {
        rl=_mul_down(i.u,i.u); ru=_mul_up(i.l,i.l);
    } else {
        rl=0.0;
        set_rounding_mode(upward);
        ru=max(i.l*i.l,i.u*i.u);
    }
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}



Interval pow(Interval i, int n) 
{
    //std::cerr<<"pow("<<i<<","<<n<<")\n";
    if(n<0) { return pow(rec(i),uint(-n)); }
    else return pow(i,uint(n));
}

Interval pow(Interval i, uint m) 
{
    //std::cerr<<"pow("<<i<<","<<m<<"u)\n";
    const Interval& nvi=i;
    if(m%2==0) { i=abs(nvi); }
    volatile double rl,ru;
    set_rounding_mode(downward);
    rl=pow_rnd(i.l,m);
    set_rounding_mode(upward);
    ru=pow_rnd(i.u,m);
    return Interval(rl,ru);
}



Interval sqrt(Interval i)
{
    rounding_mode_t rnd = get_rounding_mode();
    set_rounding_downward();
    volatile Float rl=sqrt_rnd(i.l);
    set_rounding_upward();
    volatile Float ru=sqrt_rnd(i.u);
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

Interval exp(Interval i)
{
    rounding_mode_t rnd = get_rounding_mode();
    set_rounding_downward();
    volatile Float rl=exp_rnd(i.l);
    set_rounding_upward();
    volatile Float ru=exp_rnd(i.u);
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

Interval log(Interval i)
{
    rounding_mode_t rnd = get_rounding_mode();
    set_rounding_downward();
    volatile Float rl=log_rnd(i.l);
    set_rounding_upward();
    volatile Float ru=log_rnd(i.u);
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}


template<> Interval pi<Interval>()
{
    return Interval(pi_down,pi_up);
}


Interval sin(Interval i)
{
    return cos(i-pi<Interval>()/2);
}

Interval cos(Interval i)
{
    //std::cerr<<"cos("<<i<<")"<<std::endl;
    ARIADNE_ASSERT(i.lower()<=i.upper());

    rounding_mode_t rnd = get_rounding_mode();

    static const Interval pi(pi_down,pi_up);
    if(i.radius()>2*pi_down) { return Interval(-1.0,+1.0); }

    double n=std::floor(i.lower()/(2*pi_approx)+0.5);
    i=i-2*n*pi;

    ARIADNE_ASSERT(i.lower()>=-pi_up);
    ARIADNE_ASSERT(i.lower()<=pi_up);
    
    Float rl,ru;
    if(i.lower()<=0.0) {
        if(i.upper()<=0.0) { rl=cos_down(i.lower()); ru=cos_up(i.upper()); }
        else if(i.upper()<=pi_down) { rl=cos_down(max(-i.lower(),i.upper())); ru=+1.0; }
        else { return Interval(-1.0,+1.0); }
    } else if(i.lower()<=pi_up) {
        if(i.upper()<=pi_down) { rl=cos_down(i.upper()); ru=cos_up(i.lower()); }
        else if(i.upper()<=2*pi_down) { rl=-1.0; ru=cos_up(min(i.lower(),sub_down(2*pi_down,i.upper()))); }
        else { return Interval(-1.0,+1.0); }
    } else {
        assert(false);
    }

    set_rounding_mode(rnd);

    return Interval(rl,ru);
}

Interval tan(Interval i)
{
    ARIADNE_NOT_IMPLEMENTED;
}

Interval asin(Interval i)
{
    ARIADNE_NOT_IMPLEMENTED;
}

Interval acos(Interval i)
{
    ARIADNE_NOT_IMPLEMENTED;
}

Interval atan(Interval i)
{
    ARIADNE_NOT_IMPLEMENTED;
}



#ifdef HAVE_GMPXX_H 

typedef mpq_class Rational;

Rational sqr(const Rational& q) { return q*q; }

Rational pow(const Rational& q, uint n) {
    if(n==0) { return 1; }
    Rational r=1; Rational p=q; uint m=n;
    while(m>=1) { if(m%2) { r*=p; } m/=2; p=p*p; }
    return r;
}

Rational pow(const Rational& q, int n) {
    if(n>=0) { return pow(q,uint(n)); }
    else { return pow(1/q,uint(-n)); }
}


Interval::Interval(Rational q)
    : l(down(q.get_d())), u(up(q.get_d())) { }

Interval::Interval(Rational lower, Rational upper)
    : l(down(lower.get_d())), u(up(upper.get_d())) { }

#endif // HAVE_GMPXX_H 


std::ostream& 
operator<<(std::ostream& os, const Interval& ivl)
{
    return os << '[' << ivl.l << ':' << ivl.u << ']';
}

std::istream& 
operator>>(std::istream& is, Interval& ivl)
{
    char cl,cm,cr;
    is >> cl >> ivl.l >> cm >> ivl.u >> cr;
    ARIADNE_ASSERT(is);
    ARIADNE_ASSERT(cl=='[' || cl=='(');
    ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
    ARIADNE_ASSERT(cr==']' || cr==')');
    return is;
}

} // namespace Ariadne

