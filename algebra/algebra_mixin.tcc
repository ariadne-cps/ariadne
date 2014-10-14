/***************************************************************************
 *            algebra_mixin.tcc
 *
 *  Copyright 2011  Pieter Collins
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

#include "algebra/algebra_interface.h"
#include "algebra/algebra_mixin.h"
#include "algebra/algebra.h"

#include "algebra/series.h"
#include "function/taylor_series.h"

namespace Ariadne {

namespace {
inline ApproximateFloat operator+(ApproximateFloat x1, int n2) { return x1+ApproximateFloat(n2); }
} // namespace

template<class A> EnableIfGradedAlgebra<A>
compose(const Series<typename A::NumericType>& x, const A& y)
{
    uint d=y.degree();

    A w=y - y.value();
    A r=y.create();
    r+=x[d];
    for(uint n=1; n<=d; ++n) {
        r=r*w;
        r+=x[d-n];
    }
    return r;
}


class TaylorSeries;

template<class A> EnableIfNormedAlgebra<A>
_compose(const TaylorSeries& ts, const A& tv, double eps)
{
    //std::cerr<<"_compose(TaylorSeries,A,Error)\n";
    //std::cerr<<"\n  ts="<<ts<<"\n  tv="<<tv<<"\n";
    ExactFloat& vref=const_cast<ExactFloat&>(tv.value());
    ExactFloat vtmp=vref;
    vref=0;
    A r(tv.argument_size());
    r+=ts.expansion[ts.expansion.size()-1];
    for(uint i=1; i!=ts.expansion.size(); ++i) {
        //std::cerr<<"    r="<<r<<std::endl;
        r=r*tv;
        r+=ts.expansion[ts.expansion.size()-i-1];
        r.sweep(eps);
    }
    //std::cerr<<"    r="<<r<<std::endl;
    r+=ts.error;
    //std::cerr<<"    r="<<r<<std::endl;
    vref=vtmp;
    return r;
}

template<class A> EnableIfNormedAlgebra<A>
compose(const TaylorSeries& ts, const A& tm)
{
    return _compose(ts,tm,tm.tolerance());
}


// Compose using the Taylor formula directly. The final term is the Taylor series computed
// over the range of the series. This method tends to suffer from blow-up of the
// truncation error
template<class A> EnableIfNormedAlgebra<A>
_compose1(const ValidatedSeriesFunctionPointer& fn, const A& tm, double eps)
{
    static const uint DEGREE=18;
    static const double TRUNCATION_ERROR=1e-8;
    uint d=DEGREE;
    ExactFloat c=tm.value();
    ValidatedFloat r=tm.range();
    Series<ValidatedFloat> centre_series=fn(d,ValidatedFloat(c));
    Series<ValidatedFloat> range_series=fn(d,r);

    Float truncation_error_estimate=mag(range_series[d])*pow(mag(r-c),d);
    if(truncation_error_estimate>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error_estimate
                     <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    A x=tm-c;
    A res(tm.argument_size(),tm.accuracy_ptr());
    res+=range_series[d];
    for(uint i=0; i!=d; ++i) {
        //std::cerr<<"i="<<i<<" r="<<res<<"\n";
        res=centre_series[d-i-1]+x*res;
        res.sweep(eps);
    }
    //std::cerr<<"i="<<d<<" r="<<res<<"\n";
    return res;
}

// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error. The radius of convergence of this method is still quite low,
// typically only half of the radius of convergence of the power series itself
template<class A> EnableIfNormedAlgebra<A>
_compose2(const ValidatedSeriesFunctionPointer& fn, const A& tm, double eps)
{
    static const uint DEGREE=20;
    static const Float TRUNCATION_ERROR=1e-8;
    uint d=DEGREE;
    ExactFloat c=tm.value();
    ValidatedFloat r=tm.range();
    Series<ValidatedFloat> centre_series=fn(d,ValidatedFloat(c));
    Series<ValidatedFloat> range_series=fn(d,r);

    //std::cerr<<"c="<<c<<" r="<<r<<" r-c="<<r-c<<" e="<<mag(r-c)<<"\n";
    //std::cerr<<"cs[d]="<<centre_series[d]<<" rs[d]="<<range_series[d]<<"\n";
    //std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";
    ErrorFloat truncation_error=mag(range_series[d]-centre_series[d])*pow(mag(r-c),d);
    //std::cerr<<"te="<<truncation_error<<"\n";
    if(truncation_error.raw()>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    A x=tm-c;
    A res(tm.argument_size(),tm.accuracy_ptr());
    res+=centre_series[d];
    for(uint i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        res.sweep(eps);
    }
    res+=ValidatedFloat(-truncation_error,+truncation_error);
    return res;
}


// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error. This method is better than _compose2 since the truncation error is
// assumed at the ends of the intervals
template<class A> EnableIfNormedAlgebra<A>
_compose3(const ValidatedSeriesFunctionPointer& fn, const A& tm, Float eps)
{
    static const uint DEGREE=20;
    static const Float TRUNCATION_ERROR=1e-8;
    uint d=DEGREE;
    ExactFloat c=tm.value();
    ValidatedFloat r=tm.range();
    Series<ValidatedFloat> centre_series=fn(d,ValidatedFloat(c));
    Series<ValidatedFloat> range_series=fn(d,r);

    //std::cerr<<"c="<<c<<" r="<<r<<" r-c="<<r-c<<" e="<<mag(r-c)<<"\n";
    //std::cerr<<"cs[d]="<<centre_series[d]<<" rs[d]="<<range_series[d]<<"\n";
    //std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";
    ValidatedFloat se=range_series[d]-centre_series[d];
    ValidatedFloat e=r-c;
    ValidatedFloat p=pow(e,d-1);
    p=ValidatedFloat(-(-p.lower()*(-e.lower())),p.upper()*e.upper());
    //std::cerr<<"se="<<se<<" e="<<e<<" p="<<p<<std::endl;
    // FIXME: Here we assume the dth derivative of f is monotone increasing
    Float truncation_error=max(-(-se.lower()*(-p.lower())),se.upper()*p.upper()).raw();
    //std::cerr<<"te="<<truncation_error<<"\n";
    if(truncation_error>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    A x=tm;
    A res(tm.argument_size(),tm.accuracy_ptr());
    res+=centre_series[d];
    for(uint i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        //res.sweep(eps);
    }
    res+=ValidatedFloat(-truncation_error,+truncation_error);
    return res;
}


template<class A> EnableIfNormedAlgebra<A>
_compose(const ValidatedSeriesFunctionPointer& fn, const A& tm, Float eps)
{
    return _compose3(fn,tm,eps);
}



///////////////////////////////////////////////////////////////////////////////

// Algebraic and trancendental functions
//   bounded domain (rec,sqrt,log,tan)
//   unbounded domain (exp,sin,cos)

namespace {
inline int pow2(uint k) { return 1<<k; }
inline int powm1(uint k) { return (k%2) ? -1 : +1; }
double rec_fac_up(uint n) { set_rounding_upward(); double r=1.0; for(uint i=1; i<=n; ++i) { r/=i; } return r; }
}

template<class A> EnableIfNormedAlgebra<A>
sqrt(const A& x)
{
    typedef typename A::NumericType X;

    //std::cerr<<"rec(A)\n";
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    auto avg=x.average();
    auto rad=x.radius();
    ValidatedFloat rng=avg+ValidatedFloat(-rad,+rad);

    if(rng.lower()<=0) {
        ARIADNE_THROW(DomainException,"sqrt",x);
    }

    set_rounding_upward();
    UpperFloat eps=rng.radius()/rng.midpoint();
    set_rounding_to_nearest();
    assert(eps<1);
    uint d=integer_cast<int>((log((1-eps)*ExactFloat(x.tolerance()))/log(eps)+1));
    //std::cerr<<"x="<<x<<std::endl;
    //std::cerr<<"x/a="<<x/a<<" a="<<a<<std::endl;
    A y=(x/avg)-X(1.0);
    //std::cerr<<"y="<<y<<std::endl;
    A z=x.create();
    Series<X> sqrt_series=Series<X>::sqrt(d,X(1));
    //std::cerr<<"sqrt_series="<<sqrt_series<<std::endl;
    //std::cerr<<"y="<<y<<std::endl;
    z+=sqrt_series[d-1];
    for(uint i=0; i!=d; ++i) {
        z=sqrt_series[d-i-1] + z * y;
        //std::cerr<<"z="<<z<<std::endl;
    }
    auto trunc_err=pow(eps,d)/(1-eps)*mag(sqrt_series[d]);
    //std::cerr<<"te="<<trunc_err<<" te*[-1,+1]="<<trunc_err*ValidatedFloat(-1,1)<<std::endl;
    z+=z.create_ball(trunc_err);
    //std::cerr<<"z="<<z<<std::endl;
    X sqrta=sqrt(numeric_cast<X>(avg));
    //std::cerr<<"sqrt(a)="<<sqrta<<std::endl;
    z*=sqrt(numeric_cast<X>(avg));
    //std::cerr<<"z="<<z<<std::endl;
    return z;
}
struct Foo { };
template<class A> EnableIfNormedAlgebra<A>
rec(const A& x)
{
    typedef typename A::NumericType X;
    //std::cerr<<"rec(A)\n";
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    ExactFloat avg=x.average();
    ErrorFloat rad=x.radius();
    ValidatedFloat rng=avg+ValidatedFloat(-rad,+rad);
    if(rng.upper()>=0 && rng.lower()<=0) {
        ARIADNE_THROW(DivideByZeroException,"rec(A x)","x="<<x<<", x.range()="<<rng);
    }
    set_rounding_upward();
    ErrorFloat eps=rng.radius()/abs(rng.midpoint());
    set_rounding_to_nearest();
    assert(eps<1);

    uint d=integer_cast<uint>((log((1-eps)*ExactFloat(x.tolerance()))/log(eps))+1);

    A y=1-(x/numeric_cast<X>(avg));
    A z=x.create();
    z+=Float(d%2?-1:+1);
    //std::cerr<<"  y="<<y<<"\n";
    //std::cerr<<"  z="<<z<<"\n";
    for(uint i=0; i!=d; ++i) {
        z=1.0 + z * y;
        //std::cerr<<"  z="<<z<<"\n";
    }

    // Compute the truncation error
    ErrorFloat te=pow(eps,d)/(1-eps);
    z+=z.create_ball(te);
    //std::cerr<<"  z="<<z<<"\n";
    z/=avg;
    //std::cerr<<"  z="<<z<<"\n";
    return z;
}

template<class A> EnableIfNormedAlgebra<A>
log(const A& x)
{
    typedef typename A::NumericType X;
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    auto avg=x.average();
    auto rad=x.radius();
    ValidatedFloat rng=avg+ValidatedFloat(-rad,+rad);
    if(rng.lower()<=0) {
        ARIADNE_THROW(DomainException,"sqrt",rng);
    }
    set_rounding_upward();
    ErrorFloat eps=rng.radius()/avg;
    set_rounding_to_nearest();
    assert(eps<1);
    uint d=integer_cast<uint>((log((1-eps)*ExactFloat(x.tolerance()))/log(eps)+1));
    A y=x/avg-X(1);
    A z=x.create();
    z+=Float(d%2?-1:+1)/d;
    for(uint i=1; i!=d; ++i) {
        z=Float((d-i)%2?+1:-1)/(d-i) + z * y;
    }
    z=z*y;
    ErrorFloat trunc_err=pow(eps,d)/(1-eps)/d;
    return z+log(numeric_cast<X>(avg))+ValidatedFloat(-trunc_err,+trunc_err);
}

// Use special code to utilise exp(ax+b)=exp(x)^a*exp(b)
template<class A> EnableIfNormedAlgebra<A> exp(const A& x)
{
    typedef typename A::NumericType X;
    // FIXME: Truncation error may be incorrect

    // Scale to unit interval
    ExactFloat xavg=x.average();
    A y = x-xavg;

    ErrorFloat xrad=x.radius();
    RawFloat xtol = x.tolerance();

    uint sfp=0; // A number such that 2^sfp>rad(x.range())
    while(Float(1<<sfp)<xrad.raw()) { ++sfp; }
    ExactFloat sf=ExactFloat(RawFloat(1)/(1<<sfp));
    y*=sf;
    ErrorFloat yrad=xrad*sf;

    static const uint degree = 7;
    A res=x.create();
    ErrorFloat truncation_error = (pow(yrad,degree+1));
    res += numeric_cast<X>(ValidatedFloat(-truncation_error,+truncation_error));
    for(uint i=0; i!=degree; ++i) {
        res/=(degree-i);
        res=y*res+1.0;
    }

    // Square r a total of sfp times
    A square=x.create();
    for(uint i=0; i!=sfp; ++i) {
        res=res*res;
        square.clear();
    }

    // Multiply by exp(xv)
    res*=Ariadne::exp(ValidatedFloat(xavg));

    return res;
    //return _compose(&Series<ValidatedFloat>::exp,x,x.sweep_threshold());
}

// Use special code to utilise sin(x+2pi)=sin(x)
// and that the power series is of the form x*f(x^2)
template<class A> EnableIfNormedAlgebra<A>
sin(const A& x)
{
    typedef typename A::NumericType X;
    // FIXME: Truncation error may be incorrect
    A z=x.create();
    A y=z;
    A s=z;
    A r=z;

    Float two_pi_approx=2*pi_approx;
    Float xavg=x.average().value();
    int n=integer_cast<int>(floor(RawFloat(xavg/two_pi_approx + 0.5)));

    Real const& pi=Ariadne::pi;
    auto two_pi=2*pi;
    y=x-n*two_pi;

    s=sqr(y);

    int d=8; // TODO: Change number of terms to be dependent on tolerance
    ErrorFloat srad=s.radius();
    ErrorFloat truncation_error=pow_up(srad.value(),d+1)*rec_fac_up((d+1)*2);

    // Compute x(1-y/6+y^2/120-y^3/5040+... = x(1-y/6*(1-y/20*(1-y/42*...)
    r=1;
    for(int i=0; i!=d; ++i) {
        r/=X(-2*(d-i)*(2*(d-i)+1));
        r*=s;
        r+=X(1.0);
    }
    r=y*r;

    r+=r.create_ball(truncation_error);

    return r;
}

// Use special code to utilise sin(x+2pi)=sin(x)
// and that the power series is of the form f(x^2)
template<class A> EnableIfNormedAlgebra<A>
cos(const A& x)
{
    typedef typename A::NumericType X;

    // FIXME: Truncation error may be incorrect
    A z=x.create();
    A y=z;
    A s=z;
    A r=z;

    Float two_pi_approx=2*pi_approx;
    int n=integer_cast<int>(floor(Float(x.average())/two_pi_approx + 0.5));

    auto two_pi=2*pi;
    y=x-n*two_pi;

    //if(y.error()>two_pi/2) {
    //    r.error()=1.0;
    //} else
    {
        s=sqr(y);

        int d=8; // TODO: Change number of terms to be dependent on tolerance
        ErrorFloat srad=s.radius();
        ErrorFloat truncation_error=pow_up(srad.value(),d+1)*rec_fac_up((d+1)*2);

        // Compute 1-y/2+y^2/24-y^3/720+... = (1-y/2*(1-y/12*(1-y/30*...)
        r=1.0;
        for(int i=0; i!=d; ++i) {
            r/=X(-2*(d-i)*(2*(d-i)-1));
            r*=s;
            r+=1.0;
        }

        r+=r.create_ball(truncation_error);
    }

    return r;
}

template<class A> EnableIfNormedAlgebra<A>
tan(const A& x)
{
    return sin(x)*rec(cos(x));
}

template<class A> EnableIfNormedAlgebra<A>
asin(const A& x)
{
    ARIADNE_NOT_IMPLEMENTED;
/*
    static const uint DEG=18;
    typedef typename A::NumericType X;
    Float xavg = x.average();
    Float xrad = x.radius();
    ValidatedFloat xrng = xavg + ValidatedFloat(-xrad,+xrad);
    return compose(TaylorSeries(DEG,&Series<X>::asin,xavg,xrng),x);
*/
}

template<class A> EnableIfNormedAlgebra<A>
acos(const A& x)
{
    ARIADNE_NOT_IMPLEMENTED;
/*
    static const uint DEG=18;
    typedef typename A::NumericType X;
    Float xavg = x.average();
    Float xrad = x.radius();
    ValidatedFloat xrng = xavg + ValidatedFloat(-xrad,+xrad);
    return compose(TaylorSeries(DEG,&Series<X>::acos,xavg,xrng),x);
*/
}

template<class A> EnableIfNormedAlgebra<A>
atan(const A& x)
{
    ARIADNE_NOT_IMPLEMENTED;
/*
    static const uint DEG=18;
    typedef typename A::NumericType X;
    Float xavg = x.average();
    Float xrad = x.radius();
    ValidatedFloat xrng = xavg + ValidatedFloat(-xrad,+xrad);
    return compose(TaylorSeries(DEG,&Series<X>::atan,xavg,xrng),x);
*/
}



} // namespace Ariadne
