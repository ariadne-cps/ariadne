/***************************************************************************
 *            algebra_operations.tcc
 *
 *  Copyright 2011-15  Pieter Collins
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

#include "utility/exceptions.h"
#include "algebra/algebra_operations.h"

#include "algebra/series.h"
#include "function/taylor_series.h"

namespace Ariadne {

struct Factorial {
    Nat _n;
    Factorial(Nat n) : _n(n) { }
    operator Float64Bounds() { Float64Bounds r=1; for(Nat i=1; i<=_n; ++i) { r*=i; } return r; }
    friend Float64Bounds rec(Factorial x) { return rec(Float64Bounds(x)); }
};

template<class A> EnableIfGradedAlgebra<A>
compose(const Series<typename A::NumericType>& x, const A& y)
{
    Nat d=y.degree();

    A w=y - y.value();
    A r=y.create();
    r+=x[d];
    for(Nat n=1; n<=d; ++n) {
        r=r*w;
        r+=x[d-n];
    }
    return r;
}


template<class X> class TaylorSeries;

template<class A> EnableIfNormedAlgebra<A>
_compose(const TaylorSeries<Float64Bounds>& ts, const A& tv, double eps)
{
    //std::cerr<<"_compose(TaylorSeries,A,ErrorTag)\n";
    //std::cerr<<"\n  ts="<<ts<<"\n  tv="<<tv<<"\n";
    Float64Value& vref=const_cast<Float64Value&>(tv.value());
    Float64Value vtmp=vref;
    vref=0;
    A r(tv.argument_size());
    r+=ts[ts.degree()];
    for(Nat i=1; i<=ts.degree(); ++i) {
        //std::cerr<<"    r="<<r<<std::endl;
        r=r*tv;
        r+=ts[ts.degree()-i];
        r.sweep(eps);
    }
    //std::cerr<<"    r="<<r<<std::endl;
    r+=ts.error;
    //std::cerr<<"    r="<<r<<std::endl;
    vref=vtmp;
    return r;
}

template<class A> EnableIfNormedAlgebra<A>
compose(const TaylorSeries<Float64Bounds>& ts, const A& tm)
{
    return _compose(ts,tm,tm.tolerance());
}


// Compose using the Taylor formula directly. The final term is the Taylor series computed
// over the range of the series. This method tends to suffer from blow-up of the
// truncation error
template<class A> EnableIfNormedAlgebra<A>
_compose1(const AnalyticFunction& fn, const A& tm, double eps)
{
    static const Nat DEGREE=18;
    static const double TRUNCATION_ERROR=1e-8;
    Nat d=DEGREE;
    Float64Value c=tm.value();
    Float64Bounds r=tm.range();
    Series<Float64Bounds> centre_series=fn.series(Float64Bounds(c));
    Series<Float64Bounds> range_series=fn.series(r);

    Float64 truncation_error_estimate=mag(range_series[d])*pow(mag(r-c),d);
    if(truncation_error_estimate>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error_estimate
                     <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    A x=tm-c;
    A res(tm.argument_size(),tm.accuracy_ptr());
    res+=range_series[d];
    for(Nat i=0; i!=d; ++i) {
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
_compose2(const AnalyticFunction& fn, const A& tm, double eps)
{
    static const Nat DEGREE=20;
    static const Float64 TRUNCATION_ERROR=1e-8;
    Nat d=DEGREE;
    Float64Value c=tm.value();
    Float64Bounds r=tm.range();
    Series<Float64Bounds> centre_series=fn.series(Float64Bounds(c));
    Series<Float64Bounds> range_series=fn.series(r);

    //std::cerr<<"c="<<c<<" r="<<r<<" r-c="<<r-c<<" e="<<mag(r-c)<<"\n";
    //std::cerr<<"cs[d]="<<centre_series[d]<<" rs[d]="<<range_series[d]<<"\n";
    //std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";
    Float64Error truncation_error=mag(range_series[d]-centre_series[d])*pow(mag(r-c),d);
    //std::cerr<<"te="<<truncation_error<<"\n";
    if(truncation_error.raw()>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    A x=tm-c;
    A res(tm.argument_size(),tm.accuracy_ptr());
    res+=centre_series[d];
    for(Nat i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        res.sweep(eps);
    }
    res+=Float64Bounds(-truncation_error,+truncation_error);
    return res;
}


// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error. This method is better than _compose2 since the truncation error is
// assumed at the ends of the intervals
template<class A> EnableIfNormedAlgebra<A>
_compose3(const AnalyticFunction& fn, const A& tm, Float64 eps)
{
    static const Nat DEGREE=20;
    static const Float64 TRUNCATION_ERROR=1e-8;
    Nat d=DEGREE;
    Float64Value c=tm.value();
    Float64Bounds r=tm.range();
    Series<Float64Bounds> centre_series=fn.series(Float64Bounds(c));
    Series<Float64Bounds> range_series=fn.series(r);

    //std::cerr<<"c="<<c<<" r="<<r<<" r-c="<<r-c<<" e="<<mag(r-c)<<"\n";
    //std::cerr<<"cs[d]="<<centre_series[d]<<" rs[d]="<<range_series[d]<<"\n";
    //std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";
    Float64Bounds se=range_series[d]-centre_series[d];
    Float64Bounds e=r-c;
    Float64Bounds p=pow(e,d-1);
    p=Float64Bounds(-(-p.lower()*(-e.lower())),p.upper()*e.upper());
    //std::cerr<<"se="<<se<<" e="<<e<<" p="<<p<<std::endl;
    // FIXME: Here we assume the dth derivative of f is monotone increasing
    Float64 truncation_error=max(-(-se.lower()*(-p.lower())),se.upper()*p.upper()).raw();
    //std::cerr<<"te="<<truncation_error<<"\n";
    if(truncation_error>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    A x=tm;
    A res(tm.argument_size(),tm.accuracy_ptr());
    res+=centre_series[d];
    for(Nat i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        //res.sweep(eps);
    }
    res+=Float64Bounds(-truncation_error,+truncation_error);
    return res;
}


template<class A> EnableIfNormedAlgebra<A>
_compose(const AnalyticFunction& fn, const A& tm, Float64 eps)
{
    return _compose3(fn,tm,eps);
}



///////////////////////////////////////////////////////////////////////////////

// Algebraic and trancendental functions
//   singleton domain (rec,sqrt,log,tan)
//   unbounded domain (exp,sin,cos)

namespace {
inline Int pow2(Nat k) { return 1<<k; }
inline Int powm1(Nat k) { return (k%2) ? -1 : +1; }
double rec_fac_up(Nat n) { Float64::set_rounding_upward(); double r=1; for(Nat i=1; i<=n; ++i) { r/=i; } return r; }
}

template<class A> EnableIfNormedAlgebra<A>
sqrt(const A& x)
{
    typedef typename A::NumericType X;

    //std::cerr<<"rec(A)\n";
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    auto tol=cast_exact(x.tolerance());
    auto avg=cast_exact(x.average());
    auto rad=cast_exact(x.radius());

    if(avg<=rad) {
        ARIADNE_THROW(DomainException,"log(A x)","x="<<x<<"\n");
    }

    Float64Error eps=mag(rad/avg);
    ARIADNE_DEBUG_ASSERT(eps<1);

    Series<X> sqrt_series=Series<X>::sqrt(X(1));
    Nat d=integer_cast<Int>((log((1-eps)*tol)/log(eps)+1));
    auto trunc_err=pow(eps,d)/cast_positive(1-eps)*mag(sqrt_series[d]);

    A y=x/avg-1;
    A z=x.create();
    z+=sqrt_series[d-1];
    for(Nat i=0; i!=d; ++i) {
        z=sqrt_series[d-i-1] + z * y;
    }
    z+=z.create_ball(trunc_err);
    z*=sqrt(avg);
    return z;
}

template<class A> EnableIfNormedAlgebra<A>
rec(const A& x)
{
    typedef typename A::NumericType X;
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    auto tol=cast_exact(x.tolerance());
    auto avg=cast_exact(x.average());
    auto rad=cast_exact(x.radius());

    if(decide(rad>=abs(avg))) {
        ARIADNE_THROW(DivideByZeroException,"rec(A x)","x="<<x<<", avg="<<avg<<", rad="<<rad<<"\n");
    }

    auto eps=mag(rad/avg);
    ARIADNE_DEBUG_ASSERT(eps<1);

    // Compute the degree and truncation error
    Nat d=integer_cast<Nat>((log((1-eps)*tol)/log(eps))+1);
    auto te=pow(eps,d)/cast_positive(1-eps);

    A y=1-x/avg;
    A z=x.create();
    z+=(d%2?-1:+1);
    for(Nat i=0; i!=d; ++i) {
        z=1 + z * y;
    }

    z+=z.create_ball(te);
    z/=avg;
    return z;
}

template<class A> EnableIfNormedAlgebra<A>
log(const A& x)
{
    typedef typename A::NumericType X;
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    auto tol=cast_exact(x.tolerance());
    auto avg=cast_exact(x.average());
    auto rad=cast_exact(x.radius());

    if(avg<=rad) {
        ARIADNE_THROW(DomainException,"log(A x)","x="<<x<<"\n");
    }

    auto eps=mag(rad/avg);
    ARIADNE_DEBUG_ASSERT(eps<1);

    Nat d=integer_cast<Nat>((log((1-eps)*tol)/log(eps)+1));
    auto trunc_err=pow(eps,d)/cast_positive(1-eps)/d;

    A y=x/avg-X(1);
    A z=x.create();
    z+=X(d%2?-1:+1)/d;
    for(Nat i=1; i!=d; ++i) {
        z=X((d-i)%2?+1:-1)/(d-i) + z * y;
    }
    z=z*y;

    z+=z.create_ball(trunc_err);
    z+=log(numeric_cast<X>(avg));
    return z;
}

// Use special code to utilise exp(ax+b)=exp(x)^a*exp(b)
template<class A> EnableIfNormedAlgebra<A> exp(const A& x)
{
    typedef typename A::NumericType X;

    auto avg=x.average();
    auto rad=x.radius();
    auto tol = cast_exact(x.tolerance());

    // Scale to unit interval
    Nat sfp=0; // A number such that 2^sfp>rad(x.range())
    while(decide(Float64Value(two_exp(sfp))<rad)) { ++sfp; }
    Float64Value sf=two_exp(sfp);
    A y = (x-avg)/sf;
    auto yrad=rad*mag(sf);

    // Find the required degree
    Nat deg = 0;
    auto trunc_err=pow(yrad,0u);
    trunc_err*=2u;
    do {
        ++deg;
        trunc_err=pow(yrad,deg)*mag(rec(Factorial(deg))*(deg+1u)/deg);
    } while(decide(trunc_err>tol));

    A z=x.create_constant(1);
    for(Nat i=0; i!=deg; ++i) {
        z/=(deg-i);
        z=y*z+1;
    }
    z+=z.create_ball(trunc_err);

    // Square r a total of sfp times
    for(Nat i=0; i!=sfp; ++i) {
        z=z*z;
    }

    // Multiply by exp(xv)
    z*=exp(avg);

    return z;
}

// Use special code to utilise sin(x+2pi)=sin(x)
// and that the power series is of the form x*f(x^2)
template<class A> EnableIfNormedAlgebra<A>
sin(const A& x)
{
    typedef typename A::NumericType X;
    Real const& pi=Ariadne::pi;
    // FIXME: Truncation error may be incorrect

    auto tol = cast_exact(x.tolerance());
    auto avg=x.average();
    auto rad=x.radius();
    auto rng=avg.pm(rad);
    Int n=integer_cast<Int>( round(avg/pi) );

    A y=x-(2*n)*X(pi);

    A s=sqr(y);

    // Find the required degree
    Nat deg = 0;
    auto trunc_err=pow(rad,0u)*2u;
    do {
        ++deg;
        trunc_err=pow(rad,deg)*mag(rec(Factorial(deg))*(deg+1u)/deg);
    } while(decide(trunc_err>tol));

    // Compute x(1-y/6+y^2/120-y^3/5040+... = x(1-y/6*(1-y/20*(1-y/42*...)
    A z=x.create_constant(1);
    for(Int i=0; i!=deg; ++i) {
        z/=X(-Int(2*(deg-i)*(2*(deg-i)+1)));
        z*=s;
        z+=1;
    }
    z=y*z;

    z+=z.create_ball(trunc_err);

    return z;
}

// Use special code to utilise sin(x+2pi)=sin(x)
// and that the power series is of the form f(x^2)
template<class A> EnableIfNormedAlgebra<A>
cos(const A& x)
{
    typedef typename A::NumericType X;

    auto tol = cast_exact(x.tolerance());
    auto avg=x.average();
    auto rad=x.radius();

    Float64 two_pi_approx=2*pi_approx();
    Int n=integer_cast<Int>( round(avg/pi) );

    A y=x-(2*n)*X(pi);

    A s=sqr(y);

    // Find the required degree
    Nat deg = 0;
    auto trunc_err=pow(rad,0u)*2u;
    do {
        ++deg;
        trunc_err=pow(rad,deg)*mag(rec(Factorial(deg))*(deg+1)/deg);
    } while(decide(trunc_err>tol));

    // Compute 1-y/2+y^2/24-y^3/720+... = (1-y/2*(1-y/12*(1-y/30*...)
    A z=x.create_constant(1);
    for(Int i=0; i!=deg; ++i) {
        z/=X(-Int(2*(deg-i)*(2*(deg-i)-1)));
        z*=s;
        z+=1;
    }

    z+=z.create_ball(trunc_err);

    return z;
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
    static const Nat DEG=18;
    typedef typename A::NumericType X;
    Float64 xavg = x.average();
    Float64 xrad = x.radius();
    Float64Bounds xrng = xavg + Float64Bounds(-xrad,+xrad);
    return compose(TaylorSeries(DEG,&Series<X>::asin,xavg,xrng),x);
*/
}

template<class A> EnableIfNormedAlgebra<A>
acos(const A& x)
{
    ARIADNE_NOT_IMPLEMENTED;
/*
    static const Nat DEG=18;
    typedef typename A::NumericType X;
    Float64 xavg = x.average();
    Float64 xrad = x.radius();
    Float64Bounds xrng = xavg + Float64Bounds(-xrad,+xrad);
    return compose(TaylorSeries(DEG,&Series<X>::acos,xavg,xrng),x);
*/
}

template<class A> EnableIfNormedAlgebra<A>
atan(const A& x)
{
    ARIADNE_NOT_IMPLEMENTED;
/*
    static const Nat DEG=18;
    typedef typename A::NumericType X;
    Float64 xavg = x.average();
    Float64 xrad = x.radius();
    Float64Bounds xrng = xavg + Float64Bounds(-xrad,+xrad);
    return compose(TaylorSeries(DEG,&Series<X>::atan,xavg,xrng),x);
*/
}


template<class A> int instantiate_transcendental() {
    auto rec_ptr  = (A(*)(A const&)) &rec;
    auto sqrt_ptr = (A(*)(A const&)) &sqrt;
    auto exp_ptr  = (A(*)(A const&)) &exp;
    auto log_ptr  = (A(*)(A const&)) &log;
    auto sin_ptr  = (A(*)(A const&)) &sin;
    auto cos_ptr  = (A(*)(A const&)) &cos;
    auto tan_ptr  = (A(*)(A const&)) &tan;
    auto atan_ptr  = (A(*)(A const&)) &atan;

    typedef std::size_t size_t;
    return (size_t)rec_ptr + (size_t)sqrt_ptr + (size_t)exp_ptr + (size_t)log_ptr
                + (size_t)sin_ptr + (size_t)cos_ptr + (size_t)tan_ptr + (size_t)atan_ptr;
}

template<class A> int instantiate_ordered() {
    auto abs_ptr  = (A(*)(A const&)) &abs;
    auto max_ptr = (A(*)(A const&, A const&)) &max;
    auto min_ptr = (A(*)(A const&, A const&)) &min;
    typedef std::size_t size_t;
    return (size_t)abs_ptr + (size_t)max_ptr + (size_t)min_ptr;
}

} // namespace Ariadne
