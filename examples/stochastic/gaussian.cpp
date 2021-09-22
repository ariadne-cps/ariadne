/***************************************************************************
 *            gaussian.cpp
 *
 *  Copyright  2021  Pieter Collins
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

#include "ariadne_main.hpp"

#include "numeric/real.hpp"
#include "function/function.hpp"
//#include "function/function_patch.hpp"
#include "geometry/set.hpp"
#include "geometry/valuation.hpp"
#include "geometry/union_of_intervals.hpp"

#include "function/function_mixin.tpl.hpp"


namespace Ariadne {

void foo() {
    using P=ValidatedTag;
    RealVariable x("x");
    Function<P,Real(Real)> f=make_function(x,exp(-x));
    SweeperDP swp;
    auto tx = ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate({IntervalDomainType(-1,+1)},0u,swp);
    auto mx = ValidatedScalarMultivariateFunctionModelDP(tx);
    std::cout<<"tx="<<tx<<"\n";
    std::cout<<compose(f,tx)<<"\n";
    std::cout<<"mx="<<mx<<"\n";
    std::cout<<compose(f,mx)<<"\n";
}



template<class SIG, class PR, class PRE=PR> using ValidatedFunctionModel = FunctionModel<ValidatedTag,SIG,PR,PRE>;

template<class P> class Function<P, Real(Real,Real)> {
    Function<P, Real(RealVector)> _f;
  public:
    Function(EffectiveFunction<Real(Real,Real)> const& f) : _f(f._f) { }
    Function<P,Real(RealVector)> multivariate_version() const { return _f; }
    template<class A> A operator() (A const& a1, A const& a2) const {
        // FIXME
        if constexpr (Same<A,ValidatedScalarMultivariateFunctionModel<DP>>) {
            ValidatedVectorMultivariateFunctionModel<DP> va(2,a1); va[1]=a2; return compose(_f,va);
        } else if constexpr (Same<A,ValidatedScalarMultivariateTaylorFunctionModelDP>) {
            ValidatedVectorMultivariateTaylorFunctionModelDP va(2,a1); va[1]=a2; return compose(_f,va);
        } else {
            return _f(Vector<A>({a1,a2}));
        }
    }
    friend OutputStream& operator<<(OutputStream& os, Function<P,Real(Real,Real)> const& f) { return os << f._f; }
  private:
    Function(RealVariable x, RealVariable y, RealExpression const& e)
        : _f(make_function({x,y},e)) { }
    friend class Function<ValidatedTag,Real(Real,Real)>;
    friend EffectiveFunction<Real(Real,Real)> make_function(RealVariable, RealVariable, RealExpression const&);
};

template<class G> G compose(ValidatedFunction<Real(Real,Real)> const& f, Tuple<G,G> const& g) {
    return f(std::get<0>(g), std::get<1>(g));
}

EffectiveFunction<Real(Real,Real)> make_function(RealVariable x, RealVariable y, RealExpression const& e) {
    return EffectiveFunction<Real(Real,Real)>(x,y,e); }

template<class P> class Function<P, RealValuation(Real)> {
    Function<P, Real(Real,Real)> _f;
  public:
    Function(Function<P,Real(Real,Real)> f)
        : _f(f) { }
    RealValuation operator() (RealValuation mu) const;
  public:
    Function<P,Real(Real,Real)> density() const { return this->_f; }
    template<class FLT> Bounds<FLT> _image_measure(Bounds<FLT> x, UnionOfIntervals<Value<FLT>> V) {
        Function<P,Real(Real)> curried_function; }
    friend OutputStream& operator<<(OutputStream& os, Function<P,RealValuation(Real)> const& f) { return os << f._f; }
};


template<class P, class PR, class PRE> FunctionModel<P,Real(Real),PR,PRE> cast_univariate(FunctionModel<P,Real(RealVector),PR,PRE>);

template<class PR, class PRE> class FunctionModel<ValidatedTag, Real(Real), PR, PRE> {
    using P = ValidatedTag; using SIG=Real(Real);
    FunctionModel<ValidatedTag, Real(RealVector), PR, PRE> _mvf;
  private:
    friend FunctionModel<P,Real(Real),PR,PRE> cast_univariate<P,PR,PRE>(FunctionModel<P,Real(RealVector),PR,PRE>);
    FunctionModel(FunctionModel<ValidatedTag, Real(RealVector), PR, PRE> mvf) : _mvf(mvf) { }
  public:
    auto multivariate_version() const { return _mvf; }
    IntervalDomainType domain() const { return _mvf.domain()[0]; }
    template<class A> decltype(auto) operator() (A const& a) { Vector<A> v({a}); return _mvf(v); }
    friend OutputStream& operator<<(OutputStream& os, FunctionModel<P,SIG,PR,PRE> const& f) { return os << f._mvf; }
  public:
    friend FunctionModel<P,SIG,PR,PRE> antiderivative(FunctionModel<P,SIG,PR,PRE> const& fm) {
        return FunctionModel<P,SIG,PR,PRE>(antiderivative(fm._mvf,0u)); }
};


template<class P, class PR, class PRE> FunctionModel<P,Real(Real),PR,PRE> cast_univariate(FunctionModel<P,Real(RealVector),PR,PRE> f) {
    return FunctionModel<P,Real(Real),PR,PRE>(f); }


TaylorModel<ValidatedTag,FloatDP> my_compose(ValidatedFunctionModel<Real(RealVector),DP> fm, Vector<TaylorModel<ValidatedTag,FloatDP>> vt) {
    return fm.evaluate(vt); }

template<class PR, class PRE, class T> T unchecked_evaluate(ValidatedFunctionModel<Real(Real),PR,PRE> fm, T t) {
    if constexpr (Same<T,ValidatedScalarMultivariateTaylorFunctionModelDP>) {
        return compose(fm,t);
    } else if constexpr (Same<T,TaylorModel<ValidatedTag,FloatDP>>) {
        return my_compose(fm.multivariate_version(),Vector<T>(1u,t));
    } else {
        ARIADNE_NOT_IMPLEMENTED;
    }
}

template<class P, class SIG, class PR, class PRE=PR> class UncheckedFunctionModel
    : public FunctionMixin<UncheckedFunctionModel<P,SIG,PR,PRE>,P,SIG>
{
    using D = typename Function<P,SIG>::DomainType;
    FunctionModel<P,SIG,PR,PRE> _fm;

  public:
    explicit UncheckedFunctionModel(FunctionModel<P,SIG,PR,PRE> fm) : _fm(fm) { }
    virtual SizeOne argument_size() const final override { return SizeOne(); }
    virtual SizeOne result_size() const final override { return SizeOne(); }
    virtual FunctionInterface<P,SIG>* _derivative(ElementIndexType<D>) const final override { ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& _write(OutputStream& os) const final override { return os << this->_fm; }
    template<class T> Void _compute(T& r, T const& a) const { r=unchecked_evaluate(this->_fm,a); }

};
template<class P, class SIG, class PR, class PRE>
UncheckedFunctionModel<P,SIG,PR,PRE> cast_unchecked(FunctionModel<P,SIG,PR,PRE> fm) {
    return UncheckedFunctionModel(fm);
}


template<class M> class ScalarUnivariateScaledFunctionPatch {
    using P=ValidatedTag; using SIG=Real(Real);
    using PR = typename M::PrecisionType; using PRE = typename M::ErrorPrecisionType;
    ScalarScaledFunctionPatch<M> _mvfp;
  public:
    typedef M ModelType;
    operator FunctionModel<P,SIG,PR,PRE> () const;
    ScalarUnivariateScaledFunctionPatch<M>(ScalarScaledFunctionPatch<M> mvfp) : _mvfp(mvfp) { }
    IntervalDomainType domain() const { return _mvfp.domain()[0]; }
    template<class A> decltype(auto) operator() (A const& a) { Vector<A> v({a}); return _mvfp(v); }
    friend OutputStream& operator<<(OutputStream& os, ScalarUnivariateScaledFunctionPatch<M> const& f) { return os << f._mvfp; }
};

template<class P, class F> using ScalarUnivariateTaylorFunctionModel = ScalarUnivariateScaledFunctionPatch<TaylorModel<P,F>>;


template<class P, class T> using Measure = Valuation<P,T>;
template<class T> using EffectiveMeasure = Measure<EffectiveTag,T>;
template<class T> using ValidatedMeasure = Measure<ValidatedTag,T>;
using RealMeasure = Measure<EffectiveTag,Real>;
using EffectiveRealMeasure = Measure<EffectiveTag,Real>;
using ValidatedRealMeasure = Measure<ValidatedTag,Real>;


using EffectivePositiveLowerNumber = PositiveLowerNumber<EffectiveTag>;
using ValidatedPositiveLowerNumber = PositiveLowerNumber<ValidatedTag>;


template<class PR> UnionOfIntervals<FloatLowerBound<PR>> under_approximation(ValidatedOpenSet<Real> const& ops, PR pr) {
    // FIXME: Arbitrarily choose initial interval
    typedef RawFloatType<PR> FLT;

    Dyadic tolerance(1,10u);
    Interval<Value<FLT>> dom(-1,+1,pr);

    List<Interval<LowerBound<FLT>>> ivls;
    List<Interval<Value<FLT>>> pending;
    pending.append(dom);

    while (!pending.empty()) {
        Interval<Value<FLT>> ivl = pending.back(); pending.pop();
        ValidatedKleenean sbst = ops.covers(ivl);
        if (definitely(sbst)) { ivls.append(ivl); }
        else if (ivl.width() < tolerance) {
        } else {
            auto ivlp = ivl.split();
            pending.append(ivlp.first());
            pending.append(ivlp.second());
        }
    }
}

template<class T> class LebesgueMeasure {
  public:
    EffectivePositiveLowerNumber measure(EffectiveOpenSet<T> const& ops) const;
    ValidatedPositiveLowerNumber measure(ValidatedOpenSet<T> const& ops) const;
};

template<> class LebesgueMeasure<Real> {
  public:
    EffectivePositiveLowerNumber measure(EffectiveOpenSet<Real> const& ops) const;
    ValidatedPositiveLowerNumber measure(ValidatedOpenSet<Real> const& ops) const;
    template<class U> Positive<U> measure(Interval<U> const& ivl);
    template<class U> Positive<U> measure(UnionOfIntervals<U> const& ivls);
};

template<class P, class T> class AbsolutelyContinuousMeasure;
template<class T> using EffectiveAbsolutelyContinuousMeasure = AbsolutelyContinuousMeasure<EffectiveTag,T>;
template<class T> using ValidatedAbsolutelyContinuousMeasure = AbsolutelyContinuousMeasure<ValidatedTag,T>;


template<> class AbsolutelyContinuousMeasure<ValidatedTag,Real> {
    using P=ValidatedTag; using T=Real;
    using ARG=Real;
    ValidatedFunction<Real(T)> _f;
  public:
    AbsolutelyContinuousMeasure(ValidatedFunction<Real(T)> f) : _f(f) { };
    AbsolutelyContinuousMeasure(ValidatedFunctionModel<Real(T),DP> fm);
    ValidatedFunction<Real(T)> density() const { return _f; }
    ValidatedPositiveLowerNumber measure(ValidatedOpenSet<T> const& ops) const;
    PositiveFloatDPLowerBound measure(UnionOfIntervals<FloatDPLowerBound> const& ivls) {
        return this->_measure<FloatDPLowerBound>(ivls); }
//    PositiveFloatMPLowerBound measure(UnionOfIntervals<FloatMPLowerBound> const& ivls) {
//        return this->_measure<FloatMPLowerBound>(ivls); }
    friend OutputStream& operator<<(OutputStream& os, AbsolutelyContinuousMeasure<P,T> const& mu) {
        return os << "AbsolutelyContinuousMeasure(density="<<mu._f<<")"; }
  private:
    template<class U> Positive<U> _measure(UnionOfIntervals<U> const& ivls) const;

};

template<class U> Positive<U> AbsolutelyContinuousMeasure<ValidatedTag,Real>::_measure(UnionOfIntervals<U> const& ivls) const {
    assert(ivls.size()!=0);
    Positive<U> r=nul(ivls.begin()->width());
    for (auto ivl : ivls) { r=r+cast_positive(integral(_f,ivl)); }
    return r;
}

template<class U> Positive<U> LebesgueMeasure<Real>::measure(Interval<U> const& ivl) {
    return ivl.width();
}

template<class U> Positive<U> LebesgueMeasure<Real>::measure(UnionOfIntervals<U> const& ivls) {
    assert(ivls.size()!=0);
    Positive<U> r=nul(ivls.begin()->width());
    for (auto ivl : ivls) { r=r+measure(ivl); }
    return r;
}


template<class P, class M> ScaledFunctionPatch<M> compose(const ScalarUnivariateFunction<P>& g, const ScaledFunctionPatch<M>& f) {
    return ScaledFunctionPatch<M>(f.domain(),g.evaluate(f.model()));
}




template<> class Function<ValidatedTag,ValidatedRealMeasure(Real)> {
};

template<class FLT> Bounds<FLT> integral(ValidatedFunction<Real(RealVector)> const& f, Box<Interval<Value<FLT>>> const& dom) {
    using AT = Vector<Bounds<FLT>>;
    using ScalarFunctionModelType = ValidatedScalarMultivariateTaylorFunctionModel<FLT>;
    using VectorFunctionModelType = ValidatedVectorMultivariateTaylorFunctionModel<FLT>;

    Sweeper<FLT> swp;
    assert(dom.dimension()==2);
    VectorFunctionModelType id = VectorFunctionModelType::identity(dom,swp);
    ScalarFunctionModelType fm = compose(f,id);

    auto adfm1 = antiderivative(fm,1u);
    auto rfm = partial_evaluate(adfm1,1u,dom[1].upper_bound())-partial_evaluate(adfm1,1u,dom[1].lower_bound());
    auto adfm0 = antiderivative(rfm,0u);
    return evaluate(adfm0,AT({dom[0].upper_bound()}))-evaluate(adfm0,AT({dom[0].lower_bound()}));


    auto afm0 = antiderivative(fm,0u);
    auto afm1 = antiderivative(afm0,1u);
    auto& afm = afm1;
    return afm(AT{dom[0].upper_bound(),dom[1].upper_bound()}) - afm(AT{dom[0].upper_bound(),dom[1].lower_bound()})
         - afm(AT{dom[0].lower_bound(),dom[1].upper_bound()}) + afm(AT{dom[0].lower_bound(),dom[1].lower_bound()}) ;

}


template<class FLT> Bounds<FLT> integral(ValidatedFunctionModel<Real(Real),typename FLT::PrecisionType> const& fm, Interval<Value<FLT>> const& dom) {
    assert(subset(dom,fm.domain()));
    Bounds<FLT> lb=dom.lower_bound();
    Bounds<FLT> ub=dom.upper_bound();

    auto afm = antiderivative(fm);
    return afm(ub)-afm(lb);
}

template<class FLT> Bounds<FLT> integral(ValidatedFunction<Real(Real)> const& f, Interval<Value<FLT>> const& dom) {
    using ScalarFunctionModelType = ValidatedScalarMultivariateTaylorFunctionModel<FLT>;

    Sweeper<FLT> swp;
    ScalarFunctionModelType id = ScalarFunctionModelType::coordinate(Box<Interval<Value<FLT>>>({dom}),0u,swp);
    ScalarFunctionModelType fm = compose(f,id);
    Vector<Bounds<FLT>> lb={dom.lower_bound()};
    Vector<Bounds<FLT>> ub={dom.upper_bound()};

    auto afm = antiderivative(fm,0u);
    return afm(ub)-afm(lb);
}

template<class FLT> LowerBound<FLT> integral(ValidatedFunction<Real(Real)> const& f, Interval<LowerBound<FLT>> const& ivl) {
    return integral(f,cast_exact(ivl));
}

template<class FLT> UpperBound<FLT> integral(ValidatedFunction<Real(Real)> const& f, Interval<UpperBound<FLT>> const& ivl) {
    return integral(f,cast_exact(ivl));
}

// Compute function g(x) = integral(f(x,y),y=a to b);
template<class FLT> ValidatedFunctionModel<Real(RealVector),typename FLT::PrecisionType> integral(IntervalDomainType const& dom, ValidatedFunction<Real(Real,Real)> const& f, Interval<Value<FLT>> const& rng) {
    using PR = typename FLT::PrecisionType;
    using ScalarTaylorFunctionModelType = ValidatedScalarMultivariateTaylorFunctionModel<FLT>;
    using ScalarFunctionModelType = ValidatedScalarMultivariateFunctionModel<PR>;

    Sweeper<FLT> swp;
    ScalarFunctionModelType x = ScalarTaylorFunctionModelType::coordinate(Box<Interval<Value<FLT>>>({dom,rng}),0u,swp);
    ScalarFunctionModelType y = ScalarTaylorFunctionModelType::coordinate(Box<Interval<Value<FLT>>>({dom,rng}),1u,swp);
    ScalarFunctionModelType fm = f(x,y);
    Bounds<FLT> lb=rng.lower_bound();
    Bounds<FLT> ub=rng.upper_bound();

    std::cerr<<"fm="<<fm<<"\n";
    ScalarFunctionModelType afm = antiderivative(fm,1u);
    std::cerr<<"afm="<<afm<<"\n";
    return ValidatedScalarMultivariateFunctionModel<PR>(partial_evaluate(afm,1u,ub)-partial_evaluate(afm,1u,lb));
}

template Bounds<FloatDP> integral(ValidatedFunctionModel<Real(Real),DP> const&, Interval<FloatDPValue> const& rng);


template ValidatedFunctionModel<Real(RealVector),DP> integral(IntervalDomainType const& dom, ValidatedFunction<Real(Real,Real)> const&, Interval<FloatDPValue> const& rng);

template<class PR> ValidatedFunctionModel<Real(RealVector),PR> integral(ValidatedFunctionModel<Real(RealVector),PR> const& f, SizeType j, Interval<Value<RawFloatType<PR>>> const& rng) {
    using FLT=RawFloatType<PR>;
    auto dom = f.domain();
    ARIADNE_PRECONDITION(subset(rng,dom[j]));

    assert(dom.dimension()==2);
    auto adf = antiderivative(f,j);
    Bounds<FLT> lb=rng.lower_bound();
    Bounds<FLT> ub=rng.upper_bound();
    return partial_evaluate(adf,j,ub)-partial_evaluate(adf,j,lb);
}

template ValidatedFunctionModel<Real(RealVector),DP> integral(ValidatedFunctionModel<Real(RealVector),DP> const& f, SizeType j, Interval<FloatDPValue> const& rng);

// Given a density p on R, and a function taking x to a density f_x on R, then the push-forward is given
// by the density y -> integral(p(x)*f(x,y) dx)
ValidatedAbsolutelyContinuousMeasure<Real> push_forward(ValidatedAbsolutelyContinuousMeasure<Real> const& P, ValidatedFunction<RealValuation(Real)> const& F) {
    using FLT = FloatDP; using PR = FLT::PrecisionType;
    std::cerr<<"  P="<<P<<", F="<<F<<"\n";

    //auto suppx = P.support();
    auto suppx = IntervalDomainType(0,5);
    auto p=P.density();
    auto f=F.density();

    // FIXME: Specify support of F(P)
    auto domy = IntervalDomainType(-1,+1);

    using ScalarTaylorFunctionModelType = ValidatedScalarMultivariateTaylorFunctionModel<FLT>;
    using ScalarFunctionModelType = ValidatedScalarMultivariateFunctionModel<PR>;

    Sweeper<FLT> swp=ThresholdSweeper<FLT>(dp,1e-9);
    ScalarFunctionModelType x = ScalarTaylorFunctionModelType::coordinate({suppx,domy},0u,swp);
    ScalarFunctionModelType y = ScalarTaylorFunctionModelType::coordinate({suppx,domy},1u,swp);
    Tuple<ScalarFunctionModelType,ScalarFunctionModelType> xy={x,y};
    std::cerr<<"  x ="<<x<<"\n";
    std::cerr<<"  p ="<<p<<"\n";

    auto tp=compose(p,x);
    ScalarFunctionModelType h = compose(p,x) * compose(f,xy);

    auto q=integral(h, 0u, suppx);
    return ValidatedAbsolutelyContinuousMeasure<Real>(cast_unchecked(cast_univariate(q)));
}


/*
Interval<FloatDPValue> support(ValidatedMeasure<Real> P, FloatDPValue err) {
    FloatDPValue target = cast_exact(1-FloatDPLowerBound(err));
    FloatDPLowerBound measure(0,dp);
    Int i=0;
    while (not definitely(measure > target)) {
        ValidatedOpenSet<Real> pivl=Interval<FloatDPLowerBound>(+i,+i+1,dp);
        ValidatedOpenSet<Real> nivl=Interval<FloatDPLowerBound>(-i-1,-i,dp);
        measure += P.measure(pivl) + P.measure(nivl);
        ++i;
    }
    return Interval<FloatDPValue>(-i,+i,dp);
};
*/

template<class P> class PushForwardMeasure {
    EffectiveMeasure<Real> _P;
    EffectiveFunction<RealValuation(Real)> _F;
  public:
};

//static_assert(AValuation<PushForwardMeasure<EffectiveTag>,EffectiveTag,Real>);
//EffectiveMeasure<Real> push_forward(EffectiveMeasure<Real> const& P, ValidatedFunction<RealValuation(Real)> const& F);

} // namespace Ariadne


void ariadne_main() {
    RealVariable x("x");
    RealVariable y("y");

    EffectiveScalarUnivariateFunction g=make_function(x,exp(-sqr(x)/2));

    Interval<Value<FloatDP>> xivl(0,1,dp);
    Interval<LowerBound<FloatDP>> livl(0,1,dp);
    Interval<UpperBound<FloatDP>> uivl(0,1,dp);

    Bounds<FloatDP> bp=integral(g,xivl);
    auto lp=integral(g,livl);
    auto up=integral(g,uivl);

    Bounds<FloatDP>::set_output_places(16);
    std::cout << "bp=" << bp << "; lp=" << lp << "; up=" << up << "\n";

    ValidatedFunction<Real(Real,Real)> f=make_function(x,y,  exp(-(x+2)*sqr(y-x/2)/2)*sqrt((x+2)/(2*pi)) );
    ValidatedFunction<RealValuation(Real)> F(f);
    std::cout << "F=" << F << "\n";

    Interval<FloatDPValue> V(0,1,dp);
    auto PfV = integral(IntervalDomainType(-1,+1), f,V);
    std::cout << "PfV=" << PfV << "\n";

    ValidatedFunction<Real(Real)> p=make_function(x,exp(-x));
    ValidatedAbsolutelyContinuousMeasure<Real> P(p);
    std::cout << "P=" << P << "\n";


    Sweeper<FloatDP> swp = ThresholdSweeper<FloatDP>(dp,1e-8);

    FloatDPLowerBound lprob (0,dp);
    for (int i=0; i!=8; ++i) {
        IntervalDomainType xdom(i,i+1);
        IntervalDomainType ydom(-1,0);
        BoxDomainType xydom={xdom,ydom};
        auto tid=ValidatedVectorMultivariateTaylorFunctionModelDP::identity(xydom,swp);
        auto tx=ValidatedScalarMultivariateTaylorFunctionModelDP(tid[0]);
        FloatDPBounds prob = integral(compose(p,tx)*compose(f.multivariate_version(),tid), xydom);
        lprob += prob;
        std::cout << "prob=" << prob << ", lprob=" << lprob << "\n";
    }
    std::cout << "lprob=" << lprob << "\n";




//    ValidatedAbsolutelyContinuousMeasure<Real> Q=push_forward(P,F);
//    std::cout << "Q=" << Q << "\n";



/*
    UnionOfIntervals<FloatDPLowerBound> U({{-2,-1,dp}});

    ValidatedMeasure<Real> Pr=P;
    FloatDPValue err(Dyadic(1,10u),dp);
    Interval<FloatDPLowerBound> supp = support(P,err);
    std::cerr<<"supp="<<supp<<"\n";
*/
}


