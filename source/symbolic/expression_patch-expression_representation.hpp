/***************************************************************************
 *            symbolic/expression_patch.hpp
 *
 *  Copyright  2024  Pieter Collins
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

/*! \file symbolic/expression_patch.hpp
 *  \brief Symbolic expressions in named variables over bounded domains
 */

#ifndef ARIADNE_EXPRESSION_PATCH_HPP
#define ARIADNE_EXPRESSION_PATCH_HPP

#include "symbolic/expression_set.hpp"

#include "function/taylor_function.hpp"

namespace Ariadne {

Vector<RealExpression> call(EffectiveVectorMultivariateFunction const& f, Vector<RealExpression> const & ev) {
    Vector<ElementaryAlgebra<Real>> av=ev;
    Vector<ElementaryAlgebra<Real>> rav=f(av);
    return Vector<RealExpression>(rav.size(),[&rav](SizeType i){return rav[i].extract<RealExpression>();});
}



template<class P, class T> class ExpressionPatch;
typedef VariableInterval<typename IntervalDomainType::UpperBoundType> VariableIntervalDomainType;
typedef VariablesBox<IntervalDomainType> VariablesBoxDomainType;
typedef VectorVariableBox<IntervalDomainType> VectorVariableBoxDomainType;

template<class X, class F> decltype(auto) vapply(F const& f, Vector<X> const& v) {
    using Y = decltype(f(declval<X>()));
    return Vector<Y>(v.size(),[&f,&v](SizeType i){return f(v[i]);});
}

template<class T> class ConstantOrVariable;

template<> class ConstantOrVariable<Real> : public Variant<Constant<Real>,Variable<Real>> {
};

template<> class ConstantOrVariable<Vector<Real>> : public Vector<ConstantOrVariable<Real>> {
  public:
    using Vector<ConstantOrVariable<Real>>::Vector;
    template<class... AS> ConstantOrVariable(AS... as) : ConstantOrVariable(make_list(as...)) { }
  private:
    template<class... AS> List<ConstantOrVariable<Real>> _make_list(AS... as) {
        List<ConstantOrVariable<Real>> lst; this->_append_list(lst,as...); return lst; }
    Void _append_list(List<ConstantOrVariable<Real>>& lst) { }
    template<class A0, class... AS> Void _append_list(List<ConstantOrVariable<Real>>& lst, A0 a0, AS... as) {
        lst.append(a0); this->_append_list(lst,as...); }
};


struct ConstantOrVariablePatch {
  public:
    ConstantOrVariablePatch(Dyadic);
    ConstantOrVariablePatch(RealConstant);
    ConstantOrVariablePatch(VariableIntervalDomainType);
    ConstantOrVariablePatch(VariablesBoxDomainType);
    ConstantOrVariablePatch(VectorVariableBoxDomainType);
};

template<class P> VectorMultivariateFunction<P> coordinates(SizeType n, Array<SizeType> is) {
    if (is.size()==0) { return VectorMultivariateFunction<P>(0,n); }
    return Vector<ScalarMultivariateFunction<P>>(is.size(), [n,&is](SizeType j){return ScalarMultivariateFunction<P>::coordinate(n,is[j]);});
}
template<class P> VectorMultivariateFunction<P> coordinates(SizeType n, Range rng) {
    if (rng.size()==0) { return VectorMultivariateFunction<P>(0,n); }
    return Vector<ScalarMultivariateFunction<P>>(rng.size(), [n,rng](SizeType j){return ScalarMultivariateFunction<P>::coordinate(n,rng[j]);});
}


template<> class SpacePatch<Real> {
    using T=Real;
    RealSpace _spc;
    BoxDomainType _dom;
  public:
    SpacePatch(InitializerList<SpacePatch<T>> lst) : SpacePatch(*lst.begin()) {
        assert (lst.size()!=0); bool first=true;
        for (auto sp : lst) {
            if (first) { first=false; } else { _spc=join(_spc,sp._spc); _dom=product(_dom,sp._dom); } } }
    SpacePatch(VariableIntervalDomainType vivl) : _spc({vivl.variable()}), _dom({vivl.interval()}) { }
    SpacePatch(VectorVariableBoxDomainType vvbx) : _spc(vvbx.variable()), _dom(vvbx.box()) { }
    //SpacePatch(VariablesBoxDomainType);

    DimensionType dimension() const { return this->_spc.dimension(); }
    Space<T> const& space() const { return this->_spc; }
    explicit operator Space<T> const& () const { return this->_spc; }
    BoxDomainType domain() const { return this->_dom; }
};




template<class P> VectorMultivariateFunction<P> combine(VectorMultivariateFunction<P> f1, VectorMultivariateFunction<P> f2) {
//    return make_combined_function(f1,f2);
    auto n1=f1.argument_size(); auto n2=f2.argument_size(); auto n=n1+n2;
    auto p2=coordinates<P>(n,range(n1,n1+n2));
    auto p1=coordinates<P>(n,range(0,n1));
    return join(compose(f1,p1),compose(f2,p2));
}

template<> class ExpressionPatch<ValidatedTag, RealVector> {
    using P=ValidatedTag;
    using T=RealVector;
    using SIG=T(RealVector);
  private: public:
    Map<RealVariable,IntervalDomainType> _vd;
    Vector<RealExpression> _e;
  public:
    BoxRangeType range();
    SizeType size() const { return _e.size(); }

    ExpressionPatch(RealSpacePatch dom, Vector<RealExpression> expr);
    ExpressionPatch(Map<RealVariable,IntervalDomainType> dom, Vector<RealExpression> expr);
    ExpressionPatch(Vector<RealVariable> vars, BoxDomainType fdom, Function<P,SIG> f);
/*
    ExpressionPatch(InitializerList<ExpressionPatch<P,T>> lst) : ExpressionPatch(*lst.begin()) {
        assert (lst.size()!=0);
        bool first=true;
        for (auto ep : lst) {
            if (first) { first=false; }
            else { *this=_join(*this,ep); }
        }
    }
*/
    ExpressionPatch(Dyadic w) : ExpressionPatch(RealConstant(w)) { }
    ExpressionPatch(RealConstant c)
        : _vd(), _e({c}) { }
    ExpressionPatch(VariableIntervalDomainType vivl)
        : ExpressionPatch(RealSpacePatch({vivl}),{vivl.variable()}) { }
    ExpressionPatch(VectorVariableBoxDomainType vvbx)
        : ExpressionPatch(RealSpacePatch({vvbx}),Vector<RealExpression>(Expression<RealVector>(vvbx.variable()))) { }
    ExpressionPatch(RealVariable v)
        : ExpressionPatch(v|IntervalDomainType(-inf,+inf)) { }
    ExpressionPatch(RealVectorVariable vv)
        : ExpressionPatch(vv|BoxDomainType(vv.size(),IntervalDomainType(-inf,+inf))) { }

    friend auto call(FunctionPatch<P,RealVector(RealVector)>, ExpressionPatch<P,RealVector>) -> ExpressionPatch<P,RealVector>;
    friend auto join(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_join(ep1,ep2); }

    friend auto add(ExpressionPatch<P,T>, ExpressionPatch<P,T>) -> ExpressionPatch<P,T>;

    friend OutputStream& operator<<(OutputStream& os, ExpressionPatch<P,T> const& ep) {
        return os << "ExpressionPatch(spc="<<ep._vd<<", e="<<ep._e<<")"; }
  private:
    static ExpressionPatch<P,T> _join(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2);
/*
     {
        Vector<RealVariable> vars=join(ep1._vars,ep2._vars);
        BoxDomainType fdom=product(ep1._fdom,ep2._fdom);
        if (ep2._f.argument_size()==0) {
            Function<P,SIG> fc=join(ep1._f, Function<P,SIG>::constant(ep1._f.argument_size(),ep2._c));
            return ExpressionPatch<P,T>(vars,fdom,fc);
        }
        Function<P,SIG> f=combine(ep1._f,ep2._f);
        return ExpressionPatch<P,T>(vars,fdom,f);
    }
*/
};

auto call(FunctionPatch<ValidatedTag,RealVector(RealVector)>, ExpressionPatch<ValidatedTag,RealVector>)\
    -> ExpressionPatch<ValidatedTag,RealVector>;


template<class P, class... ARGS> FunctionPatch<P,RealVector(ARGS...)>::
FunctionPatch(RealSpacePatch const& spc, ExpressionPatch<P,RES> const& ep) {
    BoxDomainType dom=spc.domain();

    Array<SizeType> coords(ep._vars.size());
    for(SizeType i=0; i!=coords.size(); ++i) {
        coords[i]=spc.space()[ep._vars[i]];
    }
    VectorMultivariateFunction<P> pr=coordinates<P>(spc.dimension(),coords);

    auto swp=ThresholdSweeper<FloatDP>(dp,1e-8);

    ValidatedVectorMultivariateTaylorFunctionModel<FloatDP> proj(dom,pr,swp);
    std::cerr<<"\nep.f[R"<<ep._f.argument_size()<<"->R"<<ep._f.result_size()<<"]="<<ep._f<<"\nproj="<<proj<<"\n\n";

    *this = compose(ep._f,proj);
}

template<class P, class... ARGS> auto FunctionPatch<P,RealVector(ARGS...)>::
operator() (ConstantOrVariable<ARGS...> const& cv) const -> ExpressionPatch<P,RES> {
    static_assert (Same<Tuple<ARGS...>,Tuple<RealVector>>);

    assert(this->argument_size()==cv.size());

    BoxDomainType const& dom=this->domain();
    Map<RealVariable,IntervalDomainType> edom;

    Vector<RealExpression> p(this->argument_size());

    for (SizeType i=0; i!=cv.size(); ++i) {
        IntervalDomainType const& vdom=dom[i];
        if (std::holds_alternative<RealVariable>(cv[i])) {
            RealVariable v=std::get<RealVariable>(cv[i]);
            if (edom.has_key(v)) { edom[v]=intersection(edom[v],vdom); }
            else { edom[v]=vdom; }
            p[i]=v;
        } else {
            RealConstant c=std::get<RealConstant>(cv[i]);
            p[i]=c;
#warning FIXME: Throw exception
//            assert (contains(vdom,c));
        }
    }

    auto f=Function<P,SIG>(*this);
    return ExpressionPatch<P,RES>(edom,f(p));
}

} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_PATCH_HPP */
