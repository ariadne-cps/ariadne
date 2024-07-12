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

template<class IVL> class VectorVariableBox {
    //Variable<RealVector> _v;
    Variable<RealVector> _v;
    Box<IVL> _bx;
  public:
    VectorVariableBox(Variable<RealVector> v, Box<IVL> bx) : _v(v.name(),v.size()), _bx(bx) {
        assert(v.size()==bx.dimension()); }
    //VectorVariableBox(Variables<Real> v, Box<IVL> bx) : _v(v), _bx(bx) { }
    Variable<RealVector> variable() const { return this->_v; }
    DimensionType dimension() const { return _bx.dimension(); }
    Box<IVL> box() const { return this->_bx; }
    decltype(auto) operator[](SizeType i) const { return VariableInterval(_v[i],_bx[i]); }
};

template<class UB> VariableInterval<UB> operator|(RealVariable v, Interval<UB> ivl) { return VariableInterval<UB>(v,ivl); }
//template<class UB> VariablesBox<UB> operator|(RealVariables v, Box<UB> ivl);
template<class IVL> VectorVariableBox<IVL> operator|(RealVectorVariable v, Box<IVL> bx) { return VectorVariableBox<IVL>(v,bx); }


template<class P, class T> class ExpressionPatch;
typedef VariableInterval<typename IntervalDomainType::UpperBoundType> VariableIntervalDomainType;
typedef VariablesBox<IntervalDomainType> VariablesBoxDomainType;
typedef VectorVariableBox<IntervalDomainType> VectorVariableBoxDomainType;

template<class X, class F> decltype(auto) vapply(F const& f, Vector<X> const& v) {
    using Y = decltype(f(declval<X>()));
    return Vector<Y>(v.size(),[&f,&v](SizeType i){return f(v[i]);});
}

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
    Vector<RealVariable> _vars;
    BoxDomainType _fdom;
    Function<P,RealVector(RealVector)> _f;
    // FIXME: Avoid special handling of scalar constant
    Vector<Real> _c;
  public:
    static ExpressionPatch<P,T> _join(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) {
        Vector<RealVariable> vars=join(ep1._vars,ep2._vars);
        BoxDomainType fdom=product(ep1._fdom,ep2._fdom);
        if (ep2._f.argument_size()==0) {
            Function<P,SIG> fc=join(ep1._f, Function<P,SIG>::constant(ep1._f.argument_size(),ep2._c));
            return ExpressionPatch<P,T>(vars,fdom,fc);
        }
        Function<P,SIG> f=combine(ep1._f,ep2._f);
        return ExpressionPatch<P,T>(vars,fdom,f);
    }
  public:
    BoxRangeType range();
    SizeType size() const { return _f.result_size(); }

    ExpressionPatch(Map<RealVariable,IntervalDomainType> dom, Vector<RealExpression> expr);
    ExpressionPatch(Vector<RealVariable> vars, BoxDomainType fdom, Function<P,SIG> f)
        : _vars(vars), _fdom(fdom), _f(f)
    {
        assert(fdom.dimension()==f.argument_size());
        assert(f.argument_size()==vars.size());
    }

    ExpressionPatch(InitializerList<ExpressionPatch<P,T>> lst) : ExpressionPatch(*lst.begin()) {
        assert (lst.size()!=0);
        bool first=true;
        for (auto ep : lst) {
            if (first) { first=false; }
            else { *this=_join(*this,ep); }
        }
    }
    ExpressionPatch(Dyadic w) : ExpressionPatch(RealConstant(w)) { }
    ExpressionPatch(RealConstant c)
        : _vars(), _fdom(), _f(Function<P,SIG>::constant(0u,RealVector({c}))), _c({c}) { }
    ExpressionPatch(VariableIntervalDomainType vivl)
        : ExpressionPatch({vivl.variable()},{vivl.interval()},Function<P,SIG>::coordinates(1)) { }
    ExpressionPatch(VectorVariableBoxDomainType vvbx)
        : ExpressionPatch(Vector<RealVariable>(RealVariables(vvbx.variable())), vvbx.box(), Function<P,SIG>::identity(vvbx.dimension())) { }
    ExpressionPatch(RealVariable v)
        : ExpressionPatch({v},{IntervalDomainType(-inf,+inf)},Function<P,SIG>::coordinates(1)) { }
    ExpressionPatch(RealVectorVariable vv)
        : ExpressionPatch(Vector<RealVariable>(RealVariables(vv)), Box(vv.size(),IntervalDomainType(-inf,+inf)), Function<P,SIG>::identity(vv.size())) { }

    friend auto call(FunctionPatch<P,RealVector(RealVector)>, ExpressionPatch<P,RealVector>) -> ExpressionPatch<P,RealVector>;
    friend auto join(ExpressionPatch<P,T>, ExpressionPatch<P,T>) -> ExpressionPatch<P,T>;

    friend OutputStream& operator<<(OutputStream& os, ExpressionPatch<P,T> const& ep) {
        return os << "ExpressionPatch(vars="<<ep._vars<<", fdom="<<ep._fdom<<", fn={R"<<ep._f.argument_size()<<"->R"<<ep._f.result_size()<<"}"<<ep._f<<")"; }
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
operator() (ExpressionPatch<P,ARGS...> const& ep) const -> ExpressionPatch<P,RES> {
    assert(this->argument_size()==ep.size());
    auto f=Function<P,SIG>(*this);
    Function<P,SIG> g=ep._f;
    return ExpressionPatch<P,RES>(ep._vars,ep._fdom,compose(f,g));
}

;

} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_PATCH_HPP */
