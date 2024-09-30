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

#include <typeinfo>

#include "symbolic/variable.hpp"
#include "symbolic/expression.hpp"
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
    using Variant<Constant<Real>,Variable<Real>>::Variant;
};

using RealVectorConstant = Constant<Vector<Real>>;


template<class T> concept IsConstantOrVariable = OneOf<T,Dyadic,RealConstant,RealVariable,RealVectorConstant,RealVectorVariable>;

template<class... FS> struct AllConstantOrVariable;
template<> struct AllConstantOrVariable<> : True { };
template<class F0, class... FS> struct AllConstantOrVariable<F0,FS...> {
    static const bool value = IsConstantOrVariable<F0> && AllConstantOrVariable<FS...>::value;
};

template<class... FS> concept AllAreConstantOrVariable = AllConstantOrVariable<FS...>::value;



template<> class ConstantOrVariable<Vector<Real>> : public Vector<ConstantOrVariable<Real>> {
  public:
    using Vector<ConstantOrVariable<Real>>::Vector;
    template<class... AS> requires AllAreConstantOrVariable<AS...>
        ConstantOrVariable(AS... as) : ConstantOrVariable(_make_list(as...)) { }
  private:
    ConstantOrVariable(Constant<Vector<Real>>);
    ConstantOrVariable(Variable<Vector<Real>>);

    template<class... AS> static List<ConstantOrVariable<Real>> _make_list(AS... as) {
        List<ConstantOrVariable<Real>> lst; _append_list(lst,as...); return lst; }
    static Void _append(List<ConstantOrVariable<Real>>& lst, Dyadic a) {
        lst.append(Constant<Real>(a)); }
    static Void _append(List<ConstantOrVariable<Real>>& lst, Constant<Real> a) {
        lst.append(a); }
    static Void _append(List<ConstantOrVariable<Real>>& lst, Variable<Real> a) {
        lst.append(a); }
    static Void _append(List<ConstantOrVariable<Real>>& lst, ConstantOrVariable<Real> a) {
        lst.append(a); }
    static Void _append(List<ConstantOrVariable<Real>>& lst, Variable<Vector<Real>> a) {
        for (SizeType i=0; i!=a.size(); ++i) { lst.append(a[i]); } }
    static Void _append(List<ConstantOrVariable<Real>>& lst, Constant<Vector<Real>> a) {
        for (SizeType i=0; i!=a.size(); ++i) { lst.append(Constant<Real>(a[i])); } }
    static Void _append(List<ConstantOrVariable<Real>>& lst, ConstantOrVariable<Vector<Real>> a) {
        for (SizeType i=0; i!=a.size(); ++i) { lst.append(a[i]); } }
    static Void _append_list(List<ConstantOrVariable<Real>>& lst) { }
    template<class A0, class... AS> static Void _append_list(List<ConstantOrVariable<Real>>& lst, A0 a0, AS... as) {
        _append(lst,a0); _append_list(lst,as...); }

    template<class... AS> friend List<ConstantOrVariable<Real>> make_constant_or_variable_list(ConstantOrVariable<AS> const& ... cvs);
};

template<class... AS> List<ConstantOrVariable<Real>> make_constant_or_variable_list(ConstantOrVariable<AS> const& ... cvs) {
    return ConstantOrVariable<Vector<Real>>::_make_list(cvs...);
}




template<class T> class ConstantOrVariablePatch;

template<> class ConstantOrVariablePatch<Real> : public Variant<RealConstant,VariableIntervalDomainType> { };

template<> class ConstantOrVariablePatch<Vector<Real>>
    : public Variant<RealVectorConstant,VectorVariableBoxDomainType> { };



class AnyConstantOrVariablePatch
    : public Variant<RealConstant,RealVectorConstant,VariableIntervalDomainType,VectorVariableBoxDomainType>
{
    typedef Variant<RealConstant,RealVectorConstant,VariableIntervalDomainType,VectorVariableBoxDomainType> Base;
  public:
    AnyConstantOrVariablePatch(Dyadic h) : Base(RealConstant(h)) { }
    AnyConstantOrVariablePatch(RealConstant c) : Base(c) { }
    AnyConstantOrVariablePatch(VariableIntervalDomainType vivl) : Base(vivl) { }
    AnyConstantOrVariablePatch(VariablesBoxDomainType vsbx);
    AnyConstantOrVariablePatch(VectorVariableBoxDomainType vvbx) : Base(vvbx) { }
};

template<class... AS> List<ConstantOrVariablePatch<Real>> make_constant_or_variable_patch_list(InitializerList<AnyConstantOrVariablePatch> const& cvplst) {
    auto _append = [](List<ConstantOrVariablePatch<Real>>& lst, auto const& t) {
        using T = decltype(t);
        if constexpr (OneOf<T,Constant<RealVector>,Variable<RealVector>,Vector<RealVariable>>) {
            for (SizeType i=0; i!=t.size(); ++i) { lst.append(t[i]); }
        } else {
            lst.append(t);
        }
    };

    List<ConstantOrVariablePatch<Real>> lst;
    for(auto cvp : cvplst) { std::visit(cvp,[&](AnyConstantOrVariablePatch const& cvp){_append(lst,cvp);}); }
    return lst;
}




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
    SpacePatch(List<RealVariable> vars, BoxDomainType dom)
        : _spc(vars), _dom(dom) { }
    SpacePatch(List<RealVariable> vars, Map<RealVariable,IntervalDomainType> vivls)
        : _spc(vars), _dom(vars.size(),[&](SizeType i){return vivls[vars[i]];}) { }
    //SpacePatch(VariablesBoxDomainType);

    DimensionType dimension() const { return this->_spc.dimension(); }
    Space<T> const& space() const { return this->_spc; }
    List<RealVariable> variables() const { return this->_spc.variables(); }
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

template<> class ExpressionPatch<ValidatedTag, RealScalar> {
    using P=ValidatedTag;
    using T=RealScalar;
    using SIG=T(RealVector);
  private: public:
    Vector<RealVariable> _vars;
    FunctionPatch<P,RealScalar(RealVector)> _fp;
  public:
    IntervalRangeType range();
    SizeOne size() const { return _fp.result_size(); }

    ExpressionPatch(Map<RealVariable,IntervalDomainType> dom, RealExpression expr);

    ExpressionPatch(RealConstant);
    ExpressionPatch(VariableIntervalDomainType);

    friend auto add(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_add(ep1,ep2); }

    friend OutputStream& operator<<(OutputStream& os, ExpressionPatch<P,T> const& ep) {
        os << "ExpressionPatch({";
        for (SizeType i=0; i!=ep._vars.size(); ++i) {
            if (i!=0) { os << ","; } os << ep._vars[i]<<"|"<<ep._fp.domain()[i];
        }
        return os << "}, fn=" << cast_unchecked(ep._fp) << ")";
    }
  private:
    static ExpressionPatch<P,T> _add(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) {
#warning
        ARIADNE_NOT_IMPLEMENTED;
    }
};

template<> class ExpressionPatch<ValidatedTag, RealVector> {
    using P=ValidatedTag;
    using T=RealVector;
    using SIG=T(RealVector);
  private: public:
    Vector<RealVariable> _vars;
    FunctionPatch<P,RealVector(RealVector)> _fp;
  public:
    Map<RealVariable,IntervalDomainType> domains() const;

    BoxRangeType range();
    SizeType size() const { return _fp.result_size(); }

    ExpressionPatch(RealSpacePatch dom, Vector<RealExpression> expr);
    ExpressionPatch(RealSpacePatch dom, Function<P,SIG> f);
    ExpressionPatch(RealSpace spc, FunctionPatch<P,SIG> fp);
    ExpressionPatch(Map<RealVariable,IntervalDomainType> dom, Vector<RealExpression> expr);
    ExpressionPatch(Vector<RealVariable> vars, FunctionPatch<P,SIG> f);
    ExpressionPatch(Vector<RealVariable> vars, BoxDomainType fdom, Function<P,SIG> f);

    typedef Variant<ExpressionPatch<P,RealVector>,ExpressionPatch<P,RealScalar>,VariableIntervalDomainType,RealConstant,Dyadic> PreExpressionPatch;

    ExpressionPatch(PreExpressionPatch var);
    ExpressionPatch(InitializerList<PreExpressionPatch> lst);
    ExpressionPatch(InitializerList<ExpressionPatch<ValidatedTag,RealScalar>> lst);

    friend auto call(FunctionPatch<P,RealVector(RealVector)>, ExpressionPatch<P,RealVector>) -> ExpressionPatch<P,RealVector>;

    friend auto join(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_join(ep1,ep2); }

    friend auto add(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_add(ep1,ep2); }

    friend OutputStream& operator<<(OutputStream& os, ExpressionPatch<P,T> const& ep) {
        return ep._write(os); }

    friend ExpressionPatch<ValidatedTag,RealVector> FunctionPatch<ValidatedTag,RealVector(RealVector)>::operator() (ExpressionPatch<ValidatedTag,RealVector> const& ep) const;
  private:

    Void _adjoin(ExpressionPatch<P,RealScalar> const& sep);
    static ExpressionPatch<P,T> _join(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2);
    static ExpressionPatch<P,T> _add(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2);

    OutputStream& _write(OutputStream& os) const;
};

//typedef Variant<ExpressionPatch<P,RealVector>,ExpressionPatch<P,RealScalar>,VariableIntervalDomainType,RealConstant,Dyadic> PreExpressionPatch;

ExpressionPatch<ValidatedTag,RealScalar>::ExpressionPatch(RealConstant var)
    : _vars(), _fp(BoxDomainType({}), Function<ValidatedTag,RealScalar(RealVector)>::constant(0u,var)) { }

ExpressionPatch<ValidatedTag,RealScalar>::ExpressionPatch(VariableIntervalDomainType vardom)
    : _vars(vardom.variable()), _fp(BoxDomainType({vardom.interval()}), Function<ValidatedTag,RealScalar(RealVector)>::coordinate(1u,0u)) { }

ExpressionPatch<ValidatedTag,RealVector> make_expression_patch(ExpressionPatch<ValidatedTag,RealVector>::PreExpressionPatch var) {
    auto maker = [](auto cve) {
        using CVE=decltype(cve);
        std::cerr<<"  cve="<<cve<<"\n";
        if constexpr (Same<CVE, Dyadic>) {
            auto sep=ExpressionPatch<ValidatedTag,RealScalar>(RealConstant(cve)); return ExpressionPatch<ValidatedTag,RealVector>({sep}); }
        else if constexpr (OneOf<CVE, RealConstant,VariableIntervalDomainType>) {
            ExpressionPatch<ValidatedTag,RealScalar> sep(cve); return ExpressionPatch<ValidatedTag,RealVector>({sep}); }
        else if constexpr (Same<CVE,ExpressionPatch<ValidatedTag,RealScalar>>) {
            return ExpressionPatch<ValidatedTag,RealVector>(cve._vars,VectorMultivariateFunctionPatch<ValidatedTag>(List<ScalarMultivariateFunctionPatch<ValidatedTag>>({cve._fp})));
        } else {
            return cve;
        }
    };

    return std::visit( maker, var);
}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(PreExpressionPatch var)
    : ExpressionPatch(make_expression_patch(var)) { }


ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(RealSpacePatch spcdom, Vector<RealExpression> expr)
    : ExpressionPatch(spcdom,Function<P,SIG>(spcdom.variables(),expr)) {
}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(RealSpacePatch spcdom, Function<P,SIG> f)
    : _vars(spcdom.variables()), _fp(spcdom.domain(),f)
{
    ARIADNE_PRECONDITION(spcdom.dimension()==f.argument_size());
}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(RealSpace spc, FunctionPatch<P,SIG> fp)
    : _vars(spc.variables()), _fp(fp)
{
    ARIADNE_PRECONDITION(spc.dimension()==fp.argument_size());
}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(Map<RealVariable,IntervalDomainType> dom, Vector<RealExpression> expr)
    : ExpressionPatch(RealSpacePatch(List<RealVariable>(dom.keys()),dom),expr) {
}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(Vector<RealVariable> vars, BoxDomainType fdom, Function<P,SIG> f)
    : _vars(vars), _fp(fdom,f)
{
    assert(fdom.dimension()==f.argument_size());
    assert(f.argument_size()==vars.size());
}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(Vector<RealVariable> vars, FunctionPatch<P,SIG> fp)
    : _vars(vars), _fp(fp)
{
    assert(fp.argument_size()==vars.size());
}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(InitializerList<ExpressionPatch<ValidatedTag,RealScalar>> lst)
    : ExpressionPatch(*lst.begin())
{
    assert (lst.size()!=0);
    bool first=true;
    for (auto ep : lst) {
        if (first) { first=false; }
        else { this->_adjoin(ep); }
    }
}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(InitializerList<PreExpressionPatch> lst)
    : ExpressionPatch(*lst.begin())
{
    assert (lst.size()!=0);
    bool first=true;
    for (auto ep : lst) {
        if (first) { first=false; }
        else { *this=_join(*this,ep); }
        std::cerr<<"*this="<<*this<<"\n";
    }
}




auto call(FunctionPatch<ValidatedTag,RealVector(RealVector)>, ExpressionPatch<ValidatedTag,RealVector>)
    -> ExpressionPatch<ValidatedTag,RealVector>;


template<class P, class... ARGS> FunctionPatch<P,RealVector(ARGS...)>::
FunctionPatch(RealSpacePatch const& spc, ExpressionPatch<P,RES> const& ep) {
    BoxDomainType dom=spc.domain();

    Array<SizeType> coords(ep._vars.size());
    for(SizeType i=0; i!=coords.size(); ++i) {
        coords[i]=spc.space()[ep._vars[i]];
    }
    VectorMultivariateFunction<P> pr=coordinates<P>(spc.dimension(),coords);

    ValidatedVectorMultivariateFunctionPatch proj=factory(ep._fp).create(dom,pr);

    *this = compose(ep._fp,proj);
}

#warning
/*
Vector<RealVariable> _vars;
BoxDomainType _fdom;
Function<P,RealVector(RealVector)> _f;
*/

template<class P, class RES, class... ARGS> auto
make_expression_patch(FunctionPatch<P,RES(ARGS...)> const& fp, List<ConstantOrVariable<Real>> const& cvl) {
    ARIADNE_ASSERT(cvl.size()==fp.argument_size());

    List<RealVariable> variable_list;
    Map<RealVariable,SizeType> indices;
    Map<RealVariable,IntervalDomainType> domains;
    SizeType k=0;
    for (SizeType i=0; i!=cvl.size(); ++i) {
        if (std::holds_alternative<Variable<Real>>(cvl[i])) {
            RealVariable v=std::get<RealVariable>(cvl[i]);
            if (not indices.has_key(v)) {
                variable_list.append(v);
                indices[v]=k;
                domains[v]=fp.domain()[i];
                ++k;
            } else {
                domains[v]=intersection(domains[v],fp.domain()[i]);
            }
        }
    }

    Vector<RealVariable> variables(variable_list);
    BoxDomainType domain(variables.size(),[&](SizeType i){return domains[variables[i]];});

    FunctionPatch<P,RealVector(ARGS...)> projection(domain, Function<P,RealVector(ARGS...)>::zeros(fp.argument_size(),domain.dimension()));
    for (SizeType i=0; i!=fp.argument_size(); ++i) {
       if (std::holds_alternative<Variable<Real>>(cvl[i])) {
            RealVariable v=std::get<RealVariable>(cvl[i]);
            projection[i]=Function<P,RealScalar(ARGS...)>::coordinate(domain.dimension(), indices[v]);
        } else {
            RealConstant c=std::get<RealConstant>(cvl[i]);
            projection[i]=Function<P,RealScalar(ARGS...)>::constant(domain.dimension(),c);
        }
    }
    FunctionPatch<P,RealVector(ARGS...)> function_patch=compose(fp,projection);
    return ExpressionPatch<P,RealVector>(variables,function_patch);
}

template<class P,class T1, class T2>
Tuple<RealSpace,MultivariateFunctionPatch<P,T1>,MultivariateFunctionPatch<P,T2>>
make_common_variables(ExpressionPatch<P,T1> ep1, ExpressionPatch<P,T2> ep2) {
    List<RealVariable> variables;
    Map<RealVariable,SizeType> indices;
    Map<RealVariable,IntervalDomainType> domains;
    SizeType as1=ep1._vars.size();
    SizeType as2=ep2._vars.size();
    for (SizeType i=0; i!=as1; ++i) {
        auto v=ep1._vars[i];
        auto d=ep1._fp.domain()[i];
        variables.append(v);
        indices[v]=i;
        domains.insert(v,d);
    }
    SizeType asr=as1;
    for (SizeType i=0; i!=as2; ++i) {
        auto v=ep2._vars[i];
        auto d=ep2._fp.domain()[i];
        if (domains.has_key(v)) {
            domains[v]=intersection(domains[v],d);
        } else {
            variables.append(v);
            indices[v]=asr; ++asr;
            domains.insert(v,d);
        }
    }
    std::cerr<<"vars1="<<ep1._vars<<"\n";
    std::cerr<<"vars2="<<ep2._vars<<"\n";
    std::cerr<<"variables="<<variables<<"\n";
    std::cerr<<"indices="<<indices<<"\n";
    std::cerr<<"domains="<<domains<<"\n";

    BoxDomainType domain(domains.size(),[&](SizeType i){return domains[variables[i]];});
    std::cerr<<"domain="<<domain<<"\n";
    VectorMultivariateFunctionPatch<P> p1=factory(ep1._fp).create_zeros(as1,domain);
    VectorMultivariateFunctionPatch<P> p2=factory(ep2._fp).create_zeros(as2,domain);
    for (SizeType i=0; i!=as1; ++i) {
        p1[i]=ScalarMultivariateFunction<P>::coordinate(asr,i);
    }
    std::cerr<<"p1="<<p1<<"\n";
    std::cerr<<"  p2="<<p2<<"\n";
    for (SizeType i=0; i!=as2; ++i) {
        auto v=ep2._vars[i]; p2[i]=ScalarMultivariateFunction<P>::coordinate(asr,indices[v]);
    }
    std::cerr<<"p2="<<p2<<"\n";

    RealSpace rspc(variables);
    auto pf1=compose(ep1._fp,p1);
    std::cerr<<"pf1="<<pf1<<"\n";
    auto pf2=compose(ep2._fp,p2);
    std::cerr<<"pf2="<<pf2<<"\n";

    return std::make_tuple(rspc,pf1,pf2);
}

Void ExpressionPatch<ValidatedTag,RealVector>::_adjoin(ExpressionPatch<P,RealScalar> const& ep) {
    ExpressionPatch<P,RealVector> const& ep1=*this; ExpressionPatch<P,RealScalar> const& ep2=ep;
    auto sff = make_common_variables(ep1,ep2); using std::get;
    *this = ExpressionPatch<P,T>(get<0>(sff), join(get<1>(sff),get<2>(sff)));
}

auto ExpressionPatch<ValidatedTag,RealVector>::_join(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
    std::cerr<<"_join(ep1,ep2):\n  ep1="<<ep1<<"\n  ep2="<<ep2<<"\n";
    auto sff = make_common_variables(ep1,ep2); using std::get;
    std::cerr<<"    sff="<<sff<<"\n";
    return ExpressionPatch<P,T>(get<0>(sff), join(get<1>(sff),get<2>(sff)));
}


auto ExpressionPatch<ValidatedTag,RealVector>::_add(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
    auto sff = make_common_variables(ep1,ep2); using std::get;
    return ExpressionPatch<P,T>(get<0>(sff), operator+(get<1>(sff),get<2>(sff)));
}



auto ExpressionPatch<ValidatedTag,RealVector>::_write(OutputStream& os) const -> OutputStream& {
    ExpressionPatch<P,T> const& ep=*this;
    os << "ExpressionPatch({";
    for (SizeType i=0; i!=ep._vars.size(); ++i) {
        if (i!=0) { os << ","; } os << ep._vars[i]<<"|"<<ep._fp.domain()[i];
    }
    return os <<"}, f=" << static_cast<Function<P,SIG>>(ep._fp) << ")";
}



template<class T0, class... TS> inline OutputStream& operator<<(OutputStream& os, Variant<T0,TS...> const& var) {
    std::visit([&os](auto t){os<<t;},var); return os; }

template<class P, class... ARGS> auto FunctionPatch<P,RealScalar(ARGS...)>::
operator() (ConstantOrVariable<ARGS>... cvs) const -> ExpressionPatch<P,RealScalar> {
    List<ConstantOrVariable<Real>> cvl=make_constant_or_variable_list(cvs...);
    return make_expression_patch(*this,cvl);
}

template<class P, class... ARGS> auto FunctionPatch<P,RealVector(ARGS...)>::
operator() (InitializerList<AnyConstantOrVariablePatch> const& cvs) const -> ExpressionPatch<P,RealVector> {
    std::cerr<<"FunctionPatch(InitializerList<AnyConstantOrVariablePatch> cvs)\n";
    List<ConstantOrVariablePatch<Real>> cvl=make_constant_or_variable_patch_list(cvs);
    return make_expression_patch(*this,cvl);
}

template<class P, class... ARGS> auto FunctionPatch<P,RealVector(ARGS...)>::
operator() (ConstantOrVariable<ARGS>... cvs) const -> ExpressionPatch<P,RealVector> {
    std::cerr<<"FunctionPatch(ConstantOrVariable<ARGS>... cvs): ";
    List<ConstantOrVariable<Real>> cvl=make_constant_or_variable_list(cvs...);
    std::cerr<<"cvs="<<std::make_tuple(cvs...)<<", ";
    std::cerr<<"cvl="<<cvl<<"\n";
    return make_expression_patch(*this,cvl);
}

template<class P, class... ARGS> auto FunctionPatch<P,RealVector(ARGS...)>::
operator() (ExpressionPatch<P,ARGS...> const& ep) const -> ExpressionPatch<P,RealVector> {
    std::cerr<<"FunctionPatch(ExpressionPatch<P,ARGS...> ep): ep="<<ep<<"\n";
    static_assert (Same<Tuple<ARGS...>,Tuple<RealVector>>);
    return ExpressionPatch<P,RealVector>(ep._vars,compose(*this,ep._fp));
}

} // namespace Ariadne

#include "algebra/sweeper.hpp"
#include "function/function_traits.hpp"
#include "function/taylor_function.hpp"

namespace Ariadne {
/*
inline FunctionPatch<EffectiveTag,RealScalar(RealVector)> make_function_patch(const BoxDomainType& dom, const Function<EffectiveTag,RealScalar(RealVector)>& f) {
    ARIADNE_NOT_IMPLEMENTED;
}
inline FunctionPatch<EffectiveTag,RealVector(RealVector)> make_function_patch(const BoxDomainType& dom, const Function<EffectiveTag,RealVector(RealVector)>& f) {
    ARIADNE_NOT_IMPLEMENTED;
}
*/

inline FunctionPatch<ValidatedTag,RealScalar(RealVector)> make_function_patch(const BoxDomainType& dom, const Function<ValidatedTag,RealScalar(RealVector)>& f) {
    return ValidatedScalarMultivariateTaylorFunctionModel<FloatDP>(dom,f,ThresholdSweeper<FloatDP>(dp,1e-8));
}
inline FunctionPatch<ValidatedTag,RealVector(RealVector)> make_function_patch(const BoxDomainType& dom, const Function<ValidatedTag,RealVector(RealVector)>& f) {
    return ValidatedVectorMultivariateTaylorFunctionModel<FloatDP>(dom,f,ThresholdSweeper<FloatDP>(dp,1e-8));
}

template<class P, class... ARGS> FunctionPatch<P,RealScalar(ARGS...)>::FunctionPatch(const DomainType& dom, const Function<P,SIG>& f)
    : FunctionPatch(make_function_patch(dom,f)) { }
template<class P, class... ARGS> FunctionPatch<P,RealVector(ARGS...)>::FunctionPatch(const DomainType& dom, const Function<P,SIG>& f)
    : FunctionPatch(make_function_patch(dom,f)) { }

} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_PATCH_HPP */
