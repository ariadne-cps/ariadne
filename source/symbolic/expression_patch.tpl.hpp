/***************************************************************************
 *            symbolic/expression_patch.tpl.hpp
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

#ifndef ARIADNE_EXPRESSION_PATCH_TPL_HPP
#define ARIADNE_EXPRESSION_PATCH_TPL_HPP

#include "symbolic/expression.hpp"
#include "symbolic/expression_patch.hpp"

#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/restricted_function.hpp"

namespace Ariadne {


#warning Don't need ScalarOrVector
/*
template<class X> class ScalarOrVector
    : public Variant<Scalar<X>,Vector<X>>
{
  public:
    template<class S> requires Convertible<S,Scalar<X>> ScalarOrVector(S const& s)
        : Variant<Scalar<X>,Vector<X>>(Scalar<X>(s)) { }
    template<class V> requires Convertible<V,Vector<X>> ScalarOrVector(V const& v)
        : Variant<Scalar<X>,Vector<X>>(Scalar<X>(v)) { }
};


template<template<class>class T, class R> class OfScalarOrVector
    : public Variant<T<Scalar<R>>,T<Vector<R>>>
{
    using Base = Variant<T<Scalar<R>>,T<Vector<R>>>;
  public:
    template<class S> requires Convertible<S,T<Scalar<R>>> OfScalarOrVector(S const& s)
        : Base(T<Scalar<R>>(s)) { }
    template<class V> requires Convertible<V,T<Vector<R>>> OfScalarOrVector(V const& v)
        : Base(T<Vector<R>>(v)) { }
};

template<class R> Vector<R> join(List<ScalarOrVector<R>> const& lst);
template<template<class>class T, class R> T<Vector<R>> join(List<OfScalarOrVector<T,R>> const&);


template<class R> class ConstantOrVariable
    : public Variant<Constant<R>,Variable<R>>
{
    using Base = Variant<Constant<R>,Variable<R>>;
  public:
    template<class C> requires Convertible<C,R> ConstantOrVariable(C const& c)
        : Base(Constant<R>(R(c))) { }
    template<class C> requires Convertible<C,Constant<R>> and (not Convertible<C,R>) ConstantOrVariable(C const& c)
        : Base(Constant<R>(c)) { }
    template<class V> requires Convertible<V,Variable<R>> ConstantOrVariable(V const& v)
        : Base(Variable<R>(v)) { }
};

template<class R> class ConstantOrVariablePatch
    : public Variant<Constant<R>,VariablePatch<R>>
{
    using Base = Variant<Constant<R>,VariablePatch<R>>;
  public:
    template<class C> requires Convertible<C,R> ConstantOrVariablePatch(C const& c)
        : Base(Constant<R>(R(c))) { }
    template<class C> requires Convertible<C,Constant<R>> and (not Convertible<C,R>) ConstantOrVariablePatch(C const& c)
        : Base(Constant<R>(c)) { }
    template<class V> requires Convertible<V,VariablePatch<R>> ConstantOrVariablePatch(V const& v)
        : Base(Variable<R>(v)) { }
};
*/


Scalar<RealExpression> call(EffectiveScalarMultivariateFunction const& f, Vector<RealExpression> const & ev) {
    Vector<ElementaryAlgebra<Real>> av=Vector<ElementaryAlgebra<Real>>(ev);
    ElementaryAlgebra<Real> rav=f(av);
    return rav.extract<RealExpression>();
}
Vector<RealExpression> call(EffectiveVectorMultivariateFunction const& f, Vector<RealExpression> const & ev) {
    Vector<ElementaryAlgebra<Real>> av=Vector<ElementaryAlgebra<Real>>(ev);
    Vector<ElementaryAlgebra<Real>> rav=f(av);
    return Vector<RealExpression>(rav.size(),[&rav](SizeType i){return rav[i].extract<RealExpression>();});
}


template<class P> VectorMultivariateFunction<P> coordinates(SizeType n, Array<SizeType> is) {
    if (is.size()==0) { return VectorMultivariateFunction<P>(0,n); }
    return Vector<ScalarMultivariateFunction<P>>(is.size(), [n,&is](SizeType j){return ScalarMultivariateFunction<P>::coordinate(n,is[j]);});
}
template<class P> VectorMultivariateFunction<P> coordinates(SizeType n, Range rng) {
    if (rng.size()==0) { return VectorMultivariateFunction<P>(0,n); }
    return Vector<ScalarMultivariateFunction<P>>(rng.size(), [n,rng](SizeType j){return ScalarMultivariateFunction<P>::coordinate(n,rng[j]);});
}


template<class P> VectorMultivariateFunction<P> combine(VectorMultivariateFunction<P> f1, VectorMultivariateFunction<P> f2) {
//    return make_combined_function(f1,f2);
    auto n1=f1.argument_size(); auto n2=f2.argument_size(); auto n=n1+n2;
    auto p2=coordinates<P>(n,range(n1,n1+n2));
    auto p1=coordinates<P>(n,range(0,n1));
    return join(compose(f1,p1),compose(f2,p2));
}


template<class P, class T1, class T2>
Tuple<Vector<RealVariable>,MultivariateFunctionPatch<P,T1>,MultivariateFunctionPatch<P,T2>>
make_common_variables(RestrictedExpression<P,T1> ep1, RestrictedExpression<P,T2> ep2);


namespace {
Vector<RealVariable> variables(Vector<RealVariablePatch> ptchs) {
    return Vector<RealVariable>(ptchs.size(),[&ptchs](SizeType i){return ptchs[i].variable();});
}
BoxDomainType euclidean_set(Vector<RealVariablePatch> ptchs) {
    return BoxDomainType(ptchs.size(),[&ptchs](SizeType i){return ptchs[i].interval();});
}
RestrictedExpression<ValidatedTag,Real> make_restricted_expression(RealVariableDomainMap vdoms, RealExpression expr) {
    Vector<RealVariable> vars(List<RealVariable>(vdoms.keys()));
    Function<ValidatedTag,Real(Vector<Real>)> f(vars,expr);
    BoxDomainType dom(vars.size(),[&vars,&vdoms](SizeType i){return vdoms[vars[i]];});
    return RestrictedExpression<ValidatedTag,Real>(vars,dom,f);
}
RestrictedExpression<ValidatedTag,Vector<Real>> make_restricted_expression(RealVariableDomainMap vdoms, Vector<RealExpression> expr) {
    Vector<RealVariable> vars(List<RealVariable>(vdoms.keys()));
    Function<ValidatedTag,Vector<Real>(Vector<Real>)> f(vars,expr);
    BoxDomainType dom(vars.size(),[&vars,&vdoms](SizeType i){return vdoms[vars[i]];});
    return RestrictedExpression<ValidatedTag,Vector<Real>>(vars,dom,f);
}
}// namespace



RestrictedExpression<ValidatedTag,RealScalar>::RestrictedExpression(Constant<Real> const& c)
    : RestrictedExpression(RealVariableDomainMap(),Expression<Real>(c))
{
}

//RestrictedExpression<ValidatedTag,RealScalar>::RestrictedExpression(Variable<Real> const& v)
//    : RestrictedExpression(RealVariableDomainMap(),Expression<Real>(v))
//{
//}

RestrictedExpression<ValidatedTag,RealScalar>::RestrictedExpression(VariablePatchType<Real> const& vptch)
    : RestrictedExpression(RealVariableDomainMap(vptch),Expression<Real>(vptch.variable()))
{
}

RestrictedExpression<ValidatedTag,RealScalar>::RestrictedExpression(Vector<RealVariable> vars, FunctionPatch<P,SIG> fp)
    : _vars(vars), _fp(fp)
{
    assert(fp.argument_size()==vars.size());
}

RestrictedExpression<ValidatedTag,RealScalar>::RestrictedExpression(Vector<RealVariable> vars, BoxDomainType fdom, Function<P,SIG> f)
    : _vars(vars), _fp(fdom,f)
{
    assert(fdom.dimension()==f.argument_size());
    assert(f.argument_size()==vars.size());
}

RestrictedExpression<ValidatedTag,RealScalar>::RestrictedExpression(Vector<RealVariablePatch> ptchs, Function<P,SIG> f)
    : RestrictedExpression(variables(ptchs),euclidean_set(ptchs),f) {
}

RestrictedExpression<ValidatedTag,RealScalar>::RestrictedExpression(RealVariableDomainMap doms, RealExpression expr)
    : RestrictedExpression(make_restricted_expression(doms,expr)) {
}

auto RestrictedExpression<ValidatedTag,RealScalar>::domains() const -> Map<RealVariable,IntervalDomainType> {
    Map<RealVariable,IntervalDomainType> doms;
    for (SizeType i=0; i!=_vars.size(); ++i) {
        RealVariable v=_vars[i];
        IntervalDomainType dom=_fp.domain()[i];
        if (doms.has_key(v)) { doms[v]=intersection(doms[v],dom); }
        else { doms[v]=dom; }
    }
    return doms;
}




//RestrictedExpression<ValidatedTag,RealVector>::RestrictedExpression(Variable<RealVector> const& v)
//{
//    ARIADNE_NOT_IMPLEMENTED;
//}

RestrictedExpression<ValidatedTag,RealVector>::RestrictedExpression(VariablePatch<RealVector> const& vdoms)
    : _vars(List<RealVariable>(vdoms.variables()))
    , _fp(ValidatedVectorMultivariateRestrictedFunction(ValidatedVectorMultivariateFunction::identity(this->_vars.size()),vdoms.euclidean_set(RealSpace(this->_vars))))
{
}


RestrictedExpression<ValidatedTag,RealVector>::RestrictedExpression(Vector<RealVariable> vars, FunctionPatch<P,SIG> fp)
    : _vars(vars), _fp(fp)
{
    assert(fp.argument_size()==vars.size());
}

RestrictedExpression<ValidatedTag,RealVector>::RestrictedExpression(Vector<RealVariable> vars, BoxDomainType fdom, Function<P,SIG> f)
    : _vars(vars), _fp(fdom,f)
{
    assert(fdom.dimension()==f.argument_size());
    assert(f.argument_size()==vars.size());
}

RestrictedExpression<ValidatedTag,RealVector>::RestrictedExpression(Vector<RealVariablePatch> ptchs, Function<P,SIG> f)
    : _vars(variables(ptchs)), _fp(euclidean_set(ptchs),f)
{
    ARIADNE_PRECONDITION(ptchs.size()==f.argument_size());
}

RestrictedExpression<ValidatedTag,RealVector>::RestrictedExpression(RealVariableDomainMap doms, Vector<RealExpression> expr)
    : RestrictedExpression(make_restricted_expression(doms,expr))
{
}

using ValidatedRealScalarRestrictedExpression = RestrictedExpression<ValidatedTag,RealScalar>;
using ValidatedRealVectorRestrictedExpression = RestrictedExpression<ValidatedTag,RealVector>;

template<class T1, class T2> OutputStream& operator<<(OutputStream& os, Variant<T1,T2> const& var) {
    std::visit([&os](auto t){os<<t;},var); return os; }

RestrictedExpression<ValidatedTag,RealVector> join(List<ScalarOrVectorRestrictedExpression<ValidatedTag,Real>> lst) {
    std::cerr<<"join(List<ScalarOrVectorRestrictedExpression>)\n  lst="<<lst<<"\n";
    ARIADNE_ASSERT(not empty(lst));
    RestrictedExpression<ValidatedTag,RealVector> res = std::holds_alternative<ValidatedRealScalarRestrictedExpression>(lst.front()) ?
        ValidatedRealVectorRestrictedExpression({std::get<ValidatedRealScalarRestrictedExpression>(lst.front())}) : std::get<ValidatedRealVectorRestrictedExpression>(lst.front());
    std::cerr<<"res="<<res<<"\n";
    for (SizeType i=1; i!=lst.size(); ++i) { res=std::visit([&res](auto e){return join(res,e);},lst[i]); }
    std::cerr<<"res="<<res<<"\n";
    return res;
}


RestrictedExpression<ValidatedTag,RealVector>::RestrictedExpression(InitializerList<ScalarOrVectorRestrictedExpression<ValidatedTag,Real>> lst)
    : RestrictedExpression(join(lst))
{
}



auto RestrictedExpression<ValidatedTag,RealVector>::domains() const -> Map<RealVariable,IntervalDomainType> {
    Map<RealVariable,IntervalDomainType> doms;
    for (SizeType i=0; i!=_vars.size(); ++i) {
        RealVariable v=_vars[i];
        IntervalDomainType dom=_fp.domain()[i];
        if (doms.has_key(v)) { doms[v]=intersection(doms[v],dom); }
        else { doms[v]=dom; }
    }
    return doms;
}


#warning Overloads provided ad-hoc to allow certain constructions for systems
ValidatedScalarRestrictedExpression operator-(VariableIntervalDomainType xr, Real c) {
    RealVariableDomainMap doms(xr);
    ValidatedScalarRestrictedExpression ep1(doms,xr.variable());
    ValidatedScalarRestrictedExpression ep2(doms,c);
    return ep1-ep2;
}

ValidatedVectorRestrictedExpression operator-(ValidatedVectorRestrictedExpression ep1, VectorVariableBoxDomainType vp2) {
    RealVariableDomainMap doms(vp2);
    ValidatedVectorRestrictedExpression ep2(doms,Vector<Expression<Real>>(Expression<RealVector>(vp2.variable())));
    return ep1-ep2;
}


template<class R> Vector<R> join(List<ScalarOrVector<R>> const& lst) {
    ARIADNE_ASSERT(not lst.empty());
    Vector<R> r = std::holds_alternative<Scalar<R>>(lst[0]) ? Vector<R>(1u,std::get<Scalar<R>>(lst[0])) : std::get<Vector<R>>(lst[0]);
    for (SizeType i=1u; i!=lst.size(); ++i) { r = std::visit([&r](auto sv){return join(r,sv);},lst[i]); }
    return r;
}

/*
template<class CVE> Void
append_scalars(List<Variant<RealConstant,VariableIntervalDomainType,ValidatedScalarRestrictedExpression>>& lst, CVE const& cve) {
    if constexpr (Same<CVE,VariablesBoxDomainType>) {
        for (SizeType i=0; i!=cve.size(); ++i) { lst.append(cve[i]); }
    } else if constexpr (Same<CVE,ValidatedVectorRestrictedExpression>) {
        for (SizeType i=0; i!=cve.size(); ++i) { lst.append(ValidatedScalarRestrictedExpression(cve._vars,cve._fp[i])); }
    } else if constexpr (Same<CVE,Real>) {
        lst.append(RealConstant(cve));
    } else {
        lst.append(cve);
    }
}

List<Variant<RealConstant,VariableIntervalDomainType,ValidatedScalarRestrictedExpression>> make_scalar_list(InitializerList<ScalarOrVectorConstantOrVariablePatchOrRestrictedExpression<Real> lst) {
    List<Variant<RealConstant,VariableIntervalDomainType,ValidatedScalarRestrictedExpression>> res;
    for (auto cve : lst) {
        std::visit([&res](auto cve){append_scalars(res,cve);},cve);
    }
    return res;
}
*/


/*
RestrictedExpression<ValidatedTag,RealVector>
make_expression_patch(List<Variant<RealConstant,RealVariablePatchType,ValidatedScalarRestrictedExpression>> const& lst) {
    List<RealVariable> vars;
    Map<RealVariable,SizeType> inds;
    Map<RealVariable,IntervalDomainType> doms;
    for (auto cve : lst) {
        if (std::holds_alternative<VariableIntervalDomainType>(cve)) {
            VariableIntervalDomainType const& vi = std::get<VariableIntervalDomainType>(cve);
            RealVariable const& v=vi.variable();
            IntervalDomainType dom=vi.interval();
            if (doms.has_key(v)) {
                doms[v]=intersection(doms[v],dom);
            } else {
                inds[v]=vars.size();
                vars.append(v);
                doms[v]=dom;
            }
        } else if (std::holds_alternative<ValidatedScalarRestrictedExpression>(cve)) {
            ValidatedScalarRestrictedExpression const& ep = std::get<ValidatedScalarRestrictedExpression>(cve);
            Vector<RealVariable> const& argvars=ep._vars;
            for (SizeType i=0; i!=argvars.size(); ++i) {
                RealVariable var = argvars[i];
                IntervalDomainType dom = ep._fp.domain()[i];
                if (doms.has_key(var)) {
                    doms[var]=intersection(doms[var],dom);
                } else {
                    inds[var]=vars.size();
                    vars.append(var);
                    doms[var]=dom;
                }
            }
        }
    }
    BoxDomainType bxdom(vars.size(), [&](SizeType i){return doms[vars[i]]; });

    ValidatedScalarMultivariateTaylorFunctionModelDP tfm(bxdom,ThresholdSweeper<FloatDP>(dp,1e-8));
    auto fctry = factory(tfm);

    List<ValidatedScalarMultivariateFunctionPatch> fps;
    for (auto cve : lst) {
        if (std::holds_alternative<RealConstant>(cve)) {
            RealConstant const& rc = std::get<RealConstant>(cve);
            fps.append(fctry.create_constant(rc));
        } else if (std::holds_alternative<VariableIntervalDomainType>(cve)) {
            VariableIntervalDomainType const& vi = std::get<VariableIntervalDomainType>(cve);
            RealVariable const& v=vi.variable();
            SizeType ind=inds[v];
            fps.append(fctry.create_coordinate(ind));
        } else if (std::holds_alternative<ValidatedScalarRestrictedExpression>(cve)) {
            ValidatedScalarRestrictedExpression const& ep = std::get<ValidatedScalarRestrictedExpression>(cve);
            ValidatedVectorMultivariateFunctionPatch projection=fctry.create_zeros(ep._vars.size());
            Vector<RealVariable> const& argvars=ep._vars;
            for (SizeType k=0; k!=argvars.size(); ++k) {
                projection[k]=static_cast<ValidatedScalarMultivariateFunctionPatch>(fctry.create_coordinate(inds[argvars[k]]));
            }
            fps.append(compose(ep._fp,projection));
        }
    }
    return ValidatedVectorRestrictedExpression(Vector<RealVariable>(vars),fps);
}
*/



auto call(FunctionPatch<ValidatedTag,RealVector(RealVector)>, RestrictedExpression<ValidatedTag,RealVector>)
    -> RestrictedExpression<ValidatedTag,RealVector>;


template<class P, class... ARGS> FunctionPatch<P,RealScalar(ARGS...)>::
FunctionPatch(RealSpacePatch const& spc, RestrictedExpression<P,RES> const& ep) {
    BoxDomainType dom=spc.domain();

    Array<SizeType> coords(ep._vars.size());
    for(SizeType i=0; i!=coords.size(); ++i) {
        coords[i]=spc.space()[ep._vars[i]];
    }
    VectorMultivariateFunction<P> pr=coordinates<P>(spc.dimension(),coords);

    ValidatedVectorMultivariateFunctionPatch proj=factory(ep._fp).create(dom,pr);

    *this = compose(ep._fp,proj);
}

template<class P, class... ARGS> FunctionPatch<P,RealVector(ARGS...)>::
FunctionPatch(RealSpacePatch const& spc, RestrictedExpression<P,RES> const& ep) {
    BoxDomainType dom=spc.domain();

    Array<SizeType> coords(ep._vars.size());
    for(SizeType i=0; i!=coords.size(); ++i) {
        coords[i]=spc.space()[ep._vars[i]];
    }
    VectorMultivariateFunction<P> pr=coordinates<P>(spc.dimension(),coords);

    ValidatedVectorMultivariateFunctionPatch proj=factory(ep._fp).create(dom,pr);

    *this = compose(ep._fp,proj);
}

/*
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

    return RestrictedExpression<P,RES>(variables,domain,compose(fp,projection));
}
*/

template<class P,class T1, class T2>
Tuple<Vector<RealVariable>,MultivariateFunctionPatch<P,T1>,MultivariateFunctionPatch<P,T2>>
make_common_variables(RestrictedExpression<P,T1> ep1, RestrictedExpression<P,T2> ep2) {

    // For each variable used in ep1 or ep2, give it an index,
    // and find its domain, which is the intersection of the two original domains for a common variable.
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

    // Construct the ordered domain
    BoxDomainType domain(domains.size(),[&](SizeType i){return domains[variables[i]];});

    // Make projection functions from the common variables onto those used by each individual expression
    VectorMultivariateFunctionPatch<P> p1=factory(ep1._fp).create_zeros(as1,domain);
    VectorMultivariateFunctionPatch<P> p2=factory(ep2._fp).create_zeros(as2,domain);
    for (SizeType i=0; i!=as1; ++i) {
        auto v=ep1._vars[i]; p1[i]=ScalarMultivariateFunction<P>::coordinate(asr,indices[v]);
    }
    for (SizeType i=0; i!=as2; ++i) {
        auto v=ep2._vars[i]; p2[i]=ScalarMultivariateFunction<P>::coordinate(asr,indices[v]);
    }

    // Make projection functions from the common variables onto those used by each individual expression
    Vector<RealVariable> vars(variables);
    MultivariateFunctionPatch<P,T1> pf1=compose(ep1._fp,p1);
    MultivariateFunctionPatch<P,T2> pf2=compose(ep2._fp,p2);

    return std::make_tuple(vars,pf1,pf2);
}

template<> inline String class_name<Vector<Real>>() { return "RealVector"; }

auto RestrictedExpression<ValidatedTag,RealVector>::_join(RestrictedExpression<P,VT> ep1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
    auto sff = make_common_variables(ep1,ep2); using std::get;
    auto const& dom=std::get<0>(sff);
    auto const& fp1=get<1>(sff);
    auto const& fp2=get<2>(sff);
    return RestrictedExpression<P,T>(dom,join(fp1,fp2));
}

auto RestrictedExpression<ValidatedTag,RealVector>::_join(RestrictedExpression<P,VT> ep1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,VT> {
    auto sff = make_common_variables(ep1,ep2); using std::get;
    auto const& dom=std::get<0>(sff);
    auto const& fp1=get<1>(sff);
    auto const& fp2=get<2>(sff);
    join(fp1,fp2);
    return RestrictedExpression<P,T>(dom,join(fp1,fp2));
}




auto RestrictedExpression<ValidatedTag,RealScalar>::_write(OutputStream& os) const -> OutputStream& {
    RestrictedExpression<P,T> const& ep=*this;
    os << "RestrictedExpression({";
    for (SizeType i=0; i!=ep._vars.size(); ++i) {
        if (i!=0) { os << ","; } os << ep._vars[i]<<"|"<<ep._fp.domain()[i];
    } os << std::flush;
    return os <<"}, f=" << static_cast<Function<P,SIG>>(ep._fp) << ")";
}

auto RestrictedExpression<ValidatedTag,RealVector>::_write(OutputStream& os) const -> OutputStream& {
    RestrictedExpression<P,T> const& ep=*this;
    os << "RestrictedExpression(" << this->_vars << ", " << this->_fp << "\n";
    os << "RestrictedExpression({";
    for (SizeType i=0; i!=ep._vars.size(); ++i) {
        if (i!=0) { os << ","; } os << ep._vars[i]<<"|"<<ep._fp.domain()[i];
    }
//    return os <<"}, f=" << static_cast<Function<P,SIG>>(ep._fp) << ")";
    Function<P,SIG> f=cast_unrestricted(ep._fp);
    return os <<"}," << f << ")";
}




template<class P, class... ARGS> auto
FunctionPatch<P,RealScalar(ARGS...)>::operator() (RestrictedExpression<P,RealVector> const& ep) const -> RestrictedExpression<P,RES> {
    return RestrictedExpression<P,RES>(ep._vars, compose(*this,ep._fp));
}

template<class P, class... ARGS> auto
FunctionPatch<P,RealVector(ARGS...)>::operator() (RestrictedExpression<P,RealVector> const& ep) const -> RestrictedExpression<P,RES> {
    std::cerr<<"FunctionPatch::operator()\n";
    return RestrictedExpression<P,RES>(ep._vars, compose(*this,ep._fp));
}

} // namespace Ariadne


#include "function/taylor_function.hpp"
namespace Ariadne {
template<class P, class... ARGS>
FunctionPatch<P,RealScalar(ARGS...)>::FunctionPatch(const DomainType& dom, const Function<P,SIG>& f)
    : FunctionPatch(ValidatedScalarMultivariateTaylorFunctionModelDP(dom,f,ThresholdSweeper<FloatDP>(dp,1e-8))) { }
template<class P, class... ARGS>
FunctionPatch<P,RealVector(ARGS...)>::FunctionPatch(const DomainType& dom, const Function<P,SIG>& f)
    : FunctionPatch(ValidatedVectorMultivariateTaylorFunctionModelDP(dom,f,ThresholdSweeper<FloatDP>(dp,1e-8))) { }

template<class P, class... ARGS>
FunctionPatch<P,RealScalar(ARGS...)>::FunctionPatch(const RealSpacePatch& spcptch, const Expression<RealScalar>& e)
    : FunctionPatch(RestrictedFunction<P,SIG>(Function<P,SIG>(Vector<RealVariable>(spcptch.variables()),e),spcptch.domain())) { }
template<class P, class... ARGS>
FunctionPatch<P,RealVector(ARGS...)>::FunctionPatch(const RealSpacePatch& spcptch, const Expression<RealVector>& e)
    : FunctionPatch(RestrictedFunction<P,SIG>(Function<P,SIG>(Vector<RealVariable>(spcptch.variables()),Vector<RealExpression>(e)),spcptch.domain())) { }

} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_PATCH_TPL_HPP */
