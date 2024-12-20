/***************************************************************************
 *            symbolic/expression_patch.cpp
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

#include "symbolic/expression_patch.hpp"

#include "algebra/sweeper.hpp"
#include "function/taylor_function.hpp"

#include "symbolic/expression_patch.tpl.hpp"


#warning Temporarily removed RestrictedExpression - FunctionPatch methods in FunctionPatch lines ~250, 500


namespace Ariadne {

template FunctionPatch<ValidatedTag,RealScalar(RealVector)>::FunctionPatch(const DomainType&, const Function<ValidatedTag,RealScalar(RealVector)>&);
template FunctionPatch<ValidatedTag,RealVector(RealVector)>::FunctionPatch(const DomainType&, const Function<ValidatedTag,RealVector(RealVector)>&);

template FunctionPatch<ValidatedTag,RealScalar(RealVector)>::FunctionPatch(const RealSpacePatch&, const RestrictedExpression<ValidatedTag,RealScalar>&);
template FunctionPatch<ValidatedTag,RealVector(RealVector)>::FunctionPatch(const RealSpacePatch&, const RestrictedExpression<ValidatedTag,RealVector>&);

template FunctionPatch<ValidatedTag,RealScalar(RealVector)>::FunctionPatch(const RealSpacePatch&, const Expression<RealScalar>&);
template FunctionPatch<ValidatedTag,RealVector(RealVector)>::FunctionPatch(const RealSpacePatch&, const Expression<RealVector>&);

template auto FunctionPatch<ValidatedTag,RealScalar(RealVector)>::operator() (RestrictedExpression<ValidatedTag,RealVector> const&) const -> RestrictedExpression<ValidatedTag,RealScalar>;
template auto FunctionPatch<ValidatedTag,RealVector(RealVector)>::operator() (RestrictedExpression<ValidatedTag,RealVector> const&) const -> RestrictedExpression<ValidatedTag,RealVector>;


VariableDomainMap<Real>::VariableDomainMap() {
}
VariableDomainMap<Real>::VariableDomainMap(Map<RealVariable,IntervalDomainType> m)
    : Map<RealVariable, IntervalDomainType>(m) {
}
VariableDomainMap<Real>::VariableDomainMap(VariableIntervalDomainType vivl) {
    this->insert(vivl.variable(),vivl.interval());
}
VariableDomainMap<Real>::VariableDomainMap(VectorVariableBoxDomainType vbx) {
    for (SizeType i=0; i!=vbx.variable().size(); ++i) {
        this->insert(vbx.variable()[i],vbx.box()[i]);
    }
}

VariableDomainMap<Real>::VariableDomainMap<Real>::VariableDomainMap(InitializerList<VariableDomainMap<R>> lst) {
    for (auto doms : lst) {
        for (auto vdom : doms) {
            RealVariable const& v = vdom.first;
            IntervalDomainType const& dom = vdom.second;
            if (this->has_key(v)) {
                (*this)[v]=intersection((*this)[v],dom);
            } else {
                (*this)[v]=dom;
            }
        }
    }
}


#warning
ValidatedVectorMultivariateFunction combine(const ValidatedVectorMultivariateFunction& f1, const ValidatedVectorMultivariateFunction& f2) {
    SizeType as1=f1.argument_size();
    SizeType as2=f2.argument_size();
    ARIADNE_NOT_IMPLEMENTED;
//    auto p1=ValidatedVectorMultivariateFunction::projection(as1+as2,range(0,as1));
//    auto p2=ValidatedVectorMultivariateFunction::projection(as1+as2,range(as1,as1+as2));
//    return join(compose(as1,p1),compose(as2,p2));
}

//template class RestrictedExpression<ValidatedTag,Real>;
//template class RestrictedExpression<ValidatedTag,Vector<Real>>;

template Tuple<Vector<RealVariable>,ValidatedScalarMultivariateFunctionPatch,ValidatedScalarMultivariateFunctionPatch>
make_common_variables(ValidatedRestrictedExpression<RealScalar> ep1, ValidatedRestrictedExpression<RealScalar> ep2);

template Tuple<Vector<RealVariable>,ValidatedScalarMultivariateFunctionPatch,ValidatedVectorMultivariateFunctionPatch>
make_common_variables(ValidatedRestrictedExpression<RealScalar> ep1, ValidatedRestrictedExpression<RealVector> ep2);

template Tuple<Vector<RealVariable>,ValidatedVectorMultivariateFunctionPatch,ValidatedScalarMultivariateFunctionPatch>
make_common_variables(ValidatedRestrictedExpression<RealVector> ep1, ValidatedRestrictedExpression<RealScalar> ep2);

template Tuple<Vector<RealVariable>,ValidatedVectorMultivariateFunctionPatch,ValidatedVectorMultivariateFunctionPatch>
make_common_variables(ValidatedRestrictedExpression<RealVector> ep1, ValidatedRestrictedExpression<RealVector> ep2);



} // namespace Ariadne
