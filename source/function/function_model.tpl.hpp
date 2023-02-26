/***************************************************************************
 *            function_model.tpl.hpp
 *
 *  Copyright  2011-20  Pieter Collins
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

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function/function_model.hpp"

#include "function/function_interface.hpp"
#include "function/function_mixin.hpp"
#include "function/function.hpp"

#include "numeric/operators.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/operations.hpp"
#include "function/domain.hpp"

namespace Ariadne {


/*
template<class P, class F, class FE> ScalarMultivariateFunctionModel<P,F,FE> FunctionModelFactory<P,F,FE>::create_constant(const DomainType& domain, const CanonicalNumericType<P,F,FE>& value) const {
    return ScalarMultivariateFunctionModel<P,F,FE>(this->_ptr->_create(domain,EffectiveScalarMultivariateFunction::zero(domain.size())))+value; };
template<class P, class F, class FE> VectorMultivariateFunctionModel<P,F,FE> FunctionModelFactory<P,F,FE>::create_constants(const DomainType& domain, const Vector<Number<P>>& values) const {
    typename CanonicalNumericType<P,F,FE>::PrecisionType pr=this->_ptr->create_number(0).precision();
    Vector<CanonicalNumericType<P,F,FE>> concrete_values(values,pr); return this->_ptr->create_constants(domain,concrete_values); }
template<class P, class F, class FE> VectorMultivariateFunctionModel<P,F,FE> FunctionModelFactory<P,F,FE>::create_constants(const DomainType& domain, const Vector<CanonicalNumericType<P,F,FE>>& values) const {
    return VectorMultivariateFunctionModel<P,F,FE>(this->_ptr->_create(domain,EffectiveVectorMultivariateFunction::zeros(values.size(),domain.size())))+values; };
*/


/*
template<class P, class F, class FE> ScalarMultivariateFunctionModel<P,F,FE> FunctionModelFactoryInterface<P,F,FE>::create(const DomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const {
    return this->_create(domain,function);
}
template<class P, class F, class FE> VectorMultivariateFunctionModel<P,F,FE> FunctionModelFactoryInterface<P,F,FE>::create(const DomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const {
    return this->_create(domain,function);
}

template<class P, class F, class FE> ScalarMultivariateFunctionModel<P,F,FE> FunctionModelFactoryInterface<P,F,FE>::create_zero(const DomainType& domain) const {
    return this->_create(domain,EffectiveScalarMultivariateFunction::zero(domain.size())); }
template<class P, class F, class FE> VectorMultivariateFunctionModel<P,F,FE> FunctionModelFactoryInterface<P,F,FE>::create_zeros(SizeType result_size, const DomainType& domain) const {
    return this->_create(domain,EffectiveVectorMultivariateFunction::zeros(result_size,domain.size())); }
template<class P, class F, class FE> ScalarMultivariateFunctionModel<P,F,FE> FunctionModelFactoryInterface<P,F,FE>::create_constant(const DomainType& domain, const Number<P>& value) const {
    CanonicalNumericType<P,F,FE> concrete_value=this->create_number(value); return this->create_constant(domain,concrete_value); }
template<class P, class F, class FE> ScalarMultivariateFunctionModel<P,F,FE> FunctionModelFactoryInterface<P,F,FE>::create_constant(const DomainType& domain, const CanonicalNumericType<P,F,FE>& value) const {
    return ScalarMultivariateFunctionModel<P,F,FE>(this->_create(domain,EffectiveScalarMultivariateFunction::zero(domain.size())))+value; };
template<class P, class F, class FE> VectorMultivariateFunctionModel<P,F,FE> FunctionModelFactoryInterface<P,F,FE>::create_constants(const DomainType& domain, const Vector<Number<P>>& values) const {
    typename CanonicalNumericType<P,F,FE>::PrecisionType pr=this->create_number(0).precision();
    Vector<CanonicalNumericType<P,F,FE>> concrete_values(values,pr); return this->create_constants(domain,concrete_values); }
template<class P, class F, class FE> VectorMultivariateFunctionModel<P,F,FE> FunctionModelFactoryInterface<P,F,FE>::create_constants(const DomainType& domain, const Vector<CanonicalNumericType<P,F,FE>>& values) const {
    return VectorMultivariateFunctionModel<P,F,FE>(this->_create(domain,EffectiveVectorMultivariateFunction::zeros(values.size(),domain.size())))+values; };
template<class P, class F, class FE> ScalarMultivariateFunctionModel<P,F,FE> FunctionModelFactoryInterface<P,F,FE>::create_coordinate(const DomainType& domain, SizeType index) const {
    return ScalarMultivariateFunctionModel<P,F,FE>(this->_create(domain,EffectiveScalarMultivariateFunction::coordinate(domain.size(),index))); };
template<class P, class F, class FE> ScalarMultivariateFunctionModel<P,F,FE> FunctionModelFactoryInterface<P,F,FE>::create_identity(const IntervalDomainType& domain) const {
    return this->_create(DomainType(1,domain),EffectiveScalarMultivariateFunction::coordinate(1,0)); };
template<class P, class F, class FE> VectorMultivariateFunctionModel<P,F,FE> FunctionModelFactoryInterface<P,F,FE>::create_identity(const DomainType& domain) const {
    return this->_create(domain,EffectiveVectorMultivariateFunction::identity(domain.size())); };
template<class P, class F, class FE> CanonicalNumericType<P,F,FE> FunctionModelFactoryInterface<P,F,FE>::create_number(const Number<P>& number) const {
    return this->_create(number); }
*/



} // namespace Ariadne
