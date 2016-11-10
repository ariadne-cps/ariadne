/***************************************************************************
 *            function_model.tcc
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

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function/function_model.h"

#include "function/function_interface.h"
#include "function/function_mixin.h"
#include "function/function.h"

#include "numeric/operators.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/operations.h"
#include "geometry/box.h"

namespace Ariadne {

template<class P> ScalarFunctionModel<P> ScalarFunctionModel<P>::create_zero() const {
    return ScalarFunctionModel<P>(this->_ptr->_create_zero(this->domain()));
}

template<class P> ScalarFunctionModel<P> ScalarFunctionModel<P>::create_constant(CanonicalNumericType<P> const& c) const {
    return ScalarFunctionModel<P>(this->_ptr->_create_constant(this->domain(),c));
}

template<class P> inline ScalarFunctionModel<P> ScalarFunctionModel<P>::create_coordinate(SizeType j) const {
    return ScalarFunctionModel<P>(this->_ptr->_create_coordinate(this->domain(),j));
}

template<class P> inline VectorFunctionModel<P> ScalarFunctionModel<P>::create_identity() const {
    return this->_ptr->_create_identity();
}

template<class P> inline ScalarFunctionModel<P> ScalarFunctionModel<P>::create(const ScalarFunction<P>& g) const {
    ScalarFunctionModel<P> const& f=*this; return compose(g,f.create_identity());
}


template<class P> inline Vector<ScalarFunctionModel<P>> ScalarFunctionModel<P>::create_coordinates(DomainType const& dom) const {
    Vector<ScalarFunctionModel<P>> res(dom.dimension(),this->_ptr->_create_zero(dom));
    for(SizeType i=0; i!=dom.dimension(); ++i) { res[i]=this->_ptr->_create_coordinate(dom,i); }
    return res;
}



template<class P> ScalarFunctionModel<P> FunctionModelFactoryInterface<P>::create(const ExactBoxType& domain, const ScalarFunctionInterface<P>& function) const {
    return this->_create(domain,function);
}
template<class P> VectorFunctionModel<P> FunctionModelFactoryInterface<P>::create(const ExactBoxType& domain, const VectorFunctionInterface<P>& function) const {
    return this->_create(domain,function);
}

template<class P> ScalarFunctionModel<P> FunctionModelFactoryInterface<P>::create_zero(const ExactBoxType& domain) const {
    return this->_create(domain,EffectiveScalarFunction::zero(domain.size())); }
template<class P> VectorFunctionModel<P> FunctionModelFactoryInterface<P>::create_zeros(SizeType result_size, const ExactBoxType& domain) const {
    return this->_create(domain,EffectiveVectorFunction::zeros(result_size,domain.size())); }
template<class P> ScalarFunctionModel<P> FunctionModelFactoryInterface<P>::create_constant(const ExactBoxType& domain, const Number<P>& value) const {
    CanonicalNumericType<P> concrete_value=this->create_number(value); return this->create_constant(domain,concrete_value); }
template<class P> ScalarFunctionModel<P> FunctionModelFactoryInterface<P>::create_constant(const ExactBoxType& domain, const CanonicalNumericType<P>& value) const {
    return ScalarFunctionModel<P>(this->_create(domain,EffectiveScalarFunction::zero(domain.size())))+value; };
template<class P> VectorFunctionModel<P> FunctionModelFactoryInterface<P>::create_constants(const ExactBoxType& domain, const Vector<Number<P>>& values) const {
    typename CanonicalNumericType<P>::PrecisionType pr; Vector<CanonicalNumericType<P>> concrete_values(values,pr); return this->create_constants(domain,concrete_values); }
template<class P> VectorFunctionModel<P> FunctionModelFactoryInterface<P>::create_constants(const ExactBoxType& domain, const Vector<CanonicalNumericType<P>>& values) const {
    return VectorFunctionModel<P>(this->_create(domain,EffectiveVectorFunction::zeros(values.size(),domain.size())))+values; };
template<class P> ScalarFunctionModel<P> FunctionModelFactoryInterface<P>::create_coordinate(const ExactBoxType& domain, SizeType index) const {
    return ScalarFunctionModel<P>(this->_create(domain,EffectiveScalarFunction::coordinate(domain.size(),index))); };
template<class P> ScalarFunctionModel<P> FunctionModelFactoryInterface<P>::create_identity(const ExactIntervalType& domain) const {
    return this->_create(ExactBoxType(1,domain),EffectiveScalarFunction::coordinate(1,0)); };
template<class P> VectorFunctionModel<P> FunctionModelFactoryInterface<P>::create_identity(const ExactBoxType& domain) const {
    return this->_create(domain,EffectiveVectorFunction::identity(domain.size())); };

template<class P> CanonicalNumericType<P> FunctionModelFactoryInterface<P>::create_number(const Number<P>& number) const {
    return CanonicalNumericType<P>(number,typename CanonicalNumericType<P>::PrecisionType()); }


} // namespace Ariadne
