/***************************************************************************
 *            function_mixin.tcc
 *
 *  Copyright 2008-10  Pieter Collins
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

#include "numeric.h"
#include "vector.h"
#include "differential.h"
#include "taylor_model.h"
#include "formula.h"
#include "algebra.h"

#include "function_mixin.h"

namespace Ariadne {

template<class F> ApproximateNumber ScalarFunctionMixin<F,ApproximateTag>::evaluate(const Vector<ApproximateNumber>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateDifferential ScalarFunctionMixin<F,ApproximateTag>::evaluate(const Vector<ApproximateDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateTaylorModel ScalarFunctionMixin<F,ApproximateTag>::evaluate(const Vector<ApproximateTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateFormula ScalarFunctionMixin<F,ApproximateTag>::evaluate(const Vector<ApproximateFormula>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateAlgebra ScalarFunctionMixin<F,ApproximateTag>::evaluate(const Vector<ApproximateAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateScalarFunctionInterface* ScalarFunctionMixin<F,ApproximateTag>::_clone() const { return new F(static_cast<const F&>(*this)); }

template<class F> ApproximateNumber ScalarFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ApproximateNumber>& x) const { return this->_base_evaluate(x); }
template<class F> ValidatedNumber ScalarFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ValidatedNumber>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateDifferential ScalarFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ApproximateDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> ValidatedDifferential ScalarFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ValidatedDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateTaylorModel ScalarFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ApproximateTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> ValidatedTaylorModel ScalarFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ValidatedTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateFormula ScalarFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ApproximateFormula>& x) const { return this->_base_evaluate(x); }
template<class F> ValidatedFormula ScalarFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ValidatedFormula>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateAlgebra ScalarFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ApproximateAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> ValidatedAlgebra ScalarFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ValidatedAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> ValidatedScalarFunctionInterface* ScalarFunctionMixin<F,ValidatedTag>::_clone() const { return new F(static_cast<const F&>(*this)); }

template<class F> ApproximateNumber ScalarFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ApproximateNumber>& x) const { return this->_base_evaluate(x); }
template<class F> ValidatedNumber ScalarFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ValidatedNumber>& x) const { return this->_base_evaluate(x); }
template<class F> EffectiveNumber ScalarFunctionMixin<F,EffectiveTag>::evaluate(const Vector<EffectiveNumber>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateDifferential ScalarFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ApproximateDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> ValidatedDifferential ScalarFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ValidatedDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateTaylorModel ScalarFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ApproximateTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> ValidatedTaylorModel ScalarFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ValidatedTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateFormula ScalarFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ApproximateFormula>& x) const { return this->_base_evaluate(x); }
template<class F> ValidatedFormula ScalarFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ValidatedFormula>& x) const { return this->_base_evaluate(x); }
template<class F> EffectiveFormula ScalarFunctionMixin<F,EffectiveTag>::evaluate(const Vector<EffectiveFormula>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateAlgebra ScalarFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ApproximateAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> ValidatedAlgebra ScalarFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ValidatedAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> EffectiveAlgebra ScalarFunctionMixin<F,EffectiveTag>::evaluate(const Vector<EffectiveAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> EffectiveScalarFunctionInterface* ScalarFunctionMixin<F,EffectiveTag>::_clone() const { return new F(static_cast<const F&>(*this)); }

template<class F> Vector<ApproximateNumber> ScalarFunctionMixin<F,EffectiveTag>::gradient(const Vector<ApproximateNumber>& x) const {
    return this->_base_evaluate(ApproximateDifferential::variables(1u,x)).gradient(); }
template<class F> Vector<ValidatedNumber> ScalarFunctionMixin<F,EffectiveTag>::gradient(const Vector<ValidatedNumber>& x) const {
    return this->_base_evaluate(ValidatedDifferential::variables(1u,x)).gradient(); }

template<class F> Vector<ApproximateNumber> VectorFunctionMixin<F,ApproximateTag>::evaluate(const Vector<ApproximateNumber>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ApproximateDifferential> VectorFunctionMixin<F,ApproximateTag>::evaluate(const Vector<ApproximateDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ApproximateTaylorModel> VectorFunctionMixin<F,ApproximateTag>::evaluate(const Vector<ApproximateTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ApproximateFormula> VectorFunctionMixin<F,ApproximateTag>::evaluate(const Vector<ApproximateFormula>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ApproximateAlgebra> VectorFunctionMixin<F,ApproximateTag>::evaluate(const Vector<ApproximateAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> ApproximateVectorFunctionInterface* VectorFunctionMixin<F,ApproximateTag>::_clone() const { return new F(static_cast<const F&>(*this)); }

template<class F> Vector<ApproximateNumber> VectorFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ApproximateNumber>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ValidatedNumber> VectorFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ValidatedNumber>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ApproximateDifferential> VectorFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ApproximateDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ValidatedDifferential> VectorFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ValidatedDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ApproximateTaylorModel> VectorFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ApproximateTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ValidatedTaylorModel> VectorFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ValidatedTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ApproximateFormula> VectorFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ApproximateFormula>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ValidatedFormula> VectorFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ValidatedFormula>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ApproximateAlgebra> VectorFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ApproximateAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ValidatedAlgebra> VectorFunctionMixin<F,ValidatedTag>::evaluate(const Vector<ValidatedAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> ValidatedVectorFunctionInterface* VectorFunctionMixin<F,ValidatedTag>::_clone() const { return new F(static_cast<const F&>(*this)); }

template<class F> Vector<ApproximateNumber> VectorFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ApproximateNumber>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ValidatedNumber> VectorFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ValidatedNumber>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<EffectiveNumber> VectorFunctionMixin<F,EffectiveTag>::evaluate(const Vector<EffectiveNumber>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ApproximateDifferential> VectorFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ApproximateDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ValidatedDifferential> VectorFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ValidatedDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ApproximateTaylorModel> VectorFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ApproximateTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ValidatedTaylorModel> VectorFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ValidatedTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ApproximateFormula> VectorFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ApproximateFormula>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ValidatedFormula> VectorFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ValidatedFormula>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<EffectiveFormula> VectorFunctionMixin<F,EffectiveTag>::evaluate(const Vector<EffectiveFormula>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ApproximateAlgebra> VectorFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ApproximateAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<ValidatedAlgebra> VectorFunctionMixin<F,EffectiveTag>::evaluate(const Vector<ValidatedAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<EffectiveAlgebra> VectorFunctionMixin<F,EffectiveTag>::evaluate(const Vector<EffectiveAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> EffectiveVectorFunctionInterface* VectorFunctionMixin<F,EffectiveTag>::_clone() const { return new F(static_cast<const F&>(*this)); }


} // namespace Ariadne
