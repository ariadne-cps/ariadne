/***************************************************************************
 *            function/formula_function.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include "../config.hpp"

#include "../numeric/numeric.hpp"

#include "../numeric/operators.hpp"
#include "../algebra/differential.hpp"
#include "../algebra/algebra.hpp"
#include "../function/taylor_model.hpp"

#include "../function/formula.hpp"
#include "../function/formula.tpl.hpp"

#include "../function/function.decl.hpp"
#include "../function/function.hpp"
#include "../function/function_interface.hpp"
#include "../function/function_mixin.hpp"
#include "../function/function_mixin.tpl.hpp"

#include "../function/formula_function.hpp"

#include "../symbolic/expression.hpp"


namespace Ariadne {

template<class Y> template<class X> Void ScalarUnivariateFormulaFunction<Y>::_compute(Scalar<X>& r, const Scalar<X>& x) const {
   r=Ariadne::evaluate(this->_formula,Vector<X>({x})); }

template<class Y> template<class X> Void VectorUnivariateFormulaFunction<Y>::_compute(Vector<X>& r, const Scalar<X>& x) const {
   r=Ariadne::evaluate(this->_formulae,Vector<X>({x})); }

template<class Y> template<class X> Void ScalarFormulaFunction<Y>::_compute(Scalar<X>& r, const Vector<X>& x) const {
    r=Ariadne::evaluate(this->_formula,x); }

template<class Y> template<class X> Void  VectorFormulaFunction<Y>::_compute(Vector<X>& r, const Vector<X>& x) const {
    r=Ariadne::evaluate(this->_formulae,x); }


template class FunctionMixin<ScalarUnivariateFormulaFunction<ApproximateNumber>,ApproximateTag,RealScalar(RealScalar)>;
template class FunctionMixin<VectorUnivariateFormulaFunction<ApproximateNumber>,ApproximateTag,RealVector(RealScalar)>;
template class FunctionMixin<ScalarFormulaFunction<ApproximateNumber>,ApproximateTag,RealScalar(RealVector)>;
template class FunctionMixin<VectorFormulaFunction<ApproximateNumber>,ApproximateTag,RealVector(RealVector)>;

template class FunctionMixin<ScalarUnivariateFormulaFunction<ValidatedNumber>,ValidatedTag,RealScalar(RealScalar)>;
template class FunctionMixin<VectorUnivariateFormulaFunction<ValidatedNumber>,ValidatedTag,RealVector(RealScalar)>;
template class FunctionMixin<ScalarFormulaFunction<ValidatedNumber>,ValidatedTag,RealScalar(RealVector)>;
template class FunctionMixin<VectorFormulaFunction<ValidatedNumber>,ValidatedTag,RealVector(RealVector)>;

template class FunctionMixin<ScalarUnivariateFormulaFunction<EffectiveNumber>,EffectiveTag,RealScalar(RealScalar)>;
template class FunctionMixin<VectorUnivariateFormulaFunction<EffectiveNumber>,EffectiveTag,RealVector(RealScalar)>;
template class FunctionMixin<ScalarFormulaFunction<EffectiveNumber>,EffectiveTag,RealScalar(RealVector)>;
template class FunctionMixin<VectorFormulaFunction<EffectiveNumber>,EffectiveTag,RealVector(RealVector)>;


} // namespace Ariadne
