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

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/differential.h"
#include "function/taylor_model.h"
#include "expression/formula.h"
#include "algebra/algebra.h"

#include "function/function_mixin.h"

namespace Ariadne {

template<class F,class D, class C> FunctionInterface<ApproximateTag,D,C>* FunctionMixin<F,ApproximateTag,D,C>::_clone() const {
    return new F(static_cast<const F&>(*this)); }
template<class F,class D, class C> FunctionInterface<ValidatedTag,D,C>* FunctionMixin<F,ValidatedTag,D,C>::_clone() const {
    return new F(static_cast<const F&>(*this)); }
template<class F,class D, class C> FunctionInterface<EffectiveTag,D,C>* FunctionMixin<F,EffectiveTag,D,C>::_clone() const {
    return new F(static_cast<const F&>(*this)); }

template<class F,class D, class C> auto
FunctionMixin<F,ApproximateTag,D,C>::_evaluate(const Argument<ApproximateNumber>& x) const -> Result<ApproximateNumber> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ApproximateTag,D,C>::_evaluate(const Argument<ApproximateDifferential>& x) const -> Result<ApproximateDifferential> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ApproximateTag,D,C>::_evaluate(const Argument<ApproximateTaylorModel>& x) const -> Result<ApproximateTaylorModel> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ApproximateTag,D,C>::_evaluate(const Argument<ApproximateAlgebra>& x) const -> Result<ApproximateAlgebra> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ApproximateTag,D,C>::_evaluate(const Argument<ApproximateFormula>& x) const -> Result<ApproximateFormula> {
    return this->_base_evaluate(x); }

template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<ValidatedNumber>& x) const -> Result<ValidatedNumber> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<ValidatedDifferential>& x) const -> Result<ValidatedDifferential> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<ValidatedTaylorModel>& x) const -> Result<ValidatedTaylorModel> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<ValidatedAlgebra>& x) const -> Result<ValidatedAlgebra> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<ValidatedFormula>& x) const -> Result<ValidatedFormula> {
    return this->_base_evaluate(x); }

template<class F,class D, class C> auto
FunctionMixin<F,EffectiveTag,D,C>::_evaluate(const Argument<EffectiveNumber>& x) const -> Result<EffectiveNumber> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,EffectiveTag,D,C>::_evaluate(const Argument<EffectiveAlgebra>& x) const -> Result<EffectiveAlgebra> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,EffectiveTag,D,C>::_evaluate(const Argument<EffectiveFormula>& x) const -> Result<EffectiveFormula> {
    return this->_base_evaluate(x); }



} // namespace Ariadne
