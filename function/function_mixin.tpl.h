/***************************************************************************
 *            function_mixin.tpl.h
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
 
#ifndef ARIADNE_FUNCTION_MIXIN_TCC
#define ARIADNE_FUNCTION_MIXIN_TCC

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/differential.h"
#include "function/taylor_model.h"
#include "function/formula.h"
#include "algebra/algebra.h"

#include "function/function_mixin.h"

namespace Ariadne {

template<class X> Scalar<X> create_result(SizeOne n, X z) { return z; }
template<class X> Vector<X> create_result(SizeType n, X z) { return Vector<X>(n,z); }

template<class F, class D, class C> template<class X> auto
FunctionMixin<F,Void,D,C>::_base_evaluate(const ElementType<D,X>& x) const -> ElementType<C,X> {
    ElementType<C,X> r=create_result<X>(this->codomain().dimension(),zero_element(x));
    static_cast<const F*>(this)->_compute(r,x); return std::move(r);
}

template<class F, class D> template<class X> auto
FunctionMixin<F,Void,D,IntervalDomain>::_base_evaluate(const ElementType<D,X>& x) const -> X {
    ElementType<C,X> r=create_result<X>(this->codomain().dimension(),zero_element(x));
    static_cast<const F*>(this)->_compute(r,x); return std::move(r);
}

template<class F,class D, class C> FunctionInterface<ApproximateTag,D,C>* FunctionMixin<F,ApproximateTag,D,C>::_clone() const {
    return new F(static_cast<const F&>(*this)); }
template<class F,class D, class C> FunctionInterface<ValidatedTag,D,C>* FunctionMixin<F,ValidatedTag,D,C>::_clone() const {
    return new F(static_cast<const F&>(*this)); }
template<class F,class D, class C> FunctionInterface<EffectiveTag,D,C>* FunctionMixin<F,EffectiveTag,D,C>::_clone() const {
    return new F(static_cast<const F&>(*this)); }

template<class F,class D, class C> auto
FunctionMixin<F,ApproximateTag,D,C>::_evaluate(const Argument<Float64Approximation>& x) const -> Result<Float64Approximation> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ApproximateTag,D,C>::_evaluate(const Argument<FloatMPApproximation>& x) const -> Result<FloatMPApproximation> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ApproximateTag,D,C>::_evaluate(const Argument<Differential<Float64Approximation>>& x) const -> Result<Differential<Float64Approximation>> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ApproximateTag,D,C>::_evaluate(const Argument<Differential<FloatMPApproximation>>& x) const -> Result<Differential<FloatMPApproximation>> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ApproximateTag,D,C>::_evaluate(const Argument<TaylorModel<ApproximateTag,Float64>>& x) const -> Result<TaylorModel<ApproximateTag,Float64>> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ApproximateTag,D,C>::_evaluate(const Argument<TaylorModel<ApproximateTag,FloatMP>>& x) const -> Result<TaylorModel<ApproximateTag,FloatMP>> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ApproximateTag,D,C>::_evaluate(const Argument<Algebra<ApproximateNumber>>& x) const -> Result<Algebra<ApproximateNumber>> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ApproximateTag,D,C>::_evaluate(const Argument<Formula<ApproximateNumber>>& x) const -> Result<Formula<ApproximateNumber>> {
    return this->_base_evaluate(x); }

template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<Float64Bounds>& x) const -> Result<Float64Bounds> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<FloatMPBounds>& x) const -> Result<FloatMPBounds> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<Differential<Float64Bounds>>& x) const -> Result<Differential<Float64Bounds>> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<Differential<FloatMPBounds>>& x) const -> Result<Differential<FloatMPBounds>> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<TaylorModel<ValidatedTag,Float64>>& x) const -> Result<TaylorModel<ValidatedTag,Float64>> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<TaylorModel<ValidatedTag,FloatMP>>& x) const -> Result<TaylorModel<ValidatedTag,FloatMP>> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<Algebra<ValidatedNumber>>& x) const -> Result<Algebra<ValidatedNumber>> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<Formula<ValidatedNumber>>& x) const -> Result<Formula<ValidatedNumber>> {
    return this->_base_evaluate(x); }

template<class F,class D, class C> auto
FunctionMixin<F,ValidatedTag,D,C>::_evaluate(const Argument<ValidatedScalarFunction>& x) const -> Result<ValidatedScalarFunction> {
    return this->_base_evaluate(x); }

template<class F,class D, class C> auto
FunctionMixin<F,EffectiveTag,D,C>::_evaluate(const Argument<Real>& x) const -> Result<Real> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,EffectiveTag,D,C>::_evaluate(const Argument<Algebra<Real>>& x) const -> Result<Algebra<Real>> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,EffectiveTag,D,C>::_evaluate(const Argument<Formula<Real>>& x) const -> Result<Formula<Real>> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,EffectiveTag,D,C>::_evaluate(const Argument<Algebra<EffectiveNumber>>& x) const -> Result<Algebra<EffectiveNumber>> {
    return this->_base_evaluate(x); }
template<class F,class D, class C> auto
FunctionMixin<F,EffectiveTag,D,C>::_evaluate(const Argument<Formula<EffectiveNumber>>& x) const -> Result<Formula<EffectiveNumber>> {
    return this->_base_evaluate(x); }



} // namespace Ariadne

#endif
