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

template<class F> Float ScalarFunctionMixin<F,Float>::evaluate(const Vector<Float>& x) const { return this->_base_evaluate(x); }
template<class F> FloatDifferential ScalarFunctionMixin<F,Float>::evaluate(const Vector<FloatDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> FloatTaylorModel ScalarFunctionMixin<F,Float>::evaluate(const Vector<FloatTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> FloatFormula ScalarFunctionMixin<F,Float>::evaluate(const Vector<FloatFormula>& x) const { return this->_base_evaluate(x); }
template<class F> FloatAlgebra ScalarFunctionMixin<F,Float>::evaluate(const Vector<FloatAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> FloatScalarFunctionInterface* ScalarFunctionMixin<F,Float>::_clone() const { return new F(static_cast<const F&>(*this)); }

template<class F> Float ScalarFunctionMixin<F,Interval>::evaluate(const Vector<Float>& x) const { return this->_base_evaluate(x); }
template<class F> Interval ScalarFunctionMixin<F,Interval>::evaluate(const Vector<Interval>& x) const { return this->_base_evaluate(x); }
template<class F> FloatDifferential ScalarFunctionMixin<F,Interval>::evaluate(const Vector<FloatDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> IntervalDifferential ScalarFunctionMixin<F,Interval>::evaluate(const Vector<IntervalDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> FloatTaylorModel ScalarFunctionMixin<F,Interval>::evaluate(const Vector<FloatTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> IntervalTaylorModel ScalarFunctionMixin<F,Interval>::evaluate(const Vector<IntervalTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> FloatFormula ScalarFunctionMixin<F,Interval>::evaluate(const Vector<FloatFormula>& x) const { return this->_base_evaluate(x); }
template<class F> IntervalFormula ScalarFunctionMixin<F,Interval>::evaluate(const Vector<IntervalFormula>& x) const { return this->_base_evaluate(x); }
template<class F> FloatAlgebra ScalarFunctionMixin<F,Interval>::evaluate(const Vector<FloatAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> IntervalAlgebra ScalarFunctionMixin<F,Interval>::evaluate(const Vector<IntervalAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> IntervalScalarFunctionInterface* ScalarFunctionMixin<F,Interval>::_clone() const { return new F(static_cast<const F&>(*this)); }

template<class F> Float ScalarFunctionMixin<F,Real>::evaluate(const Vector<Float>& x) const { return this->_base_evaluate(x); }
template<class F> Interval ScalarFunctionMixin<F,Real>::evaluate(const Vector<Interval>& x) const { return this->_base_evaluate(x); }
template<class F> Real ScalarFunctionMixin<F,Real>::evaluate(const Vector<Real>& x) const { return this->_base_evaluate(x); }
template<class F> FloatDifferential ScalarFunctionMixin<F,Real>::evaluate(const Vector<FloatDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> IntervalDifferential ScalarFunctionMixin<F,Real>::evaluate(const Vector<IntervalDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> FloatTaylorModel ScalarFunctionMixin<F,Real>::evaluate(const Vector<FloatTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> IntervalTaylorModel ScalarFunctionMixin<F,Real>::evaluate(const Vector<IntervalTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> FloatFormula ScalarFunctionMixin<F,Real>::evaluate(const Vector<FloatFormula>& x) const { return this->_base_evaluate(x); }
template<class F> IntervalFormula ScalarFunctionMixin<F,Real>::evaluate(const Vector<IntervalFormula>& x) const { return this->_base_evaluate(x); }
template<class F> RealFormula ScalarFunctionMixin<F,Real>::evaluate(const Vector<RealFormula>& x) const { return this->_base_evaluate(x); }
template<class F> FloatAlgebra ScalarFunctionMixin<F,Real>::evaluate(const Vector<FloatAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> IntervalAlgebra ScalarFunctionMixin<F,Real>::evaluate(const Vector<IntervalAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> RealAlgebra ScalarFunctionMixin<F,Real>::evaluate(const Vector<RealAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> RealScalarFunctionInterface* ScalarFunctionMixin<F,Real>::_clone() const { return new F(static_cast<const F&>(*this)); }

template<class F> Vector<Float> ScalarFunctionMixin<F,Real>::gradient(const Vector<Float>& x) const {
    return this->_base_evaluate(FloatDifferential::variables(1u,x)).gradient(); }
template<class F> Vector<Interval> ScalarFunctionMixin<F,Real>::gradient(const Vector<Interval>& x) const {
    return this->_base_evaluate(IntervalDifferential::variables(1u,x)).gradient(); }

template<class F> Vector<Float> VectorFunctionMixin<F,Float>::evaluate(const Vector<Float>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<FloatDifferential> VectorFunctionMixin<F,Float>::evaluate(const Vector<FloatDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<FloatTaylorModel> VectorFunctionMixin<F,Float>::evaluate(const Vector<FloatTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<FloatFormula> VectorFunctionMixin<F,Float>::evaluate(const Vector<FloatFormula>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<FloatAlgebra> VectorFunctionMixin<F,Float>::evaluate(const Vector<FloatAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> FloatVectorFunctionInterface* VectorFunctionMixin<F,Float>::_clone() const { return new F(static_cast<const F&>(*this)); }

template<class F> Vector<Float> VectorFunctionMixin<F,Interval>::evaluate(const Vector<Float>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<Interval> VectorFunctionMixin<F,Interval>::evaluate(const Vector<Interval>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<FloatDifferential> VectorFunctionMixin<F,Interval>::evaluate(const Vector<FloatDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<IntervalDifferential> VectorFunctionMixin<F,Interval>::evaluate(const Vector<IntervalDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<FloatTaylorModel> VectorFunctionMixin<F,Interval>::evaluate(const Vector<FloatTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<IntervalTaylorModel> VectorFunctionMixin<F,Interval>::evaluate(const Vector<IntervalTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<FloatFormula> VectorFunctionMixin<F,Interval>::evaluate(const Vector<FloatFormula>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<IntervalFormula> VectorFunctionMixin<F,Interval>::evaluate(const Vector<IntervalFormula>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<FloatAlgebra> VectorFunctionMixin<F,Interval>::evaluate(const Vector<FloatAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<IntervalAlgebra> VectorFunctionMixin<F,Interval>::evaluate(const Vector<IntervalAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> IntervalVectorFunctionInterface* VectorFunctionMixin<F,Interval>::_clone() const { return new F(static_cast<const F&>(*this)); }

template<class F> Vector<Float> VectorFunctionMixin<F,Real>::evaluate(const Vector<Float>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<Interval> VectorFunctionMixin<F,Real>::evaluate(const Vector<Interval>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<Real> VectorFunctionMixin<F,Real>::evaluate(const Vector<Real>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<FloatDifferential> VectorFunctionMixin<F,Real>::evaluate(const Vector<FloatDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<IntervalDifferential> VectorFunctionMixin<F,Real>::evaluate(const Vector<IntervalDifferential>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<FloatTaylorModel> VectorFunctionMixin<F,Real>::evaluate(const Vector<FloatTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<IntervalTaylorModel> VectorFunctionMixin<F,Real>::evaluate(const Vector<IntervalTaylorModel>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<FloatFormula> VectorFunctionMixin<F,Real>::evaluate(const Vector<FloatFormula>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<IntervalFormula> VectorFunctionMixin<F,Real>::evaluate(const Vector<IntervalFormula>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<RealFormula> VectorFunctionMixin<F,Real>::evaluate(const Vector<RealFormula>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<FloatAlgebra> VectorFunctionMixin<F,Real>::evaluate(const Vector<FloatAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<IntervalAlgebra> VectorFunctionMixin<F,Real>::evaluate(const Vector<IntervalAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> Vector<RealAlgebra> VectorFunctionMixin<F,Real>::evaluate(const Vector<RealAlgebra>& x) const { return this->_base_evaluate(x); }
template<class F> RealVectorFunctionInterface* VectorFunctionMixin<F,Real>::_clone() const { return new F(static_cast<const F&>(*this)); }


} // namespace Ariadne
