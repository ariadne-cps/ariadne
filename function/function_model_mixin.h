/***************************************************************************
 *            function_model_mixin.h
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

/*! \file function_model_mixin.h
 *  \brief Mixin for concrete functions on bounded domains.
 */

#ifndef ARIADNE_FUNCTION_MODEL_MIXIN_H
#define ARIADNE_FUNCTION_MODEL_MIXIN_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function/function_model_interface.h"
#include "function/function_model.h"

#include "numeric/operators.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/operations.h"
#include "geometry/box.h"

#include "function/function_interface.h"
#include "function/function_mixin.h"
#include "function/function.h"

namespace Ariadne {

template<class F, class P> class ScalarFunctionModelMixin;
template<class F, class P> class VectorFunctionModelMixin;

template<class F, class P> class ScalarFunctionModelMixin
    : public virtual ScalarFunctionModelInterface<P>
    , public ScalarFunctionMixin<F,P>
{
  public:
    F apply(OperatorCode op) const;
  public:
    ScalarFunctionModelInterface<P>* _clone() const override {
        return new F(static_cast<const F&>(*this)); }
    NormType const _norm() const override {
        return norm(static_cast<const F&>(*this)); }
    ScalarFunctionModelInterface<P>* _antiderivative(SizeType j) const override {
        return new F(antiderivative(static_cast<const F&>(*this),j)); }
    ScalarFunctionModelInterface<P>* _antiderivative(SizeType j, CanonicalNumericType<P> c) const override {
        return new F(antiderivative(static_cast<const F&>(*this),j,c)); }
     ScalarFunctionModelInterface<P>* _restriction(const ExactBoxType& d) const override {
        return new F(restriction(static_cast<const F&>(*this),d)); }
    ScalarFunctionModelInterface<P>* _apply(OperatorCode op) const override {
        return new F(this->apply(op)); }
    CanonicalNumericType<P> _unchecked_evaluate(const Vector<CanonicalNumericType<P>>& x) const override {
        return unchecked_evaluate(static_cast<const F&>(*this),x); }
    ScalarFunctionModelInterface<P>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const override {
        return new F(embed(d1,static_cast<const F&>(*this),d2)); }
    Boolean _refines(const ScalarFunctionModelInterface<P>& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return refines(static_cast<const F&>(*this),dynamic_cast<const F&>(f)); }
    Boolean _inconsistent(const ScalarFunctionModelInterface<P>& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return inconsistent(static_cast<const F&>(*this),dynamic_cast<const F&>(f)); }
    ScalarFunctionModelInterface<P>* _refinement(const ScalarFunctionModelInterface<P>& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return new F(refinement(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    Void _iadd(const CanonicalNumericType<P>& c) override {
        static_cast<F&>(*this)+=c; }
    Void _imul(const CanonicalNumericType<P>& c) override {
        static_cast<F&>(*this)*=c; }
    Void _isma(const CanonicalNumericType<P>& c, const ScalarFunctionModelInterface<P>& f) override {
        static_cast<F&>(*this)+=c*dynamic_cast<const F&>(f); }
    Void _ifma(const ScalarFunctionModelInterface<P>& f1, const ScalarFunctionModelInterface<P>& f2) override {
        static_cast<F&>(*this)+=dynamic_cast<const F&>(f1)*dynamic_cast<const F&>(f2); }
};

template<class F, class P> F ScalarFunctionModelMixin<F,P>::apply(OperatorCode op) const {
    const F& f=static_cast<const F&>(*this);
    switch(op) {
        case OperatorCode::NEG: return neg(f);
        case OperatorCode::REC: return rec(f);
        case OperatorCode::EXP: return exp(f);
        default: ARIADNE_FAIL_MSG("ScalarFunctionModel<P>::apply(OperatorCode op): Operator op="<<op<<" not implemented\n");
    }
}


template<class F, class P> class VectorFunctionModelMixin
    : public virtual VectorFunctionModelInterface<P>
    , public  VectorFunctionMixin<F,P>
{
    typedef typename Element<F>::Type ScalarFunctionType;
  public:
    virtual VectorFunctionModelInterface<P>* _clone() const override { return new F(static_cast<const F&>(*this)); }
    virtual Void _set(SizeType i, const ScalarFunctionModelInterface<P>& sf) override {
        if(!dynamic_cast<const typename F::ScalarFunctionType*>(&sf)) {
            ARIADNE_FAIL_MSG("Cannot set element of VectorFunctionModel "<<*this<<" to "<<sf<<"\n"); }
        static_cast<F&>(*this).F::set(i,dynamic_cast<const ScalarFunctionType&>(sf)); }
    virtual VectorFunctionModelInterface<P>* _derivative(SizeType j) const override {
        ARIADNE_NOT_IMPLEMENTED; }
    NormType const _norm() const override {
         return norm(static_cast<const F&>(*this)); }
    VectorFunctionModelInterface<P>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const override {
        return heap_copy(embed(d1,static_cast<const F&>(*this),d2)); }
    Void _adjoin(const ScalarFunctionModelInterface<P>& f) override {
        static_cast<F&>(*this).F::adjoin(dynamic_cast<const ScalarFunctionType&>(f)); }
    VectorFunctionModelInterface<P>* _join(const VectorFunctionModelInterface<P>& f) const override {
        return heap_copy(join(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    VectorFunctionModelInterface<P>* _combine(const VectorFunctionModelInterface<P>& f) const override {
        return heap_copy(combine(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    Vector<CanonicalNumericType<P>> _unchecked_evaluate(const Vector<CanonicalNumericType<P>>& x) const override {
        return unchecked_evaluate(static_cast<const F&>(*this),x); }
    ScalarFunctionModelInterface<P>* _compose(const ScalarFunctionInterface<P>& f) const override {
        return heap_copy(compose(f,static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<P>* _compose(const VectorFunctionInterface<P>& f) const override {
        return heap_copy(compose(f,static_cast<const F&>(*this))); }
    ScalarFunctionModelInterface<P>* _unchecked_compose(const ScalarFunctionInterface<P>& f) const override {
        return heap_copy(unchecked_compose(dynamic_cast<const ScalarFunctionType&>(f),static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<P>* _unchecked_compose(const VectorFunctionInterface<P>& f) const override {
        return heap_copy(unchecked_compose(dynamic_cast<const F&>(f),static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<P>* _partial_evaluate(SizeType j, const CanonicalNumericType<P>& c) const override {
        return heap_copy(partial_evaluate(static_cast<const F&>(*this),j,c)); }
};

} // namespace Ariadne

#endif
