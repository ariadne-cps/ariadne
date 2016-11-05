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

template<class F> class ScalarFunctionModelMixin<F,ValidatedTag>
    : public virtual ScalarFunctionModelInterface<ValidatedTag>
    , public ScalarFunctionMixin<F,ValidatedTag>
{
  public:
    F apply(OperatorCode op) const;
  public:
    ScalarFunctionModelInterface<ValidatedTag>* _clone() const {
        return new F(static_cast<const F&>(*this)); }
    NormType const _norm() const {
        return norm(static_cast<const F&>(*this)); }
    ScalarFunctionModelInterface<ValidatedTag>* _antiderivative(SizeType j) const {
        return new F(antiderivative(static_cast<const F&>(*this),j)); }
    ScalarFunctionModelInterface<ValidatedTag>* _antiderivative(SizeType j, ValidatedNumericType c) const {
        return new F(antiderivative(static_cast<const F&>(*this),j,c)); }
    ScalarFunctionModelInterface<ValidatedTag>* _apply(OperatorCode op) const {
        return new F(this->apply(op)); }
    ValidatedNumericType _unchecked_evaluate(const Vector<ValidatedNumericType>& x) const {
        return unchecked_evaluate(static_cast<const F&>(*this),x); }
    ScalarFunctionModelInterface<ValidatedTag>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const {
        return new F(embed(d1,static_cast<const F&>(*this),d2)); }
    Boolean _refines(const ScalarFunctionModelInterface<ValidatedTag>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return refines(static_cast<const F&>(*this),dynamic_cast<const F&>(f)); }
    Boolean _inconsistent(const ScalarFunctionModelInterface<ValidatedTag>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return inconsistent(static_cast<const F&>(*this),dynamic_cast<const F&>(f)); }
    ScalarFunctionModelInterface<ValidatedTag>* _refinement(const ScalarFunctionModelInterface<ValidatedTag>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return new F(refinement(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    Void _iadd(const ValidatedNumericType& c) {
        static_cast<F&>(*this)+=c; }
    Void _imul(const ValidatedNumericType& c) {
        static_cast<F&>(*this)*=c; }
    Void _isma(const ValidatedNumericType& c, const ScalarFunctionModelInterface<ValidatedTag>& f) {
        static_cast<F&>(*this)+=c*dynamic_cast<const F&>(f); }
    Void _ifma(const ScalarFunctionModelInterface<ValidatedTag>& f1, const ScalarFunctionModelInterface<ValidatedTag>& f2) {
        static_cast<F&>(*this)+=dynamic_cast<const F&>(f1)*dynamic_cast<const F&>(f2); }
};

template<class F> F ScalarFunctionModelMixin<F,ValidatedTag>::apply(OperatorCode op) const {
    const F& f=static_cast<const F&>(*this);
    switch(op) {
        case OperatorCode::NEG: return neg(f);
        case OperatorCode::REC: return rec(f);
        case OperatorCode::EXP: return exp(f);
        default: ARIADNE_FAIL_MSG("ScalarFunctionModel<ValidatedTag>::apply(OperatorCode op): Operator op="<<op<<" not implemented\n");
    }
}


template<class F> class VectorFunctionModelMixin<F,ValidatedTag>
    : public virtual VectorFunctionModelInterface<ValidatedTag>
    , public  VectorFunctionMixin<F,ValidatedTag>
{
    typedef typename Element<F>::Type ScalarFunctionType;
  public:
    virtual VectorFunctionModelInterface<ValidatedTag>* _clone() const { return new F(static_cast<const F&>(*this)); }
    virtual Void _set(SizeType i, const ScalarFunctionModelInterface<ValidatedTag>& sf) {
        if(!dynamic_cast<const typename F::ScalarFunctionType*>(&sf)) {
            ARIADNE_FAIL_MSG("Cannot set element of VectorFunctionModel "<<*this<<" to "<<sf<<"\n"); }
        static_cast<F&>(*this).F::set(i,dynamic_cast<const ScalarFunctionType&>(sf)); }
    virtual VectorFunctionModelInterface<ValidatedTag>* _derivative(SizeType j) const {
        ARIADNE_NOT_IMPLEMENTED; }
    NormType const _norm() const {
         return norm(static_cast<const F&>(*this)); }
    VectorFunctionModelInterface<ValidatedTag>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const {
        return heap_copy(embed(d1,static_cast<const F&>(*this),d2)); }
    Void _adjoin(const ScalarFunctionModelInterface<ValidatedTag>& f) {
        static_cast<F&>(*this).F::adjoin(dynamic_cast<const ScalarFunctionType&>(f)); }
    VectorFunctionModelInterface<ValidatedTag>* _join(const VectorFunctionModelInterface<ValidatedTag>& f) const {
        return heap_copy(join(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    VectorFunctionModelInterface<ValidatedTag>* _combine(const VectorFunctionModelInterface<ValidatedTag>& f) const {
        return heap_copy(combine(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    Vector<ValidatedNumericType> _unchecked_evaluate(const Vector<ValidatedNumericType>& x) const {
        return unchecked_evaluate(static_cast<const F&>(*this),x); }
    ScalarFunctionModelInterface<ValidatedTag>* _compose(const ScalarFunctionInterface<ValidatedTag>& f) const {
        return heap_copy(compose(f,static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<ValidatedTag>* _compose(const VectorFunctionInterface<ValidatedTag>& f) const {
        return heap_copy(compose(f,static_cast<const F&>(*this))); }
    ScalarFunctionModelInterface<ValidatedTag>* _unchecked_compose(const ScalarFunctionInterface<ValidatedTag>& f) const {
        return heap_copy(unchecked_compose(dynamic_cast<const ScalarFunctionType&>(f),static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<ValidatedTag>* _unchecked_compose(const VectorFunctionInterface<ValidatedTag>& f) const {
        return heap_copy(unchecked_compose(dynamic_cast<const F&>(f),static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<ValidatedTag>* _partial_evaluate(SizeType j, const ValidatedNumericType& c) const {
        return heap_copy(partial_evaluate(static_cast<const F&>(*this),j,c)); }
};

} // namespace Ariadne

#endif
