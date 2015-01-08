/***************************************************************************
 *            function_mixin.h
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

#ifndef ARIADNE_FUNCTION_MIXIN_H
#define ARIADNE_FUNCTION_MIXIN_H

#include "function/function_interface.h"

// Adaptors for classes to conform to the Function interface.

namespace Ariadne {

typedef ApproximateNumber ApproximateNumber;
typedef ValidatedNumber ValidatedNumber;
typedef EffectiveNumber EffectiveNumber;
typedef Differential<ApproximateNumber> ApproximateDifferential;
typedef Differential<ValidatedNumber> ValidatedDifferential;
typedef TaylorModel<ApproximateNumber> ApproximateTaylorModel;
typedef TaylorModel<ValidatedNumber> ValidatedTaylorModel;
typedef Formula<ApproximateNumber> ApproximateFormula;
typedef Formula<ValidatedNumber> ValidatedFormula;
typedef Formula<EffectiveNumber> EffectiveFormula;
typedef Algebra<ApproximateNumber> ApproximateAlgebra;
typedef Algebra<ValidatedNumber> ValidatedAlgebra;
typedef Algebra<EffectiveNumber> EffectiveAlgebra;

template<class T, class X> class ScalarFunctionMixin { };
template<class T, class X> class VectorFunctionMixin { };

template<class T> T* heap_copy(const T& t) { return new T(t); }

template<class F>
class ScalarFunctionMixin<F,ApproximateTag>
    : public virtual ScalarFunctionInterface<ApproximateTag>
{
  private:
    template<class X> X _base_evaluate(const Vector<X>& x) const {
        X r; static_cast<const F*>(this)->_compute(r,x); return r;}
  protected:
    ScalarFunctionMixin() { }
  public:
    virtual ApproximateNumber evaluate(const Vector<ApproximateNumber>& x) const;
    virtual ApproximateDifferential evaluate(const Vector<ApproximateDifferential>& x) const;
    virtual ApproximateFormula evaluate(const Vector<ApproximateFormula>& x) const;
    virtual ApproximateTaylorModel evaluate(const Vector<ApproximateTaylorModel>& x) const;
    virtual ApproximateAlgebra evaluate(const Vector<ApproximateAlgebra>& x) const;

    virtual OutputStream& repr(OutputStream& os) const { return this->write(os); }
    virtual ScalarFunctionInterface<ApproximateTag>* _clone() const;
};

// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F>
class ScalarFunctionMixin<F,ValidatedTag>
    : public virtual ScalarFunctionInterface<ValidatedTag>
{
  private:
    template<class XX> XX _base_evaluate(const Vector<XX>& x) const {
        XX r; static_cast<const F*>(this)->_compute(r,x); return r;}
  protected:
    ScalarFunctionMixin() { }
  public:
    virtual ApproximateNumber evaluate(const Vector<ApproximateNumber>& x) const;
    virtual ValidatedNumber evaluate(const Vector<ValidatedNumber>& x) const;
    virtual ApproximateDifferential evaluate(const Vector<ApproximateDifferential>& x) const;
    virtual ValidatedDifferential evaluate(const Vector<ValidatedDifferential>& x) const;
    virtual ApproximateFormula evaluate(const Vector<ApproximateFormula>& x) const;
    virtual ApproximateAlgebra evaluate(const Vector<ApproximateAlgebra>& x) const;
    virtual ValidatedFormula evaluate(const Vector<ValidatedFormula>& x) const;
    virtual ApproximateTaylorModel evaluate(const Vector<ApproximateTaylorModel>& x) const;
    virtual ValidatedTaylorModel evaluate(const Vector<ValidatedTaylorModel>& x) const;
    virtual ValidatedAlgebra evaluate(const Vector<ValidatedAlgebra>& x) const;

    virtual OutputStream& repr(OutputStream& os) const { return this->write(os); }
    virtual ScalarFunctionInterface<ValidatedTag>* _clone() const;
};

// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F>
class ScalarFunctionMixin<F,EffectiveTag>
    : public virtual ScalarFunctionInterface<EffectiveTag>
{
  private:
    template<class X> X _base_evaluate(const Vector<X>& x) const {
        X r; static_cast<const F*>(this)->_compute(r,x); return r;}
  protected:
    ScalarFunctionMixin() { }
  public:
    virtual ApproximateNumber evaluate(const Vector<ApproximateNumber>& x) const;
    virtual ValidatedNumber evaluate(const Vector<ValidatedNumber>& x) const;
    virtual EffectiveNumber evaluate(const Vector<EffectiveNumber>& x) const;
    virtual ApproximateDifferential evaluate(const Vector<ApproximateDifferential>& x) const;
    virtual ValidatedDifferential evaluate(const Vector<ValidatedDifferential>& x) const;
    virtual ApproximateFormula evaluate(const Vector<ApproximateFormula>& x) const;
    virtual ValidatedFormula evaluate(const Vector<ValidatedFormula>& x) const;
    virtual EffectiveFormula evaluate(const Vector<EffectiveFormula>& x) const;
    virtual ApproximateTaylorModel evaluate(const Vector<ApproximateTaylorModel>& x) const;
    virtual ValidatedTaylorModel evaluate(const Vector<ValidatedTaylorModel>& x) const;
    virtual ApproximateAlgebra evaluate(const Vector<ApproximateAlgebra>& x) const;
    virtual ValidatedAlgebra evaluate(const Vector<ValidatedAlgebra>& x) const;
    virtual EffectiveAlgebra evaluate(const Vector<EffectiveAlgebra>& x) const;

    virtual OutputStream& repr(OutputStream& os) const { return this->write(os); }
    virtual ScalarFunctionInterface<EffectiveTag>* _clone() const;

    Vector<ApproximateNumber> gradient(const Vector<ApproximateNumber>& v) const;
    Vector<ValidatedNumber> gradient(const Vector<ValidatedNumber>& v) const;
};


template<class F>
class VectorFunctionMixin<F,ApproximateTag>
    : public virtual VectorFunctionInterface<ApproximateTag>
{
  private:
    template<class X> Vector<X> _base_evaluate(const Vector<X>& a) const {
        Vector<X> r(this->result_size(),a.zero_element()); static_cast<const F*>(this)->_compute(r,a); return r; }
  protected:
    VectorFunctionMixin() { }
  public:
    virtual Vector<ApproximateNumber> evaluate(const Vector<ApproximateNumber>& x) const;
    virtual Vector<ApproximateDifferential> evaluate(const Vector<ApproximateDifferential>& x) const;
    virtual Vector<ApproximateTaylorModel> evaluate(const Vector<ApproximateTaylorModel>& x) const;
    virtual Vector<ApproximateFormula> evaluate(const Vector<ApproximateFormula>& x) const;
    virtual Vector<ApproximateAlgebra> evaluate(const Vector<ApproximateAlgebra>& x) const;
    virtual VectorFunctionInterface<ApproximateTag>* _clone() const;

    virtual OutputStream& repr(OutputStream& os) const { return this->write(os); }
};

template<class F>
class VectorFunctionMixin<F,ValidatedTag>
    : public virtual VectorFunctionInterface<ValidatedTag>
{
  private:
    template<class X> Vector<X> _base_evaluate(const Vector<X>& a) const {
        Vector<X> r(this->result_size(),a.zero_element()); static_cast<const F*>(this)->_compute(r,a); return r; }
  protected:
    VectorFunctionMixin() { }
  public:
    virtual Vector<ApproximateNumber> evaluate(const Vector<ApproximateNumber>& x) const;
    virtual Vector<ValidatedNumber> evaluate(const Vector<ValidatedNumber>& x) const;
    virtual Vector<ApproximateDifferential> evaluate(const Vector<ApproximateDifferential>& x) const;
    virtual Vector<ValidatedDifferential> evaluate(const Vector<ValidatedDifferential>& x) const;
    virtual Vector<ApproximateTaylorModel> evaluate(const Vector<ApproximateTaylorModel>& x) const;
    virtual Vector<ValidatedTaylorModel> evaluate(const Vector<ValidatedTaylorModel>& x) const;
    virtual Vector<ApproximateFormula> evaluate(const Vector<ApproximateFormula>& x) const;
    virtual Vector<ValidatedFormula> evaluate(const Vector<ValidatedFormula>& x) const;
    virtual Vector<ApproximateAlgebra> evaluate(const Vector<ApproximateAlgebra>& x) const;
    virtual Vector<ValidatedAlgebra> evaluate(const Vector<ValidatedAlgebra>& x) const;
    virtual VectorFunctionInterface<ValidatedTag>* _clone() const;

    virtual OutputStream& repr(OutputStream& os) const { return this->write(os); }
};

// A wrapper for classes with non-static _compute methods
template<class F>
class VectorFunctionMixin<F,EffectiveTag>
    : public virtual VectorFunctionInterface<EffectiveTag>
{
  private:
    template<class X> Vector<X> _base_evaluate(const Vector<X>& a) const {
        Vector<X> r(this->result_size(),a.zero_element()); static_cast<const F*>(this)->_compute(r,a); return r; }
  protected:
    VectorFunctionMixin() { }
  public:
    virtual Vector<ApproximateNumber> evaluate(const Vector<ApproximateNumber>& x) const;
    virtual Vector<ValidatedNumber> evaluate(const Vector<ValidatedNumber>& x) const;
    virtual Vector<EffectiveNumber> evaluate(const Vector<EffectiveNumber>& x) const;
    virtual Vector<ApproximateDifferential> evaluate(const Vector<ApproximateDifferential>& x) const;
    virtual Vector<ValidatedDifferential> evaluate(const Vector<ValidatedDifferential>& x) const;
    virtual Vector<ApproximateTaylorModel> evaluate(const Vector<ApproximateTaylorModel>& x) const;
    virtual Vector<ValidatedTaylorModel> evaluate(const Vector<ValidatedTaylorModel>& x) const;
    virtual Vector<ApproximateFormula> evaluate(const Vector<ApproximateFormula>& x) const;
    virtual Vector<ValidatedFormula> evaluate(const Vector<ValidatedFormula>& x) const;
    virtual Vector<EffectiveFormula> evaluate(const Vector<EffectiveFormula>& x) const;
    virtual Vector<ApproximateAlgebra> evaluate(const Vector<ApproximateAlgebra>& x) const;
    virtual Vector<ValidatedAlgebra> evaluate(const Vector<ValidatedAlgebra>& x) const;
    virtual Vector<EffectiveAlgebra> evaluate(const Vector<EffectiveAlgebra>& x) const;
    virtual VectorFunctionInterface<EffectiveTag>* _clone() const;

    virtual OutputStream& repr(OutputStream& os) const { return this->write(os); }

};


} // namespace Ariadne

#endif // ARIADNE_FUNCTION_TEMPLATE_H
