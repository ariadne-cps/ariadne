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

#include "function_interface.h"

// Adaptors for classes to conform to the Function interface.

namespace Ariadne {

typedef Differential<Float> FloatDifferential;
typedef Differential<Interval> IntervalDifferential;
typedef TaylorModel<Float> FloatTaylorModel;
typedef TaylorModel<Interval> IntervalTaylorModel;
typedef Formula<Float> FloatFormula;
typedef Formula<Interval> IntervalFormula;
typedef Formula<Real> RealFormula;

template<class T, class X> class ScalarFunctionMixin;
template<class T, class X> class VectorFunctionMixin;

template<class F>
class ScalarFunctionMixin<F,Float>
    : public ScalarFunctionInterface<Float>
{
  private:
    template<class X> X _base_evaluate(const Vector<X>& x) const {
        X r; static_cast<const F*>(this)->_compute(r,x); return r;}
  protected:
    ScalarFunctionMixin() { }
  public:
    virtual Float evaluate(const Vector<Float>& x) const;
    virtual FloatDifferential evaluate(const Vector<FloatDifferential>& x) const;
    virtual FloatFormula evaluate(const Vector<FloatFormula>& x) const;
    virtual FloatTaylorModel evaluate(const Vector<FloatTaylorModel>& x) const;

    virtual std::ostream& repr(std::ostream& os) const { return this->write(os); }
};

// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F>
class ScalarFunctionMixin<F,Interval>
    : public ScalarFunctionInterface<Interval>
{
  private:
    template<class X> X _base_evaluate(const Vector<X>& x) const {
        X r; static_cast<const F*>(this)->_compute(r,x); return r;}
  protected:
    ScalarFunctionMixin() { }
  public:
    virtual Float evaluate(const Vector<Float>& x) const;
    virtual Interval evaluate(const Vector<Interval>& x) const;
    virtual FloatDifferential evaluate(const Vector<FloatDifferential>& x) const;
    virtual IntervalDifferential evaluate(const Vector<IntervalDifferential>& x) const;
    virtual FloatFormula evaluate(const Vector<FloatFormula>& x) const;
    virtual IntervalFormula evaluate(const Vector<IntervalFormula>& x) const;
    virtual FloatTaylorModel evaluate(const Vector<FloatTaylorModel>& x) const;
    virtual IntervalTaylorModel evaluate(const Vector<IntervalTaylorModel>& x) const;

    virtual std::ostream& repr(std::ostream& os) const { return this->write(os); }
};

// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F>
class ScalarFunctionMixin<F,Real>
    : public ScalarFunctionInterface<Real>
{
  private:
    template<class X> X _base_evaluate(const Vector<X>& x) const {
        X r; static_cast<const F*>(this)->_compute(r,x); return r;}
  protected:
    ScalarFunctionMixin() { }
  public:
    virtual Float evaluate(const Vector<Float>& x) const;
    virtual Interval evaluate(const Vector<Interval>& x) const;
    virtual Real evaluate(const Vector<Real>& x) const;
    virtual FloatDifferential evaluate(const Vector<FloatDifferential>& x) const;
    virtual IntervalDifferential evaluate(const Vector<IntervalDifferential>& x) const;
    virtual FloatFormula evaluate(const Vector<FloatFormula>& x) const;
    virtual IntervalFormula evaluate(const Vector<IntervalFormula>& x) const;
    virtual RealFormula evaluate(const Vector<RealFormula>& x) const;
    virtual FloatTaylorModel evaluate(const Vector<FloatTaylorModel>& x) const;
    virtual IntervalTaylorModel evaluate(const Vector<IntervalTaylorModel>& x) const;

    virtual std::ostream& repr(std::ostream& os) const { return this->write(os); }

    Vector<Float> gradient(const Vector<Float>& v) const;
    Vector<Interval> gradient(const Vector<Interval>& v) const;
};


template<class F>
class VectorFunctionMixin<F,Float>
    : public VectorFunctionInterface<Float>
{
  private:
    template<class X> Vector<X> _base_evaluate(const Vector<X>& a) const {
        Vector<X> r(this->result_size(),a[0]*0.0); static_cast<const F*>(this)->_compute(r,a); return r; }
  protected:
    VectorFunctionMixin() { }
  public:
    virtual Vector<Float> evaluate(const Vector<Float>& x) const;
    virtual Vector<FloatDifferential> evaluate(const Vector<FloatDifferential>& x) const;
    virtual Vector<FloatTaylorModel> evaluate(const Vector<FloatTaylorModel>& x) const;
    virtual Vector<FloatFormula> evaluate(const Vector<FloatFormula>& x) const;
};

template<class F>
class VectorFunctionMixin<F,Interval>
    : public VectorFunctionInterface<Interval>
{
  private:
    template<class X> Vector<X> _base_evaluate(const Vector<X>& a) const {
        Vector<X> r(this->result_size(),a[0]*0.0); static_cast<const F*>(this)->_compute(r,a); return r; }
  protected:
    VectorFunctionMixin() { }
  public:
    virtual Vector<Float> evaluate(const Vector<Float>& x) const;
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const;
    virtual Vector<FloatDifferential> evaluate(const Vector<FloatDifferential>& x) const;
    virtual Vector<IntervalDifferential> evaluate(const Vector<IntervalDifferential>& x) const;
    virtual Vector<FloatTaylorModel> evaluate(const Vector<FloatTaylorModel>& x) const;
    virtual Vector<IntervalTaylorModel> evaluate(const Vector<IntervalTaylorModel>& x) const;
    virtual Vector<FloatFormula> evaluate(const Vector<FloatFormula>& x) const;
    virtual Vector<IntervalFormula> evaluate(const Vector<IntervalFormula>& x) const;
};

// A wrapper for classes with non-static _compute methods
template<class F>
class VectorFunctionMixin<F,Real>
    : public VectorFunctionInterface<Real>
{
  private:
    template<class X> Vector<X> _base_evaluate(const Vector<X>& a) const {
        Vector<X> r(this->result_size(),a[0]*0.0); static_cast<const F*>(this)->_compute(r,a); return r; }
  protected:
    VectorFunctionMixin() { }
  public:
    virtual Vector<Float> evaluate(const Vector<Float>& x) const;
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const;
    virtual Vector<Real> evaluate(const Vector<Real>& x) const;
    virtual Vector<FloatDifferential> evaluate(const Vector<FloatDifferential>& x) const;
    virtual Vector<IntervalDifferential> evaluate(const Vector<IntervalDifferential>& x) const;
    virtual Vector<FloatTaylorModel> evaluate(const Vector<FloatTaylorModel>& x) const;
    virtual Vector<IntervalTaylorModel> evaluate(const Vector<IntervalTaylorModel>& x) const;
    virtual Vector<FloatFormula> evaluate(const Vector<FloatFormula>& x) const;
    virtual Vector<IntervalFormula> evaluate(const Vector<IntervalFormula>& x) const;
    virtual Vector<RealFormula> evaluate(const Vector<RealFormula>& x) const;
};


} // namespace Ariadne

#endif // ARIADNE_FUNCTION_TEMPLATE_H