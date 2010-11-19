/***************************************************************************
 *            function.cc
 *
 *  Copyright 2008  Pieter Collins
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

#ifndef ARIADNE_FUNCTION_TEMPLATE_H
#define ARIADNE_FUNCTION_TEMPLATE_H

#include "function.h"
#include "polynomial.h"
#include "differential.h"
#include "taylor_model.h"
#include "formula.h"

namespace Ariadne {



// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F>
class IntervalScalarFunctionTemplate
    : public ScalarFunctionInterface<Interval>
{
  private:
    template<class R, class A> void _base_compute(R& r, const A& a) const {
        static_cast<const F*>(this)->_compute(r,a); }
  protected:
    IntervalScalarFunctionTemplate() { }
  public:
    virtual Float evaluate(const Vector<Float>& x) const {
        Float r; _base_compute(r,x); return r; }
    virtual Interval evaluate(const Vector<Interval>& x) const {
        Interval r; _base_compute(r,x); return r; }

    virtual TaylorModel<Float> evaluate(const Vector< TaylorModel<Float> >& x) const {
        TaylorModel<Float> r(TaylorModel<Float>(x[0].argument_size(),x[0].accuracy_ptr()));
        _base_compute(r,x); return r; }
    virtual TaylorModel<Interval> evaluate(const Vector< TaylorModel<Interval> >& x) const {
        TaylorModel<Interval> r(TaylorModel<Interval>(x[0].argument_size(),x[0].accuracy_ptr()));
        _base_compute(r,x); return r; }

    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const {
        Differential<Float> r(Differential<Float>(x[0].argument_size(),x[0].degree()));
        _base_compute(r,x); return r; }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const {
        Differential<Interval> r(Differential<Interval>(x[0].argument_size(),x[0].degree()));
        _base_compute(r,x); return r; }

    virtual Formula<Interval> evaluate(const Vector< Formula<Interval> >& x) const {
        Formula<Interval> r; _base_compute(r,x); return r; }

    virtual std::ostream& repr(std::ostream& os) const {
        return this->write(os); }
};

// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F>
class ScalarFunctionTemplate
    : public ScalarFunctionInterface<Real>
{
  private:
    template<class R, class A> void _base_compute(R& r, const A& a) const {
        static_cast<const F*>(this)->_compute(r,a); }
  protected:
    ScalarFunctionTemplate() { }
  public:
    virtual Float evaluate(const Vector<Float>& x) const {
        Float r; _base_compute(r,x); return r; }
    virtual Interval evaluate(const Vector<Interval>& x) const {
        Interval r; _base_compute(r,x); return r; }
    virtual Real evaluate(const Vector<Real>& x) const {
        Real r; _base_compute(r,x); return r; }

    virtual TaylorModel<Float> evaluate(const Vector< TaylorModel<Float> >& x) const {
        TaylorModel<Float> r(TaylorModel<Float>(x[0].argument_size(),x[0].accuracy_ptr()));
        _base_compute(r,x); return r; }
    virtual TaylorModel<Interval> evaluate(const Vector< TaylorModel<Interval> >& x) const {
        TaylorModel<Interval> r(TaylorModel<Interval>(x[0].argument_size(),x[0].accuracy_ptr()));
        _base_compute(r,x); return r; }

    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const {
        Differential<Float> r(Differential<Float>(x[0].argument_size(),x[0].degree()));
        _base_compute(r,x); return r; }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const {
        Differential<Interval> r(Differential<Interval>(x[0].argument_size(),x[0].degree()));
        _base_compute(r,x); return r; }

    virtual Formula<Interval> evaluate(const Vector< Formula<Interval> >& x) const {
        Formula<Interval> r; _base_compute(r,x); return r; }

    virtual std::ostream& repr(std::ostream& os) const {
        return this->write(os); }
};


// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F>
class IntervalVectorFunctionTemplate
    : public VectorFunctionInterface<Interval>
{
  private:
    template<class R, class A> void _base_compute(R& r, const A& a) const {
        static_cast<const F*>(this)->_compute(r,a); }
  protected:
    IntervalVectorFunctionTemplate() { }
  public:
    virtual Vector<Float> evaluate(const Vector<Float>& x) const {
        Vector<Float> r(this->result_size()); _base_compute(r,x); return r; }
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
        Vector<Interval> r(this->result_size()); _base_compute(r,x); return r; }

    virtual Vector< TaylorModel<Float> > evaluate(const Vector< TaylorModel<Float> >& x) const {
        Vector< TaylorModel<Float> > r(this->result_size(),TaylorModel<Float>(x[0].argument_size(),x[0].accuracy_ptr()));
        _base_compute(r,x); return r; }
    virtual Vector< TaylorModel<Interval> > evaluate(const Vector< TaylorModel<Interval> >& x) const {
        Vector< TaylorModel<Interval> > r(this->result_size(),TaylorModel<Interval>(x[0].argument_size(),x[0].accuracy_ptr()));
        _base_compute(r,x); return r; }

    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const {
        Vector< Differential<Float> > r(this->result_size(),Differential<Float>(x[0].argument_size(),x[0].degree()));
        _base_compute(r,x); return r; }
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const {
        Vector< Differential<Interval> > r(this->result_size(),Differential<Interval>(x[0].argument_size(),x[0].degree()));
        _base_compute(r,x); return r; }

    virtual Vector< Formula<Interval> > evaluate(const Vector< Formula<Interval> >& x) const {
        Vector< Formula<Interval> > r(this->result_size(),Formula<Interval>()); _base_compute(r,x); return r; }

};

// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F>
class VectorFunctionTemplate
    : public VectorFunctionInterface<Real>
{
  private:
    template<class R, class A> void _base_compute(R& r, const A& a) const {
        static_cast<const F*>(this)->_compute(r,a); }
  protected:
    VectorFunctionTemplate() { }
  public:
    virtual Vector<Float> evaluate(const Vector<Float>& x) const {
        Vector<Float> r(this->result_size()); _base_compute(r,x); return r; }
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
        Vector<Interval> r(this->result_size()); _base_compute(r,x); return r; }
    virtual Vector<Real> evaluate(const Vector<Real>& x) const {
        Vector<Real> r(this->result_size()); _base_compute(r,x); return r; }

    virtual Vector< TaylorModel<Float> > evaluate(const Vector< TaylorModel<Float> >& x) const {
        Vector< TaylorModel<Float> > r(this->result_size(),TaylorModel<Float>(x[0].argument_size(),x[0].accuracy_ptr()));
        _base_compute(r,x); return r; }
    virtual Vector< TaylorModel<Interval> > evaluate(const Vector< TaylorModel<Interval> >& x) const {
        Vector< TaylorModel<Interval> > r(this->result_size(),TaylorModel<Interval>(x[0].argument_size(),x[0].accuracy_ptr()));
        _base_compute(r,x); return r; }

    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const {
        Vector< Differential<Float> > r(this->result_size(),Differential<Float>(x[0].argument_size(),x[0].degree()));
        _base_compute(r,x); return r; }
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const {
        Vector< Differential<Interval> > r(this->result_size(),Differential<Interval>(x[0].argument_size(),x[0].degree()));
        _base_compute(r,x); return r; }

    virtual Vector< Formula<Interval> > evaluate(const Vector< Formula<Interval> >& x) const {
        Vector< Formula<Interval> > r(this->result_size(),Formula<Interval>()); _base_compute(r,x); return r; }

};


} // namespace Ariadne

#endif // ARIADNE_FUNCTION_TEMPLATE_H