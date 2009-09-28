/***************************************************************************
 *            function_set.cc
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

#include "macros.h"
#include "function.h"
#include "function_set.h"
#include "geometry.h"

namespace Ariadne {


//! \brief A scaling function \f$x_i' = o_i+l_ix_i\f$.
class VectorScalingFunction
    : public VectorFunctionInterface
{
  public:
    //! \brief The scaling function \f$x_i' = o_i+l_ix_i\f$.
    explicit VectorScalingFunction() : _o(), _l() { }
    explicit VectorScalingFunction(const Vector<Float>& origin,
                             const Vector<Float>& lengths)
        : _o(origin), _l(lengths) { ARIADNE_ASSERT(origin.size()==lengths.size()); }
    //! \brief The scaling function which takes the unit interval \f$[-1,+1]^n\f$ into \a range.
    explicit VectorScalingFunction(const Vector<Interval>& range)
        : _o(midpoint(range)), _l(range.size()) { for(uint i=0; i!=_l.size(); ++i) { _l[i]=range[i].radius(); } }
    const Vector<Float>& origin() const { return _o; }
    const Vector<Float>& lengths() const { return _l; }
    VectorScalingFunction* clone() const { return new VectorScalingFunction(*this); }
    SizeType result_size() const { return _l.size(); }
    SizeType argument_size() const { return _l.size(); }
    SmoothnessType smoothness() const { return SMOOTH; }
    Vector<Float> evaluate(const Vector<Float>& x) const {
        Vector<Float> r(this->result_size()); _compute_approx(r,x); return r; }
    Vector<Interval> evaluate(const Vector<Interval>& x) const {
        Vector<Interval> r(this->result_size()); _compute(r,x); return r; }
    Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const {
        Vector<TaylorModel> r(this->result_size(),TaylorModel(x[0].argument_size(),x[0].accuracy_ptr()));
        _compute(r,x); return r; }
    Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const {
        Vector< Differential<Float> > r(this->result_size(),Differential<Float>(x[0].argument_size(),x[0].degree()));
        _compute_approx(r,x); return r; }
    Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const { 
        Vector< Differential<Interval> > r(this->result_size(),Differential<Interval>(x[0].argument_size(),x[0].degree()));
        _compute_approx(r,x); return r; }
    Matrix<Float> jacobian(const Vector<Float>& x) const { ARIADNE_NOT_IMPLEMENTED; }
    Matrix<Interval> jacobian(const Vector<Interval>& x) const { ARIADNE_NOT_IMPLEMENTED; }
    std::ostream& write(std::ostream& os) const {
        return os << "VectorScalingFunction( o=" << this->origin() << ", l=" << this->lengths() << " )"; }
  private:
    template<class R, class A> void _compute(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=_o[i]+_l[i]*x[i]; } }
    template<class R, class A> void _compute_approx(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=_o[i]+_l[i]*x[i]; } }
  private:
    Vector<Float> _o;
    Vector<Float> _l;
};


ImageSet::ImageSet()
    : _domain(), _function(new VectorScalingFunction())
{
}

ImageSet::ImageSet(const Vector<Interval>& dom)
    //  : _domain(dom), _function_ptr(new IdentityFunction(dom.size()))
    : _domain(Vector<Interval>(dom.size(),Interval(-1,1))),
      _function(new VectorScalingFunction(dom))
{
}


ImageSet::ImageSet(const Vector<Interval>& dom, const VectorFunction& fn)
    : _domain(dom), _function(fn)
{
    ARIADNE_ASSERT(dom.size()==fn.argument_size());
}


ImageSet*
ImageSet::clone() const
{
    return new ImageSet(*this);
}


uint
ImageSet::dimension() const
{
    return this->_function.result_size();
}


tribool
ImageSet::empty() const
{
    return this->_domain.empty();
}


tribool
ImageSet::disjoint(const Box& bx) const
{
    if(dynamic_cast<const VectorScalingFunction*>(this->_function.pointer())) {
        return Ariadne::disjoint(bx,this->_function.evaluate(this->_domain));
    } else {
        static const int MAX_SUBDIVISIONS=8;
        return Ariadne::disjoint(this->domain(),this->_function,bx,radius(bx)/MAX_SUBDIVISIONS);
    }
}


tribool
ImageSet::overlaps(const Box& bx) const
{
    static const int MAX_SUBDIVISIONS=8;
    return !Ariadne::disjoint(this->domain(),this->_function,bx,radius(bx)/MAX_SUBDIVISIONS);
}


tribool
ImageSet::inside(const Box& bx) const
{
    if(Ariadne::inside(this->bounding_box(),bx)) {
        return true;
    } else {
        return indeterminate;
    }
}


Box
ImageSet::bounding_box() const
{
    return this->_function.evaluate(this->_domain);
}


void
ImageSet::draw(CanvasInterface& os) const
{
    return this->bounding_box().draw(os);
}


std::ostream&
ImageSet::write(std::ostream& os) const
{
    ARIADNE_NOT_IMPLEMENTED;
}





ConstraintSet::ConstraintSet(const Vector<Interval>& codom, const VectorFunction& fn)
    : _codomain(codom), _function(fn)
{
    ARIADNE_ASSERT(codom.size()==fn.result_size());
}

ConstraintSet*
ConstraintSet::clone() const
{
    return new ConstraintSet(*this);
}


uint
ConstraintSet::dimension() const
{
    return this->_function.argument_size();
}


tribool
ConstraintSet::disjoint(const Box& bx) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


tribool
ConstraintSet::overlaps(const Box& bx) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


tribool
ConstraintSet::covers(const Box& bx) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


std::ostream&
ConstraintSet::write(std::ostream& os) const
{
    ARIADNE_NOT_IMPLEMENTED;
}



} // namespace Ariadne
