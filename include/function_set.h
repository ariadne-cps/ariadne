/***************************************************************************
 *            function_set.h
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

/*! \file function_set.h
 *  \brief Images and preimages of boxes in Euclidean space.
 */

#ifndef ARIADNE_FUNCTION_SET_H
#define ARIADNE_FUNCTION_SET_H

#include <iosfwd>

#include <boost/shared_ptr.hpp>

#include "macros.h"
#include "container.h"
#include "numeric.h"
#include "vector.h"
#include "set_interface.h"
#include "function.h"
#include "constraint.h"
#include "graphics_interface.h"

#include "box.h"

namespace Ariadne {


class RealConstraintSet;
class RealBoundedConstraintSet;
class RealConstrainedImageSet;

class IntervalConstrainedImageSet;

class IntervalAffineConstrainedImageSet;

class Grid;
class PavingInterface;


//! \ingroup GeometryModule ExactSetSubModule
//! \brief A set defined as the intersection of an exact box with preimage of an exact box (the \em codomain) under a continuous function.
//! The set is described as \f$S=D\cap g^{-1}(C) = \{ x\in D \mid g(x)\in C\}\f$ where \f$D\f$ is the domain, \f$C\f$ is the codomain and \f$g\f$ the function.
class RealConstraintSet
    : public RegularSetInterface
{
    Nat _dimension;
    List< RealConstraint > _constraints;
  public:
    //! \brief Construct the preimage of \a C under \a g.
    RealConstraintSet(const RealVectorFunction& g, const RealBoxSet& C);
    //! \brief Construct the restriction of \a D under the constraints \a c.
    RealConstraintSet(const List<RealConstraint>& c);
    //! \brief The codomain of the set.
    const RealBoxSet codomain() const { return this->constraint_bounds(); }
    //! \brief The function used to define the constraints.
    const RealVectorFunction constraint_function() const;
    //! \brief The bounds of the constraints.
    const RealBoxSet constraint_bounds() const;
    //! \brief The number of constraints.
    Nat number_of_constraints() const { return this->_constraints.size(); };
    //! \brief The constraints.
    List<RealConstraint> const& constraints() const { return this->_constraints; }
    //! \brief The \a i<sup>th</sup> constraint.
    RealConstraint const& constraint(Nat i) const { return this->_constraints[i]; }

    RealConstraintSet* clone() const;
    Nat dimension() const;
    Tribool separated(const Box&) const;
    Tribool overlaps(const Box&) const;
    Tribool covers(const Box&) const;
    OutputStream& write(OutputStream&) const;
};

//! \ingroup GeometryModule ExactSetSubModule
//! \brief A set defined as the intersection of an exact box with preimage of an exact box (the \em codomain) under a continuous function.
//! The set is described as \f$S=D\cap g^{-1}(C) = \{ x\in D \mid g(x)\in C\}\f$ where \f$D\f$ is the domain, \f$C\f$ is the codomain and \f$g\f$ the function.
class RealBoundedConstraintSet
    : public SetInterface
    , public DrawableInterface
{
    RealBoxSet _domain;
    List< RealConstraint > _constraints;
  public:
    //! \brief Construct the preimage of \a C under \a g.
    RealBoundedConstraintSet(const RealBoxSet& D, const RealVectorFunction& g, const RealBoxSet& C);
    //! \brief Construct the restriction of \a D under the constraints \a c.
    RealBoundedConstraintSet(const RealBoxSet& D, const List<RealConstraint>& c);
    //! \brief Construct the box \a D.
    RealBoundedConstraintSet(const RealBoxSet& bx);
    //! \brief The domain of the set.
    const RealBoxSet& domain() const { return this->_domain; }
    //! \brief The codomain of the set.
    const RealBoxSet codomain() const { return this->constraint_bounds(); }
    //! \brief The function used to define the constraints.
    const RealVectorFunction constraint_function() const;
    //! \brief The bounds for the constraints.
    const RealBoxSet constraint_bounds() const;
    //! \brief The number of constraints.
    Nat number_of_constraints() const { return this->_constraints.size(); };
    //! \brief The constraints.
    List<RealConstraint> const& constraints() const { return this->_constraints; }
    //! \brief The \a i<sup>th</sup> constraint.
    RealConstraint const& constraint(Nat i) const { return this->_constraints[i]; }

    RealBoundedConstraintSet* clone() const;
    Nat dimension() const;
    Tribool separated(const Box&) const;
    Tribool overlaps(const Box&) const;
    Tribool covers(const Box&) const;
    Tribool inside(const Box&) const;
    Box bounding_box() const;
    OutputStream& write(OutputStream&) const;
    Void draw(CanvasInterface&,const Projection2d&) const;
};


class RealConstrainedImageSet
    : public LocatedSetInterface, public DrawableInterface
{
    RealBoxSet _domain;
    RealVectorFunction _function;
    List< RealConstraint > _constraints;
  public:
    //! \brief Construct the set with zero-dimensional parameterisation in zero dimensions with no constraints.
    RealConstrainedImageSet() : _domain(), _function() { }
    //! \brief Construct the box \a dom.
    RealConstrainedImageSet(const RealBoxSet& dom) : _domain(dom), _function(RealVectorFunction::identity(dom.size())) { }
    //! \brief Construct the image of \a dom under \a fn.
    RealConstrainedImageSet(const RealBoxSet& dom, const RealVectorFunction& fn) : _domain(dom), _function(fn) {
        ARIADNE_ASSERT_MSG(dom.size()==fn.argument_size(),"dom="<<dom<<", fn="<<fn); }
    //! \brief Construct the image of \a dom under \a fn, using constraint \a c.
    RealConstrainedImageSet(const RealBoxSet& dom, const RealVectorFunction& fn, const RealConstraint& c) : _domain(dom), _function(fn), _constraints(1u,c) {
        ARIADNE_ASSERT_MSG(dom.size()==fn.argument_size(),"dom="<<dom<<", fn="<<fn);
        ARIADNE_ASSERT_MSG(dom.size()==c.function().argument_size(),"dom="<<dom<<", c="<<c);
    }
    //! \brief Construct the image of \a dom under \a fn, using constraints \a c.
    RealConstrainedImageSet(const RealBoxSet& dom, const RealVectorFunction& fn, const List<RealConstraint>& c) : _domain(dom), _function(fn), _constraints(c) {
        ARIADNE_ASSERT_MSG(dom.size()==fn.argument_size(),"dom="<<dom<<", fn="<<fn); }
    //! \brief Convert from a bounded constraint set.
    RealConstrainedImageSet(const RealBoundedConstraintSet& set);
    //! \brief The domain of the set.
    const RealBoxSet& domain() const { return this->_domain; }
    //! \brief The function used to define the mapping from the parameter domain to the space.
    const RealVectorFunction& function() const { return this->_function; };
    //! \brief The bounds for the constraints.
    const RealVectorFunction constraint_function() const;
    //! \brief The bounds for the constraints.
    const RealBoxSet constraint_bounds() const;
    //! \brief The function used to define the set.
    const List<RealConstraint>& constraints() const { return this->_constraints; };
    //! \brief The number of parameters used to define the set, which equals the dimension of \f$D\f$.
    Nat number_of_parameters() const { return this->_domain.size(); };
    //! \brief The number of constraints.
    Nat number_of_constraints() const { return this->_constraints.size(); };
    //! \brief The \a i<sup>th</sup> constraint.
    RealConstraint const& constraint(Nat i) const { return this->_constraints[i]; }

    //! \brief Apply the function \f$h\f$ to obtain the set \f$h\circ f(D\cap g^{-1}(C))\f$.
    Void apply(const RealVectorFunction& h) {
        this->_function=compose(h,this->_function);
    }

    //! \brief Introduce a new constraint of the form \f$g(y)\in [c_l,c_u]\f$.
    Void new_parameter_constraint(const RealConstraint& c) {
        ARIADNE_ASSERT_MSG(c.function().argument_size()==this->domain().size(),*this<<", "<<c);
        this->_constraints.append(c); }

    //! \brief Introduce a new constraint of the form \f$g(y)\in [c_l,c_u]\f$.
    Void new_space_constraint(const RealConstraint& c) {
        ARIADNE_ASSERT_MSG(c.function().argument_size()==this->_function.result_size(),*this<<", "<<c);
        this->_constraints.append(RealConstraint(c.lower_bound(),compose(c.function(),_function),c.upper_bound())); }

    RealConstrainedImageSet* clone() const { return new RealConstrainedImageSet(*this); }
    Nat dimension() const { return this->_function.result_size(); }
    Tribool inside(const Box& bx) const { return subset(this->bounding_box(),bx); }

    //! \brief A coarse over-approximation to the set. Computed by taking the interval evaluation \f$h(D)\f$.
    Box bounding_box() const;
    //! \brief Construct an affine over-approximation
    IntervalAffineConstrainedImageSet affine_over_approximation() const;
    //! \brief Construct an affine approximation, with undefined accuracy.
    IntervalAffineConstrainedImageSet affine_approximation() const;
    //! \brief Split into two pieces by subdividing along a coordinate direction.
    Pair<RealConstrainedImageSet,RealConstrainedImageSet> split() const;
    //! \brief Split into two pieces by subdividing along the \a j<sup>th</sup> coordinate direction.
    Pair<RealConstrainedImageSet,RealConstrainedImageSet> split(Nat j) const;

    //! \brief Test if the set is disjoint from a (closed) box.
    Tribool separated(const Box&) const;
    //! \brief Test if the set overlaps (intersects the interior of) a box.
    Tribool overlaps(const Box&) const;
    //! \brief Adjoin an outer approximation to a paving.
    Void adjoin_outer_approximation_to(PavingInterface& paving, Int depth) const;

    //! \brief Test if the set satisfies the state constraint at all points.
    Tribool satisfies(const RealConstraint& c) const;

    //! \brief Draw to a canvas.
    Void draw(CanvasInterface&,const Projection2d&) const;
    //! \brief Write to an output stream.
    OutputStream& write(OutputStream&) const;
  private:
    Void affine_adjoin_outer_approximation_to(PavingInterface& paving, Int depth) const;
    Void subdivision_adjoin_outer_approximation_to(PavingInterface& paving, Int depth) const;
    Void constraint_adjoin_outer_approximation_to(PavingInterface& paving, Int depth) const;
};



//! \ingroup GeometryModule ExactSetSubModule
//! \brief A set defined as the image of the intersection of a box \f$D\f$ and a constraint set \f$g^{-1}(C)\f$ under a function \f$f\f$.
//! In other words, \f$S=f(D\cap g^{-1}(C))\f$.
class IntervalConstrainedImageSet
    : public LocatedSetInterface, public DrawableInterface
{
    IntervalVector _domain;
    Box _reduced_domain;
    IntervalVectorFunction _function;
    List< IntervalConstraint > _constraints;
  public:
    //! \brief Construct the set with zero-dimensional parameterisation in zero dimensions with no constraints.
    IntervalConstrainedImageSet() : _domain(), _function() { }
    //! \brief Construct the box \a dom.
    IntervalConstrainedImageSet(const Box& dom)
        : _domain(dom), _reduced_domain(dom)
        , _function(static_cast<IntervalVectorFunctionInterface const&>(RealVectorFunction::identity(dom.size()))) { }
    //! \brief Construct the image of \a dom under \a fn.
    IntervalConstrainedImageSet(const Vector<Interval>& dom, const IntervalVectorFunction& fn) : _domain(dom), _reduced_domain(dom), _function(fn) {
        ARIADNE_ASSERT_MSG(dom.size()==fn.argument_size(),"dom="<<dom<<", fn="<<fn); }

    //! \brief The domain of the set.
    const IntervalVector& domain() const { return this->_domain; }
    //! \brief The function used to define the set.
    const IntervalVectorFunction& function() const { return this->_function; }
    //! \brief The constraints used to define the set.
    const List<IntervalConstraint> constraints() const { return this->_constraints; }
    //! \brief The \a i<sup>th</sup> constraint.
    IntervalConstraint const constraint(Nat i) const { return this->_constraints[i]; };
    //! \brief The number of parameters used to define the set, which equals the dimension of \f$D\f$.
    Nat number_of_parameters() const { return this->_domain.size(); };
    //! \brief The number of constraints.
    Nat number_of_constraints() const { return this->_constraints.size(); };

    //! \brief Apply the function \f$h\f$ to obtain the set \f$h\circ f(D\cap g^{-1}(C))\f$.
    Void apply(const IntervalVectorFunction& h) {
        this->_function=compose(h,this->_function);
    }

    //! \brief Introduce a new constraint of the form \f$g(s)\in [c_l,c_u]\f$.
    Void new_parameter_constraint(const IntervalConstraint& c) {
        ARIADNE_ASSERT_MSG(c.function().argument_size()==this->number_of_parameters(),*this<<", "<<c);
        this->_constraints.append(c);
    }

    //! \brief Introduce a new constraint of the form \f$g(x)\in [c_l,c_u]\f$.
    Void new_space_constraint(const IntervalConstraint& c) {
        ARIADNE_ASSERT_MSG(c.function().argument_size()==this->dimension(),*this<<", "<<c);
        this->_constraints.append(IntervalConstraint(c.lower_bound(),compose(c.function(),this->function()),c.upper_bound()));
    }

    IntervalConstrainedImageSet* clone() const { return new IntervalConstrainedImageSet(*this); }
    Nat dimension() const { return this->_function.result_size(); }

    IntervalVectorFunction constraint_function() const;
    IntervalVector constraint_bounds() const;

    //! \brief Reduce the size of the domain by constraint propagation, if possible.
    Void reduce();
    //! \brief A coarse over-approximation to the set. Computed by taking the interval evaluation \f$h(D)\f$.
    Box bounding_box() const;
    //! \brief Construct an affine over-approximation
    IntervalAffineConstrainedImageSet affine_over_approximation() const;
    //! \brief Construct an affine approximation, with undefined accuracy.
    IntervalAffineConstrainedImageSet affine_approximation() const;
    //! \brief Split into two pieces by subdividing along a coordinate direction.
    Pair<IntervalConstrainedImageSet,IntervalConstrainedImageSet> split() const;
    //! \brief Split into two pieces by subdividing along the \a j<sup>th</sup> coordinate direction.
    Pair<IntervalConstrainedImageSet,IntervalConstrainedImageSet> split(Nat j) const;

    //! \brief Test if the set is empty.
    Tribool empty() const;
    //! \brief Test if the set is a strict subset of a box.
    Tribool inside(const Box& bx) const;
    //! \brief Test if the set is disjoint from a box.
    Tribool separated(const Box&) const;
    //! \brief Test if the set overlaps (intersects the interior of) a box.
    Tribool overlaps(const Box&) const;
    //! \brief Adjoin an outer approximation to a paving.
    Void adjoin_outer_approximation_to(PavingInterface& paving, Int depth) const;

    //! \brief Test if the set satisfies the state constraint at all points.
    Tribool satisfies(const IntervalConstraint& c) const;

    //! \brief Draw to a canvas.
    Void draw(CanvasInterface&,const Projection2d&) const;
    //! \brief Write to an output stream.
    OutputStream& write(OutputStream&) const;
};

OutputStream& operator<<(OutputStream&, const IntervalConstrainedImageSet&);



} //namespace Ariadne


#endif
