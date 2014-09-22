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

#include <memory>

#include "macros.h"
#include "container.h"
#include "numeric.h"
#include "vector.h"
#include "set_interface.h"
#include "function.h"
#include "function_model.h"
#include "constraint.h"
#include "graphics_interface.h"

#include "box.h"

namespace Ariadne {


class ConstraintSet;
class BoundedConstraintSet;
class ConstrainedImageSet;

typedef ConstraintSet EffectiveConstraintSet;
typedef BoundedConstraintSet EffectiveBoundedConstraintSet;
typedef ConstrainedImageSet EffectiveConstrainedImageSet;

class ValidatedConstrainedImageSet;
class ValidatedAffineConstrainedImageSet;


class Grid;
class PavingInterface;


//! \ingroup GeometryModule ExactSetSubModule
//! \brief A set defined as the intersection of an exact box with preimage of an exact box (the \em codomain) under a continuous function.
//! The set is described as \f$S=D\cap g^{-1}(C) = \{ x\in D \mid g(x)\in C\}\f$ where \f$D\f$ is the domain, \f$C\f$ is the codomain and \f$g\f$ the function.
class ConstraintSet
    : public virtual RegularSetInterface
{
    Nat _dimension;
    List< EffectiveConstraint > _constraints;
  public:
    //! \brief Construct the preimage of \a C under \a g.
    ConstraintSet(const EffectiveVectorFunction& g, const BoxSet& C);
    //! \brief Construct the restriction of \a D under the constraints \a c.
    ConstraintSet(const List<EffectiveConstraint>& c);
    //! \brief The codomain of the set.
    const BoxSet codomain() const { return this->constraint_bounds(); }
    //! \brief The function used to define the constraints.
    const EffectiveVectorFunction constraint_function() const;
    //! \brief The bounds of the constraints.
    const BoxSet constraint_bounds() const;
    //! \brief The number of constraints.
    Nat number_of_constraints() const { return this->_constraints.size(); };
    //! \brief The constraints.
    List<EffectiveConstraint> const& constraints() const { return this->_constraints; }
    //! \brief The \a i<sup>th</sup> constraint.
    EffectiveConstraint const& constraint(Nat i) const { return this->_constraints[i]; }

    ConstraintSet* clone() const;
    Nat dimension() const;
    Tribool separated(const Box&) const;
    Tribool overlaps(const Box&) const;
    Tribool covers(const Box&) const;
    OutputStream& write(OutputStream&) const;
};


//! \ingroup GeometryModule ExactSetSubModule
//! \brief A set defined as the intersection of an exact box with preimage of an exact box (the \em codomain) under a continuous function.
//! The set is described as \f$S=D\cap g^{-1}(C) = \{ x\in D \mid g(x)\in C\}\f$ where \f$D\f$ is the domain, \f$C\f$ is the codomain and \f$g\f$ the function.
class BoundedConstraintSet
    : public virtual SetInterface
    , public virtual DrawableInterface
{
    BoxSet _domain;
    List< EffectiveConstraint > _constraints;
  public:
    //! \brief Construct the preimage of \a C under \a g.
    BoundedConstraintSet(const BoxSet& D, const EffectiveVectorFunction& g, const BoxSet& C);
    //! \brief Construct the restriction of \a D under the constraints \a c.
    BoundedConstraintSet(const BoxSet& D, const List<EffectiveConstraint>& c);
    //! \brief Construct the box \a D.
    BoundedConstraintSet(const BoxSet& bx);
    //! \brief The domain of the set.
    const BoxSet& domain() const { return this->_domain; }
    //! \brief The codomain of the set.
    const BoxSet codomain() const { return this->constraint_bounds(); }
    //! \brief The function used to define the constraints.
    const EffectiveVectorFunction constraint_function() const;
    //! \brief The bounds for the constraints.
    const BoxSet constraint_bounds() const;
    //! \brief The number of constraints.
    Nat number_of_constraints() const { return this->_constraints.size(); };
    //! \brief The constraints.
    List<EffectiveConstraint> const& constraints() const { return this->_constraints; }
    //! \brief The \a i<sup>th</sup> constraint.
    EffectiveConstraint const& constraint(Nat i) const { return this->_constraints[i]; }

    BoundedConstraintSet* clone() const;
    Nat dimension() const;
    Tribool separated(const Box&) const;
    Tribool overlaps(const Box&) const;
    Tribool covers(const Box&) const;
    Tribool inside(const Box&) const;
    Box bounding_box() const;
    OutputStream& write(OutputStream&) const;
    Void draw(CanvasInterface&,const Projection2d&) const;
};

BoundedConstraintSet intersection(const ConstraintSet& cs, const BoxSet& bx);


class ConstrainedImageSet
    : public virtual LocatedSetInterface, public virtual DrawableInterface
{
    BoxSet _domain;
    EffectiveVectorFunction _function;
    List< EffectiveConstraint > _constraints;
  public:
    //! \brief Construct the set with zero-dimensional parameterisation in zero dimensions with no constraints.
    ConstrainedImageSet() : _domain(), _function() { }
    //! \brief Construct the box \a dom.
    ConstrainedImageSet(const BoxSet& dom) : _domain(dom), _function(EffectiveVectorFunction::identity(dom.size())) { }
    //! \brief Construct the image of \a dom under \a fn.
    ConstrainedImageSet(const BoxSet& dom, const EffectiveVectorFunction& fn) : _domain(dom), _function(fn) {
        ARIADNE_ASSERT_MSG(dom.size()==fn.argument_size(),"dom="<<dom<<", fn="<<fn); }
    //! \brief Construct the image of \a dom under \a fn, using constraint \a c.
    ConstrainedImageSet(const BoxSet& dom, const EffectiveVectorFunction& fn, const EffectiveConstraint& c) : _domain(dom), _function(fn), _constraints(1u,c) {
        ARIADNE_ASSERT_MSG(dom.size()==fn.argument_size(),"dom="<<dom<<", fn="<<fn);
        ARIADNE_ASSERT_MSG(dom.size()==c.function().argument_size(),"dom="<<dom<<", c="<<c);
    }
    //! \brief Construct the image of \a dom under \a fn, using constraints \a c.
    ConstrainedImageSet(const BoxSet& dom, const EffectiveVectorFunction& fn, const List<EffectiveConstraint>& c) : _domain(dom), _function(fn), _constraints(c) {
        ARIADNE_ASSERT_MSG(dom.size()==fn.argument_size(),"dom="<<dom<<", fn="<<fn); }
    //! \brief Convert from a bounded constraint set.
    ConstrainedImageSet(const BoundedConstraintSet& set);
    //! \brief The domain of the set.
    const BoxSet& domain() const { return this->_domain; }
    //! \brief The function used to define the mapping from the parameter domain to the space.
    const EffectiveVectorFunction& function() const { return this->_function; };
    //! \brief The bounds for the constraints.
    const EffectiveVectorFunction constraint_function() const;
    //! \brief The bounds for the constraints.
    const BoxSet constraint_bounds() const;
    //! \brief The function used to define the set.
    const List<EffectiveConstraint>& constraints() const { return this->_constraints; };
    //! \brief The number of parameters used to define the set, which equals the dimension of \f$D\f$.
    Nat number_of_parameters() const { return this->_domain.size(); };
    //! \brief The number of constraints.
    Nat number_of_constraints() const { return this->_constraints.size(); };
    //! \brief The \a i<sup>th</sup> constraint.
    EffectiveConstraint const& constraint(Nat i) const { return this->_constraints[i]; }

    //! \brief Apply the function \f$h\f$ to obtain the set \f$h\circ f(D\cap g^{-1}(C))\f$.
    Void apply(const EffectiveVectorFunction& h) {
        this->_function=compose(h,this->_function);
    }

    //! \brief Introduce a new constraint of the form \f$g(y)\in [c_l,c_u]\f$.
    Void new_parameter_constraint(const EffectiveConstraint& c) {
        ARIADNE_ASSERT_MSG(c.function().argument_size()==this->domain().size(),*this<<", "<<c);
        this->_constraints.append(c); }

    //! \brief Introduce a new constraint of the form \f$g(y)\in [c_l,c_u]\f$.
    Void new_space_constraint(const EffectiveConstraint& c) {
        ARIADNE_ASSERT_MSG(c.function().argument_size()==this->_function.result_size(),*this<<", "<<c);
        this->_constraints.append(EffectiveConstraint(c.lower_bound(),compose(c.function(),_function),c.upper_bound())); }

    ConstrainedImageSet* clone() const { return new ConstrainedImageSet(*this); }
    Nat dimension() const { return this->_function.result_size(); }
    Tribool inside(const Box& bx) const { return subset(this->bounding_box(),bx); }

    //! \brief A coarse over-approximation to the set. Computed by taking the interval evaluation \f$h(D)\f$.
    Box bounding_box() const;
    //! \brief Construct an affine over-approximation
    ValidatedAffineConstrainedImageSet affine_over_approximation() const;
    //! \brief Construct an affine approximation, with undefined accuracy.
    ValidatedAffineConstrainedImageSet affine_approximation() const;
    //! \brief Split into two pieces by subdividing along a coordinate direction.
    Pair<ConstrainedImageSet,ConstrainedImageSet> split() const;
    //! \brief Split into two pieces by subdividing along the \a j<sup>th</sup> coordinate direction.
    Pair<ConstrainedImageSet,ConstrainedImageSet> split(Nat j) const;

    //! \brief Test if the set is disjoint from a (closed) box.
    Tribool separated(const Box&) const;
    //! \brief Test if the set overlaps (intersects the interior of) a box.
    Tribool overlaps(const Box&) const;
    //! \brief Adjoin an outer approximation to a paving.
    Void adjoin_outer_approximation_to(PavingInterface& paving, Int depth) const;

    //! \brief Test if the set satisfies the state constraint at all points.
    Tribool satisfies(const EffectiveConstraint& c) const;

    //! \brief Draw to a canvas.
    Void draw(CanvasInterface&,const Projection2d&) const;
    //! \brief Write to an output stream.
    OutputStream& write(OutputStream&) const;
};



//! \ingroup GeometryModule ExactSetSubModule
//! \brief A set defined as the image of the intersection of a box \f$D\f$ and a constraint set \f$g^{-1}(C)\f$ under a function \f$f\f$.
//! In other words, \f$S=f(D\cap g^{-1}(C))\f$.
class ValidatedConstrainedImageSet
    : public virtual LocatedSetInterface, public virtual DrawableInterface
{
    Box _domain;
    Box _reduced_domain;
    ValidatedVectorFunction _function;
    List< ValidatedConstraint > _constraints;
  public:
    //! \brief Construct the set with zero-dimensional parameterisation in zero dimensions with no constraints.
    ValidatedConstrainedImageSet() : _domain(), _function() { }
    //! \brief Construct the box \a dom.
    ValidatedConstrainedImageSet(const Box& dom)
        : _domain(dom), _reduced_domain(dom)
        , _function(static_cast<ValidatedVectorFunctionInterface const&>(EffectiveVectorFunction::identity(dom.size()))) { }
    //! \brief Construct the image of \a dom under \a fn.
    ValidatedConstrainedImageSet(const Box& dom, const ValidatedVectorFunction& fn) : _domain(dom), _reduced_domain(dom), _function(fn) {
        ARIADNE_ASSERT_MSG(dom.size()==fn.argument_size(),"dom="<<dom<<", fn="<<fn); }
    ValidatedConstrainedImageSet(const Box& dom, const ValidatedVectorFunctionModel& fn) : _domain(dom), _reduced_domain(dom), _function(fn.raw_pointer()->_clone()) {
        ARIADNE_ASSERT_MSG(dom.size()==fn.argument_size(),"dom="<<dom<<", fn="<<fn); }
    //! \brief Construct the image of \a dom under \a fn.
    ValidatedConstrainedImageSet(const Box& dom, const ValidatedVectorFunction& fn, const List<ValidatedConstraint>& cnstr) : _domain(dom), _reduced_domain(dom), _function(fn), _constraints(cnstr) {
        ARIADNE_ASSERT_MSG(dom.size()==fn.argument_size(),"dom="<<dom<<", fn="<<fn); }

    //! \brief The domain of the set.
    const Box& domain() const { return this->_domain; }
    //! \brief The function used to define the set.
    const ValidatedVectorFunction& function() const { return this->_function; }
    //! \brief The constraints used to define the set.
    const List<ValidatedConstraint> constraints() const { return this->_constraints; }
    //! \brief The \a i<sup>th</sup> constraint.
    ValidatedConstraint const constraint(Nat i) const { return this->_constraints[i]; };
    //! \brief The number of parameters used to define the set, which equals the dimension of \f$D\f$.
    Nat number_of_parameters() const { return this->_domain.size(); };
    //! \brief The number of constraints.
    Nat number_of_constraints() const { return this->_constraints.size(); };

    //! \brief Apply the function \f$h\f$ to obtain the set \f$h\circ f(D\cap g^{-1}(C))\f$.
    Void apply(const ValidatedVectorFunction& h) {
        this->_function=compose(h,this->_function);
    }

    //! \brief Introduce a new constraint of the form \f$g(s)\in [c_l,c_u]\f$.
    Void new_parameter_constraint(const ValidatedConstraint& c) {
        ARIADNE_ASSERT_MSG(c.function().argument_size()==this->number_of_parameters(),*this<<", "<<c);
        this->_constraints.append(c);
    }

    //! \brief Introduce a new constraint of the form \f$g(x)\in [c_l,c_u]\f$.
    Void new_space_constraint(const ValidatedConstraint& c) {
        ARIADNE_ASSERT_MSG(c.function().argument_size()==this->dimension(),*this<<", "<<c);
        this->_constraints.append(ValidatedConstraint(c.lower_bound(),compose(c.function(),this->function()),c.upper_bound()));
    }

    ValidatedConstrainedImageSet* clone() const { return new ValidatedConstrainedImageSet(*this); }
    Nat dimension() const { return this->_function.result_size(); }

    ValidatedVectorFunction constraint_function() const;
    Box constraint_bounds() const;

    //! \brief Reduce the size of the domain by constraint propagation, if possible.
    Void reduce();
    //! \brief A coarse over-approximation to the set. Computed by taking the interval evaluation \f$h(D)\f$.
    Box bounding_box() const;
    //! \brief Construct an affine over-approximation
    ValidatedAffineConstrainedImageSet affine_over_approximation() const;
    //! \brief Construct an affine approximation, with undefined accuracy.
    ValidatedAffineConstrainedImageSet affine_approximation() const;
    //! \brief Split into two pieces by subdividing along a coordinate direction.
    Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet> split() const;
    //! \brief Split into two pieces by subdividing along the \a j<sup>th</sup> coordinate direction.
    Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet> split(Nat j) const;

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
    Tribool satisfies(const ValidatedConstraint& c) const;

    //! \brief Draw to a canvas.
    Void draw(CanvasInterface&,const Projection2d&) const;
    //! \brief Write to an output stream.
    OutputStream& write(OutputStream&) const;
};

ValidatedConstrainedImageSet join(const ValidatedConstrainedImageSet& set1, const ValidatedConstrainedImageSet& set2);
OutputStream& operator<<(OutputStream&, const ValidatedConstrainedImageSet&);



} //namespace Ariadne


#endif
