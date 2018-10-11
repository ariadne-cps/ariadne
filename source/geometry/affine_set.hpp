/***************************************************************************
 *            affine_set.hpp
 *
 *  Copyright  2009  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file affine_set.hpp
 *  \brief Affine sets described by equality, inequality and interval constraints over a is_bounded box.
 */

#ifndef ARIADNE_AFFINE_SET_HPP
#define ARIADNE_AFFINE_SET_HPP

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "../utility/declarations.hpp"
#include "../output/logging.hpp"
#include "../utility/container.hpp"
#include "../output/graphics_interface.hpp"
#include "../function/affine.hpp"
#include "../function/affine_model.hpp"
#include "../function/constraint.hpp"

namespace Ariadne {

template<class X> struct LinearProgram;

class Grid;
class PavingInterface;
class GridCell;
class DiscreteEvent;
class Figure;
class CanvasInterface;
struct Point2d;


typedef Constraint<AffineModel<ValidatedTag,FloatDP>,FloatDPBounds> ValidatedAffineModelConstraint;
typedef Constraint<Affine<ValidatedNumericType>,FloatDPBounds> ValidatedAffineConstraint;
typedef Constraint<Affine<EffectiveNumericType>,EffectiveNumericType> EffectiveAffineConstraint;
typedef Affine<ValidatedNumericType> ValidatedAffineFunction;

EffectiveAffineConstraint operator<=(const EffectiveNumericType& l, const EffectiveAffine& am);
EffectiveAffineConstraint operator<=(const EffectiveAffine& am, const EffectiveNumericType& u);
EffectiveAffineConstraint operator==(const EffectiveAffine& am, const EffectiveNumericType& b);

ValidatedAffineConstraint operator<=(const FloatDPBounds& l, const ValidatedAffineFunction& am);
ValidatedAffineConstraint operator<=(const ValidatedAffineFunction& am, const FloatDPBounds& u);
ValidatedAffineConstraint operator<=(const ValidatedAffineConstraint& am, const FloatDPBounds& u);
ValidatedAffineConstraint operator==(const ValidatedAffineFunction& am, const FloatDPBounds& b);

ValidatedAffineModelConstraint operator<=(const FloatDPBounds& l, const ValidatedAffineModel& am);
ValidatedAffineModelConstraint operator<=(const ValidatedAffineModel& am, const FloatDPBounds& u);
ValidatedAffineModelConstraint operator<=(const ValidatedAffineModelConstraint& am, const FloatDPBounds& u);
ValidatedAffineModelConstraint operator==(const ValidatedAffineModel& am, const FloatDPBounds& b);

//! \brief A constrained image set defined by affine functions.
//!  Defines a set of the form \f$S=\{ f(x) \mid x\in D \mid g(x)\leq 0 \wedge h(x)=0 \}\f$
//!  where \f$f\f$, \f$g\f$ and \f$h\f$ are vector-valued affine functions on \f$\R^n\f$ and \f$D\f$ is a box in \f$\R^n\f$.
//!
//! Includes the class of zonotopes \f$Z=\{ f(x) \mid x\in D\}\f$, polyhedra \f$\{ f(x) \mid x\in\R^n \mid g(x)\leq 0\}\f$ and polytopes \f$\{ f(x) \mid x\in[0,\infty)^n \mid \sum_{i=1}^{n} x_i -1 = 0\}\f$.
//! \sa ValidatedConstrainedImageSet
class ValidatedAffineConstrainedImageSet
    : public virtual ValidatedCompactSetInterface
	, public virtual DrawableInterface
	, public Loggable
{
    ExactBoxType  _domain;
    Vector<ValidatedAffineModel> _space_models;
    List<ValidatedAffineModelConstraint> _constraint_models;
  public:
    //!\brief The set \f$\{ y \mid y\in D\}\f$.
    ValidatedAffineConstrainedImageSet(const ExactBoxType& D);
    //!\brief The set \f$\{ Gy+c \mid y\in D\}\f$.
    ValidatedAffineConstrainedImageSet(const ExactBoxType& D, const Matrix<FloatDPValue>& G, const Vector<FloatDPValue>& c);
    //!\brief The set \f$\{ Gy+c \mid ||y||_\infty\leq 1\}\f$. \deprecated
    ValidatedAffineConstrainedImageSet(const Matrix<FloatDPValue>& G, const Vector<FloatDPValue>& c);
    //!\brief The set \f$\{ x_i=f_i(s) \mid s\in D \}\f$.
    ValidatedAffineConstrainedImageSet(const ExactBoxType& D, const Vector<ValidatedAffine>& f);
    //!\brief The set \f$\{ x_i=f_i(s) \mid s\in D \mid c(s) \}\f$.
    ValidatedAffineConstrainedImageSet(const ExactBoxType& D, const Vector<ValidatedAffine>& f, const List<ValidatedAffineConstraint>& c);
    //!\brief The set \f$\{ x_i=f_i(s) \mid s\in [-1,+1]^n \mid c(s) \}\f$.
    explicit ValidatedAffineConstrainedImageSet(const Vector<ValidatedAffineModel>& f, const List<ValidatedAffineModelConstraint>& c);
    explicit ValidatedAffineConstrainedImageSet(const Vector<ValidatedAffineModel>& f);

    ValidatedAffineConstrainedImageSet(const ExactBoxType& D, const Vector<ValidatedAffineModel>& f, const List<ValidatedAffineModelConstraint>& c);

    //!\brief The set \f$\{ x \mid x\in B\}\f$.
    ValidatedAffineConstrainedImageSet(const RealBox& B);
    //!\brief The set \f$\{ x \mid x\in B\}\f$.
    ValidatedAffineConstrainedImageSet(const Box<Interval<FloatDPBall>>& B);

    ValidatedAffineConstrainedImageSet* clone() const;
    Void new_parameter_constraint(const EffectiveAffineConstraint& c);
    Void new_parameter_constraint(const ValidatedAffineConstraint& c);
    Void new_constraint(const ValidatedAffineModelConstraint& c);

    DimensionType dimension() const;
    SizeType number_of_parameters() const;
    SizeType number_of_constraints() const;
    ExactBoxType domain() const;

    ValidatedKleenean is_bounded() const;
    UpperBoxType bounding_box() const;
    ValidatedLowerKleenean separated(const ExactBoxType& bx) const;
    ValidatedLowerKleenean inside(const ExactBoxType& bx) const;
    ValidatedLowerKleenean is_empty() const;

    //! \brief Compute the image of \f$S\f$ under the function \f$h\f$.
    friend ValidatedAffineConstrainedImageSet image(ValidatedAffineConstrainedImageSet set, ValidatedVectorMultivariateFunction const& h);

    Void adjoin_outer_approximation_to(PavingInterface& g, Nat depth) const;
    GridTreePaving outer_approximation(const Grid& g, Nat depth) const;
    Void robust_adjoin_outer_approximation_to(PavingInterface& paving, Nat depth) const;

    List<Point2d> boundary(Nat xc, Nat yc) const;

    virtual Void draw(CanvasInterface&, const Projection2d& p) const;
    virtual OutputStream& write(OutputStream& os) const;

  private:
    Void construct(const ExactBoxType& D, const Matrix<FloatDPValue>& G, const Vector<FloatDPValue>& c);
    Void construct_linear_program(LinearProgram<FloatDP>& lp) const;
    static Void _robust_adjoin_outer_approximation_to(PavingInterface& paving, LinearProgram<FloatDP>& lp, const Vector<FloatDP>& errors, GridCell& cell, Nat depth);
    static Void _adjoin_outer_approximation_to(PavingInterface& paving, LinearProgram<FloatDP>& lp, const Vector<FloatDP>& errors, GridCell& cell, Nat depth);
};

inline OutputStream& operator<<(OutputStream& os, const ValidatedAffineConstrainedImageSet& as) {
    return as.write(os); }



} // namespace Ariadne

#endif // ARIADNE_AFFINE_SET_HPP
