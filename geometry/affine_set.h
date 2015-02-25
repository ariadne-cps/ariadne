/***************************************************************************
 *            affine_set.h
 *
 *  Copyright  2009  Pieter Collins
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

/*! \file affine_set.h
 *  \brief Affine sets described by equality, inequality and interval constraints over a bounded box.
 */

#ifndef ARIADNE_AFFINE_SET_H
#define ARIADNE_AFFINE_SET_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "utility/declarations.h"
#include "utility/logging.h"
#include "utility/container.h"
#include "output/graphics_interface.h"
#include "function/affine.h"
#include "function/affine_model.h"
#include "function/constraint.h"

namespace Ariadne {

template<class X> class LinearProgram;

class Grid;
class PavingInterface;
class GridCell;
class DiscreteEvent;
class Figure;
class CanvasInterface;
class Point2d;


typedef Constraint<AffineModel<ValidatedNumber>,ValidatedFloat64> ValidatedAffineModelConstraint;
typedef Constraint<Affine<ValidatedNumber>,ValidatedFloat64> ValidatedAffineConstraint;
typedef Constraint<Affine<EffectiveNumber>,EffectiveNumber> EffectiveAffineConstraint;
typedef Affine<ValidatedNumber> ValidatedAffineFunction;

EffectiveAffineConstraint operator<=(const EffectiveNumber& l, const EffectiveAffine& am);
EffectiveAffineConstraint operator<=(const EffectiveAffine& am, const EffectiveNumber& u);
EffectiveAffineConstraint operator<=(const EffectiveAffine& am, const EffectiveNumber& u);
EffectiveAffineConstraint operator==(const EffectiveAffine& am, const EffectiveNumber& b);

ValidatedAffineConstraint operator<=(const ValidatedFloat64& l, const ValidatedAffineFunction& am);
ValidatedAffineConstraint operator<=(const ValidatedAffineFunction& am, const ValidatedFloat64& u);
ValidatedAffineConstraint operator<=(const ValidatedAffineConstraint& am, const ValidatedFloat64& u);
ValidatedAffineConstraint operator==(const ValidatedAffineFunction& am, const ValidatedFloat64& b);

ValidatedAffineModelConstraint operator<=(const ValidatedFloat64& l, const ValidatedAffineModel& am);
ValidatedAffineModelConstraint operator<=(const ValidatedAffineModel& am, const ValidatedFloat64& u);
ValidatedAffineModelConstraint operator<=(const ValidatedAffineModelConstraint& am, const ValidatedFloat64& u);
ValidatedAffineModelConstraint operator==(const ValidatedAffineModel& am, const ValidatedFloat64& b);

//! \brief A constrained image set defined by affine functions.
//!  Defines a set of the form \f$S=\{ f(x) \mid x\in D \mid g(x)\leq 0 \wedge h(x)=0 \}\f$
//!  where \f$f\f$, \f$g\f$ and \f$h\f$ are vector-valued affine functions on \f$\R^n\f$ and \f$D\f$ is a box in \f$\R^n\f$.
//!
//! Includes the class of zonotopes \f$Z=\{ f(x) \mid x\in D\}\f$, polyhedra \f$\{ f(x) \mid x\in\R^n \mid g(x)\leq 0\}\f$ and polytopes \f$\{ f(x) \mid x\in[0,\infty)^n \mid \sum_{i=1}^{n} x_i -1 = 0\}\f$.
//! \sa ValidatedConstrainedImageSet
class ValidatedAffineConstrainedImageSet
    : public virtual CompactSetInterface
	, public virtual DrawableInterface
	, public Loggable
{
    ExactBox  _domain;
    Vector<ValidatedAffineModel> _space_models;
    List<ValidatedAffineModelConstraint> _constraint_models;
  public:
    //!\brief The set \f$\{ Gy+c \mid y\in D\}\f$.
    ValidatedAffineConstrainedImageSet(const ExactBox& D, const Matrix<ExactFloat64>& G, const Vector<ExactFloat64>& c);
    //!\brief The set \f$\{ Gy+c \mid ||y||_\infty\leq 1\}\f$. \deprecated
    ValidatedAffineConstrainedImageSet(const Matrix<ExactFloat64>& G, const Vector<ExactFloat64>& c);
    //!\brief The set \f$\{ x_i=f_i(s) \mid s\in D \}\f$.
    ValidatedAffineConstrainedImageSet(const ExactBox& D, const Vector<ValidatedAffine>& f);
    //!\brief The set \f$\{ x_i=f_i(s) \mid s\in D \mid c(s) \}\f$.
    ValidatedAffineConstrainedImageSet(const ExactBox& D, const Vector<ValidatedAffine>& f, const List<ValidatedAffineConstraint>& c);
    //!\brief The set \f$\{ x_i=f_i(s) \mid s\in [-1,+1]^n \mid c(s) \}\f$.
    explicit ValidatedAffineConstrainedImageSet(const Vector<ValidatedAffineModel>& f, const List<ValidatedAffineModelConstraint>& c);
    explicit ValidatedAffineConstrainedImageSet(const Vector<ValidatedAffineModel>& f);

    ValidatedAffineConstrainedImageSet(const ExactBox& D, const Vector<ValidatedAffineModel>& f, const List<ValidatedAffineModelConstraint>& c);

    ValidatedAffineConstrainedImageSet* clone() const;
    Void new_parameter_constraint(const EffectiveAffineConstraint& c);
    Void new_parameter_constraint(const ValidatedAffineConstraint& c);
    Void new_constraint(const ValidatedAffineModelConstraint& c);

    Nat dimension() const;
    Nat number_of_parameters() const;
    Nat number_of_constraints() const;
    ExactBox domain() const;

    Tribool bounded() const;
    UpperBox bounding_box() const;
    Tribool separated(const ExactBox& bx) const;
    Tribool inside(const ExactBox& bx) const;
    Tribool empty() const;

    Void adjoin_outer_approximation_to(PavingInterface& g, Int depth) const;
    GridTreeSet outer_approximation(const Grid& g, Int depth) const;
    Void robust_adjoin_outer_approximation_to(PavingInterface& paving, Int depth) const;

    List<Point2d> boundary(Nat xc, Nat yc) const;

    virtual Void draw(CanvasInterface&, const Projection2d& p) const;
    virtual OutputStream& write(OutputStream& os) const;

  private:
    Void construct(const ExactBox& D, const Matrix<ExactFloat64>& G, const Vector<ExactFloat64>& c);
    Void construct_linear_program(LinearProgram<Float64>& lp) const;
    static Void _robust_adjoin_outer_approximation_to(PavingInterface& paving, LinearProgram<Float64>& lp, const Vector<Float64>& errors, GridCell& cell, Int depth);
    static Void _adjoin_outer_approximation_to(PavingInterface& paving, LinearProgram<Float64>& lp, const Vector<Float64>& errors, GridCell& cell, Int depth);
};

inline OutputStream& operator<<(OutputStream& os, const ValidatedAffineConstrainedImageSet& as) {
    return as.write(os); }



} // namespace Ariadne

#endif // ARIADNE_AFFINE_SET_H
