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

#include <boost/smart_ptr.hpp>
#include "graphics_interface.h"
#include "container.h"
#include "affine.h"
#include "affine_model.h"
#include "constraint.h"
#include "logging.h"

namespace Ariadne {

class Float;
class Interval;
template<class X> class Vector;
template<class X> class LinearProgram;
class Box;
class Grid;
class PavingInterface;
class GridCell;
class DiscreteEvent;
class Figure;
class CanvasInterface;
class Point2d;

typedef Constraint<AffineModel<Interval>,Float> IntervalAffineModelConstraint;
typedef Constraint<Affine<Interval>,Float> IntervalAffineConstraint;
typedef Constraint<Affine<Real>,Real> RealAffineConstraint;
typedef Affine<Interval> IntervalAffineFunction;

RealAffineConstraint operator<=(const Real& l, const RealAffine& am);
RealAffineConstraint operator<=(const RealAffine& am, const Real& u);
RealAffineConstraint operator<=(const RealAffine& am, const Real& u);
RealAffineConstraint operator==(const RealAffine& am, const Real& b);

IntervalAffineConstraint operator<=(const Float& l, const IntervalAffineFunction& am);
IntervalAffineConstraint operator<=(const IntervalAffineFunction& am, const Float& u);
IntervalAffineConstraint operator<=(const IntervalAffineConstraint& am, const Float& u);
IntervalAffineConstraint operator==(const IntervalAffineFunction& am, const Float& b);

IntervalAffineModelConstraint operator<=(const Float& l, const IntervalAffineModel& am);
IntervalAffineModelConstraint operator<=(const IntervalAffineModel& am, const Float& u);
IntervalAffineModelConstraint operator<=(const IntervalAffineModelConstraint& am, const Float& u);
IntervalAffineModelConstraint operator==(const IntervalAffineModel& am, const Float& b);

//! \brief A constrained image set defined by affine functions.
//!  Defines a set of the form \f$S=\{ f(x) \mid x\in D \mid g(x)\leq 0 \wedge h(x)=0 \}\f$
//!  where \f$f\f$, \f$g\f$ and \f$h\f$ are vector-valued affine functions on \f$\R^n\f$ and \f$D\f$ is a box in \f$\R^n\f$.
//!
//! Includes the class of zonotopes \f$Z=\{ f(x) \mid x\in D\}\f$, polyhedra \f$\{ f(x) \mid x\in\R^n \mid g(x)\leq 0\}\f$ and polytopes \f$\{ f(x) \mid x\in[0,\infty)^n \mid \sum_{i=1}^{n} x_i -1 = 0\}\f$.
//! \sa IntervalConstrainedImageSet
class IntervalAffineConstrainedImageSet
    : public DrawableInterface, public Loggable
{
    IntervalVector _domain;
    Vector<IntervalAffineModel> _space_models;
    List<IntervalAffineModelConstraint> _constraint_models;
  public:
    //!\brief The set \f$\{ Gy+c \mid y\in D\}\f$.
    IntervalAffineConstrainedImageSet(const Vector<Interval>& D, const Matrix<Float>& G, const Vector<Float>& c);
    //!\brief The set \f$\{ Gy+c \mid ||y||_\infty\leq 1\}\f$. \deprecated
    IntervalAffineConstrainedImageSet(const Matrix<Float>& G, const Vector<Float>& c);
    //!\brief The set \f$\{ x_i=f_i(s) \mid s\in D \}\f$.
    IntervalAffineConstrainedImageSet(const IntervalVector& D, const Vector<IntervalAffine>& f);
    //!\brief The set \f$\{ x_i=f_i(s) \mid s\in D \mid c(s) \}\f$.
    IntervalAffineConstrainedImageSet(const IntervalVector& D, const Vector<IntervalAffine>& f, const List<IntervalAffineConstraint>& c);
    //!\brief The set \f$\{ x_i=f_i(s) \mid s\in [-1,+1]^n \mid c(s) \}\f$.
    explicit IntervalAffineConstrainedImageSet(const Vector<IntervalAffineModel>& f, const List<IntervalAffineModelConstraint>& c);
    explicit IntervalAffineConstrainedImageSet(const Vector<IntervalAffineModel>& f);

    IntervalAffineConstrainedImageSet(const IntervalVector& D, const Vector<IntervalAffineModel>& f, const List<IntervalAffineModelConstraint>& c);

    IntervalAffineConstrainedImageSet* clone() const;
    void new_parameter_constraint(const RealAffineConstraint& c);
    void new_parameter_constraint(const IntervalAffineConstraint& c);
    void new_constraint(const IntervalAffineModelConstraint& c);

    uint dimension() const;
    uint number_of_parameters() const;
    uint number_of_constraints() const;
    IntervalVector domain() const;

    tribool bounded() const;
    Box bounding_box() const;
    tribool separated(const Box& bx) const;
    tribool empty() const;

    void adjoin_outer_approximation_to(PavingInterface& g, int depth) const;
    GridTreeSet outer_approximation(const Grid& g, int depth) const;
    void robust_adjoin_outer_approximation_to(PavingInterface& paving, int depth) const;

    List<Point2d> boundary(uint xc, uint yc) const;

    virtual void draw(CanvasInterface&, const Projection2d& p) const;
    virtual std::ostream& write(std::ostream& os) const;

  private:
    void construct(const Vector<Interval>& D, const Matrix<Float>& G, const Vector<Float>& c);
    void construct_linear_program(LinearProgram<Float>& lp) const;
    static void _robust_adjoin_outer_approximation_to(PavingInterface& paving, LinearProgram<Float>& lp, GridCell& cell, int depth);
    static void _adjoin_outer_approximation_to(PavingInterface& paving, LinearProgram<Float>& lp, GridCell& cell, int depth);
};

inline std::ostream& operator<<(std::ostream& os, const IntervalAffineConstrainedImageSet& as) {
    return as.write(os); }



} // namespace Ariadne

#endif // ARIADNE_AFFINE_SET_H
