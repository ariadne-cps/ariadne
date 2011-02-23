/***************************************************************************
 *            taylor_set.h
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

/*! \file taylor_set.h
 *  \brief Sets based on Taylor series.
 */
#ifndef ARIADNE_TAYLOR_SET_H
#define ARIADNE_TAYLOR_SET_H

#include <iosfwd>
#include "container.h"
#include "numeric.h"
#include "vector.h"
#include "set_interface.h"
#include "graphics_interface.h"

#include "taylor_model.h"
#include "taylor_function.h"

#include "box.h"
#include <boost/concept_check.hpp>

namespace Ariadne {

//! \related TaylorConstrainedImageSet \brief The possible types of method used to draw a nonlinear set.
enum DrawingMethod { CURVE_DRAW, BOX_DRAW, AFFINE_DRAW, GRID_DRAW };
//! \related TaylorConstrainedImageSet \brief The type of method currently used to draw a set.
//! HACK: May be replaced by more advanced functionality in the future.
extern DrawingMethod DRAWING_METHOD;
//! \related TaylorConstrainedImageSet \brief The accuracy used to draw a set.
//! HACK: May be replaced by more advanced functionality in the future.
extern unsigned int DRAWING_ACCURACY;

//! \related TaylorConstrainedImageSet \brief The possible types of method used to discretise a nonlinear set.
enum DiscretisationMethod { SUBDIVISION_DISCRETISE, AFFINE_DISCRETISE, CONSTRAINT_DISCRETISE };
//! \related TaylorConstrainedImageSet \brief The type of method currently used to discretise a nonlinear set.
//! HACK: May be replaced by more advanced functionality in the future.
extern DiscretisationMethod DISCRETISATION_METHOD;

class Float;
class Interval;
template<class X> class Vector;
typedef Vector<Interval> IntervalVector;
template<class X> class Matrix;

template<class X> class ScalarFunction;
typedef ScalarFunction<Real> RealScalarFunction;
template<class X> class VectorFunction;
typedef VectorFunction<Real> RealVectorFunction;

class NonlinearConstraint;
template<class X> class TaylorModel;
class ScalarTaylorFunction;
class VectorTaylorFunction;

class AffineSet;

template<class BS> class ListSet;

class Grid;
class GridTreeSet;


class HybridEnclosure;

//! \brief A set of the form \f$x=f(s)\f$ for \f$s\in D\f$ satisfying \f$g(s)\leq0\f$ and \f$h(s)=0\f$.
class TaylorConstrainedImageSet
    : public DrawableInterface
    , public CompactSetInterface
{
    friend class HybridEnclosure;
    Box _domain;
    VectorTaylorFunction _function;
    List<ScalarTaylorFunction> _constraints;
    List<ScalarTaylorFunction> _equations;

    mutable Box _reduced_domain;
    mutable bool _is_fully_reduced;
  public:
    //! \brief Construct a set with \f$D=\emptyset\f$ in \f$\mathbb{R}^0\f$.
    explicit TaylorConstrainedImageSet();
    //! \brief Construct a representation of the box \a bx.
    explicit TaylorConstrainedImageSet(const Box& bx);
    //! \brief Construct the set with parameter domain \a d and image function \a f.
    explicit TaylorConstrainedImageSet(const IntervalVector& d, const RealVectorFunction& f);
    //! \brief Construct the set with parameter domain \a d, image function \a f and constraints \a c.
    explicit TaylorConstrainedImageSet(const IntervalVector& d, const RealVectorFunction& f, const List<NonlinearConstraint>& c);
    //! \brief Construct the set with domain \a d, image function \a f, negative constraints \a g and equality constraints \a h.
    explicit TaylorConstrainedImageSet(const IntervalVector& d, const RealVectorFunction& f, const List<RealScalarFunction>& g, List<RealScalarFunction>& h);
    //! \brief Construct a set with a single constraint \a c. \deprecated Use a list of constraints instead
    explicit TaylorConstrainedImageSet(const IntervalVector& d, const RealVectorFunction& f, const NonlinearConstraint& c);

    //! \brief Construct the set with domain equal to the natural domain of \a f.
    explicit TaylorConstrainedImageSet(const VectorTaylorFunction& f);
    //! \brief Create a dynamically-allocated copy.
    TaylorConstrainedImageSet* clone() const;

    //! \brief The sweeper used to control the accuracy.
    Sweeper sweeper() const;
    //! \brief The parameter domain \f$D\f$.
    Vector<Interval> domain() const;
    //! \brief A subset of the parameter domain containing all feasible points.
    Vector<Interval> reduced_domain() const;
    //! \brief An over-approximation to the image of \f$D\f$ under \f$f\f$.
    Vector<Interval> codomain() const;
    //! \brief The image function \f$f\f$.
    VectorTaylorFunction const& function() const;
    VectorTaylorFunction taylor_function() const;
    RealVectorFunction real_function() const;

    //! \brief Substitutes the expression \f$x_j=v(x_1,\ldots,x_{j-1},x_{j+1}\ldots,x_n)\f$ into the function and constraints.
    //! Requires that \f$v(D_1,\ldots,D_{j-1},D_{j+1}\ldots,D_n) \subset D_j\f$ where \f$D\f$ is the domain.
    void substitute(uint j, ScalarTaylorFunction v);
    //! \brief Substitutes the expression \f$x_j=c\f$ into the function and constraints.
    void substitute(uint j, Float c);
    //! \brief Apply the map \f$r\f$ to the map \f$f\f$.
    void apply_map(RealVectorFunction r);
    //! \brief Apply the map \f$r\f$ to the map \f$f\f$.
    void apply_map(const IntervalVectorFunctionInterface& r);
    //! \brief Apply the map \f$r\f$ to the map \f$f\f$.
    void apply_map(VectorTaylorFunction r);
    //! \brief Apply the flow \f$\phi(x,t)\f$ to the map \f$f\f$.
    void apply_flow(RealVectorFunction phi, Interval time);
    //! \brief Apply the flow \f$\phi(x,t)\f$ to the map \f$f\f$.
    void apply_flow(VectorTaylorFunction phi, Interval time);
    //! \brief Apply the flow \f$\phi(x,h)\f$ to the map \f$f\f$.
    void apply_flow_step(VectorTaylorFunction phi, Float h);

    //! \brief Introduces the constraint \f$c\f$ applied to the state \f$x=f(s)\f$.
    void new_state_constraint(NonlinearConstraint c);
    //! \brief Introduces the constraint \f$c\f$ applied to the parameter \f$s\f$.
    void new_parameter_constraint(NonlinearConstraint c);

    //! \brief Introduces the constraint \f$g(s) \leq 0\f$.
    void new_negative_constraint(RealScalarFunction g);
    void new_negative_constraint(ScalarTaylorFunction g);
    //! \brief Introduces the constraint \f$h(s) = 0\f$.
    void new_zero_constraint(RealScalarFunction h);
    void new_zero_constraint(ScalarTaylorFunction h);
    //! \brief Introduces the constraint \f$h(s) = 0\f$. \deprecated
    void new_equality_constraint(RealScalarFunction h);

    //! \brief The functions \f$g\f$ defining the inequality constraints \f$g(x) \leq 0\f$.
    const List<ScalarTaylorFunction>& negative_constraints() const;
    //! \brief The functions \f$h\f$ defining the equality constraints \f$h(x) = 0\f$.
    const List<ScalarTaylorFunction>& zero_constraints() const;
    //! \brief All equality and inequality constraints.
    List<NonlinearConstraint> constraints() const;

    //! \brief The number of negative constraints.
    uint number_of_constraints() const;
    //! \brief The number of negative constraints.
    uint number_of_negative_constraints() const;
    //! \brief The number of zero constraints.
    uint number_of_zero_constraints() const;

    //! \brief The \a i<sup>th</sup> negative constraint.
    ScalarTaylorFunction negative_constraint(uint i) const;
    //! \brief The \a i<sup>th</sup> zero constraint.
    ScalarTaylorFunction zero_constraint(uint i) const;

    //! \brief  Returns true if \f$g(x)>0\f$ over the whole set,
    //! false \f$g(x)<0\f$ over the whole set,
    //! and indeterminate otherwise.
    tribool satisfies(RealScalarFunction g) const;
    //! \brief Tests if the set satisfies the constraint \a c. Returns \c true if all points in the set satisfy
    //! the constraint, and \c false if no points in the set satisfy the constraint.
    virtual tribool satisfies(NonlinearConstraint c) const;

    //! \brief The dimension of the set.
    uint dimension() const;
    //! \brief The number of parameters i.e. the dimension of the parameter domain.
    uint number_of_parameters() const;
    //! \brief A bounding box for the set.
    Box bounding_box() const;
    //! \brief A point in the image of the <em>unconstrained</em> parameter domain.
    Point centre() const;
    //! \brief An over-approximation to the radius of the set.
    Float radius() const;
    //! \brief Returns \c true if the set is definitely bounded.
    tribool bounded() const;
    //! \brief Returns \c true if the set is provably empty.
    //! May return \c false if the set can (easily) be proved to be nonempty.
    tribool empty() const;
    //! \brief Returns \c true if the set can be shown to be disjoint from \a bx.
    tribool disjoint(const Box& bx) const;
    //! \brief Returns \c true if the set can be shown to be a subset of \a bx..
    tribool inside(const Box& bx) const;
    //! \brief Returns \c true if the set can be shown to be a subset of \a bx..
    tribool subset(const Box& bx) const;

    //! \brief Reduces the size of the effective parameter domain
    //! by pruning away infeasible points. Does not affect the set as a mathematical entity.
    void reduce() const;
    //! \brief Reconditions the set to give an over-approximation with a simpler representation.
    void recondition();
    //! \brief Restrict the parameter domain to \a subdomain.
    //! \details May also restrict the domain of the defining function models,
    //! resulting in more accurate computations.
    void restrict(const IntervalVector& subdomain);
    //! \brief The set obtained by restricting to the \a subdomain.
    TaylorConstrainedImageSet restriction(const Vector<Interval>& subdomain) const;

    //! \brief Compute an outer approximation on the \a grid to the given \a depth.
    GridTreeSet outer_approximation(const Grid& grid, int depth) const;
    //! \brief Compute an outer approximation on the \a grid to the given \a depth
    //! by subdividing the parameter domain. Does not require constraint propagation
    //! or nonlinear programming, but may be inefficient.
    GridTreeSet subdivision_outer_approximation(const Grid& grid, int depth) const;
    //! \brief Compute an outer approximation on the \a grid to the given \a depth
    //! by first computing affine over-approximations of the set.
    GridTreeSet affine_outer_approximation(const Grid& grid, int depth) const;

    //! \brief Adjoin an outer approximation to the given \a depth to the \a paving.
    void adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const;
    //! \brief Adjoin an outer approximation to the given \a depth to the \a paving
    //! by subdividing the parameter domain. Does not require constraint propagation,
    //! but may be inefficient.
    void subdivision_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const;
    //! \brief Adjoin an outer approximation to the given \a depth to the \a paving
    //! by first computing affine over-approximations of the set.
    void affine_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const;
    //! \brief Adjoin an outer approximation to the given \a depth to the \a paving
    //! by using constraint propagation.
    void constraint_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const;
    //! \brief Adjoin an outer approximation to the given \a depth to the \a paving
    //! by using an interior point method to try to find good barrier functions
    //! and using constraint propagation to prove disjointness with cells.
    //! \details Potentially very efficient, but may be unreliable due to the
    //! use of nonlinear programming to find good Lyapounov multipliers for
    //! the constraints.
    void optimal_constraint_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const;

    //! \brief An approximation as an affine set.
    //! \details Most easily computed by dropping all nonlinear terms in the
    //! image and constraint functions. Potentially a very poor approximation.
    AffineSet affine_approximation() const;
    //! \brief An over-approximation as an affine set.
    //! \details Most easily computed by sweeping all nonlinear terms in the
    //! image and constraint function to constant error terms.
    //! Potentially a very poor approximation, but guaranteed to be an over-
    //! approximation.
    AffineSet affine_over_approximation() const;

    //! \brief A collection of parameter subdomains chosen to make the bounding boxes as small as possible.
    List<IntervalVector> splitting_subdomains_zeroth_order() const;
    //! \brief A collection of parameter subdomains chosen to make the set as close to affine as possible.
    List<IntervalVector> splitting_subdomains_first_order() const;
    //! \brief Split into subsets based on the given subdomains.
    List<TaylorConstrainedImageSet> split(const List<IntervalVector>& subdomains);

    //! \brief The direction along which the set should be split to reduce the bounding box.
    uint splitting_index_zeroth_order() const;
    //! \brief Split into two by splitting the parameter domain along
    //! the direction which reduces the size of the bounding box.
    Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet> split_zeroth_order() const;
    //! \brief Split into two by splitting the parameter domain along
    //! the direction which reduces the nonlinearity of the set.
    Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet> split_first_order() const;
    //! \brief Split into two by splitting the parameter domain along a suitably-chosed direction.
    //! Currently defaults to split_zeroth_order().
    Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet> split() const;
    //! \brief Split into two by splitting the parameter domain along
    //! the \a k<sup>th</sup> direction.
    Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet> split(uint k) const;

    //! \brief Draw to a canvas.
    void draw(CanvasInterface&) const;
    //! \brief Draw the bounding box to a canvas. Useful to obtain a quick and rough
    //! image or when all else fails.
    void box_draw(CanvasInterface&) const;
    //! \brief Draw the to a canvas by splitting into small enough pieces that
    //! affine over-approximations yield a good image.
    void affine_draw(CanvasInterface&, uint=1u) const;
    //! \brief Draw the to a canvas by over-approximating on a grid.
    void grid_draw(CanvasInterface&, uint=1u) const;

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream&) const;
  private:
    void _check() const;
    void _subdivision_adjoin_outer_approximation_to(GridTreeSet& gts, const Vector<Interval>& subdomain, uint depth, const Vector<Float>& errors) const;
    void _solve_zero_constraints();
  private:
    friend TaylorConstrainedImageSet product(const TaylorConstrainedImageSet&, const Interval&);
    friend TaylorConstrainedImageSet product(const TaylorConstrainedImageSet&, const Box&);
    friend TaylorConstrainedImageSet product(const TaylorConstrainedImageSet&, const TaylorConstrainedImageSet&);
};

//! \related TaylorConstrainedImageSet \brief Stream output operator.
inline std::ostream& operator<<(std::ostream& os, const TaylorConstrainedImageSet& s) { return s.write(os); }

//! \related TaylorConstrainedImageSet \brief The Cartesian product of a constrained image set with an interval in one dimension.
TaylorConstrainedImageSet product(const TaylorConstrainedImageSet& set, const Interval& ivl);
//! \related TaylorConstrainedImageSet \brief The Cartesian product of a constrained image set with a box.
TaylorConstrainedImageSet product(const TaylorConstrainedImageSet& set, const Box& bx);
//! \related TaylorConstrainedImageSet \brief The Cartesian product of two constrained image sets.
TaylorConstrainedImageSet product(const TaylorConstrainedImageSet& set1, const TaylorConstrainedImageSet& set2);

//! \related TaylorConstrainedImageSet \brief The image of the \a set under the \a function.
TaylorConstrainedImageSet apply(const IntervalVectorFunctionInterface& function, const TaylorConstrainedImageSet& set);
TaylorConstrainedImageSet apply(const RealVectorFunction& function, const TaylorConstrainedImageSet& set);
TaylorConstrainedImageSet apply(const VectorTaylorFunction& function, const TaylorConstrainedImageSet& set);
//! \related TaylorConstrainedImageSet \brief The image of the \a set under the \a function. Does not perform domain-checking.
TaylorConstrainedImageSet unchecked_apply(const VectorTaylorFunction& function, const TaylorConstrainedImageSet& set);

} //namespace Ariadne

#endif /* ARIADNE_TAYLOR_SET_H */
