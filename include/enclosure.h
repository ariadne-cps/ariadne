/***************************************************************************
 *            enclosure.h
 *
 *  Copyright  2011  Pieter Collins
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

/*! \file enclosure.h
 *  \brief Enclosure sets for continuous systems
 */

#ifndef ARIADNE_ENCLOSURE_H
#define ARIADNE_ENCLOSURE_H

#include <iosfwd>
#include "container.h"
#include "numeric.h"
#include "vector.h"
#include "set_interface.h"
#include "graphics_interface.h"

#include "function_model.h"

#include "box.h"
#include <boost/concept_check.hpp>

#ifndef ARIADNE_TAYLOR_SET_H

namespace Ariadne {

//! \related Enclosure \brief The possible types of method used to draw a nonlinear set.
enum DrawingMethod { CURVE_DRAW, BOX_DRAW, AFFINE_DRAW, GRID_DRAW };
//! \related Enclosure \brief The type of method currently used to draw a set.
//! HACK: May be replaced by more advanced functionality in the future.
extern DrawingMethod DRAWING_METHOD;
//! \related Enclosure \brief The accuracy used to draw a set.
//! HACK: May be replaced by more advanced functionality in the future.
extern unsigned int DRAWING_ACCURACY;

//! \related Enclosure \brief The possible types of method used to discretise a nonlinear set.
enum DiscretisationMethod { SUBDIVISION_DISCRETISE, AFFINE_DISCRETISE, CONSTRAINT_DISCRETISE };
//! \related Enclosure \brief The type of method currently used to discretise a nonlinear set.
//! HACK: May be replaced by more advanced functionality in the future.
extern DiscretisationMethod DISCRETISATION_METHOD;

} // namespace Ariadne

#endif

namespace Ariadne {

class Float;
class Interval;
template<class X> class Vector;
typedef Vector<Interval> IntervalVector;
template<class X> class Matrix;

template<class X> class ScalarFunction;
typedef ScalarFunction<Interval> IntervalScalarFunction;
template<class X> class VectorFunction;
typedef VectorFunction<Interval> IntervalVectorFunction;

template<class X> class FunctionModelFactoryInterface;
typedef FunctionModelFactoryInterface<Interval> IntervalFunctionModelFactoryInterface;

template<class X, class R> class Constraint;
typedef Constraint<RealScalarFunction,Real> RealConstraint;
typedef Constraint<IntervalScalarFunction,Float> IntervalConstraint;

class ValidatedAffineConstrainedImageSet;
class BoundedConstraintSet;

template<class BS> class ListSet;

class Grid;
class PavingInterface;

typedef Constraint<IntervalScalarFunctionModel,Float> IntervalConstraintModel;

//! \brief A set of the form \f$x=f(s)\f$ for \f$s\in D\f$ satisfying \f$g(s)\leq0\f$ and \f$h(s)=0\f$.
class Enclosure
    : public DrawableInterface
    , public CompactSetInterface
{
    Box _domain;
    IntervalVectorFunctionModel _space_function;
    IntervalScalarFunctionModel _time_function;
    IntervalScalarFunctionModel _dwell_time_function;
    List<IntervalConstraintModel> _constraints;
    IntervalFunctionModelFactoryInterface* _function_factory_ptr;
    mutable Box _reduced_domain;
    mutable bool _is_fully_reduced;
  public:
    //! \brief Construct a set with \f$D=\emptyset\f$ in \f$\mathbb{R}^0\f$.
    explicit Enclosure();
    //! \brief Construct a representation of the box \a bx.
    explicit Enclosure(const Box& bx, const IntervalFunctionModelFactoryInterface& fac);
    //! \brief Construct the set with parameter domain \a d and image function \a f.
    explicit Enclosure(const IntervalVector& d, const IntervalVectorFunction& f, const IntervalFunctionModelFactoryInterface& fac);
    //! \brief Construct the set with parameter domain \a d, image function \a f and constraints \a c.
    explicit Enclosure(const IntervalVector& d, const IntervalVectorFunction& f, const List<IntervalConstraint>& c, const IntervalFunctionModelFactoryInterface& fac);
    //! \brief Construct the set with parameter domain \a d, image function \a sf, time function \a tf and constraints \a c.
    explicit Enclosure(const IntervalVector& d, const IntervalVectorFunction& sf, const IntervalScalarFunction& tf, const List<IntervalConstraint>& c, const IntervalFunctionModelFactoryInterface& fac);
    //! \brief Construct the set with domain \a d, space function \a sf, time function \a tf, negative constraints \a g and equality constraints \a h.
    //!   (Not currently implemented.)
    explicit Enclosure(const IntervalVector& d, const IntervalVectorFunction& sf, const IntervalScalarFunction tf, const IntervalVectorFunction& g, const IntervalVectorFunction& h, const IntervalFunctionModelFactoryInterface& fac);
    //! \brief Construct from an exact bounded constraint \a set.
    explicit Enclosure(const BoundedConstraintSet& set, const IntervalFunctionModelFactoryInterface& fac);

    //! \brief Create a dynamically-allocated copy.
    Enclosure* clone() const;

    //! \brief The class used to create new function instances.
    const IntervalFunctionModelFactoryInterface& function_factory() const;
    //! \brief The parameter domain \f$D\f$.
    Vector<Interval> domain() const;
    Vector<Interval> parameter_domain() const;
    //! \brief A subset of the parameter domain containing all feasible points.
    Vector<Interval> reduced_domain() const;
    //! \brief An over-approximation to the image of \f$D\f$ under \f$f\f$.
    Vector<Interval> codomain() const;
    //! \brief The image function \f$f\f$.
    IntervalVectorFunctionModel const& function() const;
    IntervalVectorFunctionModel const& space_function() const;
    IntervalScalarFunctionModel const& time_function() const;
    IntervalScalarFunctionModel const& dwell_time_function() const;
    IntervalVectorFunctionModel const constraint_function() const;
    IntervalVector const constraint_bounds() const;

    //! \brief Introduces a new parameter with values in the interval \a ivl. The set itself does not change.
    void new_parameter(Interval ivl);
    //! \brief Introduces a new independent variable with values in the interval \a ivl.
    //! Equivalent to constructing the set \f$S\times I\f$.
    void new_variable(Interval ivl);
    //! \brief Substitutes the expression \f$x_j=v(x_1,\ldots,x_{j-1},x_{j+1}\ldots,x_n)\f$ into the function and constraints.
    //! Requires that \f$v(D_1,\ldots,D_{j-1},D_{j+1}\ldots,D_n) \subset D_j\f$ where \f$D\f$ is the domain.
    void substitute(uint j, IntervalScalarFunctionModel v);
    //! \brief Substitutes the expression \f$x_j=c\f$ into the function and constraints.
    void substitute(uint j, Float c);

    //! \brief Apply the map \f$r\f$ to the map \f$f\f$.
    void apply_map(IntervalVectorFunction r);
    //! \brief Apply the flow \f$\phi(x,h)\f$ to the map \f$f\f$.
    void apply_fixed_evolve_step(IntervalVectorFunction phi, Float h);
    //! \brief Apply the flow \f$xi'(s)=\phi(\xi(s),\epsilon(\xi(s)))\f$, \f$\tau'(s)=\tau(s)+\epsilon(\xi(s))\f$.
    void apply_space_evolve_step(IntervalVectorFunction phi, IntervalScalarFunction elps);
    //! \brief Apply the flow \f$xi'(s)=\phi(\xi(s),\epsilon(\xi(s),\tau(s)))\f$, \f$\tau'(s)=\tau(s)+\epsilon(\xi(s),\tau(s))\f$.
    void apply_spacetime_evolve_step(IntervalVectorFunction phi, IntervalScalarFunction elps);
    //! \brief Set \f$\xi'(s)=\phi(\xi(s),\epsilon(s))\f$ and \f$\tau'(s)=\tau(s)+\epsilon(s)\f$.
    void apply_parameter_evolve_step(IntervalVectorFunction phi, IntervalScalarFunction elps);
    //! \brief Set \f$\xi'(s)=\phi(\xi(s),\omega(s)-\tau(s))\f$ and \f$\tau'(s)=\omega(s)\f$.
    void apply_finishing_parameter_evolve_step(IntervalVectorFunction phi, IntervalScalarFunction omega);
    //! \brief Set \f$\xi'(s,r)=\phi(\xi(s),r)\f$ and \f$\tau'(s,r)=\tau(s)+r\f$ for $r\leq h.
    void apply_full_reach_step(IntervalVectorFunctionModel phi);
    //! \brief Apply the flow \f$xi'(s,r)=\phi(\xi(s),r)\f$, \f$\tau'(s,r)=\tau(s)+r\f$, \f$r\leq\elps\f$
    void apply_spacetime_reach_step(IntervalVectorFunctionModel phi, IntervalScalarFunction elps);
    //! \brief Set \f$\xi'(s,r)=\phi(\xi(s),r)\f$ and \f$\tau'(s,r)=\tau(s)+r\f$ for $r-elps(s)\leq 0$.
    void apply_parameter_reach_step(IntervalVectorFunctionModel phi, IntervalScalarFunction elps);
/*
    //! \brief Apply the flow \f$\phi(x,t)\f$ for \f$t\in[0,h]\f$
    void apply_reach_step(IntervalVectorFunction phi, Float h);
    //! \brief Apply the flow \f$\phi(x,t)\f$ for \f$t\in[0,\max(h,\epsilon(x))]\f$
    void apply_reach_step(IntervalVectorFunction phi, IntervalScalarFunction elps);
*/

    //! \brief Introduces the constraint \f$c\f$ applied to the state \f$x=f(s)\f$.
    void new_state_constraint(IntervalConstraint c);
    //! \brief Introduces the constraint \f$c\f$ applied to the parameter \f$s\f$.
    void new_parameter_constraint(IntervalConstraint c);

    //! \brief Introduces the constraint \f$-g(\xi(s)) \leq 0\f$.
    void new_positive_state_constraint(IntervalScalarFunction g);
    //! \brief Introduces the constraint \f$g(\xi(s)) \leq 0\f$.
    void new_negative_state_constraint(IntervalScalarFunction g);
    //! \brief Introduces the constraint \f$h(\xi(s)) = 0\f$.
    void new_zero_state_constraint(IntervalScalarFunction h);
    //! \brief Introduces the constraint \f$g(s) \leq 0\f$.
    void new_negative_parameter_constraint(IntervalScalarFunction g);
    //! \brief Introduces the constraint \f$h(s) = 0\f$.
    void new_zero_parameter_constraint(IntervalScalarFunction h);

    //! \brief The number of negative constraints.
    uint number_of_constraints() const;
    //! \brief All equality and inequality constraints.
    List<IntervalConstraintModel> const& constraint_models() const;
    //! \brief All equality and inequality constraints.
    List<IntervalConstraint> const constraints() const;
    //! \brief The \a i<sup>th</sup> constraint.
    IntervalConstraintModel const& constraint(uint i) const;

    //! \brief  Returns true if \f$g(x)>0\f$ over the whole set,
    //! false \f$g(x)<0\f$ over the whole set,
    //! and indeterminate otherwise.
    tribool satisfies(IntervalScalarFunction g) const;
    //! \brief Tests if the set satisfies the constraint \a c. Returns \c true if all points in the set satisfy
    //! the constraint, and \c false if no points in the set satisfy the constraint.
    virtual tribool satisfies(IntervalConstraint c) const;

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
    tribool separated(const Box& bx) const;
    //! \brief Returns \c true if the set can be shown to be a subset of \a bx..
    tribool inside(const Box& bx) const;
    //! \brief Returns \c true if the set can be shown to be a subset of \a bx..
    tribool subset(const Box& bx) const;

    //! \brief Reduces the size of the effective parameter domain
    //! by pruning away infeasible points. Does not affect the set as a mathematical entity.
    void reduce() const;
    //! \brief Reconditions the set to give an over-approximation with a simpler representation.
    void recondition();
    //! \brief Simplifies the representation by changing all uniform errors into independent variables.
    void uniform_error_recondition();
    //! \brief Simplifies the representation by choosing most significant independent variables to keep, and merging the rest into a single error for each component.
    void kuhn_recondition();
    //! \brief Restrict the parameter domain to \a subdomain.
    //! \details May also restrict the domain of the defining function models,
    //! resulting in more accurate computations.
    void restrict(const IntervalVector& subdomain);
    //! \brief The set obtained by restricting to the \a subdomain.
    Enclosure restriction(const Vector<Interval>& subdomain) const;

    //! \brief Compute an outer approximation on the \a grid to the given \a depth.
    GridTreeSet outer_approximation(const Grid& grid, int depth) const;

    //! \brief Adjoin an outer approximation to the given \a depth to the \a paving.
    void adjoin_outer_approximation_to(PavingInterface& paving, int depth) const;
    //! \brief Adjoin an outer approximation to the given \a depth to the \a paving
    //! by subdividing the parameter domain. Does not require constraint propagation,
    //! but may be inefficient.
    void subdivision_adjoin_outer_approximation_to(PavingInterface& paving, int depth) const;
    //! \brief Adjoin an outer approximation to the given \a depth to the \a paving
    //! by first computing affine over-approximations of the set.
    void affine_adjoin_outer_approximation_to(PavingInterface& paving, int depth) const;
    //! \brief Adjoin an outer approximation to the given \a depth to the \a paving
    //! by using constraint propagation.
    void constraint_adjoin_outer_approximation_to(PavingInterface& paving, int depth) const;
    //! \brief Adjoin an outer approximation to the given \a depth to the \a paving
    //! by using an interior point method to try to find good barrier functions
    //! and using constraint propagation to prove disjointness with cells.
    //! \details Potentially very efficient, but may be unreliable due to the
    //! use of nonlinear programming to find good Lyapounov multipliers for
    //! the constraints.
    void optimal_constraint_adjoin_outer_approximation_to(PavingInterface& paving, int depth) const;

    //! \brief An approximation as an affine set.
    //! \details Most easily computed by dropping all nonlinear terms in the
    //! image and constraint functions. Potentially a very poor approximation.
    ValidatedAffineConstrainedImageSet affine_approximation() const;
    //! \brief An over-approximation as an affine set.
    //! \details Most easily computed by sweeping all nonlinear terms in the
    //! image and constraint function to constant error terms.
    //! Potentially a very poor approximation, but guaranteed to be an over-
    //! approximation.
    ValidatedAffineConstrainedImageSet affine_over_approximation() const;

    //! \brief A collection of parameter subdomains chosen to make the bounding boxes as small as possible.
    List<IntervalVector> splitting_subdomains_zeroth_order() const;
    //! \brief A collection of parameter subdomains chosen to make the set as close to affine as possible.
    List<IntervalVector> splitting_subdomains_first_order() const;
    //! \brief Split into subsets based on the given subdomains.
    List<Enclosure> split(const List<IntervalVector>& subdomains);

    //! \brief The direction along which the set should be split to reduce the bounding box.
    uint splitting_index_zeroth_order() const;
    //! \brief Split into two by splitting the parameter domain along
    //! the direction which reduces the size of the bounding box.
    Pair<Enclosure,Enclosure> split_zeroth_order() const;
    //! \brief Split into two by splitting the parameter domain along
    //! the direction which reduces the nonlinearity of the set.
    Pair<Enclosure,Enclosure> split_first_order() const;
    //! \brief Split into two by splitting the parameter domain along a suitably-chosed direction.
    //! Currently defaults to split_zeroth_order().
    Pair<Enclosure,Enclosure> split() const;
    //! \brief Split into two by splitting the parameter domain along
    //! the \a k<sup>th</sup> direction.
    Pair<Enclosure,Enclosure> split(uint k) const;


    //! \brief Draw to a canvas.
    void draw(CanvasInterface& c, const Projection2d& p) const;
    //! \brief Draw the bounding box to a canvas. Useful to obtain a quick and rough
    //! image or when all else fails.
    void box_draw(CanvasInterface&, const Projection2d& p) const;
    //! \brief Draw the to a canvas by splitting into small enough pieces that
    //! affine over-approximations yield a good image.
    void affine_draw(CanvasInterface&, const Projection2d& p, uint=1u) const;
    //! \brief Draw the to a canvas by over-approximating on a grid.
    void grid_draw(CanvasInterface&, const Projection2d& p, uint=1u) const;

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream&) const;
  private:
    void _check() const;
    void _solve_zero_constraints();
    RealVectorFunction real_function() const;
  private:
    friend Enclosure product(const Enclosure&, const Interval&);
    friend Enclosure product(const Enclosure&, const Box&);
    friend Enclosure product(const Enclosure&, const Enclosure&);
};

//! \related Enclosure \brief Stream output operator.
inline std::ostream& operator<<(std::ostream& os, const Enclosure& s) { return s.write(os); }

//! \related Enclosure \brief The Cartesian product of a constrained image set with an interval in one dimension.
Enclosure product(const Enclosure& set, const Interval& ivl);
//! \related Enclosure \brief The Cartesian product of a constrained image set with a box.
Enclosure product(const Enclosure& set, const Box& bx);
//! \related Enclosure \brief The Cartesian product of two constrained image sets.
//! \precondition The time function of each set is constant with the same value.
Enclosure product(const Enclosure& set1, const Enclosure& set2);

//! \related Enclosure \brief The image of the \a set under the \a function.
Enclosure apply(const IntervalVectorFunction& function, const Enclosure& set);
//! \related Enclosure \brief The image of the \a set under the \a function. Does not perform domain-checking.
Enclosure unchecked_apply(const IntervalVectorFunctionModel& function, const Enclosure& set);

} //namespace Ariadne

#endif /* ARIADNE_TAYLOR_SET_H */
