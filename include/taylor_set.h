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

namespace Ariadne {

//! \related TaylorConstrainedImageSet \brief The possible types of method used to draw a nonlinear set.
enum DrawingMethod { CURVE_DRAW, BOX_DRAW, AFFINE_DRAW, GRID_DRAW };
//! \related TaylorConstrainedImageSet \brief The type of method currently used to draw a set.
//! HACK: May be replaced by more advanced functionality in the future.
extern DrawingMethod DRAWING_METHOD;
//! \related TaylorConstrainedImageSet \brief The accuracy used to draw a set.
//! HACK: May be replaced by more advanced functionality in the future.
extern unsigned int DRAWING_ACCURACY;

class Float;
class Interval;
template<class X> class Vector;
typedef Vector<Interval> IntervalVector;
template<class X> class Matrix;

class ScalarFunction;
class VectorFunction;
class NonlinearConstraint;
class TaylorModel;
class ScalarTaylorFunction;
class VectorTaylorFunction;
class TaylorImageSet;

class AffineSet;

template<class BS> class ListSet;

class Zonotope;
class Grid;
class GridTreeSet;

/*! \brief Sets expressed as the image of a box under a polynomial with error bounds.
 *
 *  See also TaylorModel, ScalarTaylorFunction, VectorTaylorFunction.
 */
class TaylorImageSet
    : public LocatedSetInterface
    , public DrawableInterface
{
  private:
    Vector<TaylorModel> _models;
  public:
    //! \brief Construct the origin in dimension \a d with \a ng generators.
    TaylorImageSet(uint d=0, uint ng=0);
    //! \brief Construct the image of the box \a d under the function \a f.
    TaylorImageSet(const VectorFunction& f, const Vector<Interval>& d);
    //! \brief Construct from a list of models giving set as the image of a unit box.
    TaylorImageSet(const Vector<TaylorModel>& tv);
    //! \brief The box \a bx.
    TaylorImageSet(const Vector<Interval>& bx);

    //! \brief Construct from raw data in the form of dense polynomial expansions with errors.
    TaylorImageSet(uint rs, uint as, uint deg, double x0, ...);
    //! \brief Construct from raw data in the form of a polynomial expansion with errors.
    TaylorImageSet(const Vector< Expansion<Float> >& f, const Vector<Float>& e);

    template<class E> TaylorImageSet(const ublas::vector_expression<E>& e) {
        *this = TaylorImageSet(Vector<TaylorModel>(e())); }

    //! \brief Equality operator.
    friend bool operator==(const TaylorImageSet& ts1, const TaylorImageSet& ts2);

    //! \brief Set the accuracy parameters.
    void set_accuracy(shared_ptr<TaylorModel::Accuracy> acc_ptr);
    //! \brief Get the accuracy parameters.
    shared_ptr<TaylorModel::Accuracy> accuracy_ptr() const;

    //! \brief The dimension the space lies in.
    uint dimension() const { return this->_models.size(); }
    //! \brief The number of generators of the set.
    uint generators_size() const { assert(this->_models.size()>0); return this->_models[0].argument_size(); }

    //! \brief The Taylor models used to define the set.
    const Vector<TaylorModel>& models() const { return this->_models; }
    //! \brief The domain of which the set is an image.
    Vector<Interval> domain() const { return Vector<Interval>(this->generators_size(),Interval(-1,+1)); }
    //! \brief A box bounding the range of the generating function.
    Vector<Interval> range() const {
        Vector<Interval> result(this->_models.size());
        for(uint i=0; i!=this->_models.size(); ++i) {
            result[i]=this->_models[i].range(); }
        return result; }

    uint size() const { return this->_models.size(); }
    uint result_size() const { return this->_models.size(); }
    uint argument_size() const { return this->_models[0].argument_size(); }
    const TaylorModel& operator[](uint i) const { return this->_models[i]; }
    TaylorModel& operator[](uint i) { return this->_models[i]; }


    //! \brief Create a dynamically-allocated copy.
    virtual TaylorImageSet* clone() const { return new TaylorImageSet(*this); }
    //! \brief Tests if the set is empty.
    virtual tribool empty() const;
    //! \brief Tests if the set is disjoint from a box.
    virtual tribool disjoint(const Box&) const;
    //! \brief Tests if the set overlaps a box.
    virtual tribool overlaps(const Box&) const;
    //! \brief Tests if the set is a subset of a box.
    virtual tribool inside(const Box&) const;
    //! \brief A bounding box for the set.
    virtual Box bounding_box() const;
    //! \brief Draw on a two-dimensional canvas.
    virtual void draw(CanvasInterface& c) const;
    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const;

    //! \brief The centre of the set.
    Vector<Float> centre() const;
    //! \brief The radius of the set.
    Float radius() const;
    //! \brief An over-approximation in the form of a zonotope.
    TaylorImageSet linearise() const;
    //! \brief An over-approximation in the form of a list of boxes.
    ListSet<Box> discretise(const Float& eps) const;
    //! \brief An outer-approximation on a grid.
    GridTreeSet discretise(const Grid& grid, uint depth) const;
    //! \brief Adjoin an outer-approximation on a grid to an existing set.
    GridTreeSet& discretise(GridTreeSet& grid_set, uint depth) const;
    //! \brief An over-approximation with better numerical conditioning;
    //! currently implemented as the orthogonal part of a QR factorisation.
    TaylorImageSet recondition() const;
    //! \brief Subsume the constant error terms in new generators.
    TaylorImageSet subsume() const;
    //! \brief Subsume constant error terms of magnitude greater than \a e in new generators.
    TaylorImageSet subsume(double e) const;
    //! \brief Split by subdividing along generator \a g.
    pair<TaylorImageSet,TaylorImageSet> split(uint g) const;
    //! \brief Split by subdividing along a judiciously-chosen generator.
    pair<TaylorImageSet,TaylorImageSet> split() const;
    //! \brief Subdivide into sets of maximum radius \a rad.
    ListSet<TaylorImageSet> subdivide(Float rad) const;

    //! \brief Compute an over-approximation to the image under a function.
    friend TaylorImageSet apply(const VectorFunction& f, const TaylorImageSet& s);
    //! \brief Compute an over-approximation to the image under a Taylor function approximation.
    friend TaylorImageSet apply(const VectorTaylorFunction& f, const TaylorImageSet& s);
  private:
    Matrix<Float> jacobian() const;
};

TaylorModel apply(const ScalarTaylorFunction& f, const TaylorImageSet& s);
TaylorImageSet apply(const VectorTaylorFunction& f, const TaylorImageSet& s);
TaylorModel apply(const ScalarFunction& f, const TaylorImageSet& s);
TaylorImageSet apply(const VectorFunction& f, const TaylorImageSet& s);

TaylorModel unchecked_apply(const ScalarTaylorFunction& f, const TaylorImageSet& s);
TaylorImageSet unchecked_apply(const VectorTaylorFunction& f, const TaylorImageSet& s);

GridTreeSet outer_approximation(const TaylorImageSet& set, const Grid& grid, uint depth);
void adjoin_outer_approximation(GridTreeSet& grid_set, const TaylorImageSet& set, uint depth);
Zonotope zonotope(const TaylorImageSet& ts);

TaylorImageSet product(const TaylorImageSet& set, const Interval& ivl);
TaylorImageSet product(const TaylorImageSet& set, const Box& bx);
TaylorImageSet product(const TaylorImageSet& set1, const TaylorImageSet& set2);

void standard_draw(CanvasInterface& g, const TaylorImageSet& ts);
void box_draw(CanvasInterface& g, const TaylorImageSet& ts);
void affine_draw(CanvasInterface& g, const TaylorImageSet& ts);
void curve_draw(CanvasInterface& g, const TaylorImageSet& ts);
void grid_draw(CanvasInterface& g, const TaylorImageSet& ts);


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
  public:
    //! \brief Construct a set with \f$D=\emptyset\f$ in \f$\mathbb{R}^0\f$.
    TaylorConstrainedImageSet();
    //! \brief Construct a representation of the box \a bx.
    TaylorConstrainedImageSet(const Box& bx);
    //! \brief Construct the set with parameter domain \a d and image function \a f.
    TaylorConstrainedImageSet(const IntervalVector& d, const VectorFunction& f);
    //! \brief Construct the set with parameter domain \a d, image function \a f and constraints \a c.
    TaylorConstrainedImageSet(const IntervalVector& d, const VectorFunction& f, const List<NonlinearConstraint>& c);
    //! \brief Construct the set with domain \a d, image function \a f, negative constraints \a g and equality constraints \a h.
    TaylorConstrainedImageSet(const IntervalVector& d, const VectorFunction& f, const List<ScalarFunction>& g, List<ScalarFunction>& h);
    //! \brief Construct a set with a single constraint \a c. \deprecated Use a list of constraints instead
    TaylorConstrainedImageSet(const IntervalVector& d, const VectorFunction& f, const NonlinearConstraint& c);

    //! \brief Construct the set with domain equal to the natural domain of \a f.
    explicit TaylorConstrainedImageSet(const VectorTaylorFunction& f);
    //! \brief Create a dynamically-allocated copy.
    TaylorConstrainedImageSet* clone() const;

    //! \brief The parameter domain \f$D\f$.
    Vector<Interval> domain() const;
    //! \brief An over-approximation to the image of \f$D\f$ under \f$f\f$.
    Vector<Interval> codomain() const;
    //! \brief The image function \f$f\f$.
    VectorTaylorFunction const& function() const;
    VectorTaylorFunction taylor_function() const;
    VectorFunction real_function() const;

    //! \brief Substitutes the expression \f$x_j=v(x_1,\ldots,x_{j-1},x_{j+1}\ldots,x_n)\f$ into the function and constraints.
    //! Requires that \f$v(D_1,\ldots,D_{j-1},D_{j+1}\ldots,D_n) \subset D_j\f$ where \f$D\f$ is the domain.
    void substitute(uint j, ScalarTaylorFunction v);
    //! \brief Substitutes the expression \f$x_j=c\f$ into the function and constraints.
    void substitute(uint j, Float c);
    //! \brief Apply the map \f$r\f$ to the map \f$f\f$.
    void apply_map(VectorFunction r);
    //! \brief Apply the map \f$r\f$ to the map \f$f\f$.
    void apply_map(VectorTaylorFunction r);
    //! \brief Apply the flow \f$\phi(x,t)\f$ to the map \f$f\f$.
    void apply_flow(VectorFunction phi, Interval time);
    //! \brief Apply the flow \f$\phi(x,t)\f$ to the map \f$f\f$.
    void apply_flow(VectorTaylorFunction phi, Interval time);
    //! \brief Apply the flow \f$\phi(x,h)\f$ to the map \f$f\f$.
    void apply_flow_step(VectorTaylorFunction phi, Float h);

    //! \brief Introduces the constraint \f$c\f$ applied to \f$x=f(s)\f$.
    void new_state_constraint(NonlinearConstraint c);

    //! \brief Introduces the constraint \f$g(s) \leq 0\f$.
    void new_negative_constraint(ScalarFunction g);
    void new_negative_constraint(ScalarTaylorFunction g);
    //! \brief Introduces the constraint \f$h(s) = 0\f$.
    void new_zero_constraint(ScalarFunction h);
    void new_zero_constraint(ScalarTaylorFunction h);
    //! \brief Introduces the constraint \f$h(s) = 0\f$. \deprecated
    void new_equality_constraint(ScalarFunction h);

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
    tribool satisfies(ScalarFunction g) const;
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
    //! by pruning away infeasible points.
    void reduce() const;

    //! \brief Compute an outer approximation on the \a grid to the given \a depth.
    GridTreeSet outer_approximation(const Grid& grid, int depth) const;
    //! \brief Compute an outer approximation on the \a grid to the given \a depth
    //! by subdividing the parameter domain. Does not require constraint propagation,
    //! but may be inefficient.
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

    //! \brief Split into two by splitting the parameter domain along
    //! the \a k<sup>th</sup> direction.
    Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet> split(uint k) const;
    //! \brief Restrict the parameter domain to \a subdomain.
    //! \details May also restrict the domain of the defining function models,
    //! resulting in more accurate computations.
    TaylorConstrainedImageSet restriction(const Vector<Interval>& subdomain) const;

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


} //namespace Ariadne

#endif /* ARIADNE_TAYLOR_SET_H */
