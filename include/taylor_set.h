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

namespace Ariadne {

typedef double Float;
class Interval;
template<class X> class Vector;
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

void standard_draw(CanvasInterface& g, const TaylorImageSet& ts);
void box_draw(CanvasInterface& g, const TaylorImageSet& ts);
void affine_draw(CanvasInterface& g, const TaylorImageSet& ts);
void curve_draw(CanvasInterface& g, const TaylorImageSet& ts);
void grid_draw(CanvasInterface& g, const TaylorImageSet& ts);

void plot(const char* fn, const Box& bbx, const TaylorImageSet& ts);

class HybridEnclosure;

//! \brief A set of the form \f$x=f(s)\f$ for \f$s\in D\f$ satisfying \f$g(s)\leq0\f$ and \f$h(s)=0\f$.
class TaylorConstrainedImageSet
    : public DrawableInterface
{

    Vector<Interval> _domain;
    VectorTaylorFunction _function;
    List<ScalarTaylorFunction> _constraints;
    List<ScalarTaylorFunction> _equations;
  public:
    TaylorConstrainedImageSet();
    TaylorConstrainedImageSet(Box);
    TaylorConstrainedImageSet(Box, VectorFunction);
    TaylorConstrainedImageSet(Box, VectorFunction, List<NonlinearConstraint>);
    TaylorConstrainedImageSet(Box, VectorFunction, NonlinearConstraint);
    TaylorConstrainedImageSet(Box, const List<ScalarFunction>&);

    TaylorConstrainedImageSet(const VectorTaylorFunction&);
    TaylorConstrainedImageSet* clone() const;

    Vector<Interval> domain() const;
    Vector<Interval> codomain() const;
    VectorFunction function() const;
    VectorTaylorFunction const& taylor_function() const;

    //! \brief Substitutes the expression \f$x_j=v(x_1,\ldots,x_{j-1},x_{j+1}\ldots,x_n)\f$ into the function and constraints.
    //! Requires that \f$v(D_1,\ldots,D_{j-1},D_{j+1}\ldots,D_n) \subset D_j\f$ where \f$D\f$ is the domain.
    void substitute(uint j, ScalarTaylorFunction v);
    //! \brief Apply the map \f$r\f$ to the map \f$f\f$.
    void apply_map(VectorFunction r);
    //! \brief Apply the flow \f$\phi(x,t)\f$ to the map \f$f\f$.
    void apply_flow(VectorFunction phi, Interval time);
    //! \brief Apply the flow \f$\phi(x,t)\f$ to the map \f$f\f$.
    void apply_flow(VectorTaylorFunction phi, Interval time);

    //! \brief Introduces the constraint \f$g(s) \leq 0\f$.
    void new_negative_constraint(ScalarFunction g);
    void new_negative_constraint(ScalarTaylorFunction g);
    //! \brief Introduces the constraint \f$h(s) = 0\f$.
    void new_equality_constraint(ScalarFunction h);

    //! \brief  Returns true if \f$g(x)\leq0\f$ over the whole set,
    //! false if the constraint is negative over the whole set,
    //! and indeterminate otherwise.
    tribool satisfies(ScalarFunction g) const;

    uint dimension() const;
    uint number_of_parameters() const;
    Box bounding_box() const;
    Point centre() const;
    Float radius() const;
    tribool bounded() const;
    tribool empty() const;
    tribool disjoint(Box bx) const;

    GridTreeSet outer_approximation(const Grid& grid, int depth) const;
    GridTreeSet subdivision_outer_approximation(const Grid& grid, int depth) const;
    GridTreeSet affine_outer_approximation(const Grid& grid, int depth) const;

    void adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const;
    void subdivision_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const;
    void affine_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const;
    void constraint_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const;

    AffineSet affine_approximation() const;
    AffineSet affine_over_approximation() const;

    Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet> split(uint dim) const;
    TaylorConstrainedImageSet restriction(const Vector<Interval>& subdomain) const;

    void draw(CanvasInterface&) const;
    void box_draw(CanvasInterface&) const;
    void affine_draw(CanvasInterface&, uint=1u) const;
    std::ostream& write(std::ostream&) const;
  private:
    void _check() const;
    void _subdivision_adjoin_outer_approximation_to(GridTreeSet& gts, const Vector<Interval>& subdomain, uint depth, const Vector<Float>& errors) const;
    void _solve_zero_constraints();
};

inline std::ostream& operator<<(std::ostream& os, const TaylorConstrainedImageSet& s) { return s.write(os); }


} //namespace Ariadne

#endif /* ARIADNE_TAYLOR_SET_H */
