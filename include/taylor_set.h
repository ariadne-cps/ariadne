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
#include "numeric.h"
#include "vector.h"
#include "set_interface.h"
#include "graphics_interface.h"

#include "taylor_model.h"

namespace Ariadne {

typedef double Float;
class Interval;
template<class X> class Vector;
template<class X> class Matrix;

class ScalarFunction;
class VectorFunction;
class TaylorModel;
class ScalarTaylorFunction;
class VectorTaylorFunction;
class TaylorSet;

template<class BS> class ListSet;

class Zonotope;
class Grid;
class GridTreeSet;

/*! \brief Sets expressed as the image of a box under a polynomial with error bounds.
 *
 *  See also TaylorModel, ScalarTaylorFunction, VectorTaylorFunction.
 */
class TaylorSet
    : public LocatedSetInterface
    , public DrawableInterface
{
  private:
    Vector<TaylorModel> _models;
  public:
    //! \brief Construct the origin in dimension \a d with \a ng generators.
    TaylorSet(uint d=0, uint ng=0);
    //! \brief Construct the image of the box \a d under the function \a f.
    TaylorSet(const VectorFunction& f, const Vector<Interval>& d);
    //! \brief Construct from a list of models giving set as the image of a unit box.
    TaylorSet(const Vector<TaylorModel>& tv);
    //! \brief The box \a bx.
    TaylorSet(const Vector<Interval>& bx);

    //! \brief Construct from raw data in the form of dense polynomial expansions with errors.
    TaylorSet(uint rs, uint as, uint deg, double x0, ...);
    //! \brief Construct from raw data in the form of a polynomial expansion with errors.
    TaylorSet(const Vector< Expansion<Float> >& f, const Vector<Float>& e);

    template<class E> TaylorSet(const ublas::vector_expression<E>& e) {
        *this = TaylorSet(Vector<TaylorModel>(e())); }

    //! \brief Equality operator.
    friend bool operator==(const TaylorSet& ts1, const TaylorSet& ts2);

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
    virtual TaylorSet* clone() const { return new TaylorSet(*this); }
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
    TaylorSet linearise() const;
    //! \brief An over-approximation in the form of a list of boxes.
    ListSet<Box> discretise(const Float& eps) const;
    //! \brief An outer-approximation on a grid.
    GridTreeSet discretise(const Grid& grid, uint depth) const;
    //! \brief Adjoin an outer-approximation on a grid to an existing set.
    GridTreeSet& discretise(GridTreeSet& grid_set, uint depth) const;
    //! \brief An over-approximation with better numerical conditioning;
    //! currently implemented as the orthogonal part of a QR factorisation.
    TaylorSet recondition() const;
    //! \brief Subsume the constant error terms in new generators.
    TaylorSet subsume() const;
    //! \brief Subsume constant error terms of magnitude greater than \a e in new generators.
    TaylorSet subsume(double e) const;
    //! \brief Split by subdividing along generator \a g.
    pair<TaylorSet,TaylorSet> split(uint g) const;
    //! \brief Split by subdividing along a judiciously-chosen generator.
    pair<TaylorSet,TaylorSet> split() const;
    //! \brief Subdivide into sets of maximum radius \a rad.
    ListSet<TaylorSet> subdivide(Float rad) const;

    //! \brief Compute an over-approximation to the image under a function.
    friend TaylorSet apply(const VectorFunction& f, const TaylorSet& s);
    //! \brief Compute an over-approximation to the image under a Taylor function approximation.
    friend TaylorSet apply(const VectorTaylorFunction& f, const TaylorSet& s);
  private:
    Matrix<Float> jacobian() const;
};

TaylorModel apply(const ScalarTaylorFunction& f, const TaylorSet& s);
TaylorSet apply(const VectorTaylorFunction& f, const TaylorSet& s);
TaylorModel apply(const ScalarFunction& f, const TaylorSet& s);
TaylorSet apply(const VectorFunction& f, const TaylorSet& s);

TaylorModel unchecked_apply(const ScalarTaylorFunction& f, const TaylorSet& s);
TaylorSet unchecked_apply(const VectorTaylorFunction& f, const TaylorSet& s);

GridTreeSet outer_approximation(const TaylorSet& set, const Grid& grid, uint depth);
void adjoin_outer_approximation(GridTreeSet& grid_set, const TaylorSet& set, uint depth);
Zonotope zonotope(const TaylorSet& ts);

void standard_draw(CanvasInterface& g, const TaylorSet& ts);
void box_draw(CanvasInterface& g, const TaylorSet& ts);
void affine_draw(CanvasInterface& g, const TaylorSet& ts);
void curve_draw(CanvasInterface& g, const TaylorSet& ts);
void grid_draw(CanvasInterface& g, const TaylorSet& ts);

void plot(const char* fn, const Box& bbx, const TaylorSet& ts);

class DiscreteEvent;

//! \brief A set defined as the image of a subset of a box (the \em domain) under a continuous function.
//! The domain is of the form \f$(x_1,\ldots,x_n,t_0,\ldots,t_k)\f$ where each \f$\delta_i\f$ is positive.
//! The constraints are of the form \f$a(x,t)\geq0\f$, \f$i(x,\ldots,[0,t_k])\leq 0\f$ or \f$g(x,\ldots,[0,t_k])\leq 0 \wedge g(x,\ldots,t_k)=0\f$.
//! We abbreviate the latter as \f$g(x,t,[0,t_k]) <\!\&\!= 0\f$.
class TaylorConstrainedFlowSet
    : public LocatedSetInterface
{
    Vector<Interval> _domain;
    Vector<TaylorModel> _models;
    // Constraints of the form c(x,[0,t])<=0
    std::map<DiscreteEvent,TaylorModel> _invariants;
    // Constraints of the form a(x,t)>=0
    std::map<DiscreteEvent,TaylorModel> _activations;
    // Constraints of the form g(x,t)<&=0
    std::map<DiscreteEvent,TaylorModel> _guards;
  public:
    //! \brief Construct the preimage of \a codom under \a fn.
    TaylorConstrainedFlowSet(const Vector<Interval>& dom, const VectorFunction& fn);
    //! \brief Construct the preimage of \a codom under \a fn.
    TaylorConstrainedFlowSet(const VectorTaylorFunction& fn);
    //! \brief Construct the preimage of \a codom under \a fn.
    TaylorConstrainedFlowSet(const Box& bx);
    //! \brief A dynamically-allocated copy.
    TaylorConstrainedFlowSet* clone() const;

    //! \brief Add a new constraint of the form \f$i(x,\ldots,[0,t_k])\leq0\f$.
    void new_invariant(const DiscreteEvent& e, const ScalarFunction& c);
    //! \brief Add a new constraint of the form \f$g(x,\ldots,[0,t_k])<\!\&\!=0\f$.
    void new_guard(const DiscreteEvent& e, const ScalarFunction& c);
    //! \brief Add a new constraint of the form \f$a(x,\ldots,t_k)\geq0\f$.
    void new_activation(const DiscreteEvent& e, const ScalarFunction& c);

    //! \brief Add a new constraint of the form \f$i(x,\ldots,[0,t_k])\leq0\f$.
    void new_invariant(const DiscreteEvent& e, const ScalarTaylorFunction& c);
    //! \brief Add a new constraint of the form \f$g(x,\ldots,[0,t_k])<\!\&\!=0\f$.
    void new_guard(const DiscreteEvent& e, const ScalarTaylorFunction& c);
    //! \brief Add a new constraint of the form \f$a(x,\ldots,t_k)\geq0\f$.
    void new_activation(const DiscreteEvent& e, const ScalarTaylorFunction& c);

    //! \brief The domain of the set.
    Vector<Interval> domain() const;
    //! \brief The function used to define the set.
    VectorTaylorFunction function() const;

    uint dimension() const;
    tribool empty() const;
    tribool disjoint(const Box&) const;
    tribool overlaps(const Box&) const;
    tribool inside(const Box&) const;
    Box bounding_box() const;
    GridTreeSet outer_approximation(const Grid& grid, uint depth) const;
    std::ostream& write(std::ostream&) const;
  private:
    friend TaylorConstrainedFlowSet apply(const VectorFunction& f, const TaylorConstrainedFlowSet& s);
    friend TaylorConstrainedFlowSet apply(const VectorTaylorFunction& f, const TaylorConstrainedFlowSet& s);
    friend TaylorConstrainedFlowSet apply_flow(const VectorTaylorFunction& f, const TaylorConstrainedFlowSet& s);
  private:
    void _adjoin_outer_approximation_to(GridTreeSet& gts, const Vector<Interval>& subdomain, uint depth) const;
    tribool _empty(const Vector<Interval>& subdomain, uint depth) const;
};


TaylorConstrainedFlowSet apply(const VectorFunction& f, const TaylorConstrainedFlowSet& s);
TaylorConstrainedFlowSet apply(const VectorTaylorFunction& f, const TaylorConstrainedFlowSet& s);
TaylorConstrainedFlowSet apply_flow(const VectorTaylorFunction& f, const TaylorConstrainedFlowSet& s);



} //namespace Ariadne

#endif /* ARIADNE_TAYLOR_SET_H */
