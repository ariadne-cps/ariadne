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

#include "list_set.h"
#include "taylor_model.h"

namespace Ariadne {

typedef double Float;
class Interval;
template<class X> class Vector;
template<class X> class Matrix;

class ExpressionInterface;
class FunctionInterface;
class TaylorModel;
class TaylorExpression;
class TaylorFunction;
class TaylorSet;

class Zonotope;
class Grid;
class GridTreeSet;
class GraphicsInterface;

/*! \brief Sets expressed as the image of a box under a polynomial with error bounds.
 *
 *  See also TaylorModel, TaylorExpression, TaylorFunction.
 */
class TaylorSet
    : public LocatedSetInterface
{
  private:
    Vector<TaylorModel> _models;
  public:
    //! \brief Construct the origin in dimension \a d with \a ng generators.
    TaylorSet(uint d=0, uint ng=0);
    //! \brief Construct the image of the box \a d under the function \a f.
    TaylorSet(const FunctionInterface& f, const Vector<Interval>& d);
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
    friend TaylorSet apply(const FunctionInterface& f, const TaylorSet& s);
    //! \brief Compute an over-approximation to the image under a Taylor function approximation.
    friend TaylorSet apply(const TaylorFunction& f, const TaylorSet& s);
  private:
    Matrix<Float> jacobian() const;
};

TaylorModel apply(const TaylorExpression& f, const TaylorSet& s);
TaylorSet apply(const TaylorFunction& f, const TaylorSet& s);
TaylorModel apply(const ExpressionInterface& f, const TaylorSet& s);
TaylorSet apply(const FunctionInterface& f, const TaylorSet& s);

TaylorModel unchecked_apply(const TaylorExpression& f, const TaylorSet& s);
TaylorSet unchecked_apply(const TaylorFunction& f, const TaylorSet& s);

GridTreeSet outer_approximation(const TaylorSet& set, const Grid& grid, uint depth);
void adjoin_outer_approximation(GridTreeSet& grid_set, const TaylorSet& set, uint depth);
Zonotope zonotope(const TaylorSet& ts);

void draw(GraphicsInterface& g, const TaylorSet& ts);
void box_draw(GraphicsInterface& g, const TaylorSet& ts);
void affine_draw(GraphicsInterface& g, const TaylorSet& ts);
void curve_draw(GraphicsInterface& g, const TaylorSet& ts);
void grid_draw(GraphicsInterface& g, const TaylorSet& ts);

void plot(const char* fn, const Box& bbx, const TaylorSet& ts);



} //namespace Ariadne

#endif /* ARIADNE_TAYLOR_SET_H */
