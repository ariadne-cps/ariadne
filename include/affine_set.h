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
#include "logging.h"

namespace Ariadne {

class Interval;
template<class X> class Vector;
template<class X> class LinearProgram;
class ScalarFunction;
class VectorFunction;
class ScalarTaylorFunction;
class VectorTaylorFunction;
class Box;
class Grid;
class GridTreeSet;
class DiscreteEvent;
class Figure;
class CanvasInterface;

class AffineSet
    : public DrawableInterface, public Loggable
{
    Vector<Interval> _domain;
    List< Affine<Float> > _function;
    List< Affine<Float> > _constraints;
    List< Affine<Float> > _equations;
  public:
    //!\brief The set \f$\{ Gy+c \mid y\in D\}\f$.
    AffineSet(const Vector<Interval>& D, const Matrix<Float>& G, const Vector<Float>& c);
    //!\brief The set \f$\{ Gy+c \mid ||y||_\infty\leq 1\}\f$. \deprecated
    AffineSet(const Matrix<Float>& G, const Vector<Float>& c);
    AffineSet* clone() const;
    //!\brief The constraint \f$a\cdot y\leq b\f$.
    void new_inequality_constraint(const Vector<Float>& a, const Float& b);
    //!\brief The constraint \f$a\cdot y = b\f$.
    void new_equality_constraint(const Vector<Float>& a, const Float& b);

    uint dimension() const;
    uint number_of_parameters() const;
    uint number_of_constraints() const;
    Vector<Interval> domain() const;

    tribool bounded() const;
    Box bounding_box() const;
    tribool disjoint(const Box& bx) const;
    tribool empty() const;

    void adjoin_outer_approximation_to(GridTreeSet& g, int depth) const;
    GridTreeSet outer_approximation(const Grid& g, int depth) const;
    void robust_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const;

    List<Point2d> boundary(uint xc, uint yc) const;

    virtual void draw(CanvasInterface&) const;
    virtual std::ostream& write(std::ostream& os) const;

  private:
    void construct(const Vector<Interval>& D, const Matrix<Float>& G, const Vector<Float>& c);
    void construct_linear_program(LinearProgram<Float>& lp) const;
};

inline std::ostream& operator<<(std::ostream& os, const AffineSet& as) {
    return as.write(os); }



} // namespace Ariadne

#endif // ARIADNE_AFFINE_SET_H
