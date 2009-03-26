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

class FunctionInterface;
class TaylorModel;
class TaylorFunction;
class TaylorSet;

class Zonotope;
class Grid;
class GridCell;
class GridTreeSet;
class Figure;

/*! \brief Sets expressed as the image of a box under a polynomial with error bounds. */
class TaylorSet
    : public LocatedSetInterface
{
  private:
    Vector<TaylorModel> _models;
  public:
    TaylorSet(uint d=0, uint ng=0);
    template<class XE, class XP> TaylorSet(uint rs, uint as, uint d, const XE* eps, const XP* ptr);
    TaylorSet(uint rs, uint as, uint deg, double x0, ...);
    TaylorSet(const FunctionInterface& f, const Vector<Interval>& d);
    TaylorSet(const Vector<TaylorModel>& tv);
    TaylorSet(const Vector< Expansion<Float> >& f, const Vector<Float>& e);
    TaylorSet(const Vector<Interval>& bx);

    friend bool operator==(const TaylorSet& ts1, const TaylorSet& ts2);

    uint dimension() const { return this->_models.size(); }
    uint generators_size() const { assert(this->_models.size()>0); return this->_models[0].argument_size(); }
    uint argument_size() const { return this->_models[0].argument_size(); }
    const TaylorModel& operator[](uint i) const { return this->_models[i]; }
    TaylorModel& operator[](uint i) { return this->_models[i]; }

    Vector<Interval> domain() const { return Vector<Interval>(this->generators_size(),Interval(-1,+1)); }
    Vector<Interval> range() const {
        Vector<Interval> result(this->_models.size());
        for(uint i=0; i!=this->_models.size(); ++i) {
            result[i]=this->_models[i].range(); }
        return result; }
    Vector<TaylorModel> models() const { return this->_models; }

    TaylorSet* clone() const { return new TaylorSet(*this); }
    Float radius() const;
    tribool disjoint(const Box&) const;
    tribool overlaps(const Box&) const;
    tribool inside(const Box&) const;
    Box bounding_box() const;
    std::ostream& write(std::ostream& os) const;

    TaylorSet linearise() const;
    GridTreeSet discretise(const Grid& grid, uint depth) const;
    GridTreeSet& discretise(GridTreeSet& grid_set, uint depth) const;
    TaylorSet recondition() const;
    TaylorSet subsume() const;
    pair<TaylorSet,TaylorSet> split(uint dim) const;
    pair<TaylorSet,TaylorSet> split() const;
    ListSet<TaylorSet> subdivide(Float rad) const;
  private:
    Matrix<Float> jacobian() const;
};

TaylorSet apply(const TaylorFunction& f, const TaylorSet& s);
TaylorSet apply(const FunctionInterface& f, const TaylorSet& s);

GridTreeSet outer_approximation(const TaylorSet& set, const Grid& grid, uint depth);
void adjoin_outer_approximation(GridTreeSet& grid_set, const TaylorSet& set, uint depth);
Zonotope zonotope(const TaylorSet& ts);
void draw(Figure& g, const TaylorSet& ts);
void box_draw(Figure& g, const TaylorSet& ts);
void affine_draw(Figure& g, const TaylorSet& ts);
void curve_draw(Figure& g, const TaylorSet& ts);
void grid_draw(Figure& g, const TaylorSet& ts);

void plot(const char* fn, const Box& bbx, const TaylorSet& ts);

template<class XE, class XP>
TaylorSet::TaylorSet(uint rs, uint as, uint d,
                     const XE* eps, const XP* ptr)
    : _models(rs)
{
    for(uint i=0; i!=rs; ++i) {
        _models[i]=TaylorModel(as,d,eps,ptr);
    }
}


} //namespace Ariadne

#endif /* ARIADNE_TAYLOR_SET_H */
