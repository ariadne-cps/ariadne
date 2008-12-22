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

#include "taylor_variable.h"

namespace Ariadne {

class Interval;
template<class X> class Vector;

class TaylorModel;
class ApproximateTaylorModel;

class Zonotope;
class Grid;
class GridTreeSet;
class Figure;

class TaylorSet 
    : public LocatedSetInterface
{
  private:
    Vector<TaylorVariable> _variables;
  public:
    TaylorSet(uint d=0);
    template<class XE, class XP> TaylorSet(uint rs, uint as, uint d, const XE* eps, const XP* ptr);
    TaylorSet(const TaylorModel& tm);
    TaylorSet(const ApproximateTaylorModel& atm);
    
    uint dimension() const { return this->_variables.size(); }
    uint number_of_generators() const { assert(this->_variables.size()>0); return this->_variables[0].expansion().argument_size(); }
    const TaylorVariable& operator[](uint i) const { return this->_variables[i]; }
    TaylorVariable& operator[](uint i) { return this->_variables[i]; }

    Vector<Interval> domain() const;

    TaylorSet* clone() const { return new TaylorSet(*this); } 
    tribool disjoint(const Box&) const;
    tribool overlaps(const Box&) const;
    tribool inside(const Box&) const;
    Box bounding_box() const;
    std::ostream& write(std::ostream& os) const;
    
    pair<TaylorSet,TaylorSet> split(uint dim) const;
};

GridTreeSet outer_approximation(const TaylorSet& set, uint depth);
GridTreeSet outer_approximation(const TaylorSet& set, const Grid& grid, uint depth);
Zonotope zonotope(const TaylorSet& ts);
void draw(Figure& g, const TaylorSet& ts);



template<class XE, class XP> 
TaylorSet::TaylorSet(uint rs, uint as, uint d, 
                     const XE* eps, const XP* ptr)
    : _variables(rs)
{
    for(uint i=0; i!=rs; ++i) {
        _variables[i]=TaylorVariable(as,d,eps,ptr);
    }
}



} //namespace Ariadne

#endif /* ARIADNE_APPROXIMATE_TAYLOR_MODEL_H */
