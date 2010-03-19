/***************************************************************************
 *            hybrid_enclosure.h
 *
 *  Copyright  2009-10  Pieter Collins
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

/*! \file hybrid_enclosure.h
 *  \brief Enclosure sets for hybrid systems
 */

#ifndef ARIADNE_HYBRID_ENCLOSURE_H
#define ARIADNE_CONSTRAINED_IMAGE_SET_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>
#include "taylor_function.h"
#include "graphics_interface.h"
#include "container.h"
#include "logging.h"

namespace Ariadne {

class Interval;
template<class X> class Vector;
template<class X> class LinearProgram;
class ScalarFunction;
class VectorFunction;
class ScalarTaylorFunction;
class VectorTaylorFunction;
class TaylorConstrainedImageSet;
class HybridEnclosure;
class Box;
class Grid;
class GridTreeSet;
class AffineSet;
class DiscreteEvent;
class Figure;
class CanvasInterface;


template<class BS> class HybridBasicSet;
typedef HybridBasicSet<Box> HybridBox;
class HybridGridTreeSet;

class HybridEnclosure
    : public DrawableInterface
{
    friend class ConstrainedImageSetHybridEvolver;
    typedef Vector<Interval> IntervalVector;
  private:
    DiscreteLocation _location;
    List<DiscreteEvent> _events;
    List< List<DiscreteEvent> > _constraint_events;
    TaylorConstrainedImageSet _set;
    ScalarTaylorFunction _time;
  public:
    typedef TaylorConstrainedImageSet ContinuousStateSetType;

    HybridEnclosure(const DiscreteLocation&, const Box&);
    ~HybridEnclosure();
    HybridEnclosure* clone() const;

    const DiscreteLocation& location() const;
    const TaylorConstrainedImageSet& continuous_state_set() const;

    void apply_reset(DiscreteEvent, DiscreteLocation, VectorFunction);
    void apply_flow(VectorFunction, Interval);
    void apply_flow(VectorTaylorFunction);

    void set_maximum_time(DiscreteEvent,Float);
    void set_time(DiscreteEvent,Float);
    void set_dwell_time(DiscreteEvent,Float);
    void set_dwell_time(DiscreteEvent,ScalarFunction);
    void set_dwell_time(DiscreteEvent,ScalarTaylorFunction);

    void new_invariant(DiscreteEvent,ScalarFunction,ScalarFunction);
    void new_activation(DiscreteEvent,ScalarFunction,ScalarFunction);
    void new_time_bound(DiscreteEvent,ScalarFunction);

    uint dimension() const;
    tribool empty() const;

    HybridBox bounding_box() const;
    tribool disjoint(const HybridBox& bx) const;
    void adjoin_outer_approximation_to(HybridGridTreeSet& paving, int depth) const;
    Pair<HybridEnclosure,HybridEnclosure> split(uint dim) const;

    void draw(CanvasInterface&) const;
    std::ostream& write(std::ostream&) const;
};

inline std::ostream& operator<<(std::ostream& os, const HybridEnclosure& s) { return s.write(os); }


} // namespace Ariadne

#endif // ARIADNE_CONSTRAINED_IMAGE_SET_H
