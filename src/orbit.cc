/***************************************************************************
 *            orbit.cc
 *
 *  Copyright 2007-8  Pieter Collins
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

#include <utility>

#include "orbit.h"

#include "box.h"
#include "point.h"
#include "curve.h"
#include "taylor_set.h"
#include "list_set.h"
#include "hybrid_set.h"

#include "hybrid_time.h"

namespace Ariadne {


Orbit<Point>::Orbit(const Point& pt)
    : _curve(new InterpolatedCurve(0.0,pt))
{ }

void
Orbit<Point>::insert(Time t, const Point& pt)
{
    this->_curve->insert(t,pt);
}


Orbit<HybridPoint>::Orbit(const HybridPoint& pt)
    : _curves(new std::vector<HybridInterpolatedCurve>(1u,make_pair(pt.first,InterpolatedCurve(pt.second))))
{ }

uint
Orbit<HybridPoint>::size() const
{
    return this->_curves->size();
}

const InterpolatedCurve&
Orbit<HybridPoint>::curve(uint m) const
{
    return (*this->_curves)[m].second;
}

void
Orbit<HybridPoint>::insert(HybridTime ht, HybridPoint& hpt)
{
    ARIADNE_ASSERT((uint)ht.discrete_time()<=this->size());
    if(this->size()==ht.discrete_time()) {
        this->_curves->push_back(make_pair(hpt.location(),InterpolatedCurve(hpt.continuous_state_set())));
    } else {
        (*this->_curves)[ht.discrete_time()].second.insert(ht.continuous_time(),hpt.continuous_state_set());
    }
}


template<>
std::ostream&
operator<<(std::ostream& os, const Orbit< HybridPoint >& orb)
{
    return os << orb.curves();
}


struct Orbit<TaylorImageSet>::Data {
    Data(const TaylorImageSet& initial_set)
        : initial(initial_set) { }
    TaylorImageSet initial;
    TaylorImageSetList reach;
    TaylorImageSetList intermediate;
    TaylorImageSetList final;
};

Orbit<TaylorImageSet>::
Orbit(const TaylorImageSet& initial_set)
    : _data(new Data(initial_set))
{
}

void
Orbit<TaylorImageSet>::
adjoin_reach(const TaylorImageSet& set)
{
    this->_data->reach.adjoin(set);
}

void
Orbit<TaylorImageSet>::
adjoin_intermediate(const TaylorImageSet& set)
{
    this->_data->intermediate.adjoin(set);
}

void
Orbit<TaylorImageSet>::
adjoin_final(const TaylorImageSet& set)
{
    this->_data->final.adjoin(set);
}


void
Orbit<TaylorImageSet>::
adjoin_reach(const TaylorImageSetList& list_set)
{
    this->_data->reach.adjoin(list_set);
}

void
Orbit<TaylorImageSet>::
adjoin_intermediate(const TaylorImageSetList& list_set)
{
    this->_data->intermediate.adjoin(list_set);
}

void
Orbit<TaylorImageSet>::
adjoin_final(const TaylorImageSetList& list_set)
{
    this->_data->final.adjoin(list_set);
}


TaylorImageSet const&
Orbit<TaylorImageSet>::
initial() const
{
    return this->_data->initial;
}

TaylorImageSetList const&
Orbit<TaylorImageSet>::
reach() const
{
    return this->_data->reach;
}

TaylorImageSetList const&
Orbit<TaylorImageSet>::
intermediate() const
{
    return this->_data->intermediate;
}

TaylorImageSetList const&
Orbit<TaylorImageSet>::
final() const
{
    return this->_data->final;
}




struct Orbit<HybridTaylorImageSet>::Data {
    Data(const HybridTaylorImageSet& initial_set)
        : initial(initial_set) { }
    HybridTaylorImageSet initial;
    HybridTaylorImageSetList reach;
    HybridTaylorImageSetList intermediate;
    HybridTaylorImageSetList final;
};

Orbit<HybridTaylorImageSet>::
Orbit(const HybridTaylorImageSet& initial_set)
    : _data(new Data(initial_set))
{
}

void
Orbit<HybridTaylorImageSet>::
adjoin_reach(const HybridTaylorImageSet& set)
{
    this->_data->reach.adjoin(set);
}

void
Orbit<HybridTaylorImageSet>::
adjoin_intermediate(const HybridTaylorImageSet& set)
{
    this->_data->intermediate.adjoin(set);
}

void
Orbit<HybridTaylorImageSet>::
adjoin_final(const HybridTaylorImageSet& set)
{
    this->_data->final.adjoin(set);
}


void
Orbit<HybridTaylorImageSet>::
adjoin_reach(const HybridTaylorImageSetList& list_set)
{
    this->_data->reach.adjoin(list_set);
}

void
Orbit<HybridTaylorImageSet>::
adjoin_intermediate(const HybridTaylorImageSetList& list_set)
{
    this->_data->intermediate.adjoin(list_set);
}

void
Orbit<HybridTaylorImageSet>::
adjoin_final(const HybridTaylorImageSetList& list_set)
{
    this->_data->final.adjoin(list_set);
}


HybridTaylorImageSet const&
Orbit<HybridTaylorImageSet>::
initial() const
{
    return this->_data->initial;
}

HybridTaylorImageSetList const&
Orbit<HybridTaylorImageSet>::
reach() const
{
    return this->_data->reach;
}

HybridTaylorImageSetList const&
Orbit<HybridTaylorImageSet>::
intermediate() const
{
    return this->_data->intermediate;
}

HybridTaylorImageSetList const&
Orbit<HybridTaylorImageSet>::
final() const
{
    return this->_data->final;
}





struct Orbit<GridCell>::Data {
    Data(const Grid& grid) : initial(grid), reach(grid), intermediate(grid), final(grid) { }
    GridTreeSet initial;
    GridTreeSet reach;
    GridTreeSet intermediate;
    GridTreeSet final;
};

Orbit<GridCell>::
Orbit(const GridTreeSet& initial_set)
    : _data(new Data(initial_set.grid()))
{
    this->_data->initial=initial_set;
}


Orbit<GridCell>::
Orbit(const GridTreeSet& initial_set,
      const GridTreeSet& reach_set,
      const GridTreeSet& intermediate_set,
      const GridTreeSet& final_set)
    : _data(new Data(initial_set.grid()))
{
    this->_data->initial=initial_set;
    this->_data->reach=reach_set;
    this->_data->intermediate=intermediate_set;
    this->_data->final=final_set;
}


GridTreeSet const&
Orbit<GridCell>::
initial() const
{
    return this->_data->initial;
}

GridTreeSet const&
Orbit<GridCell>::
reach() const
{
    return this->_data->reach;
}

GridTreeSet const&
Orbit<GridCell>::
intermediate() const
{
    return this->_data->intermediate;
}

GridTreeSet const&
Orbit<GridCell>::
final() const
{
    return this->_data->final;
}



struct Orbit<HybridGridCell>::Data {
    Data(const HybridGrid& grid)
        : initial(grid), reach(grid), intermediate(grid), final(grid) { }
    HybridGridTreeSet initial;
    HybridGridTreeSet reach;
    HybridGridTreeSet intermediate;
    HybridGridTreeSet final;
};

Orbit<HybridGridCell>::
Orbit(const HybridGridTreeSet& initial_set)
    : _data(new Data(initial_set.grid()))
{
    this->_data->initial=initial_set;
}

Orbit<HybridGridCell>::
Orbit(const HybridGridTreeSet& initial_set,
      const HybridGridTreeSet& reach_set,
      const HybridGridTreeSet& intermediate_set,
      const HybridGridTreeSet& final_set)
    : _data(new Data(initial_set.grid()))
{
    this->_data->initial=initial_set;
    this->_data->reach=reach_set;
    this->_data->intermediate=intermediate_set;
    this->_data->final=final_set;
}

HybridGridTreeSet const&
Orbit<HybridGridCell>::
initial() const
{
    return this->_data->initial;
}

HybridGridTreeSet const&
Orbit<HybridGridCell>::
reach() const
{
    return this->_data->reach;
}

HybridGridTreeSet const&
Orbit<HybridGridCell>::
intermediate() const
{
    return this->_data->intermediate;
}

HybridGridTreeSet const&
Orbit<HybridGridCell>::
final() const
{
    return this->_data->final;
}


template<>
std::ostream&
operator<<(std::ostream& os, const Orbit<TaylorImageSet>& orb)
{
    os << "Orbit(\n  initial=" << Orbit<TaylorImageSet>::EnclosureListType(orb.initial()).bounding_boxes()
       << "\n  intermediate=" << orb.intermediate().bounding_boxes()
       << "\n  reach=" << orb.reach().bounding_boxes()
       << "\n  final=" << orb.final().bounding_boxes()
       << ")\n";
    return os;
}


template<>
std::ostream&
operator<<(std::ostream& os, const Orbit<HybridTaylorImageSet>& orb)
{
    os << "Orbit(\n  initial=" << Orbit<HybridTaylorImageSet>::EnclosureListType(orb.initial()).bounding_boxes()
       << "\n  intermediate=" << orb.intermediate().bounding_boxes()
       << "\n  reach=" << orb.reach().bounding_boxes()
       << "\n  final=" << orb.final().bounding_boxes()
       << ")\n";
    return os;
}

void draw(CanvasInterface& graphic, const Orbit<TaylorImageSet>& orbit)
{
    orbit.reach().draw(graphic);
    orbit.initial().draw(graphic);
    orbit.final().draw(graphic);
}



void draw(CanvasInterface& graphic, const Orbit<HybridPoint>& orbit)
{
    for(uint i=0; i<=orbit.size(); ++i) {
        orbit.curve(i).draw(graphic);
    }
}


} // namespace Ariadne
