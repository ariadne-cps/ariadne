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
    : _curves(new std::vector< std::pair<DiscreteState,InterpolatedCurve> >(1u,make_pair(pt.first,InterpolatedCurve(pt.second))))
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
    ARIADNE_ASSERT((uint)ht.discrete_time<=this->size());
    if(this->size()==ht.discrete_time) {
        this->_curves->push_back(make_pair(hpt.first,InterpolatedCurve(hpt.second)));
    } else {
        (*this->_curves)[ht.discrete_time].second.insert(ht.continuous_time,hpt.second);
    }
}


template<> 
std::ostream& 
operator<<(std::ostream& os, const Orbit< HybridPoint >& orb)
{
    return os << orb.curves();
}


struct Orbit<TaylorSet>::Data {
    Data(const TaylorSet& initial_set) 
        : initial(initial_set) { }
    TaylorSet initial;
    TaylorSetList reach;
    TaylorSetList intermediate;
    TaylorSetList final;
};

Orbit<TaylorSet>::
Orbit(const TaylorSet& initial_set)
    : _data(new Data(initial_set))
{
}

void
Orbit<TaylorSet>::
adjoin_reach(const TaylorSet& set)
{
    this->_data->reach.adjoin(set);
}

void
Orbit<TaylorSet>::
adjoin_intermediate(const TaylorSet& set)
{
    this->_data->intermediate.adjoin(set);
}

void
Orbit<TaylorSet>::
adjoin_final(const TaylorSet& set)
{
    this->_data->final.adjoin(set);
}


void
Orbit<TaylorSet>::
adjoin_reach(const TaylorSetList& list_set)
{
    this->_data->reach.adjoin(list_set);
}

void
Orbit<TaylorSet>::
adjoin_intermediate(const TaylorSetList& list_set)
{
    this->_data->intermediate.adjoin(list_set);
}

void
Orbit<TaylorSet>::
adjoin_final(const TaylorSetList& list_set)
{
    this->_data->final.adjoin(list_set);
}


TaylorSet const&
Orbit<TaylorSet>::
initial() const
{
    return this->_data->initial;
}

TaylorSetList const&
Orbit<TaylorSet>::
reach() const
{
    return this->_data->reach;
}

TaylorSetList const&
Orbit<TaylorSet>::
intermediate() const
{
    return this->_data->intermediate;
}

TaylorSetList const&
Orbit<TaylorSet>::
final() const
{
    return this->_data->final;
}




struct Orbit<HybridTaylorSet>::Data {
    Data(const HybridTaylorSet& initial_set) 
        : initial(initial_set) { }
    HybridTaylorSet initial;
    HybridTaylorSetList reach;
    HybridTaylorSetList intermediate;
    HybridTaylorSetList final;
};

Orbit<HybridTaylorSet>::
Orbit(const HybridTaylorSet& initial_set)
    : _data(new Data(initial_set))
{
}

void
Orbit<HybridTaylorSet>::
adjoin_reach(const HybridTaylorSet& set)
{
    this->_data->reach.adjoin(set);
}

void
Orbit<HybridTaylorSet>::
adjoin_intermediate(const HybridTaylorSet& set)
{
    this->_data->intermediate.adjoin(set);
}

void
Orbit<HybridTaylorSet>::
adjoin_final(const HybridTaylorSet& set)
{
    this->_data->final.adjoin(set);
}


void
Orbit<HybridTaylorSet>::
adjoin_reach(const HybridTaylorSetList& list_set)
{
    this->_data->reach.adjoin(list_set);
}

void
Orbit<HybridTaylorSet>::
adjoin_intermediate(const HybridTaylorSetList& list_set)
{
    this->_data->intermediate.adjoin(list_set);
}

void
Orbit<HybridTaylorSet>::
adjoin_final(const HybridTaylorSetList& list_set)
{
    this->_data->final.adjoin(list_set);
}


HybridTaylorSet const&
Orbit<HybridTaylorSet>::
initial() const
{
    return this->_data->initial;
}

HybridTaylorSetList const&
Orbit<HybridTaylorSet>::
reach() const
{
    return this->_data->reach;
}

HybridTaylorSetList const&
Orbit<HybridTaylorSet>::
intermediate() const
{
    return this->_data->intermediate;
}

HybridTaylorSetList const&
Orbit<HybridTaylorSet>::
final() const
{
    return this->_data->final;
}





struct Orbit<GridCell>::Data {
    Data(const GridCell& initial_set) 
        : initial(initial_set) { }
    GridCell initial;
    GridTreeSet reach;
    GridTreeSet intermediate;
    GridTreeSet final;
};

Orbit<GridCell>::
Orbit(const GridCell& initial_set)
    : _data(new Data(initial_set))
{
}

Orbit<GridCell>::
Orbit(const GridCell& initial_set,
      const GridTreeSet& reach_set,
      const GridTreeSet& intermediate_set,
      const GridTreeSet& final_set)
    : _data(new Data(initial_set))
{
    this->_data->reach=reach_set;
    this->_data->intermediate=intermediate_set;
    this->_data->final=final_set;
    std::cerr<<"Orbit<GridCell>::Orbit(...)\n  initial="<<initial_set<<"\n  reach="<<reach_set<<"\n  intermediate="<<intermediate_set<<"\n  final="<<final_set<<"\n\n"<<std::flush;
    std::cerr<<"orbit="<<*this<<"\n\n"<<std::flush;
}


GridCell const&
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
    Data(const HybridGridCell& initial_set) 
        : initial(initial_set) { }
    HybridGridCell initial;
    HybridGridTreeSet reach;
    HybridGridTreeSet intermediate;
    HybridGridTreeSet final;
};

Orbit<HybridGridCell>::
Orbit(const HybridGridCell& initial_set)
    : _data(new Data(initial_set))
{
}

Orbit<HybridGridCell>::
Orbit(const HybridGridCell& initial_set,
      const HybridGridTreeSet& reach_set,
      const HybridGridTreeSet& intermediate_set,
      const HybridGridTreeSet& final_set)
    : _data(new Data(initial_set))
{
    this->_data->reach=reach_set;
    this->_data->intermediate=intermediate_set;
    this->_data->final=final_set;
}

HybridGridCell const&
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
operator<<(std::ostream& os, const Orbit<TaylorSet>& orb)
{
    os << "Orbit(\n  initial=" << Orbit<TaylorSet>::EnclosureListType(orb.initial()).bounding_boxes()
       << "\n  intermediate=" << orb.intermediate().bounding_boxes()
       << "\n  reach=" << orb.reach().bounding_boxes()
       << "\n  final=" << orb.final().bounding_boxes()
       << ")\n";
    return os;
}


template<> 
std::ostream& 
operator<<(std::ostream& os, const Orbit<HybridTaylorSet>& orb)
{
    os << "Orbit(\n  initial=" << Orbit<HybridTaylorSet>::EnclosureListType(orb.initial()).bounding_boxes()
       << "\n  intermediate=" << orb.intermediate().bounding_boxes()
       << "\n  reach=" << orb.reach().bounding_boxes()
       << "\n  final=" << orb.final().bounding_boxes()
       << ")\n";
    return os;
}


} // namespace Ariadne
