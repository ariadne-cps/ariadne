/***************************************************************************
 *            hybrid_paving.cpp
 *
 *  Copyright 2008--17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "../function/functional.hpp"
#include "../config.hpp"

#include "../hybrid/hybrid_expression_set.hpp"
#include "../hybrid/hybrid_set.hpp"
#include "../hybrid/hybrid_paving.hpp"

#include "../numeric/real.hpp"

#include "../symbolic/expression_set.hpp"
#include "../geometry/function_set.hpp"

#include "../hybrid/hybrid_space.hpp"
#include "../hybrid/hybrid_time.hpp"
#include "../hybrid/hybrid_orbit.hpp"
#include "../hybrid/hybrid_automaton_interface.hpp"
#include "../output/graphics.hpp"
#include "../hybrid/hybrid_graphics.hpp"
#include "../numeric/rounding.hpp"
#include "../symbolic/assignment.hpp"
#include "../output/graphics_interface.hpp"
#include "../geometry/function_set.hpp"

namespace Ariadne {


HybridGridTreePaving::ConstIterator HybridGridTreePaving::begin() const {
    return ConstIterator(this->_map,this->_hgrid.space(),false);
}

HybridGridTreePaving::ConstIterator HybridGridTreePaving::end() const {
    return ConstIterator(this->_map,this->_hgrid.space(),true);
}


Void HybridGridTreePaving::adjoin(const HybridGridCell& hgc) {
    this->_provide_location(hgc.location()).adjoin(hgc.euclidean_set());
}

Void HybridGridTreePaving::adjoin(const ListSet<HybridGridCell>& hgcls) {
    for(ListSet<HybridGridCell>::ConstIterator iter=hgcls.begin(); iter!=hgcls.end(); ++iter) {
        this ->adjoin(*iter);
    }
}

Void HybridGridTreePaving::adjoin(const HybridGridTreePaving& hgts) {
    for(HybridGridTreePaving::LocationsConstIterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
        this->_provide_location(_loc_iter->first).adjoin(_loc_iter->second);
    }
}

Void HybridGridTreePaving::remove(const HybridGridTreePaving& hgts) {
    for(HybridGridTreePaving::LocationsConstIterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
        this->_provide_location(_loc_iter->first).remove(_loc_iter->second);
    }
}

Void HybridGridTreePaving::restrict(const HybridGridTreePaving& hgts) {
    for(HybridGridTreePaving::LocationsConstIterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
        this->_provide_location(_loc_iter->first).restrict(_loc_iter->second);
    }
}

Void HybridGridTreePaving::restrict_to_height(Nat h) {
    for(LocationsIterator _loc_iter=locations_begin(); _loc_iter!=locations_end(); ++_loc_iter) {
        _loc_iter->second.restrict_to_height(h);
    }
}

Void HybridGridTreePaving::adjoin_inner_approximation(const HybridSetInterface& hset, const Nat depth) {
    Set<DiscreteLocation> locations=hset.locations();
    for(auto location : locations) {
        RealSpace space = this->space(location);
        this->_provide_location(location).adjoin_inner_approximation(hset.euclidean_set(location,space),depth);
    }
}

Void HybridGridTreePaving::adjoin_inner_approximation(const HybridExactBoxes& hbxs, const Nat depth) {
    for(HybridExactBoxes::ConstIterator _loc_iter=hbxs.begin();
            _loc_iter!=hbxs.end(); ++_loc_iter) {
        DiscreteLocation const& loc=_loc_iter->first;
        LabelledSet<ExactBoxType> const& vbx=_loc_iter->second;
        ARIADNE_ASSERT(vbx.space() == this->space(loc));
        this->_provide_location(loc).adjoin_inner_approximation(vbx.euclidean_set(),depth);
    }
}

Void HybridGridTreePaving::adjoin_lower_approximation(const HybridOvertSetInterface& hs, const Nat height, const Nat depth) {
    Set<DiscreteLocation> hlocs=dynamic_cast<const HybridBoundedSetInterface&>(hs).locations();
    for(Set<DiscreteLocation>::ConstIterator _loc_iter=hlocs.begin();
            _loc_iter!=hlocs.end(); ++_loc_iter) {
        DiscreteLocation loc=*_loc_iter;
        RealSpace spc=this->space(loc);
        this->_provide_location(loc).adjoin_lower_approximation(hs.euclidean_set(loc,spc),height,depth);
    }
}

Void HybridGridTreePaving::adjoin_outer_approximation(const HybridCompactSetInterface& hs, const Nat depth) {
    Set<DiscreteLocation> hlocs=hs.locations();
    for(Set<DiscreteLocation>::ConstIterator _loc_iter=hlocs.begin();
            _loc_iter!=hlocs.end(); ++_loc_iter) {
        DiscreteLocation loc=*_loc_iter;
        RealSpace spc=this->space(loc);
        this->_provide_location(loc).adjoin_outer_approximation(hs.euclidean_set(loc,spc),depth);
    }
}

Void HybridGridTreePaving::adjoin_outer_approximation(const HybridExactBoxes& hbxs, const Nat depth) {
    for(HybridExactBoxes::ConstIterator _loc_iter=hbxs.begin();
            _loc_iter!=hbxs.end(); ++_loc_iter) {
        DiscreteLocation const& loc=_loc_iter->first;
        LabelledSet<ExactBoxType> const& vbx=_loc_iter->second;
        ARIADNE_ASSERT(vbx.space() == this->space(loc));
        this->_provide_location(_loc_iter->first).adjoin_outer_approximation(vbx.euclidean_set(),depth);
    }
}


GridTreePaving& HybridGridTreePaving::operator[](DiscreteLocation q) {
    return this->_provide_location(q);
}

const GridTreePaving& HybridGridTreePaving::operator[](DiscreteLocation q) const {
    ARIADNE_ASSERT_MSG(this->has_location(q),"q="<<q);
    return const_cast<HybridGridTreePaving*>(this)->_provide_location(q);
}

Bool HybridGridTreePaving::is_empty() const {
    for(LocationsConstIterator _loc_iter=this->locations_begin();
        _loc_iter!=this->locations_end(); ++_loc_iter) {
        if(!_loc_iter->second.is_empty()) { return false; }
    }
    return true;
}

SizeType HybridGridTreePaving::size() const {
    SizeType result=0;
    for(LocationsConstIterator _loc_iter=this->locations_begin(); _loc_iter!=this->locations_end(); ++_loc_iter) {
        result+=_loc_iter->second.size();
    }
    return result;
}

HybridListSet<ExactBoxType> HybridGridTreePaving::boxes() const {
    HybridListSet<ExactBoxType> result;
    for(ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        result.adjoin(iter->location(),iter->euclidean_set().box());
    }
    return result;
}

Void HybridGridTreePaving::mince(Nat depth) {
    for(LocationsIterator _loc_iter=this->locations_begin();
        _loc_iter!=this->locations_end(); ++_loc_iter) {
        _loc_iter->second.mince(depth);
    }
}

Void HybridGridTreePaving::recombine() {
    for(LocationsIterator _loc_iter=this->locations_begin();
        _loc_iter!=this->locations_end(); ++_loc_iter) {
        _loc_iter->second.recombine();
    }
}

ValidatedLowerKleenean HybridGridTreePaving::separated(const HybridExactBox& hbx) const {
    LocationsConstIterator _loc_iter = this->_map.find( hbx.location() );
    return _loc_iter != this->locations_end() || _loc_iter->second.separated( hbx.euclidean_set() );
}

ValidatedLowerKleenean HybridGridTreePaving::overlaps(const HybridExactBox& hbx) const {
    LocationsConstIterator _loc_iter = this->_map.find( hbx.location() );
    return _loc_iter != this->locations_end() && _loc_iter->second.overlaps( hbx.euclidean_set() );
}

ValidatedLowerKleenean HybridGridTreePaving::covers(const HybridExactBox& hbx) const {
    LocationsConstIterator _loc_iter=this->_map.find(hbx.location());
    return _loc_iter!=this->locations_end() && _loc_iter->second.covers( hbx.euclidean_set() );
}

ValidatedLowerKleenean HybridGridTreePaving::inside(const HybridExactBoxes& hbx) const  {
    for( LocationsConstIterator _loc_iter = this->locations_begin(); _loc_iter != this->locations_end(); ++_loc_iter ) {
        if( !_loc_iter->second.is_empty() ) {
            DiscreteLocation const& loc = _loc_iter->first;
            RealSpace spc=this->space(loc);
            return this->euclidean_set(loc).inside(hbx.euclidean_set(loc,spc));
        }
    }
    return true;
}

HybridUpperBoxes HybridGridTreePaving::bounding_box() const {
    HybridExactBoxes result;
    for( LocationsConstIterator _loc_iter = this->locations_begin(); _loc_iter != this->locations_end(); ++_loc_iter ) {
        if( !_loc_iter->second.is_empty() ) {
            DiscreteLocation const& loc = _loc_iter->first;
            RealSpace const& spc=this->space(loc);
            result.insert(loc,spc,cast_exact_box(_loc_iter->second.bounding_box()));
        }
    }
    return result;
}

Bool subset(const HybridGridTreePaving& hgts1, const HybridGridTreePaving& hgts2) {
    for( typename HybridGridTreePaving::LocationsConstIterator _loc_iter = hgts1.locations_begin(); _loc_iter != hgts1.locations_end(); ++_loc_iter ) {
        if( !_loc_iter->second.is_empty() ) {
            DiscreteLocation const& loc = _loc_iter->first;
            RealSpace spc=hgts1.space(loc);
            if(not subset(hgts1.euclidean_set(loc),hgts2.euclidean_set(loc,spc))) {
                return false;
            }
        }
    }
    return true;
}

OutputStream& HybridGridTreePaving::write(OutputStream& os) const {
    return os << this->_map;
}

Void HybridGridTreePaving::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axis_variables) const {
    for(LocationsConstIterator loc_iter=this->locations_begin(); loc_iter!=this->locations_end(); ++loc_iter) {
        if(locations.empty() || locations.contains(loc_iter->first)) {
            RealSpace const& space=this->space(loc_iter->first);
            Projection2d projection(space.dimension(),space.index(axis_variables.x_variable()),space.index(axis_variables.y_variable()));
            loc_iter->second.draw(canvas,projection);
        }
    }
}


GridTreePaving& HybridGridTreePaving::_provide_location(const DiscreteLocation& q) {
    std::map<DiscreteLocation,GridTreePaving>::iterator iter=this->_map.find(q);
    if(iter==this->_map.end()) {
        this->_map.insert(std::make_pair(q,GridTreePaving(this->_hgrid[q])));
        iter=this->_map.find(q);
    }
    return iter->second; }

} // namespace Ariadne
