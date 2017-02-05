/***************************************************************************
 *            hybrid_set.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "function/functional.hpp"
#include "config.h"

#include "hybrid/hybrid_expression_set.hpp"
#include "hybrid/hybrid_set.hpp"

#include "numeric/real.hpp"

#include "expression/expression_set.hpp"
#include "geometry/function_set.hpp"

#include "hybrid/hybrid_space.hpp"
#include "hybrid/hybrid_time.hpp"
#include "hybrid/hybrid_orbit.hpp"
#include "hybrid/hybrid_automaton_interface.hpp"
#include "output/graphics.hpp"
#include "hybrid/hybrid_graphics.hpp"
#include <boost/concept_check.hpp>
#include "numeric/rounding.hpp"
#include "expression/assignment.hpp"
#include "output/graphics_interface.hpp"
#include "geometry/function_set.hpp"

namespace Ariadne {

template<> inline Float64Value numeric_cast<Float64Value>(Real const& r) {
    return cast_exact(r.get(Precision64()));
}

ExactBoxType over_approximation(RealBox const&);
ExactBoxType under_approximation(RealBox const&);
ExactBoxType approximation(RealBox const&);

HybridExactBox under_approximation(HybridRealBox const& hbx) {
    return HybridExactBox(hbx.location(),hbx.space(),under_approximation(hbx.euclidean_set()));
}

HybridExactBox over_approximation(HybridRealBox const& hbx) {
    return HybridExactBox(hbx.location(),hbx.space(),over_approximation(hbx.euclidean_set()));
}


Orbit<HybridExactPoint>::Orbit(const HybridExactPoint& hpt)
    : _curves_ptr(new std::vector<HybridInterpolatedCurve>(1u,HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(0,hpt.point()))))
{ }

Nat
Orbit<HybridExactPoint>::size() const
{
    return this->_curves_ptr->size();
}

const InterpolatedCurve&
Orbit<HybridExactPoint>::curve(Nat m) const
{
    return (*this->_curves_ptr)[m].euclidean_set();
}

Void
Orbit<HybridExactPoint>::insert(HybridTime ht, const HybridExactPoint& hpt)
{
    ARIADNE_ASSERT(ht.discrete_time()<=this->size());
    // FIXME: Should allow non-exact times
    Real time=ht.continuous_time();
    Float64Value flt_time=cast_exact(time.get(Precision64()));
    ARIADNE_ASSERT(decide(Real(flt_time)==time));
    if(this->size()==ht.discrete_time()) {
        this->_curves_ptr->push_back(HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(flt_time,hpt.point())));
    } else {
        (*this->_curves_ptr)[ht.discrete_time().get_si()].euclidean_set().insert(flt_time,hpt.point());
    }
}

Void Orbit<HybridExactPoint>::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axes) const {
    const Orbit<HybridExactPoint>& orbit=*this;
    for(Nat i=0; i!=orbit._curves_ptr->size(); ++i) {
        HybridInterpolatedCurve const& hcurve=this->_curves_ptr->at(i);
        if(locations.empty() || locations.contains(hcurve.location())) {
            RealSpace const& space=hcurve.space();
            if(valid_axis_variables(space,axes)) {
                hcurve.euclidean_set().draw(canvas,projection(space,axes));
            }
        }
    }
}

template<>
OutputStream&
operator<<(OutputStream& os, const Orbit< HybridExactPoint >& orb)
{
    return os << orb.curves();
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



Void Orbit<HybridEnclosure>::draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d& v) const {
    this->reach().draw(c,l,v);
}

template<>
OutputStream&
operator<<(OutputStream& os, const Orbit< HybridEnclosure >& orb)
{
    os << "Orbit(\n  initial=" << orb.initial()
       << "\n  intermediate=" << orb.intermediate()
       << "\n  reach=" << orb.reach()
       << "\n  final=" << orb.final()
       << ")\n";
    return os;
}






Map<RealVariable,RealInterval> make_map(const List<RealVariableInterval>& b) {
    Map<RealVariable,RealInterval> res;
    for(Nat i=0; i!=b.size(); ++i) {
        res.insert(b[i].variable(),RealInterval(b[i].lower(),b[i].upper()));
    }
    return res;
}

template<class X> HybridPoint<X>::HybridPoint(const DiscreteLocation& q, const Map<RealVariable,X>& x)
    : HybridBasicSet<Point<X>>(q,make_list(x.keys()),Point<X>(x.size()))
{
    SizeType i=0;
    for(auto iter=x.begin(); iter!=x.end(); ++iter, ++i) {
        this->point()[i]=iter->second;
    }
}

template<class X> HybridPoint<X>::HybridPoint(const DiscreteLocation& q, const Map<RealVariable,Real>& x)
    : HybridBasicSet<Point<X>>(q,make_list(x.keys()),Point<X>(x.size()))
{
    SizeType i=0;
    for(auto iter=x.begin(); iter!=x.end(); ++iter, ++i) {
        this->point()[i]=numeric_cast<X>(iter->second);
    }
}

template<class X> HybridPoint<X>::HybridPoint(const DiscreteLocation& q, const List<Assignment<RealVariable,Real>>& x)
    : HybridBasicSet<Point<X>>(q,left_hand_sides(x),Point<X>(x.size()))
{
    SizeType i=0;
    for(auto iter=x.begin(); iter!=x.end(); ++iter, ++i) {
        this->point()[i]=numeric_cast<X>(iter->right_hand_side());
    }
}

template<class X> HybridPoint<X>::HybridPoint(const DiscreteLocation& q, const InitializerList<Assignment<RealVariable,Real>>& x)
    : HybridPoint<X>(q,List<Assignment<RealVariable,Real>>(x))
{
}

template<class X> Map<RealVariable,X> HybridPoint<X>::values() const {
    Map<RealVariable,X> r;
    for(SizeType i=0; i!=this->space().dimension(); ++i) {
        r.insert(this->space()[i],this->point()[i]);
    }
    return r;
}



Void RealHybridVariablesBox::draw(CanvasInterface& c, const Set<DiscreteLocation>& qs, const Variables2d& vs) const {
    if(qs.empty() || qs.contains(this->location())) {
        RealSpace spc(List<RealVariable>(this->variables()));
        Projection2d prj(spc.dimension(),spc.index(vs.x_variable()),spc.index(vs.y_variable()));
        this->euclidean_set(spc).draw(c,prj);
    }
}

HybridConstraintSet::HybridConstraintSet()
    : _sets()
{
}

HybridConstraintSet::HybridConstraintSet(const DiscreteLocation& loc,
                                                 const List<ContinuousPredicate>& cnstr)
    : _sets()
{
    _sets.insert(loc,RealExpressionConstraintSet(cnstr));
}

HybridConstraintSet* HybridConstraintSet::clone() const {
    return new HybridConstraintSet(*this);
}

Set<RealVariable> HybridConstraintSet::variables(DiscreteLocation loc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return _sets[loc].variables();
}

ConstraintSet const HybridConstraintSet::euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return ConstraintSet(this->_sets[loc].euclidean_set(spc));
}

RegularSetInterface* HybridConstraintSet::_euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    return new ConstraintSet(this->euclidean_set(loc,spc));
}

ValidatedSierpinskian HybridConstraintSet::overlaps(const HybridExactBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).overlaps(bx.euclidean_set());
    } else {
        return false;
    }
}

ValidatedSierpinskian HybridConstraintSet::covers(const HybridExactBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).covers(bx.euclidean_set());
    } else {
        return bx.euclidean_set().is_empty();
    }
}

ValidatedSierpinskian HybridConstraintSet::separated(const HybridExactBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).separated(bx.euclidean_set());
    } else {
        return true;
    }
}

OutputStream& HybridConstraintSet::write(OutputStream& os) const {
    return os << "HybridConstraintSet( "<< this->_sets << " )";
}



HybridBoundedConstraintSet::HybridBoundedConstraintSet()
    : _sets()
{
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                               const InitializerList<RealVariableInterval>& bnd)
    : _sets()
{
    _sets.insert(loc,RealExpressionBoundedConstraintSet(bnd));
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                               const InitializerList<RealVariableInterval>& bnd,
                                                               const InitializerList<ContinuousPredicate>& cnstr)
    : _sets()
{
    _sets.insert(loc,RealExpressionBoundedConstraintSet(bnd,cnstr));
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                               const RealVariablesBox& bx)
    : _sets()
{
    _sets.insert(loc,RealExpressionBoundedConstraintSet(bx));
}

HybridBoundedConstraintSet* HybridBoundedConstraintSet::clone() const {
    return new HybridBoundedConstraintSet(*this);
}

Set<RealVariable> HybridBoundedConstraintSet::variables(DiscreteLocation loc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return _sets[loc].variables();
}

BoundedConstraintSet const HybridBoundedConstraintSet::euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    // FIXME: Should be no need to cache Euclidean sets.
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return BoundedConstraintSet(this->_sets[loc].euclidean_set(spc));
}

BoundedConstraintSet* HybridBoundedConstraintSet::_euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    return new BoundedConstraintSet(this->euclidean_set(loc,spc));
}

ValidatedSierpinskian HybridBoundedConstraintSet::overlaps(const HybridExactBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).overlaps(bx.euclidean_set());
    } else {
        return false;
    }
}

ValidatedSierpinskian HybridBoundedConstraintSet::covers(const HybridExactBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).covers(bx.euclidean_set());
    } else {
        return bx.euclidean_set().is_empty();
    }
}

ValidatedSierpinskian HybridBoundedConstraintSet::separated(const HybridExactBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).separated(bx.euclidean_set());
    } else {
        return true;
    }
}

ValidatedSierpinskian HybridBoundedConstraintSet::inside(const HybridExactBoxes& bxs) const {
    ValidatedSierpinskian result=true;
    for(Map<DiscreteLocation,RealExpressionBoundedConstraintSet>::ConstIterator iter=this->_sets.begin(); iter!=this->_sets.end(); ++iter) {
        DiscreteLocation const& loc=iter->first;
        RealExpressionBoundedConstraintSet const& set = iter->second;
        Set<RealVariable> vars=set.variables();
        RealSpace const& spc = bxs[loc].space();
        ExactBoxType const& bx = bxs[loc].euclidean_set();
        result = result && set.euclidean_set(spc).inside(bx);
    }
    return result;
}

DiscreteLocation HybridBoundedConstraintSet::location() const {
    ARIADNE_ASSERT(this->_sets.size()==1);
    return this->_sets.begin()->first;
}

Set<DiscreteLocation> HybridBoundedConstraintSet::locations() const {
    return this->_sets.keys();
}

HybridExactBoxes HybridBoundedConstraintSet::bounding_box() const {
    HybridExactBoxes result;
    for(Map<DiscreteLocation,RealExpressionBoundedConstraintSet>::ConstIterator iter=this->_sets.begin(); iter!=this->_sets.end(); ++iter) {
        RealSpace spc=make_list(iter->second.variables());
        RealVariablesBox bnds=iter->second.bounds();
        RealBox rbx=bnds.euclidean_set(spc);
        ExactBoxType ebx=over_approximation(rbx);
        result.insert(iter->first,spc,ebx);
    }
    return result;
}

OutputStream& HybridBoundedConstraintSet::write(OutputStream& os) const {
    return os << "HybridBoundedConstraintSet( "<< this->_sets << " )";
}

Void HybridBoundedConstraintSet::draw(CanvasInterface& c, const Set<DiscreteLocation>& q, const Variables2d& p) const {
    if(q.empty() || q.contains(this->location())) {
        Set<RealVariable> variables=this->variables(this->location());
        RealSpace space(List<RealVariable>(variables.begin(),variables.end()));
        this->euclidean_set(this->location(),space).draw(c,projection(space,p));
    }
}

template<class EBS> Void HybridBasicSet<EBS>::adjoin_outer_approximation_to(HybridGridTreeSet& paving, Int depth) const {
    if(this->space()==paving.space(this->location())) {
        paving[this->location()].adjoin_outer_approximation(this->euclidean_set(),depth);
    } else {
        ARIADNE_FAIL_MSG("HybridSet's state variables "<<this->space()<<
                         " do not match variables "<<paving.space()<<" of paving in location "<<this->location());
    }
}

template Void HybridBasicSet<Enclosure>::adjoin_outer_approximation_to(HybridGridTreeSet& paving, Int depth) const;


template<class BS> Void draw_hybrid_basic_set(CanvasInterface& canvas, const DiscreteLocation& location, const Variables2d& axes, const HybridBasicSet<BS>& set)
{
    if(set.location()==location) {
        Projection2d projection(set.euclidean_set().dimension(),set.space().index(axes.x_variable()),set.space().index(axes.y_variable()));
        set.euclidean_set().draw(canvas,projection);
    }
}




template<class IVL> Void HybridBox<IVL>::draw(CanvasInterface& c, const Set<DiscreteLocation>& qs, const Variables2d& p) const {
    DiscreteLocation const& q=this->location();
    if(qs.empty() || qs.contains(q)) {
        Ariadne::draw_hybrid_basic_set(c,q,p,static_cast<HybridBasicSet<Box<IVL>>const&>(*this));
    }
}

template<class IVL> Void HybridBoxes<IVL>::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axis_variables) const {
    for(auto loc_iter=this->begin(); loc_iter!=this->end(); ++loc_iter) {
        if(locations.empty() || locations.contains(loc_iter->first)) {
            RealSpace const& space=this->space(loc_iter->first);
            Projection2d projection(space.dimension(),space.index(axis_variables.x_variable()),space.index(axis_variables.y_variable()));
            loc_iter->second.euclidean_set().draw(canvas,projection);
        }
    }
}



template<class DS, class HBS>
HybridSpaceSetConstIterator<DS,HBS>::
HybridSpaceSetConstIterator(const std::map<DiscreteLocation,DS>& map, const HybridSpace& spc, Bool end)
    : _loc_begin(map.begin()),
      _loc_end(map.end()),
      _loc_iter(end?_loc_end:_loc_begin),
      hspc(spc)
{
    if(_loc_iter!=_loc_end) {
        _bs_iter=_loc_iter->second.begin();
        this->increment_loc();
    }
}

template<class DS, class HBS> inline
Void
HybridSpaceSetConstIterator<DS,HBS>::increment_loc()
{
    while(_bs_iter==_loc_iter->second.end()) {
        ++_loc_iter;
        if(_loc_iter==_loc_end) { return; }
        _bs_iter=_loc_iter->second.begin();
    }
}



HybridGridTreeSet::ConstIterator HybridGridTreeSet::begin() const {
    return ConstIterator(this->_map,this->_hgrid.space(),false);
}

HybridGridTreeSet::ConstIterator HybridGridTreeSet::end() const {
    return ConstIterator(this->_map,this->_hgrid.space(),true);
}


Void HybridGridTreeSet::adjoin(const HybridGridCell& hgc) {
    this->_provide_location(hgc.location()).adjoin(hgc.euclidean_set());
}

Void HybridGridTreeSet::adjoin(const ListSet<HybridGridCell>& hgcls) {
    for(ListSet<HybridGridCell>::ConstIterator iter=hgcls.begin(); iter!=hgcls.end(); ++iter) {
        this ->adjoin(*iter);
    }
}

Void HybridGridTreeSet::adjoin(const HybridGridTreeSet& hgts) {
    for(HybridGridTreeSet::LocationsConstIterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
        this->_provide_location(_loc_iter->first).adjoin(_loc_iter->second);
    }
}

Void HybridGridTreeSet::remove(const HybridGridTreeSet& hgts) {
    for(HybridGridTreeSet::LocationsConstIterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
        this->_provide_location(_loc_iter->first).remove(_loc_iter->second);
    }
}

Void HybridGridTreeSet::restrict(const HybridGridTreeSet& hgts) {
    for(HybridGridTreeSet::LocationsConstIterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
        this->_provide_location(_loc_iter->first).restrict(_loc_iter->second);
    }
}

Void HybridGridTreeSet::restrict_to_height(Nat h) {
    for(LocationsIterator _loc_iter=locations_begin(); _loc_iter!=locations_end(); ++_loc_iter) {
        _loc_iter->second.restrict_to_height(h);
    }
}

Void HybridGridTreeSet::adjoin_inner_approximation(const HybridExactBoxes& hbxs, const Int depth) {
    for(HybridExactBoxes::ConstIterator _loc_iter=hbxs.begin();
            _loc_iter!=hbxs.end(); ++_loc_iter) {
        DiscreteLocation const& loc=_loc_iter->first;
        ExpressionSet<ExactBoxType> const& vbx=_loc_iter->second;
        ARIADNE_ASSERT(vbx.space() == this->space(loc));
        this->_provide_location(loc).adjoin_inner_approximation(vbx.euclidean_set(),depth);
    }
}

Void HybridGridTreeSet::adjoin_lower_approximation(const HybridOvertSetInterface& hs, const Int height, const Int depth) {
    Set<DiscreteLocation> hlocs=dynamic_cast<const HybridBoundedSetInterface&>(hs).locations();
    for(Set<DiscreteLocation>::ConstIterator _loc_iter=hlocs.begin();
            _loc_iter!=hlocs.end(); ++_loc_iter) {
        DiscreteLocation loc=*_loc_iter;
        RealSpace spc=this->space(loc);
        this->_provide_location(loc).adjoin_lower_approximation(hs.euclidean_set(loc,spc),height,depth);
    }
}

Void HybridGridTreeSet::adjoin_outer_approximation(const HybridCompactSetInterface& hs, const Int depth) {
    Set<DiscreteLocation> hlocs=hs.locations();
    for(Set<DiscreteLocation>::ConstIterator _loc_iter=hlocs.begin();
            _loc_iter!=hlocs.end(); ++_loc_iter) {
        DiscreteLocation loc=*_loc_iter;
        RealSpace spc=this->space(loc);
        this->_provide_location(loc).adjoin_outer_approximation(hs.euclidean_set(loc,spc),depth);
    }
}

Void HybridGridTreeSet::adjoin_outer_approximation(const HybridExactBoxes& hbxs, const Int depth) {
    for(HybridExactBoxes::ConstIterator _loc_iter=hbxs.begin();
            _loc_iter!=hbxs.end(); ++_loc_iter) {
        DiscreteLocation const& loc=_loc_iter->first;
        ExpressionSet<ExactBoxType> const& vbx=_loc_iter->second;
        ARIADNE_ASSERT(vbx.space() == this->space(loc));
        this->_provide_location(_loc_iter->first).adjoin_outer_approximation(vbx.euclidean_set(),depth);
    }
}


GridTreeSet& HybridGridTreeSet::operator[](DiscreteLocation q) {
    return this->_provide_location(q);
}

const GridTreeSet& HybridGridTreeSet::operator[](DiscreteLocation q) const {
    ARIADNE_ASSERT_MSG(this->has_location(q),"q="<<q);
    return const_cast<HybridGridTreeSet*>(this)->_provide_location(q);
}

Bool HybridGridTreeSet::is_empty() const {
    for(LocationsConstIterator _loc_iter=this->locations_begin();
        _loc_iter!=this->locations_end(); ++_loc_iter) {
        if(!_loc_iter->second.is_empty()) { return false; }
    }
    return true;
}

SizeType HybridGridTreeSet::size() const {
    SizeType result=0;
    for(LocationsConstIterator _loc_iter=this->locations_begin(); _loc_iter!=this->locations_end(); ++_loc_iter) {
        result+=_loc_iter->second.size();
    }
    return result;
}

HybridListSet<ExactBoxType> HybridGridTreeSet::boxes() const {
    HybridListSet<ExactBoxType> result;
    for(ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        result.adjoin(iter->location(),iter->euclidean_set().box());
    }
    return result;
}

Void HybridGridTreeSet::mince(Int depth) {
    for(LocationsIterator _loc_iter=this->locations_begin();
        _loc_iter!=this->locations_end(); ++_loc_iter) {
        _loc_iter->second.mince(depth);
    }
}

Void HybridGridTreeSet::recombine() {
    for(LocationsIterator _loc_iter=this->locations_begin();
        _loc_iter!=this->locations_end(); ++_loc_iter) {
        _loc_iter->second.recombine();
    }
}

ValidatedSierpinskian HybridGridTreeSet::separated(const HybridExactBox& hbx) const {
    LocationsConstIterator _loc_iter = this->_map.find( hbx.location() );
    return _loc_iter != this->locations_end() || _loc_iter->second.separated( hbx.euclidean_set() );
}

ValidatedSierpinskian HybridGridTreeSet::overlaps(const HybridExactBox& hbx) const {
    LocationsConstIterator _loc_iter = this->_map.find( hbx.location() );
    return _loc_iter != this->locations_end() && _loc_iter->second.overlaps( hbx.euclidean_set() );
}

ValidatedSierpinskian HybridGridTreeSet::covers(const HybridExactBox& hbx) const {
    LocationsConstIterator _loc_iter=this->_map.find(hbx.location());
    return _loc_iter!=this->locations_end() && _loc_iter->second.covers( hbx.euclidean_set() );
}

ValidatedSierpinskian HybridGridTreeSet::inside(const HybridExactBoxes& hbx) const  {
    for( LocationsConstIterator _loc_iter = this->locations_begin(); _loc_iter != this->locations_end(); ++_loc_iter ) {
        if( !_loc_iter->second.is_empty() ) {
            DiscreteLocation const& loc = _loc_iter->first;
            RealSpace spc=this->space(loc);
            return this->euclidean_set(loc).inside(hbx.euclidean_set(loc,spc));
        }
    }
    return true;
}

HybridUpperBoxes HybridGridTreeSet::bounding_box() const {
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

OutputStream& HybridGridTreeSet::write(OutputStream& os) const {
    return os << this->_map;
}

Void HybridGridTreeSet::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axis_variables) const {
    for(LocationsConstIterator loc_iter=this->locations_begin(); loc_iter!=this->locations_end(); ++loc_iter) {
        if(locations.empty() || locations.contains(loc_iter->first)) {
            RealSpace const& space=this->space(loc_iter->first);
            Projection2d projection(space.dimension(),space.index(axis_variables.x_variable()),space.index(axis_variables.y_variable()));
            loc_iter->second.draw(canvas,projection);
        }
    }
}


GridTreeSet& HybridGridTreeSet::_provide_location(const DiscreteLocation& q) {
    std::map<DiscreteLocation,GridTreeSet>::iterator iter=this->_map.find(q);
    if(iter==this->_map.end()) {
        this->_map.insert(std::make_pair(q,GridTreeSet(this->_hgrid[q])));
        iter=this->_map.find(q);
    }
    return iter->second; }


template<> String class_name<RealInterval>() { return "RealInterval"; }
template<> String class_name<InterpolatedCurve>() { return "InterpolatedCurve"; }
template<> String class_name<Box<RealInterval>>() { return "RealBox"; }
template<> String class_name<ExactBoxType>() { return "ExactFloat64Box"; }

template<> String class_name<GridCell>() { return "GridCell"; }


// Instantiations
template class HybridBasicSet<ExactBoxType>;

template class HybridPoint<ExactNumericType>;
template class HybridBox<ExactIntervalType>;
template class HybridBoxes<ExactIntervalType>;

template class HybridBox<RealInterval>;


template class HybridSpaceSetConstIterator<GridTreeSet, HybridGridCell>;

} // namespace Ariadne
