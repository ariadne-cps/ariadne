/***************************************************************************
 *            hybrid_set.h
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

/*! \file hybrid_set.h
 *  \brief Sets in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_SET_H
#define ARIADNE_HYBRID_SET_H

#include <map>

#include <boost/iterator.hpp>
#include <boost/iterator_adaptors.hpp>

#include <memory>

#include "utility/macros.h"
#include "utility/stlio.h"
#include "utility/declarations.h"
#include "utility/container.h"
#include "geometry/function_set.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/curve.h"

#include "expression/expression_set.h"

#include "hybrid/hybrid_set_interface.h"
#include "hybrid/hybrid_space.h"
#include "hybrid/hybrid_grid.h"
#include "geometry/point.h"
#include "geometry/box.h"

#ifdef ARIADNE_ENABLE_SERIALIZATION
#include "output/serialization.h"
#endif /* ARIADNE_ENABLE_SERIALIZATION */

#include "hybrid/hybrid_graphics_interface.h"

namespace Ariadne {

// Classes defined in this file
class HybridPoint;
class HybridBoxType;
typedef HybridBoxType HybridUpperBoxType;
class HybridBoxes;
typedef HybridBoxes HybridUpperBoxes;
class HybridGridTreeSet;
template<class ES> class HybridBasicSet;
template<class ES> class ListSet< HybridBasicSet<ES> >;

class HybridBoundedConstraintSet;
typedef HybridBoundedConstraintSet HybridSet;

template<class DS, class HBS> class HybridSetConstIterator;



//! \ingroup ExpressionSetSubModule
//! \ingroup HybridSetSubModule
//! \brief A hybrid set defined by a box in a single location.
// FIXME: Merge with HybridBox
class HybridBoxSet
    : public Pair<DiscreteLocation,RealVariablesBox>
    , public virtual HybridDrawableInterface
{
  public:
    HybridBoxSet(const DiscreteLocation& loc, const RealVariablesBox& bx)
        : Pair<DiscreteLocation,RealVariablesBox>(loc,bx) { }
    HybridBoxSet(const DiscreteLocation& loc, const RealSpace& spc, const RealBox& bx)
        : Pair<DiscreteLocation,RealVariablesBox>(loc,RealVariablesBox(spc,bx)) { }
    //! \brief The location in which the box is defined.
    DiscreteLocation location() const { return this->first; }
    //! \brief The active variables in the location \a loc.
    virtual Set<RealVariable> variables() const { return this->second.variables(); }
    //! \brief The subset of \f$\mathbb{R}^n\f$ obtained by ordering the variables as defined by \a spc.
    RealBox euclidean_set(const RealSpace& spc) const { return this->second.euclidean_set(spc); }

    virtual Void draw(CanvasInterface&, const Set<DiscreteLocation>&, const Variables2d&) const override;
};

//! \ingroup ExpressionSetSubModule
//! \ingroup HybridSetSubModule
//! \brief A hybrid set defined by a constraint system in each location.
class HybridConstraintSet
    : public virtual HybridRegularSetInterface
{
    Map<DiscreteLocation, RealExpressionConstraintSet> _sets;
  public:
    HybridConstraintSet();
    //! \brief Construct a set in a single \a location with a list of \a bounds on the variables and nonlinear \a constraints.
    HybridConstraintSet(const DiscreteLocation& location,
                            const List<ContinuousPredicate>& constraints);

    virtual HybridConstraintSet* clone() const override;

    //! \brief The active variables in the location \a loc.
    virtual Set<RealVariable> variables(DiscreteLocation loc) const override;
    //! \brief The subset of \f$\mathbb{R}^n\f$ obtained by restricting to location \a loc and ordering the variables as defined by \a spc.
    ConstraintSet const euclidean_set(DiscreteLocation loc, RealSpace spc) const;

    virtual ValidatedSierpinskian overlaps(const HybridBoxType& bx) const override;
    virtual ValidatedSierpinskian separated(const HybridBoxType& bx) const override;
    virtual ValidatedSierpinskian covers(const HybridBoxType& bx) const override;

    virtual OutputStream& write(OutputStream& os) const override;
  protected:
    virtual RegularSetInterface* _euclidean_set(DiscreteLocation loc, RealSpace spc) const override;
};

//! \ingroup ExpressionSetSubModule
//! \ingroup HybridSetSubModule
//! \brief A hybrid set defined by the intersection of a box and a constraint system in each location.
class HybridBoundedConstraintSet
    : public virtual HybridSetInterface
    , public virtual HybridDrawableInterface
{
    Map<DiscreteLocation, RealExpressionBoundedConstraintSet> _sets;
  public:
    HybridBoundedConstraintSet();
    HybridBoundedConstraintSet(const HybridBoxSet& bx);
    HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                   const RealVariablesBox& bx);
    HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                   const InitializerList<RealVariableInterval>& bnd);
    //! \brief Construct a set in a single \a location with a list of \a bounds on the variables and nonlinear \a constraints.
    HybridBoundedConstraintSet(const DiscreteLocation& location,
                                   const InitializerList<RealVariableInterval>& bounds,
                                   const InitializerList<ContinuousPredicate>& constraints);

    DiscreteLocation location() const;

    virtual HybridBoundedConstraintSet* clone() const override;

    //! \brief The set of discrete locations in which the set is nontrivial.
    virtual Set<DiscreteLocation> locations() const override;
    //! \brief The active variables in the location \a loc.
    virtual Set<RealVariable> variables(DiscreteLocation loc) const override;
    //! \brief The subset of \f$\mathbb{R}^n\f$ obtained by restricting to location \a loc and ordering the variables as defined by \a spc.
    BoundedConstraintSet const euclidean_set(DiscreteLocation loc, RealSpace spc) const;

    virtual ValidatedSierpinskian overlaps(const HybridBoxType& bx) const override;
    virtual ValidatedSierpinskian inside(const HybridBoxes& bx) const override;

    virtual ValidatedSierpinskian separated(const HybridBoxType& bx) const override;
    virtual ValidatedSierpinskian covers(const HybridBoxType& bx) const override;
    virtual HybridUpperBoxes bounding_box() const override;

    virtual OutputStream& write(OutputStream& os) const override;
    virtual Void draw(CanvasInterface&, const Set<DiscreteLocation>&, const Variables2d&) const override;
  protected:
    virtual BoundedConstraintSet* _euclidean_set(DiscreteLocation loc, RealSpace spc) const override;
};



//! \ingroup HybridSetSubModule
//! \brief A hybrid set defined in a single location obtained from a Euclidean set by naming variables.
template<class EBS>
class HybridBasicSet
{
    Tuple<DiscreteLocation,RealSpace,EBS> _tuple;
  public:
    //! \brief The type of the Euclidean set used to describe the hybrid set.
    typedef EBS ContinuousSetType;
    HybridBasicSet()
        : _tuple(DiscreteLocation(),RealSpace(),EBS()) { }
    //! \brief Construct a set in location \a loc, with variables ordered by \a spc, defined by Euclidean set \a ebs.
    HybridBasicSet(const DiscreteLocation& loc, const RealSpace& spc, const ContinuousSetType& ebs)
        : _tuple(loc,spc,ebs) { }
    //! \brief The location the set is contained in.
    const DiscreteLocation& location() const { return get_first(this->_tuple); }
    //! \brief The ordering of variables used to define the set.
    const RealSpace& space() const { return get_second(this->_tuple); }
    //! \brief The continuous Euclidean subset.
    const ContinuousSetType& continuous_set() const { return get_third(this->_tuple); }
    ContinuousSetType& continuous_set() { return get_third(this->_tuple); }

    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const HybridBasicSet<EBS>& hbs) {
        return os << "HybridBasicSet{ " << hbs.location() << ", " << hbs.space() << ", " << hbs.continuous_set() << " )";
    }

    //! \brief Test for equality
    friend Bool operator==(const HybridBasicSet<EBS>& hset1, const HybridBasicSet<EBS>& hset2) {
        if(hset1.location()==hset2.location()) {
            ARIADNE_ASSERT(hset1.space()==hset2.space());
            return hset1.continuous_set() == hset2.continuous_set();
        } else {
            return false;
        }
    }
};

class HybridPoint
    : public HybridBasicSet<ExactPoint>
{
    typedef ExactPoint::ValueType ValueType;
  public:
    HybridPoint() : HybridBasicSet<ExactPoint>() { }
    HybridPoint(const DiscreteLocation& q, const RealSpace& spc, const ExactPoint& pt) : HybridBasicSet<ExactPoint>(q,spc,pt) { }
    HybridPoint(const DiscreteLocation& q, const Map<RealVariable,ValueType>& val);
    HybridPoint(const DiscreteLocation& q, const Map<RealVariable,Real>& val);
    HybridPoint(const DiscreteLocation& q, const List< Assignment<RealVariable,Real> >& val);
    HybridPoint(const DiscreteLocation& q, const InitializerList< Assignment<RealVariable,Real> >& val);
    ExactPoint& point() { return this->continuous_set(); }
    const ExactPoint& point() const { return this->continuous_set(); }
    const ExactPoint& real_values() const { return this->continuous_set(); }
    Map<RealVariable,ValueType> values() const;
};

//! \ingroup HybridSetSubModule
//! \brief A box in a location of a hybrid space.
//! \details Primarily used as a basic set against which abstract set properties can be tested.
// FIXME: Merge with HybridBoxSet
class HybridBoxType
    : public HybridBasicSet<ExactBoxType>
    , public virtual HybridDrawableInterface
{
  public:
    HybridBoxType(const DiscreteLocation& loc, const ExactVariablesBoxType& bx)
        : HybridBasicSet<ExactBoxType>(loc,bx.space(),bx.continuous_set()) { }
    HybridBoxType(const DiscreteLocation& loc, const List<RealVariableInterval>& bnds)
        : HybridBoxType(loc,ExactVariablesBoxType(bnds)) { }
    HybridBoxType(const DiscreteLocation& loc, const RealSpace& spc, const ExactBoxType& bx)
        : HybridBasicSet<ExactBoxType>(loc,spc,bx) { }

    //! \brief The subset of \f$\mathbb{R}^n\f$ obtained by restricting to location \a loc and ordering the variables as defined by \a spc.
    ExactBoxType euclidean_set(const RealSpace& spc) const {
        if(spc==this->space()) { return this->continuous_set(); }
        else { return ExactVariablesBoxType(this->space(),this->continuous_set() ).euclidean_set(spc); }
    }

    virtual Void draw(CanvasInterface& c, const Set<DiscreteLocation>& q, const Variables2d& v) const override;
};

//! \ingroup HybridSetSubModule
//! \brief A collection of boxes, one in each location of a hybrid space.
//! \details Primarily used to represent bounds for a compact hybrid set.
class HybridBoxes
    : public Map<DiscreteLocation,ExactVariablesBoxType>
{
  public:
    //! \brief Set the continuous state set in location \a loc to \a vbx.
    Void insert(const DiscreteLocation& loc, const ExactVariablesBoxType& vbx) {
        this->Map<DiscreteLocation,ExactVariablesBoxType>::insert(loc,vbx); }
    //! \brief Set the continuous state set in location \a loc to box \a bx using \a spc to order the variables.
    Void insert(const DiscreteLocation& loc, const RealSpace& spc, const ExactBoxType& bx) {
        this->insert(loc,ExactVariablesBoxType(spc,bx)); }

    //! \brief The set of discrete locations in which the set is nontrivial.
    Set<DiscreteLocation> locations() const { return this->keys(); }
    //! \brief The ordering of variables used to define the Euclidean box in location \a loc.
    RealSpace const& space(const DiscreteLocation& loc) const { return this->operator[](loc).space(); }
    //! \brief The Euclidean box in location \a loc.
    ExactBoxType const& continuous_set(const DiscreteLocation& loc) const { return this->operator[](loc).continuous_set(); }

    //! \brief The subset of \f$\mathbb{R}^n\f$ obtained by restricting to location \a loc and ordering the variables as defined by \a spc.
    ExactBoxType euclidean_set(const DiscreteLocation& loc, const RealSpace& spc) const {
        return this->operator[](loc).euclidean_set(spc); }
};


template<class ES> class HybridListSetConstIterator;

template<class ES>
class HybridListSetConstIterator
    : public boost::iterator_facade<
                HybridListSetConstIterator<ES>,
                HybridBasicSet<ES>,
                boost::forward_traversal_tag,
                HybridBasicSet<ES> const&
             >

{
  public:
    typedef HybridBasicSet<ES> const& Reference;
  public:
    HybridListSetConstIterator(const Map<DiscreteLocation,Pair<RealSpace,ListSet<ES> > >& map, Bool);
    Bool equal(const HybridListSetConstIterator<ES>&) const;
    const HybridBasicSet<ES>& dereference() const;
    Void increment();
  private:
    Void _increment_loc();
  private:
    typedef typename Map< DiscreteLocation, Pair<RealSpace,ListSet<ES> > >::ConstIterator LocationsIterator;
    typedef typename ListSet<ES>::ConstIterator BasicSetIterator;
    LocationsIterator _loc_begin;
    LocationsIterator _loc_end;
    LocationsIterator _loc_iter;
    BasicSetIterator _bs_iter;
    mutable HybridBasicSet<ES> hybrid_set;
};


template<class ES> inline
HybridListSetConstIterator<ES>::
HybridListSetConstIterator(const Map<DiscreteLocation,Pair<RealSpace,ListSet<ES> > >& map, Bool end)
    : _loc_begin(map.begin()),
      _loc_end(map.end()),
      _loc_iter(end?_loc_end:_loc_begin)
{
    if(_loc_iter!=_loc_end) {
        _bs_iter=_loc_iter->second.second.begin();
        this->_increment_loc();
    }
}


template<class ES> inline
Bool
HybridListSetConstIterator<ES>::equal(const HybridListSetConstIterator<ES>& other) const
{
    return this->_loc_iter==other._loc_iter && (this->_loc_iter==this->_loc_end || this->_bs_iter==other._bs_iter);
}


template<class ES> inline
HybridBasicSet<ES> const&
HybridListSetConstIterator<ES>::dereference() const
{
    this->hybrid_set=HybridBasicSet<ES>(_loc_iter->first,_loc_iter->second.first,*this->_bs_iter);
    return this->hybrid_set;
}


template<class ES> inline
Void
HybridListSetConstIterator<ES>::increment()
{
    ++this->_bs_iter;
    this->_increment_loc();
}

template<class ES> inline
Void
HybridListSetConstIterator<ES>::_increment_loc()
{
    while(_bs_iter==_loc_iter->second.second.end()) {
        ++_loc_iter;
        if(_loc_iter==_loc_end) { return; }
        _bs_iter=_loc_iter->second.second.begin();
    }
}


template<class ES> class HybridListSet;
template<class ES> OutputStream& operator<<(OutputStream& os, const HybridListSet<ES>& hls);

//! \ingroup HybridSetSubModule
//! A set comprising a %ListSet in each location.
template<class ES>
class HybridListSet
{
    Map<DiscreteLocation, Pair< RealSpace, ListSet<ES> > > _locations;
  public:
    typedef typename Map<DiscreteLocation,Pair< RealSpace, ListSet<ES> > >::ConstIterator LocationsConstIterator;
    typedef typename Map<DiscreteLocation,Pair< RealSpace, ListSet<ES> > >::Iterator LocationsIterator;
    typedef HybridListSetConstIterator<ES> ConstIterator;

    HybridListSet() { }
    HybridListSet(const HybridBasicSet<ES>& hes) { this->adjoin(hes); }

    ConstIterator begin() const {
        return ConstIterator(this->_locations,false); }
    ConstIterator end() const {
        return ConstIterator(this->_locations,true); }

    LocationsConstIterator locations_begin() const {
        return this->_locations.begin(); }
    LocationsConstIterator locations_end() const {
        return this->_locations.end(); }

    //! \brief Returns the number of basic hybrid sets forming this object.
    SizeType size() const {
        SizeType s = 0;
        for(LocationsConstIterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            s += _loc_iter->second.second.size();
        }
        return s;
    }

    //using std::map< DiscreteLocation,ListSet<ES> >::insert;

    //using std::map<DiscreteLocation,ListSet<ES> >::operator[];
    //ListSet<ES>& operator[](const DiscreteLocation& q) {
    //    return this->_locations[q].second; }
    const ListSet<ES>& operator[](const DiscreteLocation& q) const {
        ARIADNE_ASSERT_MSG(this->find(q)!=this->locations_end(),(*this)<<" has no location "<<q);
        return this->find(q)->second->second; }

    Void insert(const DiscreteLocation& loc, const RealSpace& spc) {
        this->_locations.insert(loc, make_pair(spc,ListSet<ES>())); }
    Void adjoin(const DiscreteLocation& loc, const RealSpace& spc, const ES& es) {
        LocationsIterator loc_iter=this->_locations.find(loc);
        if(loc_iter!=this->_locations.end()) {
            ARIADNE_ASSERT_MSG(loc_iter->second.first==spc,
                               "Space "<<loc_iter->second.first<<" of location "<<loc_iter->first<<" differs from "<<spc); }
        else { this->insert(loc,spc); loc_iter=this->_locations.find(loc); }
        loc_iter->second.second.adjoin(es); }
    Void adjoin(const DiscreteLocation& loc, const ES& es) {
        ARIADNE_ASSERT_MSG(this->_locations.find(loc)!=this->_locations.end(),(*this)<<" has no location "<<loc);
        this->_locations[loc].second.adjoin(es); }
    Void adjoin(const HybridBasicSet<ES>& hes) {
        this->adjoin(hes.location(),hes.space(),hes.continuous_set()); }
    Void adjoin(const HybridListSet<ES>& hls) {
        for(LocationsConstIterator _loc_iter=hls.locations_begin();
            _loc_iter!=hls.locations_end(); ++_loc_iter) {
            (*this)[_loc_iter->first].adjoin(_loc_iter->second); } }

    HybridListSet<UpperBoxType> bounding_boxes() const {
        HybridListSet<UpperBoxType> result;
        for(LocationsConstIterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            result[_loc_iter->first]=_loc_iter->second.bounding_boxes(); }
        return result; }

    friend
    OutputStream& operator<<(OutputStream& os, const HybridListSet<ES>& hls) {
        return os << hls._locations; }
};

template<class ES>
class ListSet< HybridBasicSet<ES> >
    : public HybridListSet<ES>
{
  public:
    ListSet() { }
    ListSet(const HybridBasicSet<ES>& hes) { this->adjoin(hes); }
};


template<class ES>
OutputStream&
operator<<(OutputStream& os,
           const ListSet< HybridBasicSet<ES> >& ls) {
    const std::map< DiscreteLocation, ListSet<ES> >& hls=ls;
    return os << "HybridListSet" << hls;
}




//! \ingroup HybridSetModule
//! \brief A cell associated with a grid in a hybrid space.
class HybridGridCell
    : public HybridBasicSet<GridCell>
{
  public:
    HybridGridCell()
        : HybridBasicSet<GridCell>() { }
    HybridGridCell(DiscreteLocation q,const RealSpace& s,const GridCell& gc)
        : HybridBasicSet<GridCell>(q,s,gc) { }
    HybridBoxType box() const { return HybridBoxType(this->location(),this->space(),this->continuous_set().box()); }
};

class HybridGridTreeSet;



template<class HDS1, class HDS2>
Void adjoin_denotable_set(HDS1& hds1, const HDS2 hds2) {
    for(typename HDS2::LocationsConstIterator loc2_iter=
            hds2.locations_begin(); loc2_iter!=hds2.locations_end(); ++loc2_iter)
        {
            typename HDS1::LocationsIterator loc1_iter=hds1.find(loc2_iter->first);
            if(loc1_iter==hds1.locations_end()) {
                hds1.insert(make_pair(loc2_iter->first,
                                      typename HDS2::value_type::second_type(loc2_iter->second)));
            } else {
                loc1_iter->second.adjoin(loc2_iter->second);
            }
        }
}


template<class DS, class HBS>
class HybridSetConstIterator
    : public boost::iterator_facade<HybridSetConstIterator<DS,HBS>,
                                    HBS,
                                    boost::forward_traversal_tag,
                                    HBS const&
                                    >

{
  public:
    typedef HBS const& Reference;
  public:
    HybridSetConstIterator(const std::map<DiscreteLocation,DS>&, const HybridSpace& hspc, Bool);
    Bool equal(const HybridSetConstIterator<DS,HBS>&) const;
    const HBS& dereference() const;
    Void increment();
  private:
    Void increment_loc();
  private:
    typename std::map< DiscreteLocation,DS>::const_iterator _loc_begin;
    typename std::map< DiscreteLocation,DS>::const_iterator _loc_end;
    typename std::map< DiscreteLocation,DS>::const_iterator _loc_iter;
    typename DS::ConstIterator _bs_iter;
    HybridSpace hspc;
    mutable HBS hybrid_set;
};



//! \ingroup HybridSetSubModule
//! A set comprising a %GridTreeSet in each location.
class HybridGridTreeSet
    : public HybridDrawableInterface
{
  public:
    HybridGrid _hgrid;
    Map<DiscreteLocation,GridTreeSet> _map;
  public:
    typedef Map<DiscreteLocation,GridTreeSet>::Iterator LocationsIterator;
    typedef Map<DiscreteLocation,GridTreeSet>::ConstIterator LocationsConstIterator;
    typedef HybridSetConstIterator<GridTreeSet,HybridGridCell> ConstIterator;
  public:
    //!
    LocationsIterator locations_begin() { return this->_map.begin(); }
    //!
    LocationsIterator locations_end() { return this->_map.end(); }
    //!
    LocationsConstIterator locations_begin() const { return this->_map.begin(); }
    //!
    LocationsConstIterator locations_end() const { return this->_map.end(); }
    //!
    ConstIterator begin() const { return ConstIterator(this->_map,this->_hgrid.space(),false); }
    //!
    ConstIterator end() const { return ConstIterator(this->_map,this->_hgrid.space(),true); }
  public:
    //! Construct from a grid.
    HybridGridTreeSet(const HybridGrid& hgrid) : _hgrid(hgrid), _map() { }

    //!
    HybridGrid grid() const { return this->_hgrid; }

    //! Test if \a q is a location of the set i.e. corresponds to a valid location of the underlying hybrid space.
    Bool has_location(DiscreteLocation q) const { return _hgrid.has_location(q); }
    //! Test if \a q is a nontrivial location of the set i.e. contained in the map of <DiscreteLocation,GridTreeSet> pairs.
    Bool nontrivial_location(DiscreteLocation q) const { return _map.has_key(q); }
    //! The continuous state space corresponding to location \a q.
    RealSpace space(DiscreteLocation q) const { return _hgrid.space(q); }
    //! The continuous state space corresponding to location \a q.
    const GridTreeSet& continuous_set(DiscreteLocation q) const { return _map[q]; }
    //! The continuous state space corresponding to location \a q.
    const GridTreeSet& euclidean_set(DiscreteLocation q, const RealSpace& s) const {
        ARIADNE_ASSERT_MSG(s==this->space(q),"Variable ordering in HybridGridTreeSet location "<<q<<" is "<<this->space(q)<<", "
                                             "which does not match requested ordering "<<q);
        return _map[q]; }

    //!
    Void insert(DiscreteLocation q, const GridTreeSet& gts) {
        this->_map.insert(q,gts); }
    Void insert(Pair<DiscreteLocation,GridTreeSet>& qgts) {
        this->_map.insert(qgts); }

    //!
    Void adjoin(DiscreteLocation q, const GridCell& c) {
        this->_provide_location(q).adjoin(c); }

    //!
    Void adjoin(const HybridGridCell& hgc) {
        this->_provide_location(hgc.location()).adjoin(hgc.continuous_set()); }

    //!
    Void adjoin(const ListSet<HybridGridCell>& hgcls) {
        for(ListSet<HybridGridCell>::ConstIterator iter=hgcls.begin(); iter!=hgcls.end(); ++iter) {
            this ->adjoin(*iter); } }

    //!
    Void adjoin(const HybridGridTreeSet& hgts) {
        for(HybridGridTreeSet::LocationsConstIterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
            this->_provide_location(_loc_iter->first).adjoin(_loc_iter->second); } }

    //!
    Void remove(const HybridGridTreeSet& hgts) {
        for(HybridGridTreeSet::LocationsConstIterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
            this->_provide_location(_loc_iter->first).remove(_loc_iter->second); } }

    //!
    Void restrict(const HybridGridTreeSet& hgts) {
        for(HybridGridTreeSet::LocationsConstIterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
            this->_provide_location(_loc_iter->first).restrict(_loc_iter->second); } }

    //!
    Void restrict_to_height(Nat h) {
        for(LocationsIterator _loc_iter=locations_begin(); _loc_iter!=locations_end(); ++_loc_iter) {
            _loc_iter->second.restrict_to_height(h); } }

    //!
    Void adjoin_inner_approximation(const HybridBoxes& hbxs, const Int depth) {
        for(HybridBoxes::ConstIterator _loc_iter=hbxs.begin();
                _loc_iter!=hbxs.end(); ++_loc_iter) {
            DiscreteLocation const& loc=_loc_iter->first;
            ExactVariablesBoxType const& vbx=_loc_iter->second;
            ARIADNE_ASSERT(vbx.space() == this->space(loc));
            this->_provide_location(loc).adjoin_inner_approximation(vbx.continuous_set(),depth); } }

    //!
    Void adjoin_lower_approximation(const HybridOvertSetInterface& hs, const Int height, const Int depth) {
        Set<DiscreteLocation> hlocs=dynamic_cast<const HybridBoundedSetInterface&>(hs).locations();
        for(Set<DiscreteLocation>::ConstIterator _loc_iter=hlocs.begin();
                _loc_iter!=hlocs.end(); ++_loc_iter) {
            DiscreteLocation loc=*_loc_iter;
            RealSpace spc=this->space(loc);
            this->_provide_location(loc).adjoin_lower_approximation(hs.euclidean_set(loc,spc),height,depth); } }

    //!
    Void adjoin_outer_approximation(const HybridCompactSetInterface& hs, const Int depth) {
        Set<DiscreteLocation> hlocs=hs.locations();
        for(Set<DiscreteLocation>::ConstIterator _loc_iter=hlocs.begin();
                _loc_iter!=hlocs.end(); ++_loc_iter) {
            DiscreteLocation loc=*_loc_iter;
            RealSpace spc=this->space(loc);
            this->_provide_location(loc).adjoin_outer_approximation(hs.euclidean_set(loc,spc),depth); } }

    //!
    Void adjoin_outer_approximation(const HybridBoxes& hbxs, const Int depth) {
        for(HybridBoxes::ConstIterator _loc_iter=hbxs.begin();
                _loc_iter!=hbxs.end(); ++_loc_iter) {
            DiscreteLocation const& loc=_loc_iter->first;
            ExactVariablesBoxType const& vbx=_loc_iter->second;
            ARIADNE_ASSERT(vbx.space() == this->space(loc));
            this->_provide_location(_loc_iter->first).adjoin_outer_approximation(vbx.continuous_set(),depth); } }

    //!
    template<class S> Void adjoin_outer_approximation(DiscreteLocation q, const S& s) {
        this->_provide_location(q).adjoin_outer_approximation(s); }

    //!
    GridTreeSet& operator[](DiscreteLocation q) {
        return this->_provide_location(q);
    }

    //!
    const GridTreeSet& operator[](DiscreteLocation q) const {
        ARIADNE_ASSERT_MSG(this->has_location(q),"q="<<q);
        return const_cast<HybridGridTreeSet*>(this)->_provide_location(q);
    }

    //!
    Bool is_empty() const {
        for(LocationsConstIterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            if(!_loc_iter->second.is_empty()) { return false; } }
        return true; }

    //!
    SizeType size() const {
        SizeType result=0;
        for(LocationsConstIterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            result+=_loc_iter->second.size(); }
        return result; }

    //!
    HybridListSet<ExactBoxType> boxes() const {
        HybridListSet<ExactBoxType> result;
        for(ConstIterator iter=this->begin();
            iter!=this->end(); ++iter) {
            result.adjoin(iter->location(),iter->continuous_set().box()); }
        return result; }

    //!
    Void mince(Int depth) {
        for(LocationsIterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            _loc_iter->second.mince(depth); } }

    //!
    Void recombine() {
        for(LocationsIterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            _loc_iter->second.recombine(); } }
  public:
    //@{ \name HybridSetInterface methods

    //!
    HybridGridTreeSet* clone() const { return new HybridGridTreeSet(*this); }

    //!
    HybridSpace space() const { return this->grid().space(); }

    //!
    ValidatedSierpinskian separated(const HybridBoxType& hbx) const {
        LocationsConstIterator _loc_iter = this->_map.find( hbx.location() );
        return _loc_iter != this->locations_end() || _loc_iter->second.separated( hbx.continuous_set() );
    }

    //!
    ValidatedSierpinskian overlaps(const HybridBoxType& hbx) const {
        LocationsConstIterator _loc_iter = this->_map.find( hbx.location() );
        return _loc_iter != this->locations_end() && _loc_iter->second.overlaps( hbx.continuous_set() );
    }

    //!
    ValidatedSierpinskian covers(const HybridBoxType& hbx) const {
        LocationsConstIterator _loc_iter=this->_map.find(hbx.location());
        return _loc_iter!=this->locations_end() && _loc_iter->second.covers( hbx.continuous_set() );
    }

    //!
    ValidatedSierpinskian inside(const HybridBoxes& hbx) const  {
        for( LocationsConstIterator _loc_iter = this->locations_begin(); _loc_iter != this->locations_end(); ++_loc_iter ) {
            if( !_loc_iter->second.is_empty() ) {
                DiscreteLocation const& loc = _loc_iter->first;
                RealSpace spc=this->space(loc);
                return this->continuous_set(loc).inside(hbx.euclidean_set(loc,spc));
            }
        }
        return true;
    }

    //!
    HybridUpperBoxes bounding_box() const {
        HybridBoxes result;
        for( LocationsConstIterator _loc_iter = this->locations_begin(); _loc_iter != this->locations_end(); ++_loc_iter ) {
            if( !_loc_iter->second.is_empty() ) {
                DiscreteLocation const& loc = _loc_iter->first;
                RealSpace const& spc=this->space(loc);
                result.insert(loc,spc,cast_exact_box(_loc_iter->second.bounding_box()));
            }
        }
        return result;
    }

    //!
    OutputStream& write(OutputStream& os) const {
        return os << this->_map;
    }

    //!
    Void draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d&v) const;

    //!
    friend inline OutputStream&
    operator<<(OutputStream& os, const HybridGridTreeSet& hgts) {
        return os << "HybridGridTreeSet(" << hgts._map << ")"; }

  private:
    GridTreeSet& _provide_location(const DiscreteLocation& q) {
        std::map<DiscreteLocation,GridTreeSet>::iterator iter=this->_map.find(q);
        if(iter==this->_map.end()) {
            this->_map.insert(std::make_pair(q,GridTreeSet(this->_hgrid[q])));
            iter=this->_map.find(q);
        }
        return iter->second; }
};




template<class DS, class HBS> inline
HybridSetConstIterator<DS,HBS>::
HybridSetConstIterator(const std::map<DiscreteLocation,DS>& map, const HybridSpace& spc, Bool end)
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
Bool
HybridSetConstIterator<DS,HBS>::equal(const HybridSetConstIterator<DS,HBS>& other) const
{
    return this->_loc_iter==other._loc_iter && (this->_loc_iter==this->_loc_end || this->_bs_iter==other._bs_iter);
}


template<class DS, class HBS> inline
HBS const&
HybridSetConstIterator<DS,HBS>::dereference() const
{
    this->hybrid_set=HBS(_loc_iter->first,this->hspc[_loc_iter->first],*this->_bs_iter);
    return this->hybrid_set;
}


template<class DS, class HBS> inline
Void
HybridSetConstIterator<DS,HBS>::increment()
{
    ++this->_bs_iter;
    this->increment_loc();
}

template<class DS, class HBS> inline
Void
HybridSetConstIterator<DS,HBS>::increment_loc()
{
    while(_bs_iter==_loc_iter->second.end()) {
        ++_loc_iter;
        if(_loc_iter==_loc_end) { return; }
        _bs_iter=_loc_iter->second.begin();
    }
}


} // namespace Ariadne

#ifdef ARIADNE_ENABLE_SERIALIZATION
  namespace boost { namespace serialization {
  template<class A> void serialize(A& archive, const Ariadne::HybridGridTreeSet& set, const unsigned int version);
  template<class A> void serialize(A& archive, const Ariadne::DiscreteLocation& state, const unsigned int version);
  }}
#endif /* ARIADNE_ENABLE_SERIALIZATION */

#endif // ARIADNE_HYBRID_SET_H
