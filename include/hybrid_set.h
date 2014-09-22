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

#include "macros.h"
#include "stlio.h"
#include "container.h"
#include "function_set.h"
#include "list_set.h"
#include "grid_set.h"
#include "curve.h"

#include "expression_set.h"

#include "hybrid_set_interface.h"
#include "hybrid_space.h"
#include "hybrid_grid.h"
#include "point.h"
#include "box.h"

#ifdef ARIADNE_ENABLE_SERIALIZATION
#include "serialization.h"
#endif /* ARIADNE_ENABLE_SERIALIZATION */

#include "graphics_interface.h"

namespace Ariadne {

// Classes defined in this file
class HybridPoint;
class HybridBox;
class HybridBoxes;
class HybridGridTreeSet;
template<class ES> class HybridBasicSet;
template<class ES> class ListSet< HybridBasicSet<ES> >;

class HybridBoundedConstraintSet;
typedef HybridBoundedConstraintSet HybridSet;

template<class DS, class HBS> class HybridSetConstIterator;



//! \ingroup ExpressionSetSubModule
//! \ingroup HybridSetSubModule
//! \brief A hybrid set defined by a box in a single location.
class HybridBoxSet
    : public Pair<DiscreteLocation,RealVariablesBox>
    , public virtual HybridDrawableInterface
{
  public:
    HybridBoxSet(const DiscreteLocation& loc, const RealVariablesBox& bx)
        : Pair<DiscreteLocation,RealVariablesBox>(loc,bx) { }
    HybridBoxSet(const DiscreteLocation& loc, const RealSpace& spc, const BoxSet& bx)
        : Pair<DiscreteLocation,RealVariablesBox>(loc,RealVariablesBox(spc,bx)) { }
    //! \brief The location in which the box is defined.
    DiscreteLocation location() const { return this->first; }
    //! \brief The active variables in the location \a loc.
    virtual Set<RealVariable> variables() const { return this->second.variables(); }
    //! \brief The subset of \f$\mathbb{R}^n\f$ obtained by ordering the variables as defined by \a spc.
    BoxSet euclidean_set(const RealSpace& spc) const { return this->second.euclidean_set(spc); }

    virtual void draw(CanvasInterface&, const Set<DiscreteLocation>&, const Variables2d&) const override;
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

    virtual tribool overlaps(const HybridBox& bx) const override;
    virtual tribool separated(const HybridBox& bx) const override;
    virtual tribool covers(const HybridBox& bx) const override;

    virtual std::ostream& write(std::ostream& os) const override;
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
                                   const List<RealVariableInterval>& bnd);
    //! \brief Construct a set in a single \a location with a list of \a bounds on the variables and nonlinear \a constraints.
    HybridBoundedConstraintSet(const DiscreteLocation& location,
                                   const List<RealVariableInterval>& bounds,
                                   const List<ContinuousPredicate>& constraints);

    DiscreteLocation location() const;

    virtual HybridBoundedConstraintSet* clone() const override;

    //! \brief The set of discrete locations in which the set is nontrivial.
    virtual Set<DiscreteLocation> locations() const override;
    //! \brief The active variables in the location \a loc.
    virtual Set<RealVariable> variables(DiscreteLocation loc) const override;
    //! \brief The subset of \f$\mathbb{R}^n\f$ obtained by restricting to location \a loc and ordering the variables as defined by \a spc.
    BoundedConstraintSet const euclidean_set(DiscreteLocation loc, RealSpace spc) const;

    virtual tribool overlaps(const HybridBox& bx) const override;
    virtual tribool inside(const HybridBoxes& bx) const override;

    virtual tribool separated(const HybridBox& bx) const override;
    virtual tribool covers(const HybridBox& bx) const override;
    virtual HybridBoxes bounding_box() const override;

    virtual std::ostream& write(std::ostream& os) const override;
    virtual void draw(CanvasInterface&, const Set<DiscreteLocation>&, const Variables2d&) const override;
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
    const DiscreteLocation& location() const { return this->_tuple.first; }
    //! \brief The ordering of variables used to define the set.
    const RealSpace& space() const { return this->_tuple.second; }
    //! \brief The continuous Euclidean subset.
    const ContinuousSetType& continuous_set() const { return this->_tuple.third; }
    ContinuousSetType& continuous_set() { return this->_tuple.third; }
};
//! \relates HybridBasicSet \brief Write to an output stream.
template<class EBS> OutputStream& operator<<(OutputStream& os, const HybridBasicSet<EBS>& hbs) {
    return os << "HybridBasicSet{ " << hbs.location() << ", " << hbs.space() << ", " << hbs.continuous_set() << " )";
}

template<class EBS> inline
bool operator==(const HybridBasicSet<EBS> hset1, const HybridBasicSet<EBS>& hset2) {
    if(hset1.location()==hset2.location()) {
        ARIADNE_ASSERT(hset1.space()==hset2.space());
        return hset1.continuous_set() == hset2.continuous_set();
    } else {
        return false;
    }
}

class HybridPoint
    : public HybridBasicSet<Point>
{
    typedef Point::ValueType ValueType;
  public:
    HybridPoint() : HybridBasicSet<Point>() { }
    HybridPoint(const DiscreteLocation& q, const RealSpace& spc, const Point& pt) : HybridBasicSet<Point>(q,spc,pt) { }
    HybridPoint(const DiscreteLocation& q, const Map<Identifier,ValueType>& val);
    HybridPoint(const DiscreteLocation& q, const Map<Identifier,Real>& val);
    HybridPoint(const DiscreteLocation& q, const List< Assignment<RealVariable,Real> >& val);
    Point& point() { return this->continuous_set(); }
    const Point& point() const { return this->continuous_set(); }
    const Point& real_values() const { return this->continuous_set(); }
    Map<RealVariable,ValueType> values() const;
};

//! \ingroup HybridSetSubModule
//! \brief A box in a location of a hybrid space.
//! \details Primarily used as a basic set against whabstract set properties can be tested.
class HybridBox
    : public HybridBasicSet<Box>
    , public virtual HybridDrawableInterface
{
  public:
    HybridBox(const DiscreteLocation& loc, const VariablesBox& bx)
        : HybridBasicSet<Box>(loc,bx.space(),bx.continuous_set()) { }
    HybridBox(const DiscreteLocation& loc, const List<RealVariableInterval>& bnds)
        : HybridBox(loc,VariablesBox(bnds)) { }
    HybridBox(const DiscreteLocation& loc, const RealSpace& spc, const Box& bx)
        : HybridBasicSet<Box>(loc,spc,bx) { }

   //! \brief The subset of \f$\mathbb{R}^n\f$ obtained by restricting to location \a loc and ordering the variables as defined by \a spc.
    Box euclidean_set(const RealSpace& spc) const {
        if(spc==this->space()) { return this->continuous_set(); }
        else { return VariablesBox(this->space(),this->continuous_set() ).euclidean_set(spc); }
    }

    virtual void draw(CanvasInterface& c, const Set<DiscreteLocation>& q, const Variables2d& v) const override;
};

//! \ingroup HybridSetSubModule
//! \brief A collection of boxes, one in each location of a hybrid space.
//! \details Primarily used to represent bounds for a compact hybrid set.
class HybridBoxes
    : public Map<DiscreteLocation,VariablesBox>
{
  public:
    //! \brief Set the continuous state set in location \a loc to \a vbx.
    void insert(const DiscreteLocation& loc, const VariablesBox& vbx) {
        this->Map<DiscreteLocation,VariablesBox>::insert(loc,vbx); }
    //! \brief Set the continuous state set in location \a loc to box \a bx using \a spc to order the variables.
    void insert(const DiscreteLocation& loc, const RealSpace& spc, const Box& bx) {
        this->insert(loc,VariablesBox(spc,bx)); }

    //! \brief The set of discrete locations in which the set is nontrivial.
    Set<DiscreteLocation> locations() const { return this->keys(); }
    //! \brief The ordering of variables used to define the Euclidean box in location \a loc.
    RealSpace const& space(const DiscreteLocation& loc) const { return this->operator[](loc).space(); }
    //! \brief The Euclidean box in location \a loc.
    Box const& continuous_set(const DiscreteLocation& loc) const { return this->operator[](loc).continuous_set(); }

    //! \brief The subset of \f$\mathbb{R}^n\f$ obtained by restricting to location \a loc and ordering the variables as defined by \a spc.
    Box euclidean_set(const DiscreteLocation& loc, const RealSpace& spc) const {
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
    typedef HybridBasicSet<ES> const& reference;
  public:
    HybridListSetConstIterator(const Map<DiscreteLocation,Pair<RealSpace,ListSet<ES> > >& map, bool);
    bool equal(const HybridListSetConstIterator<ES>&) const;
    const HybridBasicSet<ES>& dereference() const;
    void increment();
  private:
    void _increment_loc();
  private:
    typedef typename Map< DiscreteLocation, Pair<RealSpace,ListSet<ES> > >::const_iterator LocationsIterator;
    typedef typename ListSet<ES>::const_iterator BasicSetIterator;
    LocationsIterator _loc_begin;
    LocationsIterator _loc_end;
    LocationsIterator _loc_iter;
    BasicSetIterator _bs_iter;
    mutable HybridBasicSet<ES> hybrid_set;
};


template<class ES> inline
HybridListSetConstIterator<ES>::
HybridListSetConstIterator(const Map<DiscreteLocation,Pair<RealSpace,ListSet<ES> > >& map, bool end)
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
bool
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
void
HybridListSetConstIterator<ES>::increment()
{
    ++this->_bs_iter;
    this->_increment_loc();
}

template<class ES> inline
void
HybridListSetConstIterator<ES>::_increment_loc()
{
    while(_bs_iter==_loc_iter->second.second.end()) {
        ++_loc_iter;
        if(_loc_iter==_loc_end) { return; }
        _bs_iter=_loc_iter->second.second.begin();
    }
}


template<class ES> class HybridListSet;
template<class ES> std::ostream& operator<<(std::ostream& os, const HybridListSet<ES>& hls);

//! \ingroup HybridSetSubModule
//! A set comprising a %ListSet in each location.
template<class ES>
class HybridListSet
{
    Map<DiscreteLocation, Pair< RealSpace, ListSet<ES> > > _locations;
  public:
    typedef typename Map<DiscreteLocation,Pair< RealSpace, ListSet<ES> > >::const_iterator locations_const_iterator;
    typedef typename Map<DiscreteLocation,Pair< RealSpace, ListSet<ES> > >::iterator locations_iterator;
    typedef HybridListSetConstIterator<ES> const_iterator;

    HybridListSet() { }
    HybridListSet(const HybridBasicSet<ES>& hes) { this->adjoin(hes); }

    const_iterator begin() const {
        return const_iterator(this->_locations,false); }
    const_iterator end() const {
        return const_iterator(this->_locations,true); }

    locations_const_iterator locations_begin() const {
        return this->_locations.begin(); }
    locations_const_iterator locations_end() const {
        return this->_locations.end(); }

    //! \brief Returns the number of basic hybrid sets forming this object.
    size_t size() const {
        size_t s = 0;
        for(locations_const_iterator _loc_iter=this->locations_begin();
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

    void insert(const DiscreteLocation& loc, const RealSpace& spc) {
        this->_locations.insert(loc, make_pair(spc,ListSet<ES>())); }
    void adjoin(const DiscreteLocation& loc, const RealSpace& spc, const ES& es) {
        locations_iterator loc_iter=this->_locations.find(loc);
        if(loc_iter!=this->_locations.end()) {
            ARIADNE_ASSERT_MSG(loc_iter->second.first==spc,
                               "Space "<<loc_iter->second.first<<" of location "<<loc_iter->first<<" differs from "<<spc); }
        else { this->insert(loc,spc); loc_iter=this->_locations.find(loc); }
        loc_iter->second.second.adjoin(es); }
    void adjoin(const DiscreteLocation& loc, const ES& es) {
        ARIADNE_ASSERT_MSG(this->_locations.find(loc)!=this->_locations.end(),(*this)<<" has no location "<<loc);
        this->_locations[loc].second.adjoin(es); }
    void adjoin(const HybridBasicSet<ES>& hes) {
        this->adjoin(hes.location(),hes.space(),hes.continuous_set()); }
    void adjoin(const HybridListSet<ES>& hls) {
        for(locations_const_iterator _loc_iter=hls.locations_begin();
            _loc_iter!=hls.locations_end(); ++_loc_iter) {
            (*this)[_loc_iter->first].adjoin(_loc_iter->second); } }

    HybridListSet<Box> bounding_boxes() const {
        HybridListSet<Box> result;
        for(locations_const_iterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            result[_loc_iter->first]=_loc_iter->second.bounding_boxes(); }
        return result; }

    friend
    std::ostream& operator<<(std::ostream& os, const HybridListSet<ES>& hls) {
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
std::ostream&
operator<<(std::ostream& os,
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
    HybridBox box() const { return HybridBox(this->location(),this->space(),this->continuous_set().box()); }
};

class HybridGridTreeSet;



template<class HDS1, class HDS2>
void adjoin_denotable_set(HDS1& hds1, const HDS2 hds2) {
    for(typename HDS2::locations_const_iterator loc2_iter=
            hds2.locations_begin(); loc2_iter!=hds2.locations_end(); ++loc2_iter)
        {
            typename HDS1::locations_iterator loc1_iter=hds1.find(loc2_iter->first);
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
    typedef HBS const& reference;
  public:
    HybridSetConstIterator(const std::map<DiscreteLocation,DS>&, const HybridSpace& hspc, bool);
    bool equal(const HybridSetConstIterator<DS,HBS>&) const;
    const HBS& dereference() const;
    void increment();
  private:
    void increment_loc();
  private:
    typename std::map< DiscreteLocation,DS>::const_iterator _loc_begin;
    typename std::map< DiscreteLocation,DS>::const_iterator _loc_end;
    typename std::map< DiscreteLocation,DS>::const_iterator _loc_iter;
    typename DS::const_iterator _bs_iter;
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
    typedef Map<DiscreteLocation,GridTreeSet>::iterator locations_iterator;
    typedef Map<DiscreteLocation,GridTreeSet>::const_iterator locations_const_iterator;
    typedef HybridSetConstIterator<GridTreeSet,HybridGridCell> const_iterator;
  public:
    //!
    locations_iterator locations_begin() { return this->_map.begin(); }
    //!
    locations_iterator locations_end() { return this->_map.end(); }
    //!
    locations_const_iterator locations_begin() const { return this->_map.begin(); }
    //!
    locations_const_iterator locations_end() const { return this->_map.end(); }
    //!
    const_iterator begin() const { return const_iterator(this->_map,this->_hgrid.space(),false); }
    //!
    const_iterator end() const { return const_iterator(this->_map,this->_hgrid.space(),true); }
  public:
    //! Construct from a grid.
    HybridGridTreeSet(const HybridGrid& hgrid) : _hgrid(hgrid), _map() { }

    //!
    HybridGrid grid() const { return this->_hgrid; }

    //! Test if \a q is a location of the set i.e. corresponds to a valid location of the underlying hybrid space.
    bool has_location(DiscreteLocation q) const { return _hgrid.has_location(q); }
    //! Test if \a q is a nontrivial location of the set i.e. contained in the map of <DiscreteLocation,GridTreeSet> pairs.
    bool nontrivial_location(DiscreteLocation q) const { return _map.has_key(q); }
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
    void insert(DiscreteLocation q, const GridTreeSet& gts) {
        this->_map.insert(q,gts); }
    void insert(std::pair<DiscreteLocation,GridTreeSet>& qgts) {
        this->_map.insert(qgts); }

    //!
    void adjoin(DiscreteLocation q, const GridCell& c) {
        this->_provide_location(q).adjoin(c); }

    //!
    void adjoin(const HybridGridCell& hgc) {
        this->_provide_location(hgc.location()).adjoin(hgc.continuous_set()); }

    //!
    void adjoin(const ListSet<HybridGridCell>& hgcls) {
        for(ListSet<HybridGridCell>::const_iterator iter=hgcls.begin(); iter!=hgcls.end(); ++iter) {
            this ->adjoin(*iter); } }

    //!
    void adjoin(const HybridGridTreeSet& hgts) {
        for(HybridGridTreeSet::locations_const_iterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
            this->_provide_location(_loc_iter->first).adjoin(_loc_iter->second); } }

    //!
    void remove(const HybridGridTreeSet& hgts) {
        for(HybridGridTreeSet::locations_const_iterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
            this->_provide_location(_loc_iter->first).remove(_loc_iter->second); } }

    //!
    void restrict(const HybridGridTreeSet& hgts) {
        for(HybridGridTreeSet::locations_const_iterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
            this->_provide_location(_loc_iter->first).restrict(_loc_iter->second); } }

    //!
    void restrict_to_height(uint h) {
        for(locations_iterator _loc_iter=locations_begin(); _loc_iter!=locations_end(); ++_loc_iter) {
            _loc_iter->second.restrict_to_height(h); } }

    //!
    void adjoin_inner_approximation(const HybridBoxes& hbxs, const int depth) {
        for(HybridBoxes::const_iterator _loc_iter=hbxs.begin();
                _loc_iter!=hbxs.end(); ++_loc_iter) {
            DiscreteLocation const& loc=_loc_iter->first;
            VariablesBox const& vbx=_loc_iter->second;
            ARIADNE_ASSERT(vbx.space() == this->space(loc));
            this->_provide_location(loc).adjoin_inner_approximation(vbx.continuous_set(),vbx.continuous_set(),depth); } }

    //!
    void adjoin_lower_approximation(const HybridOvertSetInterface& hs, const int height, const int depth) {
        Set<DiscreteLocation> hlocs=dynamic_cast<const HybridBoundedSetInterface&>(hs).locations();
        for(Set<DiscreteLocation>::const_iterator _loc_iter=hlocs.begin();
                _loc_iter!=hlocs.end(); ++_loc_iter) {
            DiscreteLocation loc=*_loc_iter;
            RealSpace spc=this->space(loc);
            this->_provide_location(loc).adjoin_lower_approximation(hs.euclidean_set(loc,spc),height,depth); } }

    //!
    void adjoin_outer_approximation(const HybridCompactSetInterface& hs, const int depth) {
        Set<DiscreteLocation> hlocs=hs.locations();
        for(Set<DiscreteLocation>::const_iterator _loc_iter=hlocs.begin();
                _loc_iter!=hlocs.end(); ++_loc_iter) {
            DiscreteLocation loc=*_loc_iter;
            RealSpace spc=this->space(loc);
            this->_provide_location(loc).adjoin_outer_approximation(hs.euclidean_set(loc,spc),depth); } }

    //!
    void adjoin_outer_approximation(const HybridBoxes& hbxs, const int depth) {
        for(HybridBoxes::const_iterator _loc_iter=hbxs.begin();
                _loc_iter!=hbxs.end(); ++_loc_iter) {
            DiscreteLocation const& loc=_loc_iter->first;
            VariablesBox const& vbx=_loc_iter->second;
            ARIADNE_ASSERT(vbx.space() == this->space(loc));
            this->_provide_location(_loc_iter->first).adjoin_outer_approximation(vbx.continuous_set(),depth); } }

    //!
    template<class S> void adjoin_outer_approximation(DiscreteLocation q, const S& s) {
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
    bool empty() const {
        for(locations_const_iterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            if(!_loc_iter->second.empty()) { return false; } }
        return true; }

    //!
    size_t size() const {
        size_t result=0;
        for(locations_const_iterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            result+=_loc_iter->second.size(); }
        return result; }

    //!
    HybridListSet<Box> boxes() const {
        HybridListSet<Box> result;
        for(const_iterator iter=this->begin();
            iter!=this->end(); ++iter) {
            result.adjoin(iter->location(),iter->continuous_set().box()); }
        return result; }

    //!
    void mince(int depth) {
        for(locations_iterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            _loc_iter->second.mince(depth); } }

    //!
    void recombine() {
        for(locations_iterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            _loc_iter->second.recombine(); } }
  public:
    //@{ \name HybridSetInterface methods

    //!
    HybridGridTreeSet* clone() const { return new HybridGridTreeSet(*this); }

    //!
    HybridSpace space() const { return this->grid().space(); }

    //!
    tribool separated(const HybridBox& hbx) const {
        locations_const_iterator _loc_iter = this->_map.find( hbx.location() );
        return _loc_iter != this->locations_end() || _loc_iter->second.separated( hbx.continuous_set() );
    }

    //!
    tribool overlaps(const HybridBox& hbx) const {
        locations_const_iterator _loc_iter = this->_map.find( hbx.location() );
        return _loc_iter != this->locations_end() && _loc_iter->second.overlaps( hbx.continuous_set() );
    }

    //!
    tribool covers(const HybridBox& hbx) const {
        locations_const_iterator _loc_iter=this->_map.find(hbx.location());
        return _loc_iter!=this->locations_end() && _loc_iter->second.covers( hbx.continuous_set() );
    }

    //!
    tribool inside(const HybridBoxes& hbx) const  {
        for( locations_const_iterator _loc_iter = this->locations_begin(); _loc_iter != this->locations_end(); ++_loc_iter ) {
            if( !_loc_iter->second.empty() ) {
                DiscreteLocation const& loc = _loc_iter->first;
                RealSpace spc=this->space(loc);
                return this->continuous_set(loc).inside(hbx.euclidean_set(loc,spc));
            }
        }
        return true;
    }

    //!
    HybridBoxes bounding_box() const {
        HybridBoxes result;
        for( locations_const_iterator _loc_iter = this->locations_begin(); _loc_iter != this->locations_end(); ++_loc_iter ) {
            if( !_loc_iter->second.empty() ) {
                DiscreteLocation const& loc = _loc_iter->first;
                RealSpace const& spc=this->space(loc);
                result.insert(loc,spc,_loc_iter->second.bounding_box());
            }
        }
        return result;
    }

    //!
    std::ostream& write(std::ostream& os) const {
        return os << this->_map;
    }

    //!
    void draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d&v) const;

    //!
    friend inline std::ostream&
    operator<<(std::ostream& os, const HybridGridTreeSet& hgts) {
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
HybridSetConstIterator(const std::map<DiscreteLocation,DS>& map, const HybridSpace& spc, bool end)
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
bool
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
void
HybridSetConstIterator<DS,HBS>::increment()
{
    ++this->_bs_iter;
    this->increment_loc();
}

template<class DS, class HBS> inline
void
HybridSetConstIterator<DS,HBS>::increment_loc()
{
    while(_bs_iter==_loc_iter->second.end()) {
        ++_loc_iter;
        if(_loc_iter==_loc_end) { return; }
        _bs_iter=_loc_iter->second.begin();
    }
}

#ifdef ARIADNE_ENABLE_SERIALIZATION
  template<class A> void serialize(A& archive, HybridGridTreeSet& set, const unsigned int version) {
      archive & static_cast<std::map<DiscreteLocation,GridTreeSet>&>(set); }
#endif /* ARIADNE_ENABLE_SERIALIZATION */


} // namespace Ariadne

#ifdef ARIADNE_ENABLE_SERIALIZATION
  namespace boost { namespace serialization {
  template<class A> void serialize(A& archive, const Ariadne::HybridGridTreeSet& set, const uint version);
  template<class A> void serialize(A& archive, const Ariadne::DiscreteLocation& state, const uint version);
  }}
#endif /* ARIADNE_ENABLE_SERIALIZATION */

#endif // ARIADNE_HYBRID_SET_H
