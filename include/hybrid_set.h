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

#include <boost/shared_ptr.hpp>

#include "macros.h"
#include "stlio.h"
#include "container.h"
#include "function_set.h"
#include "list_set.h"
#include "grid_set.h"
#include "curve.h"

#include "hybrid_set_interface.h"
#include "point.h"
#include "box.h"

#include "serialization.h"

#include "graphics_interface.h"

namespace Ariadne {

class HybridGridTreeSet;
class HybridImageSet;
class HybridConstraintSet;


template<class HBS> class HybridBasicSetExpression { };
template<class HDS> class HybridDenotableSetExpression { };

//! \brief A hybrid space \f$\bigsqcup_{q\in Q} \R^{d_q}\f$ with discrete states \f$Q\f$.
class HybridSpace
    : public std::map<DiscreteLocation,uint>
{
  public:
    //! \brief The interface satisified by bounded sets in the space.
    typedef HybridBoundedSetInterface BoundedSetInterfaceType;
    //! \brief The interface satisified by overt sets in the space.
    typedef HybridOvertSetInterface OvertSetInterfaceType;
    //! \brief The interface satisified by over sets in the space.
    typedef HybridOpenSetInterface OpenSetInterfaceType;
    //! \brief The interface satisified by closed sets in the space.
    typedef HybridClosedSetInterface ClosedSetInterfaceType;
    //! \brief The interface satisified by compact sets in the space.
    typedef HybridCompactSetInterface CompactSetInterfaceType;
    //! \brief The interface satisified by regular sets in the space.
    typedef HybridRegularSetInterface RegularSetInterfaceType;
    //! \brief The interface satisified by located sets in the space.
    typedef HybridLocatedSetInterface LocatedSetInterfaceType;
    //! \brief The type of approximations to sets in the space.
    typedef HybridGridTreeSet SetApproximationType;

    typedef std::map<DiscreteLocation,uint>::const_iterator
    locations_const_iterator;

    HybridSpace() : std::map<DiscreteLocation,uint>() { }
    template<class SET> HybridSpace(const std::map<DiscreteLocation,SET>& qsmap) {
        for(typename std::map<DiscreteLocation,SET>::const_iterator loc_iter
                =qsmap.begin(); loc_iter!=qsmap.end(); ++loc_iter) {
            this->insert(make_pair(loc_iter->first,loc_iter->second.dimension())); }
    }
    template<class HSET> HybridSpace(const HSET& set) {
        for(typename HSET::locations_const_iterator loc_iter
                =set.locations_begin(); loc_iter!=set.locations_end(); ++loc_iter) {
            this->insert(make_pair(loc_iter->first,loc_iter->second.dimension())); }
    }

    locations_const_iterator locations_begin() const {
        return this->std::map<DiscreteLocation,uint>::begin(); }
    locations_const_iterator locations_end() const {
        return this->std::map<DiscreteLocation,uint>::end(); }
};

inline
HybridBoxes
bounding_boxes(const std::map<DiscreteLocation,uint> space, Interval bound)
{
    HybridBoxes result;
    for(std::map<DiscreteLocation,uint>::const_iterator loc_iter=space.begin();
        loc_iter!=space.end(); ++loc_iter)
        {
            result.insert(make_pair(loc_iter->first,Box(loc_iter->second, bound)));
        }
    return result;
}


template<class BS>
class HybridBasicSet
    : public std::pair<DiscreteLocation,BS>
{
  public:
    typedef BS ContinuousStateSetType;
    HybridBasicSet(const DiscreteLocation& q, const BS& s) : std::pair<DiscreteLocation,BS>(q,s) { }
    HybridBasicSet(const std::pair<DiscreteLocation,BS>& p) : std::pair<DiscreteLocation,BS>(p) { }
    const DiscreteLocation& location() const { return this->first; }
    const ContinuousStateSetType& continuous_state_set() const { return this->second; }
};

template< class DS, class HBS >
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
    HybridSetConstIterator(const std::map<DiscreteLocation,DS>&, bool);
    bool equal(const HybridSetConstIterator<DS,HBS>&) const;
    const HBS& dereference() const;
    void increment();
  private:
    void increment_loc();
  private:
    typename std::map< DiscreteLocation,DS>::const_iterator loc_begin;
    typename std::map< DiscreteLocation,DS>::const_iterator loc_end;
    typename std::map< DiscreteLocation,DS>::const_iterator loc_iter;
    typename DS::const_iterator bs_iter;
    mutable HBS hybrid_set;
};


//! A set comprising an ImageSet in each location.
class HybridImageSet
    : public std::map<DiscreteLocation,ImageSet>
    , public HybridLocatedSetInterface
{
  public:
    typedef std::map<DiscreteLocation,ImageSet>::iterator locations_iterator;
    typedef std::map<DiscreteLocation,ImageSet>::const_iterator locations_const_iterator;
    locations_iterator locations_begin() {
        return this->std::map<DiscreteLocation,ImageSet>::begin(); }
    locations_iterator locations_end() {
        return this->std::map<DiscreteLocation,ImageSet>::end(); }
    locations_const_iterator locations_begin() const {
        return this->std::map<DiscreteLocation,ImageSet>::begin(); }
    locations_const_iterator locations_end() const {
        return this->std::map<DiscreteLocation,ImageSet>::end(); }

    using std::map<DiscreteLocation,ImageSet>::insert;

    virtual HybridImageSet* clone() const { return new HybridImageSet(*this); }
    virtual HybridSpace space() const { return HybridSpace(*this); }
    virtual ImageSet& operator[](DiscreteLocation q) {
        return this->std::map<DiscreteLocation,ImageSet>::operator[](q); }
    virtual ImageSet const& operator[](DiscreteLocation q) const {
        ARIADNE_ASSERT(this->find(q)!=this->locations_end());
        return this->find(q)->second; }
    virtual tribool overlaps(const HybridBox& hbx) const {
        locations_const_iterator loc_iter=this->find(hbx.first);
        if(loc_iter==this->locations_end()) { return false; }
        return loc_iter->second.overlaps(hbx.second); }
    virtual tribool disjoint(const HybridBox& hbx) const {
        locations_const_iterator loc_iter=this->find(hbx.first);
        if(loc_iter==this->locations_end()) { return true; }
        return loc_iter->second.disjoint(hbx.second); }
    virtual tribool inside(const HybridBoxes& hbx) const  {
        tribool result = true;
        for(locations_const_iterator loc_iter=this->locations_begin(); loc_iter!=this->locations_end(); ++loc_iter) {
            if(!loc_iter->second.empty()) {
                HybridBoxes::const_iterator hbx_loc_iter=hbx.find(loc_iter->first);
                if(hbx_loc_iter==hbx.end()) { result=false; }
                else { result = result && loc_iter->second.inside(hbx_loc_iter->second); }
                if(result==false) { return result; }
            }
        }
        return result; }
    virtual HybridBoxes bounding_box() const {
        HybridBoxes result;
        for(locations_const_iterator loc_iter=this->locations_begin(); loc_iter!=this->locations_end(); ++loc_iter) {
            if(!loc_iter->second.empty()) { result.insert(std::make_pair(loc_iter->first,loc_iter->second.bounding_box())); } }
        return result; }
    virtual std::ostream& write(std::ostream& os) const { return os << "HybridImageSet(...)"; }
};


//! A set comprising a ConstraintSet in each location.
class HybridConstraintSet
    : public std::map<DiscreteLocation,ConstraintSet>
{
};

//! A set comprising a %ListSet in each location.
template<class ES>
class HybridListSet
    : public std::map<DiscreteLocation,ListSet<ES> >
{
  public:
    typedef typename std::map< DiscreteLocation,ListSet<ES> >::iterator locations_iterator;
    typedef typename std::map< DiscreteLocation,ListSet<ES> >::const_iterator locations_const_iterator;
    typedef HybridSetConstIterator< ListSet<ES>, std::pair<DiscreteLocation,ES> > const_iterator;

    HybridListSet() { }
    HybridListSet(const std::pair<DiscreteLocation,ES>& hes) { this->adjoin(hes); }

    locations_iterator locations_begin() {
        return this->std::map<DiscreteLocation,ListSet<ES> >::begin(); }
    locations_iterator locations_end() {
        return this->std::map<DiscreteLocation,ListSet<ES> >::end(); }
    locations_const_iterator locations_begin() const {
        return this->std::map<DiscreteLocation,ListSet<ES> >::begin(); }
    locations_const_iterator locations_end() const {
        return this->std::map<DiscreteLocation,ListSet<ES> >::end(); }
    const_iterator begin() const {
        return const_iterator(*this,false); }
    const_iterator end() const {
        return const_iterator(*this,true); }

    /*! \brief Returns the number of basic hybrid sets forming this object. */
    size_t size() const {
        size_t s = 0;
        for(locations_const_iterator loc_iter=this->locations_begin();
            loc_iter!=this->locations_end(); ++loc_iter) {
            s += loc_iter->second.size();
        }
        return s;
    }

    //using std::map< DiscreteLocation,ListSet<ES> >::insert;

    //using std::map<DiscreteLocation,ListSet<ES> >::operator[];
    ListSet<ES>& operator[](const DiscreteLocation& q) {
        return this->std::map<DiscreteLocation,ListSet<ES> >::operator[](q); }
    const ListSet<ES>& operator[](const DiscreteLocation& q) const {
        ARIADNE_ASSERT_MSG(this->find(q)!=this->locations_end(),(*this)<<" has no location "<<q);
        return this->find(q)->second; }

    void adjoin(const DiscreteLocation& q, const ES& es) {
        (*this)[q].adjoin(es); }
    void adjoin(const std::pair<DiscreteLocation,ES>& hes) {
        (*this)[hes.first].adjoin(hes.second); }
    void adjoin(const HybridListSet<ES>& hls) {
        for(locations_const_iterator loc_iter=hls.locations_begin();
            loc_iter!=hls.locations_end(); ++loc_iter) {
            (*this)[loc_iter->first].adjoin(loc_iter->second); } }
    void adjoin(const ListSet<std::pair<DiscreteLocation,ES> >& hls) {
        for(locations_const_iterator loc_iter=hls.locations_begin();
            loc_iter!=hls.locations_end(); ++loc_iter) {
            (*this)[loc_iter->first].adjoin(loc_iter->second); } }
    void adjoin(const HybridBasicSet<ES>& hbs) {
        (*this)[hbs.location()].adjoin(hbs.continuous_state_set()); }

    HybridListSet<Box> bounding_boxes() const {
        HybridListSet<Box> result;
        for(locations_const_iterator loc_iter=this->locations_begin();
            loc_iter!=this->locations_end(); ++loc_iter) {
            result[loc_iter->first]=loc_iter->second.bounding_boxes(); }
        return result; }

    HybridSpace space() const { return HybridSpace(*this); }
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

class HybridAutomaton;

class HybridGrid
{
    Map<DiscreteLocation,Grid> _grids;
  public:
    HybridGrid();
    HybridGrid(const HybridSpace& hs, double l);
    void insert(DiscreteLocation q, const Grid& g);
    Grid operator[](const DiscreteLocation& loc) const;
    Grid& operator[](const DiscreteLocation& loc);

    bool has_location(const DiscreteLocation& q) const;

    friend inline std::ostream&
    operator<<(std::ostream& os, const HybridGrid& hgrid) {
        return os << "HybridGrid(" << hgrid._grids << ")"; }

};




class HybridGridCell
    : public std::pair<DiscreteLocation,GridCell>
{
  public:
    HybridGridCell()
        : std::pair<DiscreteLocation,GridCell>() { }
    HybridGridCell(DiscreteLocation q,const GridCell& gc)
        : std::pair<DiscreteLocation,GridCell>(q,gc) { }
    HybridGridCell(const std::pair<DiscreteLocation,GridCell>& hgc)
        : std::pair<DiscreteLocation,GridCell>(hgc) { }
    HybridGridCell(const std::pair<const DiscreteLocation,GridCell>& hgc)
        : std::pair<DiscreteLocation,GridCell>(hgc.first,hgc.second) { }
    HybridBox box() const { return HybridBox(this->first,this->second.box()); }
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


//! A set comprising a %GridTreeSet in each location.
class HybridGridTreeSet
    : public std::map<DiscreteLocation,GridTreeSet>
{
    HybridGrid _hgrid;
  public:
    typedef std::map<DiscreteLocation,GridTreeSet>::iterator locations_iterator;
    typedef std::map<DiscreteLocation,GridTreeSet>::const_iterator locations_const_iterator;
    typedef HybridSetConstIterator<GridTreeSet,HybridGridCell> const_iterator;
  public:
    locations_iterator locations_begin() {
        return this->std::map<DiscreteLocation,GridTreeSet>::begin(); }
    locations_iterator locations_end() {
        return this->std::map<DiscreteLocation,GridTreeSet>::end(); }
    locations_const_iterator locations_begin() const {
        return this->std::map<DiscreteLocation,GridTreeSet>::begin(); }
    locations_const_iterator locations_end() const {
        return this->std::map<DiscreteLocation,GridTreeSet>::end(); }
    const_iterator begin() const {
        return const_iterator(*this,false); }
    const_iterator  end() const {
        return const_iterator(*this,true); }
  public:
    HybridGridTreeSet() { }
/*
    HybridGridTreeSet(const HybridSpace& hspace) {
        for(HybridSpace::locations_const_iterator loc_iter = hspace.
                locations_begin(); loc_iter!=hspace.locations_end(); ++loc_iter) {
            this->insert(make_pair(loc_iter->first,Grid(loc_iter->second))); } }
*/
    HybridGridTreeSet(const HybridGrid& hgrid) : _hgrid(hgrid) { }

    HybridGrid grid() const { return this->_hgrid; }

    bool has_location(DiscreteLocation q) const {
        return this->find(q)!=this->std::map<DiscreteLocation,GridTreeSet>::end(); }

    void adjoin(DiscreteLocation q, const GridCell& c) {
        this->_provide_location(q).adjoin(c); }

    void adjoin(const HybridGridCell& hgc) {
        this->_provide_location(hgc.first).adjoin(hgc.second); }

    void adjoin(const HybridGridTreeSet& hgts) {
        for(HybridGridTreeSet::locations_const_iterator loc_iter=hgts.locations_begin(); loc_iter!=hgts.locations_end(); ++loc_iter) {
            this->_provide_location(loc_iter->first).adjoin(loc_iter->second); } }

    void remove(const HybridGridTreeSet& hgts) {
        for(HybridGridTreeSet::locations_const_iterator loc_iter=hgts.locations_begin(); loc_iter!=hgts.locations_end(); ++loc_iter) {
            if(this->has_location(loc_iter->first)) {
                this->find(loc_iter->first)->second.remove(loc_iter->second); } } }

    void restrict(const HybridGridTreeSet& hgts) {
        for(HybridGridTreeSet::locations_const_iterator loc_iter=hgts.locations_begin(); loc_iter!=hgts.locations_end(); ++loc_iter) {
            if(this->has_location(loc_iter->first)) {
                this->find(loc_iter->first)->second.restrict(loc_iter->second); } } }

    void restrict_to_height(uint h) {
        for(locations_iterator loc_iter=locations_begin(); loc_iter!=locations_end(); ++loc_iter) {
            loc_iter->second.restrict_to_height(h); } }

    void adjoin_inner_approximation(const HybridBoxes& hbxs, const int depth) {
        for(HybridBoxes::const_iterator loc_iter=hbxs.begin();
                loc_iter!=hbxs.end(); ++loc_iter) {
            DiscreteLocation loc=loc_iter->first;
            this->_provide_location(loc).adjoin_inner_approximation(loc_iter->second,loc_iter->second,depth); } }

    void adjoin_lower_approximation(const HybridOvertSetInterface& hs, const int height, const int depth) {
        HybridSpace hspc=hs.space();
        for(HybridSpace::const_iterator loc_iter=hspc.begin();
                loc_iter!=hspc.end(); ++loc_iter) {
            DiscreteLocation loc=loc_iter->first;
            this->_provide_location(loc).adjoin_lower_approximation(hs[loc],height,depth); } }

    void adjoin_outer_approximation(const HybridCompactSetInterface& hs, const int depth) {
        HybridSpace hspc=hs.space();
        for(HybridSpace::const_iterator loc_iter=hspc.begin();
                loc_iter!=hspc.end(); ++loc_iter) {
            this->_provide_location(loc_iter->first).adjoin_outer_approximation(hs[loc_iter->first],depth); } }

    void adjoin_outer_approximation(const HybridBoxes& hbxs, const int depth) {
        for(HybridBoxes::const_iterator loc_iter=hbxs.begin();
                loc_iter!=hbxs.end(); ++loc_iter) {
            this->_provide_location(loc_iter->first).adjoin_outer_approximation(loc_iter->second,depth); } }

    template<class S> void adjoin_outer_approximation(DiscreteLocation q, const S& s) {
        this->_provide_location(q).adjoin_outer_approximation(s); }

    GridTreeSet& operator[](DiscreteLocation q) {
        return this->_provide_location(q);
    }

    const GridTreeSet& operator[](DiscreteLocation q) const {
        ARIADNE_ASSERT(this->has_location(q));
        return this->find(q)->second;
    }

    bool empty() const {
        for(locations_const_iterator loc_iter=this->locations_begin();
            loc_iter!=this->locations_end(); ++loc_iter) {
            if(!loc_iter->second.empty()) { return false; } }
        return true; }

    size_t size() const {
        size_t result=0;
        for(locations_const_iterator loc_iter=this->locations_begin();
            loc_iter!=this->locations_end(); ++loc_iter) {
            result+=loc_iter->second.size(); }
        return result; }

    HybridListSet<Box> boxes() const {
        HybridListSet<Box> result;
        for(const_iterator iter=this->begin();
            iter!=this->end(); ++iter) {
            result[iter->first].adjoin(iter->second.box()); }
        return result; }
    void mince(int depth) {
        for(locations_iterator loc_iter=this->locations_begin();
            loc_iter!=this->locations_end(); ++loc_iter) {
            loc_iter->second.mince(depth); } }
    void recombine() {
        for(locations_iterator loc_iter=this->locations_begin();
            loc_iter!=this->locations_end(); ++loc_iter) {
            loc_iter->second.recombine(); } }
  public:
    // HybridSetInterface methods
    HybridGridTreeSet* clone() const { return new HybridGridTreeSet(*this); }
    HybridSpace space() const { return HybridSpace(*this); }

    bool disjoint(const HybridBox& hbx) const {
        locations_const_iterator loc_iter = this->find( hbx.first );
        return loc_iter != this->locations_end() || loc_iter->second.disjoint( hbx.second );
    }

    bool overlaps(const HybridBox& hbx) const {
        locations_const_iterator loc_iter = this->find( hbx.first );
        return loc_iter != this->locations_end() && loc_iter->second.overlaps( hbx.second );
    }

    bool superset(const HybridBox& hbx) const {
        locations_const_iterator loc_iter=this->find(hbx.first);
        return loc_iter!=this->locations_end() && loc_iter->second.superset( hbx.second );
    }

    bool subset(const HybridBoxes& hbx) const  {
        for( locations_const_iterator loc_iter = this->locations_begin(); loc_iter != this->locations_end(); ++loc_iter ) {
            if( !loc_iter->second.empty() ) {
                HybridBoxes::const_iterator hbx_loc_iter = hbx.find( loc_iter->first );
                if( hbx_loc_iter != hbx.end() && ! loc_iter->second.subset( hbx_loc_iter->second ) ) {
                    return false;
                }
            }
        }
        return true;
    }

    HybridBoxes bounding_box() const {
        HybridBoxes result;
        for( locations_const_iterator loc_iter = this->locations_begin(); loc_iter != this->locations_end(); ++loc_iter ) {
            if( !loc_iter->second.empty() ) {
                result.insert( std::make_pair( loc_iter->first, loc_iter->second.bounding_box() ) );
            }
        }
        return result;
    }

    std::ostream& write(std::ostream& os) const {
        return os << static_cast<const std::map<DiscreteLocation,GridTreeSet>&>(*this);
    }

  private:
    GridTreeSet& _provide_location(const DiscreteLocation& q) {
        std::map<DiscreteLocation,GridTreeSet>::iterator iter=this->find(q);
        if(iter==this->std::map<DiscreteLocation,GridTreeSet>::end()) {
            this->insert(std::make_pair(q,this->_hgrid[q]));
            iter=this->find(q);
        }
        return iter->second; }
    /*
      friend class boost::serialization::access;
      template<class Archive> void serialize(Archive & ar, const unsigned int version) {
      ar & static_cast<std::map<int,GridTreeSet>&>(*this); }
    */
};



template<class DS, class HBS> inline
HybridSetConstIterator<DS,HBS>::
HybridSetConstIterator(const std::map<DiscreteLocation,DS>& map, bool end)
    : loc_begin(map.begin()),
      loc_end(map.end()),
      loc_iter(end?loc_end:loc_begin)
{
    if(loc_iter!=loc_end) {
        bs_iter=loc_iter->second.begin();
        this->increment_loc();
    }
}


template<class DS, class HBS> inline
bool
HybridSetConstIterator<DS,HBS>::equal(const HybridSetConstIterator<DS,HBS>& other) const
{
    return this->loc_iter==other.loc_iter && (this->loc_iter==this->loc_end || this->bs_iter==other.bs_iter);
}


template<class DS, class HBS> inline
HBS const&
HybridSetConstIterator<DS,HBS>::dereference() const
{
    this->hybrid_set=HBS(loc_iter->first,*this->bs_iter);
    return this->hybrid_set;
}


template<class DS, class HBS> inline
void
HybridSetConstIterator<DS,HBS>::increment()
{
    ++this->bs_iter;
    this->increment_loc();
}

template<class DS, class HBS> inline
void
HybridSetConstIterator<DS,HBS>::increment_loc()
{
    while(bs_iter==loc_iter->second.end()) {
        ++loc_iter;
        if(loc_iter==loc_end) { return; }
        bs_iter=loc_iter->second.begin();
    }
}


template<class A> void serialize(A& archive, HybridGridTreeSet& set, const unsigned int version) {
    archive & static_cast<std::map<DiscreteLocation,GridTreeSet>&>(set); }

template<class BS> inline
void
draw(FigureInterface& figure, const HybridBasicSet<BS>& hs) {
    draw(figure,hs.continuous_state_set());
}



template<class DS> inline
void
draw(FigureInterface& figure, const std::map<DiscreteLocation,DS>& hds) {
    for(typename std::map<DiscreteLocation,DS>::const_iterator loc_iter=hds.begin();
        loc_iter!=hds.end(); ++loc_iter) {
        draw(figure,loc_iter->second);
        //figure.draw(loc_iter->second);
        //figure << loc_iter->second;
    }
}

template<class BS> inline FigureInterface& operator<<(FigureInterface& figure, const HybridBasicSet<BS>& hs) {
    draw(figure,hs); return figure;
}

template<class DS> inline FigureInterface& operator<<(FigureInterface& figure, const std::map<DiscreteLocation,DS>& hs) {
    draw(figure,hs); return figure;
}


} // namespace Ariadne

namespace boost { namespace serialization {
template<class A> void serialize(A& archive, const Ariadne::HybridGridTreeSet& set, const uint version);
template<class A> void serialize(A& archive, const Ariadne::DiscreteLocation& state, const uint version);
}}

#endif // ARIADNE_HYBRID_SET_H
