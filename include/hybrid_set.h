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
#include "function_set.h"
#include "list_set.h"
#include "grid_set.h"

#include "hybrid_set_interface.h"

#include "serialization.h"

namespace Ariadne {

typedef int DiscreteState;

class HybridSpace 
  : public std::map<DiscreteState,uint>
{
 public:
  typedef std::map<DiscreteState,uint>::const_iterator 
    locations_const_iterator;
    
  HybridSpace() : std::map<DiscreteState,uint>() { }
  template<class HSET> HybridSpace(const HSET& set) {
    for(typename HSET::locations_const_iterator loc_iter
          =set.locations_begin(); loc_iter!=set.locations_end(); ++loc_iter) {
      this->insert(make_pair(loc_iter->first,loc_iter->second.dimension())); }
  }

  locations_const_iterator locations_begin() const { 
    return this->std::map<DiscreteState,uint>::begin(); }
  locations_const_iterator locations_end() const { 
    return this->std::map<DiscreteState,uint>::end(); }
};

inline
HybridBoxes 
bounding_boxes(const std::map<DiscreteState,uint> space, Interval bound)
{
  HybridBoxes result; 
  for(std::map<DiscreteState,uint>::const_iterator loc_iter=space.begin();
      loc_iter!=space.end(); ++loc_iter)
  {
    result.insert(make_pair(loc_iter->first,Box(loc_iter->second, bound)));
  }
  return result;
}


template< class DS, class HBS >
class HybridSetConstIterator
  : public boost::iterator_facade<HybridSetConstIterator<DS,HBS>,
                                  HBS,
                                  boost::forward_traversal_tag,
                                  HBS
                                  >

{
 public:
  HybridSetConstIterator(const std::map<DiscreteState,DS>&, bool);
  bool equal(const HybridSetConstIterator<DS,HBS>&) const;
  HBS dereference() const;
  void increment();
 private:
  void increment_loc();
 private:
  typename std::map< DiscreteState,DS>::const_iterator loc_begin;
  typename std::map< DiscreteState,DS>::const_iterator loc_end;
  typename std::map< DiscreteState,DS>::const_iterator loc_iter;
  typename DS::const_iterator bs_iter;
};


//! A set comprising a %ListSet in each location.
class HybridImageSet
  : public std::map<DiscreteState,ImageSet> 
  , public HybridLocatedSetInterface
{
 public:
  typedef std::map<DiscreteState,ImageSet>::iterator locations_iterator;
  typedef std::map<DiscreteState,ImageSet>::const_iterator locations_const_iterator;
  locations_iterator locations_begin() { 
    return this->std::map<DiscreteState,ImageSet>::begin(); }
  locations_iterator locations_end() { 
    return this->std::map<DiscreteState,ImageSet>::end(); }
  locations_const_iterator locations_begin() const { 
    return this->std::map<DiscreteState,ImageSet>::begin(); }
  locations_const_iterator locations_end() const { 
    return this->std::map<DiscreteState,ImageSet>::end(); }

  using std::map<DiscreteState,ImageSet>::insert;

  virtual HybridImageSet* clone() const { return new HybridImageSet(*this); }
  virtual HybridSpace space() const { return HybridSpace(*this); }
  virtual ImageSet& operator[](DiscreteState q) {
    return this->std::map<DiscreteState,ImageSet>::operator[](q); }
  virtual ImageSet const& operator[](DiscreteState q) const {
    ARIADNE_ASSERT(this->find(q)!=this->locations_end());
    return this->find(q)->second; }
  virtual tribool intersects(const HybridBox& hbx) const { 
    locations_const_iterator loc_iter=this->find(hbx.first);
    return loc_iter!=this->locations_end()
      && loc_iter->second.intersects(hbx.second); }
  virtual tribool disjoint(const HybridBox& hbx) const { 
    locations_const_iterator loc_iter=this->find(hbx.first);
    return loc_iter!=this->locations_end()
      || loc_iter->second.disjoint(hbx.second); }
  virtual tribool subset(const HybridBoxes& hbx) const  { 
    for(locations_const_iterator loc_iter=this->begin(); loc_iter!=this->locations_end(); ++loc_iter) {
      if(!loc_iter->second.empty()) { 
        HybridBoxes::const_iterator hbx_loc_iter=hbx.find(loc_iter->first); 
        if(hbx_loc_iter!=hbx.end() && !loc_iter->second.subset(hbx_loc_iter->second)) { return false; }
      } } return true; }
  virtual HybridBoxes bounding_box() const {  
    HybridBoxes result;
    for(locations_const_iterator loc_iter=this->begin(); loc_iter!=this->locations_end(); ++loc_iter) {
      if(!loc_iter->second.empty()) { result.insert(std::make_pair(loc_iter->first,loc_iter->second.bounding_box())); } } 
    return result; }
  virtual std::ostream& write(std::ostream& os) const { return os << "HybridImageSet(...)"; }
};


//! A set comprising a %ListSet in each location.
template<class ES>
class HybridListSet 
  : public std::map<DiscreteState,ListSet<ES> > 
{
 public:
  typedef typename std::map< DiscreteState,ListSet<ES> >::iterator locations_iterator;
  typedef typename std::map< DiscreteState,ListSet<ES> >::const_iterator locations_const_iterator;
  typedef HybridSetConstIterator< ListSet<ES>, std::pair<DiscreteState,ES> > const_iterator;

  HybridListSet() { }
  HybridListSet(const std::pair<DiscreteState,ES>& hes) { this->adjoin(hes); }

  locations_iterator locations_begin() { 
    return this->std::map<DiscreteState,ListSet<ES> >::begin(); }
  locations_iterator locations_end() { 
    return this->std::map<DiscreteState,ListSet<ES> >::end(); }
  locations_const_iterator locations_begin() const { 
    return this->std::map<DiscreteState,ListSet<ES> >::begin(); }
  locations_const_iterator locations_end() const { 
    return this->std::map<DiscreteState,ListSet<ES> >::end(); }
  const_iterator begin() const { 
    return const_iterator(*this,false); }
  const_iterator end() const { 
    return const_iterator(*this,true); }

  using std::map< DiscreteState,ListSet<ES> >::insert;

  //using std::map<DiscreteState,ListSet<ES> >::operator[];
  ListSet<ES>& operator[](const DiscreteState& q) {
    return this->std::map<DiscreteState,ListSet<ES> >::operator[](q); }
  const ListSet<ES>& operator[](const DiscreteState& q) const {
    ARIADNE_ASSERT_MSG(this->find(q)!=this->locations_end(),(*this)<<" has no location "<<q);
    return this->find(q)->second; }
  
  void adjoin(const DiscreteState& q, const ES& es) {
    (*this)[q].adjoin(es); }
  void adjoin(const std::pair<DiscreteState,ES>& hes) {
    (*this)[hes.first].adjoin(hes.second); }
  void adjoin(const HybridListSet<ES>& hls) {
    for(locations_const_iterator loc_iter=hls.locations_begin();
        loc_iter!=hls.locations_end(); ++loc_iter) {  
      (*this)[loc_iter->first].adjoin(loc_iter->second); } }

  HybridSpace space() const { return HybridSpace(*this); }
};



template<class ES>
class ListSet< std::pair<DiscreteState,ES> >
  : public HybridListSet<ES>
{
 public:
  ListSet() { }
  ListSet(const std::pair<DiscreteState,ES>& hes) { this->adjoin(hes); }
};

template<class ES>
std::ostream& 
operator<<(std::ostream& os, 
           const ListSet< std::pair<DiscreteState,ES> >& ls) {
  const std::map< DiscreteState, ListSet<ES> >& hls=ls;
  return os << "HybridListSet" << hls; 
}
    

class HybridGrid
  : public std::map<DiscreteState,Grid>
{
 public:
  typedef std::map<DiscreteState,Grid>::const_iterator 
    locations_const_iterator;
  HybridGrid(const HybridSpace& hspc, const Float l=1.0) { 
    for(HybridSpace::locations_const_iterator loc_iter=hspc.begin();
        loc_iter!=hspc.end(); ++loc_iter) {
      this->insert(make_pair(loc_iter->first,Grid(loc_iter->second,l))); } }
  template<class HGSET> HybridGrid(const HGSET& set) {
    for(typename HGSET::locations_const_iterator loc_iter=set.
          locations_begin(); loc_iter!=set.locations_end(); ++loc_iter) {
      this->insert(make_pair(loc_iter->first,loc_iter->second.grid())); }
  }

  locations_const_iterator locations_begin() const { 
    return this->std::map<DiscreteState,Grid>::begin(); }
  locations_const_iterator locations_end() const { 
    return this->std::map<DiscreteState,Grid>::end(); }
  const Grid& operator[](DiscreteState q) const {
    ARIADNE_ASSERT(this->find(q)!=this->locations_end());
    return this->find(q)->second; }
};


class HybridGridCell 
  : public std::pair<DiscreteState,GridCell> 
{
 public:
  HybridGridCell(DiscreteState q,const GridCell& gc)
    : std::pair<DiscreteState,GridCell>(q,gc) { }
  HybridGridCell(const std::pair<DiscreteState,GridCell>& hgc) 
    : std::pair<DiscreteState,GridCell>(hgc) { }
  HybridGridCell(const std::pair<const DiscreteState,GridCell>& hgc) 
    : std::pair<DiscreteState,GridCell>(hgc.first,hgc.second) { }
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
  : public std::map<DiscreteState,GridTreeSet>
{
 public:
  typedef std::map<DiscreteState,GridTreeSet>::iterator locations_iterator;
  typedef std::map<DiscreteState,GridTreeSet>::const_iterator locations_const_iterator;
  typedef HybridSetConstIterator<GridTreeSet,HybridGridCell> const_iterator;
 public:
  locations_iterator locations_begin() { 
    return this->std::map<DiscreteState,GridTreeSet>::begin(); }
  locations_iterator locations_end() { 
    return this->std::map<DiscreteState,GridTreeSet>::end(); }
  locations_const_iterator locations_begin() const { 
    return this->std::map<DiscreteState,GridTreeSet>::begin(); }
  locations_const_iterator locations_end() const { 
    return this->std::map<DiscreteState,GridTreeSet>::end(); }
  const_iterator begin() const { 
    return const_iterator(*this,false); }
  const_iterator  end() const { 
    return const_iterator(*this,true); }
 public:
  HybridGridTreeSet() { }
  HybridGridTreeSet(const HybridSpace& hspace) {
    for(HybridSpace::locations_const_iterator loc_iter = hspace.
          locations_begin(); loc_iter!=hspace.locations_end(); ++loc_iter) {
      this->insert(make_pair(loc_iter->first,Grid(loc_iter->second))); } }
  HybridGridTreeSet(const HybridGrid& hgrid) { 
    for(HybridGrid::locations_const_iterator loc_iter = hgrid.
          locations_begin(); loc_iter!=hgrid.locations_end(); ++loc_iter) {
      this->insert(make_pair(loc_iter->first,Grid(loc_iter->second))); } }
  HybridGridTreeSet(const HybridGridCell& hgc) { 
    this->adjoin(hgc); }

  HybridGrid grid() const { return HybridGrid(*this); }
  
  bool has_location(DiscreteState q) const {
    return this->find(q)!=this->std::map<DiscreteState,GridTreeSet>::end(); }

  void adjoin(DiscreteState q, const GridCell& c) {
    this->find(q)->second.adjoin(c); }

  void adjoin(const HybridGridCell& hgc) {
    this->find(hgc.first)->second.adjoin(hgc.second); }

  void adjoin(const HybridGridTreeSet& hgts) {
    for(HybridGridTreeSet::locations_const_iterator loc_iter=hgts.locations_begin(); loc_iter!=hgts.locations_end(); ++loc_iter) {
      if(!this->has_location(loc_iter->first)) {
        GridTreeSet new_gts(loc_iter->second.grid());
        this->insert(make_pair(loc_iter->first,new_gts)); }
      this->find(loc_iter->first)->second.adjoin(loc_iter->second); } }

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
        loc_iter!=hbxs.end(); ++loc_iter)
    {
      DiscreteState loc=loc_iter->first;
      if(!this->has_location(loc)) {
        this->insert(make_pair(loc,GridTreeSet(loc_iter->second.dimension()))); }
      (*this)[loc].adjoin_inner_approximation(loc_iter->second,depth); } }

  void adjoin_lower_approximation(const HybridOvertSetInterface& hs, const int depth) { 
    HybridSpace hspc=hs.space();
    for(HybridSpace::const_iterator loc_iter=hspc.begin();
        loc_iter!=hspc.end(); ++loc_iter)
    {
      DiscreteState loc=loc_iter->first;
      if(!this->has_location(loc)) {
        this->insert(make_pair(loc,GridTreeSet(loc_iter->second))); }
      (*this)[loc].adjoin_lower_approximation(hs[loc],depth); } }

  void adjoin_outer_approximation(const HybridCompactSetInterface& hs, const int depth) { 
    HybridSpace hspc=hs.space();
    for(HybridSpace::const_iterator loc_iter=hspc.begin();
        loc_iter!=hspc.end(); ++loc_iter)
    {
      DiscreteState loc=loc_iter->first;
      if(!this->has_location(loc)) {
        this->insert(make_pair(loc,GridTreeSet(loc_iter->second))); }
      (*this)[loc].adjoin_outer_approximation(hs[loc],depth); } }

  void adjoin_outer_approximation(const HybridBoxes& hbxs, const int depth) { 
    for(HybridBoxes::const_iterator loc_iter=hbxs.begin();
        loc_iter!=hbxs.end(); ++loc_iter)
    {
      DiscreteState loc=loc_iter->first;
      if(!this->has_location(loc)) {
        this->insert(make_pair(loc,GridTreeSet(loc_iter->second.dimension()))); }
      (*this)[loc].adjoin_outer_approximation(loc_iter->second,depth); } }

  template<class S> void adjoin_outer_approximation(DiscreteState q, const S& s) {
    this->operator[](q).adjoin_outer_approximation(s); }

  GridTreeSet& operator[](DiscreteState q) {
    ARIADNE_ASSERT(this->has_location(q));
    return this->find(q)->second;
  }

  const GridTreeSet& operator[](DiscreteState q) const {
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
  tribool disjoint(const HybridBox& hbx) const { 
    locations_const_iterator loc_iter=this->find(hbx.first);
    return loc_iter!=this->locations_end()
      || Ariadne::disjoint(loc_iter->second,hbx.second); }
  tribool intersects(const HybridBox& hbx) const { 
    locations_const_iterator loc_iter=this->find(hbx.first);
    return loc_iter!=this->locations_end()
      && Ariadne::intersects(loc_iter->second,hbx.second); }
  tribool superset(const HybridBox& hbx) const { 
    locations_const_iterator loc_iter=this->find(hbx.first);
    return loc_iter!=this->locations_end()
      && Ariadne::superset(loc_iter->second,hbx.second); }
  tribool subset(const HybridBoxes& hbx) const  { 
    for(locations_const_iterator loc_iter=this->locations_begin(); loc_iter!=this->locations_end(); ++loc_iter) {
      if(!loc_iter->second.empty()) { 
        HybridBoxes::const_iterator hbx_loc_iter=hbx.find(loc_iter->first); 
        if(hbx_loc_iter!=hbx.end() && !Ariadne::subset(loc_iter->second,hbx_loc_iter->second)) { return false; }
      } } return true; }
  HybridBoxes bounding_box() const {  
    HybridBoxes result;
    for(locations_const_iterator loc_iter=this->locations_begin(); loc_iter!=this->locations_end(); ++loc_iter) {
      if(!loc_iter->second.empty()) { result.insert(std::make_pair(loc_iter->first,Ariadne::bounding_box(loc_iter->second))); } } 
    return result; }
  std::ostream& write(std::ostream& os) const { 
    return os << static_cast<const std::map<DiscreteState,GridTreeSet>&>(*this); }

 private:
  /*
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version) {
      ar & static_cast<std::map<int,GridTreeSet>&>(*this); }
  */
};


template<class DS, class HBS> inline
HybridSetConstIterator<DS,HBS>::
HybridSetConstIterator(const std::map<DiscreteState,DS>& map, bool end)
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
HBS
HybridSetConstIterator<DS,HBS>::dereference() const
{
  return HBS(loc_iter->first,*this->bs_iter);
  //this->hybrid_cell=HBS(loc_iter->first,*this->bs_iter);
  //return this->hybrid_cell;
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
  archive & static_cast<std::map<int,GridTreeSet>&>(set); }


template<class G, class BS> inline 
void 
draw(G& graphic, const std::pair<int,BS>& hs) { 
  draw(graphic,hs.second); 
}

template<class G, class DS> inline 
void 
draw(G& graphic, const std::map<int,DS>& hds) { 
  for(typename std::map<int,DS>::const_iterator loc_iter=hds.begin();
      loc_iter!=hds.end(); ++loc_iter) {
    draw(graphic,loc_iter->second); 
  }
}


} // namespace Ariadne

namespace boost { namespace serialization {
template<class A> void serialize(A& archive, const Ariadne::HybridGridTreeSet& set, const uint version);
}}

#endif // ARIADNE_HYBRID_SET_H
