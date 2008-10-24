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


#include "macros.h"
#include "stlio.h"
#include "set.h"
#include "list_set.h"
#include "grid_set.h"

#include "hybrid_set_interface.h"


namespace Ariadne {

typedef int DiscreteState;

class HybridSpace 
  : public std::map<DiscreteState,uint>
{
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


// FIXME: This class doesn't work

template< class DS, class HBS=std::pair<DiscreteState,typename DS::value_type> >
class HybridSetIterator
  : public boost::iterator_facade<HybridSetIterator<DS>,
                                  HBS,
                                  boost::forward_traversal_tag,
                                  HBS
                                  >

{
 public:
  HybridSetIterator(const std::map<DiscreteState,DS>&, bool);
  bool equal(const HybridSetIterator<DS>&) const;
  HBS dereference() const;
  void increment();
  const DiscreteState& discrete_state() const;
  const typename DS::value_type& continuous_state_set() const;
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
  virtual HybridSpace space() const { ARIADNE_NOT_IMPLEMENTED; }
  virtual tribool intersects(const HybridBox& hbx) const { 
    locations_const_iterator loc_iter=this->find(hbx.first);
    return loc_iter!=this->locations_end()
      && loc_iter->second.intersects(hbx.second); }
  virtual tribool disjoint(const HybridBox& hbx) const { 
    locations_const_iterator loc_iter=this->find(hbx.first);
    return loc_iter!=this->locations_end()
      || loc_iter->second.disjoint(hbx.second); }
  virtual tribool subset(const HybridBoxes& hbx) const { ARIADNE_NOT_IMPLEMENTED; } 
  virtual HybridBoxes bounding_box() const { ARIADNE_NOT_IMPLEMENTED; } 
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
  typedef HybridSetIterator< ListSet<ES> > const_iterator;
  locations_iterator locations_begin() { 
    return this->std::map<DiscreteState,ListSet<ES> >::begin(); }
  locations_iterator locations_end() { 
    return this->std::map<DiscreteState,ListSet<ES> >::end(); }
  locations_const_iterator locations_begin() const { 
    return this->std::map<DiscreteState,ListSet<ES> >::begin(); }
  locations_const_iterator locations_end() const { 
    return this->std::map<DiscreteState,ListSet<ES> >::end(); }
  const_iterator begin() const { 
    return HybridSetIterator< ListSet<ES> >(*this,false); }
  const_iterator end() const { 
    return HybridSetIterator< ListSet<ES> >(*this,true); }

  using std::map< DiscreteState,ListSet<ES> >::insert;

  void adjoin(const DiscreteState& q, const ES& es) {
    (*this)[q].adjoin(es); }
  void adjoin(const std::pair<DiscreteState,ES>& hes) {
    (*this)[hes.first].adjoin(hes.second); }
  void adjoin(const HybridListSet<ES>& hls) {
    for(locations_const_iterator loc_iter=hls.locations_begin();
        loc_iter!=hls.locations_end(); ++loc_iter)
    {  
      (*this)[loc_iter->first].adjoin(loc_iter->second); 
    }
  }
};



template<class ES>
class ListSet< std::pair<DiscreteState,ES> >
  : public HybridListSet<ES>
{
};


class HybridGrid
  : public std::map<DiscreteState,Grid>
{
 public:
  HybridGrid(const HybridSpace& hspc, const Float& l) { ARIADNE_NOT_IMPLEMENTED; }
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

class HybridGridCellListSet
  : public std::vector<HybridGridCell>
{
 public:
  HybridGridCellListSet() { ARIADNE_NOT_IMPLEMENTED; };
  HybridGridCellListSet(const HybridGridTreeSet&) { ARIADNE_NOT_IMPLEMENTED; };
  void adjoin(const HybridGridCell&) { ARIADNE_NOT_IMPLEMENTED; };
  void adjoin(const HybridGridCellListSet&) { ARIADNE_NOT_IMPLEMENTED; };
  void remove(const HybridGridCellListSet&) { ARIADNE_NOT_IMPLEMENTED; };
  void restrict(const HybridBoxes&) { ARIADNE_NOT_IMPLEMENTED; };
  void unique_sort() { ARIADNE_NOT_IMPLEMENTED; };
};




//! A set comprising a %GridTreeSet in each location.
class HybridGridTreeSet 
  : public virtual HybridSetInterface,
    public std::map<DiscreteState,GridTreeSet>
{
  class Iterator { 
   public:
    const HybridGridCell& operator*() const { ARIADNE_NOT_IMPLEMENTED; }
    const HybridGridCell* operator->() const { ARIADNE_NOT_IMPLEMENTED; }
    void operator++() { ARIADNE_NOT_IMPLEMENTED; }
    bool operator==(const Iterator&) const { ARIADNE_NOT_IMPLEMENTED; }
    bool operator!=(const Iterator&) const { ARIADNE_NOT_IMPLEMENTED; }
  };
 public:
  typedef Iterator iterator;
  typedef Iterator const_iterator;
 public:
    Iterator begin() const { ARIADNE_NOT_IMPLEMENTED; }
    Iterator end() const { ARIADNE_NOT_IMPLEMENTED; }
    Iterator begin() { ARIADNE_NOT_IMPLEMENTED; }
    Iterator end() { ARIADNE_NOT_IMPLEMENTED; }
 public:
  HybridGridTreeSet() { ARIADNE_NOT_IMPLEMENTED; }
  HybridGridTreeSet(const HybridSpace& hspace) { ARIADNE_NOT_IMPLEMENTED; }
  HybridGridTreeSet(const HybridGrid& hgrid) { ARIADNE_NOT_IMPLEMENTED; }
  HybridGridTreeSet(const HybridGridCellListSet& hgcls) { ARIADNE_NOT_IMPLEMENTED; }

  HybridGrid grid() const { ARIADNE_NOT_IMPLEMENTED; }
  
  bool has_location(DiscreteState q) const {
    return this->find(q)!=this->std::map<DiscreteState,GridTreeSet>::end(); }

  void adjoin(DiscreteState q, const GridCell& c) {
    this->find(q)->second.adjoin(c); }

  void adjoin(const HybridGridCell& hgc) {
    this->find(hgc.first)->second.adjoin(hgc.second); }

  void adjoin(const HybridGridCellListSet& hgcls) {
    for(HybridGridCellListSet::const_iterator iter=hgcls.begin(); iter!=hgcls.end(); ++iter) {
      this->find(iter->first)->second.adjoin(iter->second); } }

  void adjoin(const HybridGridTreeSet& hgts) {
    for(HybridGridTreeSet::const_iterator iter=hgts.begin(); iter!=hgts.end(); ++iter) {
      this->find(iter->first)->second.adjoin(iter->second); } }

  void adjoin_lower_approximation(const HybridOvertSetInterface& hs) { ARIADNE_NOT_IMPLEMENTED; };
  void adjoin_outer_approximation(const HybridCompactSetInterface& hs) { ARIADNE_NOT_IMPLEMENTED; };

  template<class S> void adjoin_outer_approximation(DiscreteState q, const S& s) {
    this->operator[](q).adjoin_outer_approximation(s); }

  const GridTreeSet& operator[](DiscreteState q) const {
    ARIADNE_ASSERT(this->has_location(q));
    return const_cast<HybridGridTreeSet*>(this)->operator[](q);
  }

  HybridGridCellListSet cells() const { ARIADNE_NOT_IMPLEMENTED; }
  HybridListSet<Box> boxes() const { ARIADNE_NOT_IMPLEMENTED; }
 public:
  // HybridSetInterface methods
  HybridGridTreeSet* clone() const { return new HybridGridTreeSet(*this); }
  HybridSpace space() const { ARIADNE_NOT_IMPLEMENTED; }
  tribool disjoint(const HybridBox& hbx) const { ARIADNE_NOT_IMPLEMENTED; }
  tribool intersects(const HybridBox& hbx) const { ARIADNE_NOT_IMPLEMENTED; }
  tribool superset(const HybridBox& hbx) const { ARIADNE_NOT_IMPLEMENTED; }
  tribool subset(const HybridBoxes& hbxs) const { ARIADNE_NOT_IMPLEMENTED; }
  HybridBoxes bounding_box() const { ARIADNE_NOT_IMPLEMENTED; }
  std::ostream& write(std::ostream& os) const { 
    return os << static_cast<const std::map<DiscreteState,GridTreeSet>&>(*this); }
};




template<class DS, class HBS> inline
HybridSetIterator<DS,HBS>::HybridSetIterator(const std::map<DiscreteState,DS>& map, bool end)
  : loc_begin(map.begin()),
    loc_end(map.end()),
    loc_iter(end?loc_end:loc_begin),
    bs_iter()
{
  if(loc_iter!=loc_end) {
    bs_iter=loc_iter->second.begin();
    this->increment_loc();
  }
}


template<class DS, class HBS> inline
bool
HybridSetIterator<DS,HBS>::equal(const HybridSetIterator<DS>& other) const
{
  return this->loc_iter==other.loc_iter && (this->loc_iter==this->loc_end || this->bs_iter==other.bs_iter);
}


template<class DS, class HBS> inline
HBS
HybridSetIterator<DS,HBS>::dereference() const
{
  return HBS(loc_iter->first,*this->bs_iter);
}


template<class DS, class HBS> inline
void
HybridSetIterator<DS,HBS>::increment() 
{
  ++this->bs_iter;
  this->increment_loc();
}

template<class DS, class HBS> inline
void
HybridSetIterator<DS,HBS>::increment_loc() 
{
  while(bs_iter==loc_iter->second.end()) {
    ++loc_iter;
    if(loc_iter==loc_end) { return; } 
    bs_iter=loc_iter->second.begin();
  }
}

template<class DS, class HBS> inline
const DiscreteState&
HybridSetIterator<DS,HBS>::discrete_state() const 
{ 
  return loc_iter->first; 
}

template<class DS, class HBS> inline
const typename DS::value_type&
HybridSetIterator<DS,HBS>::continuous_state_set() const { 
  return *bs_iter; 
}

template<class G, class BS> inline 
void 
plot(G& graphic, const std::pair<int,BS>& hs) { 
  plot(graphic,hs.second); 
}


} // namespace Ariadne

#endif // ARIADNE_HYBRID_SET_H
