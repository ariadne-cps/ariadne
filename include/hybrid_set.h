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

//! A set comprising a GridTreeSet in each location.
template<class ES>
class HybridListSet 
  : public ListSet< std::pair<DiscreteState,ES> >
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




//! A set comprising a GridTreeSet in each location.
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


} // namespace Ariadne

#endif // ARIADNE_HYBRID_SET_H
