/***************************************************************************
 *            hybrid_set.h
 *
 *  Copyright  2006  Alberto Casagrande,  Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
#ifndef _ARIADNE_HYBRID_SET_H
#define _ARIADNE_HYBRID_SET_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "../declarations.h"
#include "../geometry/grid_set.h"
#include "../system/hybrid_automaton.h"

namespace Ariadne {  
namespace Geometry {

typedef size_type location_type;
extern int verbosity;
  
class HybridSystemError : public std::runtime_error {
 public:
  HybridSystemError(const std::string& what) : std::runtime_error(what) { }
};

/*! \ingroup HybridSet
 *  \brief A base class for representing subsets of hybrid state spaces.
 */
template< class Set >
class HybridSet
{
 public:
  typedef typename Set::real_type real_type;
  typedef typename Set::state_type state_type;
 
 public:
  /*! \brief Construct a set for \a nq discrete modes. */
  HybridSet(location_type nq) 
    : _component_sets(nq) { }
    
  /*! \brief Construct a set for \a nq discrete modes, each based on the set \a s. */
  HybridSet(location_type nq, const Set& s) 
    : _component_sets(nq,s) { }
    
  /*! \brief Construct a set for \a n discrete modes, each based on the same cells in the same grid. */
  template<class S> HybridSet(const HybridSet<S>& hs) : _component_sets() {
    for(location_type q=0; q!=hs.number_of_discrete_locations(); ++q) {
      this->_component_sets.push_back(Set(hs[q]));
    }
  }
  
  /*! \brief The number of discrete locations or components comprising the set. */
  location_type number_of_discrete_locations() const { return _component_sets.size(); }
  /*! \brief A reference to the state set corresponding to discrete location \a q. */
  Set& operator[] (const location_type& q) { return _component_sets[q]; }
  /*! \brief The state set corresponding to discrete location \a q. */
  const Set& operator[] (const location_type& q) const { return _component_sets[q]; }
  /*! \brief Clear all discrete locations. */
  void clear() { 
    for(location_type q=0; q!=this->number_of_discrete_locations(); ++q) { 
      (*this)[q].clear(); } 
  }
  /*! \brief Returns true if every component is empty. */
  tribool empty() const { 
    tribool result=true; 
    for(size_type q=0; q!=this->number_of_discrete_locations(); ++q) {
      result=result && (*this)[q].empty();
      if(!result) { return result; } }
    return result;
  }
  
  /*! \brief Adjoin another hybrid set. */
  template<class S> void adjoin(const HybridSet<S>& hs) {
    if(this->number_of_discrete_locations()!=hs.number_of_discrete_locations()) {
      throw std::runtime_error("Invalid number of discrete components");
    }
    for(location_type q=0; q!=this->number_of_discrete_locations(); ++q) {
      (*this)[q].adjoin(hs[q]);
    }
  }
  
 private:
  std::vector< Set > _component_sets;
};

/*! \ingroup HybridSet
 *  \brief A hybrid set comprising of a GridMaskSet for every component.
 */
template< class R >
class HybridGridMaskSet
  : public HybridSet< GridMaskSet<R> >
{
 public:
  /*! \brief Construct a set for \a n discrete modes, each based on the same cells in the same grid. */
  HybridGridMaskSet(location_type nq, const FiniteGrid<R>& fg) 
    : HybridSet< GridMaskSet<R> >(nq,GridMaskSet<R>(fg)) { }
    
 
  /*! \brief Construct a set for \a n discrete modes, based on a list of finite grids. */
  template<class FGC> HybridGridMaskSet(const FGC& fgs);
};



/*! \ingroup HybridSet
 *  \brief A hybrid set comprising of a GridCellListSet for every component.
 */
template< class R >
class HybridGridCellListSet
  : public HybridSet< GridCellListSet<R> >
{
 public:
  /*! \brief Construct a set for \a n discrete modes, based on a fixed grid. */
  HybridGridCellListSet(location_type nq, const Grid<R>& g) 
    : HybridSet< GridCellListSet<R> >(nq,GridCellListSet<R>(g)) { }
    
  /*! \brief Construct a set for \a n discrete modes, based on a list of finite grids. */
  HybridGridCellListSet(const HybridSet< GridMaskSet<R> >& hgms)
    : HybridSet< GridCellListSet<R> >(hgms) { }
};


template<class R> 
inline
HybridGridCellListSet<R>
regular_intersection(const HybridGridCellListSet<R>& hgcl, const HybridGridMaskSet<R>& hgms) {
  if(verbosity>5) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
  
  if(hgcl.number_of_discrete_locations()!=hgms.number_of_discrete_locations()) {
    throw HybridSystemError("Intersection of sets with different numbers of discrete locations");
  }
  
  HybridGridCellListSet<R> result=hgcl;
  result.clear();
  for(location_type q=0; q!=hgcl.number_of_discrete_locations(); ++q) {
    result[q]=regular_intersection(hgcl[q],hgms[q]);
  }
  return result;
}


template<class R> 
inline
HybridGridCellListSet<R>
regular_intersection(const HybridGridMaskSet<R>& hgms, const HybridGridCellListSet<R>& hgcl) {
  return regular_intersection(hgcl,hgms);
}

   
}
}

#endif /* _ARIADNE_HYBRID_EVOLVER_H */
