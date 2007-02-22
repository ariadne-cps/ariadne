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
#include <sstream>

#include "../declarations.h"
#include "../geometry/grid_set.h"
#include "../system/hybrid_automaton.h"

namespace Ariadne {  
  namespace Geometry {

    typedef size_type location_type;



    class HybridSystemError : public std::runtime_error {
     public:
      HybridSystemError(const std::string& what) : std::runtime_error(what) { }
    };



    /*! \ingroup HybridSet
     *  \brief A base class for representing subsets of hybrid state spaces.
     */
    template< class S >
    class HybridSet
    {
     public:
      typedef typename S::real_type real_type;
      typedef typename S::state_type state_type;
      
     public:
      /*! \brief Construct a set for \a nq discrete modes. */
      HybridSet(location_type nq);
      
      /*! \brief Construct a set for \a nq discrete modes, each based on the set \a s. */
      HybridSet(location_type nq, const S& s);
      
      /*! \brief Construct a set for \a n discrete modes, each based on the same cells in the same grid. */
      template<class S1> HybridSet(const HybridSet<S1>& hs);
      
      /*! \brief The number of discrete locations or components comprising the set. */
      location_type number_of_discrete_locations() const;
      /*! \brief A reference to the state set corresponding to discrete location \a q. */
      S& operator[](const location_type& q);
      /*! \brief The state set corresponding to discrete location \a q. */
      const S& operator[](const location_type& q) const;
      
      /*! \brief Clear all discrete locations. */
      void clear();
      /*! \brief Returns true if every component is empty. */
      tribool empty() const;
      
      /*! \brief Adjoin another hybrid set. */
      template<class S1> void adjoin(const HybridSet<S1>& hs);
      
     private:
      std::vector< S > _component_sets;
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
      HybridGridMaskSet(location_type nq, const FiniteGrid<R>& fg);
      
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
      HybridGridCellListSet(location_type nq, const Grid<R>& g);
      
      /*! \brief Construct a set for \a n discrete modes, based on a list of finite grids. */
      HybridGridCellListSet(const HybridSet< GridMaskSet<R> >& hgms);
    };

    template<class R> 
    HybridGridCellListSet<R>
    regular_intersection(const HybridGridCellListSet<R>& hgcl, const HybridGridMaskSet<R>& hgms);

    template<class R> 
    HybridGridCellListSet<R>
    regular_intersection(const HybridGridMaskSet<R>& hgms, const HybridGridCellListSet<R>& hgcl);
   
  }
}

#include "hybrid_set.inline.h"

#endif /* _ARIADNE_HYBRID_EVOLVER_H */
