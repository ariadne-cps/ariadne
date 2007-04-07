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
 
#ifndef ARIADNE_HYBRID_SET_H
#define ARIADNE_HYBRID_SET_H

#include <set>
#include <map>

#include <string>
#include <iostream>
#include <sstream>

#include "../geometry/set_reference.h"
#include "../geometry/grid_set.h"
#include "../geometry/hybrid_space.h"

namespace Ariadne {  
  namespace Geometry {

    /*! \brief The type identifying a discrete locatation of a hybrid system. */
    typedef id_type location_type;


    class HybridSystemError : public std::runtime_error {
     public:
      HybridSystemError(const std::string& what) : std::runtime_error(what) { }
    };


    /*! \ingroup HybridSet
     *  \brief A base class for representing subsets of hybrid state spaces.
     */
    template< class S >
    class HybridSetBase
    {
     public:
      typedef typename S::real_type real_type;
      typedef typename S::state_type state_type;
      typedef S set_type;
      typedef typename std::map<location_type,S>::iterator iterator;
      typedef typename std::map<location_type,S>::const_iterator const_iterator;
     public:
      /*! \brief Virtual destructor. */
      virtual ~HybridSetBase();
      /*! \brief Construct a set with no locations. */
      HybridSetBase();
      /*! \brief Construct from a dictionary of location identifyers and dimensions. */
      HybridSetBase(const HybridSpace& locations);
      /*! \brief Copy constructor. */
      HybridSetBase(const HybridSetBase<S>& hs);
      /*! \brief Conversion constructor from another hybrid set. */
      template<class S1> explicit HybridSetBase(const HybridSetBase<S1>& hs);
      
      /*! \brief Create a new location with dimension \a d. */
      S& new_location(location_type q, dimension_type d);
      /*! \brief Create a new location based on the set \a s. */
      S& new_location(location_type q, const S& s);
      /*! \brief Create a new location by constructing a set from an object of type T. */
      template<class T> S& new_location(location_type q, const T& t);

      /*! \brief The discrete locations of the set. */
      HybridSpace locations() const;
      /*! \brief The number of discrete locations or components comprising the set. */
      location_type number_of_locations() const;
      /*! \brief Check if the hybrid set has a component for discrete location \a q. */
      bool has_location(location_type q) const;
      /*! \brief A reference to the state set corresponding to discrete location \a q. */
      S& operator[](location_type q);
      /*! \brief The state set corresponding to discrete location \a q. */
      const S& operator[](location_type q) const;
      
      /*! \brief Clear all discrete locations. */
      void clear();
      /*! \brief Returns true if every component is empty. */
      tribool empty() const;
      /*! \brief Returns true if every component is bounded. */
      tribool bounded() const;
      
      /*! \brief Adjoin the set \a s to location \a q. */
      template<class S1> void adjoin(location_type q, const S1& s);
      /*! \brief Adjoin another hybrid set. */
      template<class S1> void adjoin(const HybridSetBase<S1>& hs);
      /*! \brief Restrict to hybrid set. */
      template<class S1> void restrict(const HybridSetBase<S1>& hs);
      /*! \brief Remove the set \a hs. */
      template<class S1> void remove(const HybridSetBase<S1>& hs);
      
      /*! \brief An iterator to the beginning of the component sets. */
      iterator begin();
      /*! \brief A constant iterator to the beginning of the component sets. */
      const_iterator begin() const;
      /*! \brief An iterator to the end of the component sets. */
      iterator end();
      /*! \brief A constant iterator to the end of the component sets. */
      const_iterator end() const;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
     private:
      std::map< location_type, S > _component_sets;
    };

    template<class S> 
    std::ostream& operator<<(std::ostream& os, const HybridSetBase<S>& hs);

    template<class S1, class S2 > 
    tribool
    subset(const HybridSetBase<S1>&, const HybridSetBase<S2>&);

    template<class S1, class S2 > 
    tribool
    disjoint(const HybridSetBase<S1>&, const HybridSetBase<S2>&);


    /*! \ingroup HybridSet
     *  \brief A hybrid set comprising of a GridMaskSet for every component.
     */
    template< class R >
    class HybridGridMaskSet
      : public HybridSetBase< GridMaskSet<R> >
    {
     public:
      HybridGridMaskSet()
        : HybridSetBase< GridMaskSet<R> >() { }
      HybridGridMaskSet(const HybridGridMaskSet<R>& hgms)
        : HybridSetBase< GridMaskSet<R> >(hgms) { }
      size_type size() const;
      size_type capacity() const;

      virtual std::ostream& write(std::ostream& os) const;
    };




    /*! \ingroup HybridSet
     *  \brief A hybrid set comprising of a GridCellListSet for every component.
     */
    template< class R >
    class HybridGridCellListSet
      : public HybridSetBase< GridCellListSet<R> >
    {
     public:
      HybridGridCellListSet()
        : HybridSetBase< GridCellListSet<R> >() { }
      HybridGridCellListSet(const HybridGridCellListSet<R>& hgcls)
        : HybridSetBase< GridCellListSet<R> >(hgcls) { }
      HybridGridCellListSet(const HybridGridMaskSet<R>& hgms)
        : HybridSetBase< GridCellListSet<R> >(hgms) { }
      size_type size() const;
      void unique_sort();
    };



    template<class R> 
    HybridGridCellListSet<R>
    difference(const HybridGridCellListSet<R>& hgcl, const HybridGridMaskSet<R>& hgms);

    template<class R> 
    HybridGridMaskSet<R>
    difference(const HybridGridMaskSet<R>& hgcl, const HybridGridMaskSet<R>& hgms);

    template<class R> 
    HybridGridCellListSet<R>
    regular_intersection(const HybridGridCellListSet<R>& hgcl, const HybridGridMaskSet<R>& hgms);

    template<class R> 
    HybridGridCellListSet<R>
    regular_intersection(const HybridGridMaskSet<R>& hgms, const HybridGridCellListSet<R>& hgcl);
   
    template<class R> 
    HybridGridMaskSet<R>
    regular_intersection(const HybridGridMaskSet<R>& hgms1, const HybridGridMaskSet<R>& hgms2);
   
    
    
    /*! \ingroup HybridSet
     *  \brief A hybrid set comprising of a GridMaskSet for every component.
     */
    template< class R >
    class HybridSet
      : public HybridSetBase< SetReference<R> >
    {
     public:
      HybridSet()
        : HybridSetBase< SetReference<R> >() { }
      HybridSet(const HybridSet<R>& hs)
        : HybridSetBase< SetReference<R> >(hs) { }
      template<class S1> HybridSet(const HybridSetBase<S1>& hs)
        : HybridSetBase< SetReference<R> >(hs) { }
    };


    
  }
}

#include "hybrid_set.inline.h"

#endif /* ARIADNE_HYBRID_SET_H */
