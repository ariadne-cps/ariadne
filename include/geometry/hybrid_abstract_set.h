/***************************************************************************
 *            hybrid_abstract_set.h
 *
 *  Copyright  2006-7  Alberto Casagrande,  Pieter Collins
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
 
#ifndef ARIADNE_HYBRID_ABSTRACT_SET_H
#define ARIADNE_HYBRID_ABSTRACT_SET_H

#include <set>
#include <map>

#include <string>
#include <iostream>
#include <sstream>

#include <boost/shared_ptr.hpp>

#include "../geometry/set_interface.h"
#include "../geometry/grid_cell_list_set.h"
#include "../geometry/grid_mask_set.h"
#include "../geometry/discrete_state.h"
#include "../geometry/hybrid_space.h"
#include "../geometry/hybrid_set_iterator.h"

#include "../geometry/hybrid_basic_set.h"
#include "../geometry/hybrid_denotable_set.h"

namespace Ariadne {
  namespace Geometry {

    class basic_set_tag;
    class denotable_set_tag;
    class list_set_tag;


    /*! \ingroup HybridSet
     *  \brief A class for representing subsets of hybrid state spaces.
     */
    template< class S >
    class HybridAbstractSet
    {
      typedef typename S::real_type R;
     public:
      typedef DiscreteState discrete_state_type;
      typedef boost::shared_ptr<S> pointer_type;
      typedef typename S::real_type real_type;
      typedef typename S::state_type state_type;
      typedef S set_type;
      typedef typename std::map<discrete_state_type,pointer_type>::iterator locations_iterator;
      typedef typename std::map<discrete_state_type,pointer_type>::const_iterator locations_const_iterator;
     protected:
      /*! \brief Construct a set with no locations. */
      HybridAbstractSet();
      /*! \brief Construct from a dictionary of location identifyers and dimensions. */
      HybridAbstractSet(const HybridSpace& locations);
      /*! \brief Copy constructor. */
      HybridAbstractSet(const HybridAbstractSet<S>& hs);
      /*! \brief Copy assignment operator. */
      HybridAbstractSet<S>& operator=(const HybridAbstractSet<S>& hs);
      /*! \brief Conversion constructor from another hybrid set. */
      template<class DS> explicit HybridAbstractSet(const HybridDenotableSet<DS>& hs);
      /*! \brief Conversion constructor from another hybrid set. */
      template<class S1> explicit HybridAbstractSet(const HybridAbstractSet<S1>& hs);
     public:
      /*! \brief Virtual destructor. */
      virtual ~HybridAbstractSet();
      
      /*! \brief Create a new location with dimension \a d. */
      S& new_location(discrete_state_type q, dimension_type d);
      /*! \brief Create a new location based on the set \a s. */
      S& new_location(discrete_state_type q, const S& s);
      /*! \brief Create a new location based on the parameter \a t. */
      template<class T> S& new_location(discrete_state_type q, const T& t);
      /*! \brief Create a new location based on the parameters \a t1 and \a t2. */
      template<class T1, class T2> S& new_location(discrete_state_type q, const T1& t1, const T2& t2);
      /*! \brief Create a new location based on the rectangle \a r. */
      S& new_location(discrete_state_type q, const Rectangle<R>& r);
      /*! \brief Create a new location based on the polyhedron \a p. */
      S& new_location(discrete_state_type q, const Polyhedron<R>& p);

      /*! \brief The discrete locations of the set. */
      HybridSpace locations() const;
      /*! \brief The number of discrete locations or components comprising the set. */
      size_type number_of_locations() const;
      /*! \brief Check if the hybrid set has a component for discrete location \a q. */
      bool has_location(discrete_state_type q) const;
      /*! \brief A reference to the state set corresponding to discrete location \a q. */
      S& operator[](discrete_state_type q);
      /*! \brief The state set corresponding to discrete location \a q. */
      const S& operator[](discrete_state_type q) const;
      
      /*! \brief Clear all discrete locations. */
      void clear();
      /*! \brief Returns true if every component is empty. */
      tribool empty() const;
      /*! \brief Returns true if every component is bounded. */
      tribool bounded() const;
      
      /*! \brief An iterator to the beginning of the component sets. */
      locations_iterator locations_begin();
      /*! \brief A constant iterator to the beginning of the component sets. */
      locations_const_iterator locations_begin() const;
      /*! \brief An iterator to the end of the component sets. */
      locations_iterator locations_end();
      /*! \brief A constant iterator to the end of the component sets. */
      locations_const_iterator locations_end() const;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
     private:
      std::map< discrete_state_type, pointer_type > _component_sets;
    };


    template<class S> 
    std::ostream& operator<<(std::ostream& os, const HybridAbstractSet<S>& hs);

    template<class S1, class S2 > 
    tribool
    subset(const HybridAbstractSet<S1>&, const HybridAbstractSet<S2>&);

    template<class S1, class S2 > 
    tribool
    disjoint(const HybridAbstractSet<S1>&, const HybridAbstractSet<S2>&);


    /*! \ingroup HybridSet
     *  \brief A map on cells of a hybrid grid set.
     */
    template< class R >
    class HybridSet
      : public HybridAbstractSet< SetInterface<R> >
    {
     public:
      HybridSet() : HybridAbstractSet< SetInterface<R> >() { }
      template<class T> HybridSet(const T& t) : HybridAbstractSet< SetInterface<R> >(t) { }
    };
  


    /*! \ingroup HybridSet
     *  \brief A map on cells of a hybrid grid set.
     */
    template< class R >
    class HybridGridMultiMap 
    { 
      typedef HybridGridCell<R> hybrid_basic_set_type;
      typedef HybridGridCellListSet<R> hybrid_denotable_set_type;
     public:
      HybridGridMultiMap(const HybridGrid<R>&, const HybridGrid<R>&) { assert(false); }
      bool has_key(const hybrid_basic_set_type&) const { assert(false); return false; }
      void set_image(const hybrid_basic_set_type&,const hybrid_denotable_set_type&) { assert(false); }
      const hybrid_denotable_set_type& image(const hybrid_basic_set_type&) const { assert(false); }
    };


   
  }
}

#include "hybrid_abstract_set.inline.h"

#endif /* ARIADNE_HYBRID_ABSTRACT_SET_H */
