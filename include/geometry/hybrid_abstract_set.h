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

#include "geometry/set_interface.h"
#include "geometry/grid_cell_list_set.h"
#include "geometry/grid_mask_set.h"
#include "geometry/discrete_state.h"
#include "geometry/hybrid_space.h"
#include "geometry/hybrid_set_iterator.h"

#include "geometry/hybrid_basic_set.h"
#include "geometry/hybrid_denotable_set.h"

namespace Ariadne {
  

    class basic_set_tag;
    class denotable_set_tag;
    class list_set_tag;

 
    /*! \ingroup HybridSet \ingroup SetInterface
     *  \brief Interface for hybrid sets. 
     */
    template< class R >
    class HybridSetInterface
      : public SetInterface< HybridBox<R> >
    { };

    /*! \ingroup HybridSet
     *  \brief A class for representing subsets of hybrid state spaces.
     */
    template< class R >
    class HybridSet
      : public SetInterface< HybridBox<R> >
    {
      typedef SetInterface< Box<R> > S;
      typedef typename S::basic_set_type BS;
      typedef boost::shared_ptr<S> SetPointer;
    public:
      typedef abstract_set_tag set_category;
      typedef R real_type;
      typedef HybridBox<R> basic_set_type;
      typedef DiscreteState discrete_state_type;
      typedef typename std::map<DiscreteState,SetPointer>::iterator locations_iterator;
      typedef typename std::map<DiscreteState,SetPointer>::const_iterator locations_const_iterator;
     public:
      /*! \brief Construct a set with no locations. */
      HybridSet();
      /*! \brief Construct from a dictionary of location identifyers and dimensions. */
      HybridSet(const HybridSpace& locations);
      /*! \brief Copy constructor. */
      HybridSet(const HybridSet<R>& hs);
      /*! \brief Copy assignment operator. */
      HybridSet<R>& operator=(const HybridSet<R>& hs);
      /*! \brief Conversion constructor from another hybrid set. */
      template<class DS> explicit HybridSet(const HybridDenotableSet<DS>& hs);
     public:
      /*! \brief Virtual destructor. */
      virtual ~HybridSet();
      
      virtual dimension_type dimension() const { throw NotImplemented(__PRETTY_FUNCTION__); }
      virtual HybridSet<R>* clone() const { throw NotImplemented(__PRETTY_FUNCTION__); }
   
      /*! \brief Create a new location with dimension \a d. */
      S& new_location(DiscreteState q, dimension_type d);
      /*! \brief Create a new location based on the set \a s. */
      S& new_location(DiscreteState q, const S& s);
      /*! \brief Create a new location based on the parameter \a t. */
      template<class T> S& new_location(DiscreteState q, const T& t);
      /*! \brief Create a new location based on the parameters \a t1 and \a t2. */
      template<class T1, class T2> S& new_location(DiscreteState q, const T1& t1, const T2& t2);
      /*! \brief Create a new location based on the rectangle \a r. */
      S& new_location(DiscreteState q, const Rectangle<R>& r);
      /*! \brief Create a new location based on the polyhedron \a p. */
      S& new_location(DiscreteState q, const Polyhedron<R>& p);

      /*! \brief The discrete locations of the set. */
      HybridSpace locations() const;
      /*! \brief The number of discrete locations or components comprising the set. */
      size_type number_of_locations() const;
      /*! \brief Check if the hybrid set has a component for discrete location \a q. */
      bool has_location(DiscreteState q) const;
      /*! \brief A reference to the state set corresponding to discrete location \a q. */
      S& operator[](DiscreteState q);
      /*! \brief The state set corresponding to discrete location \a q. */
      const S& operator[](DiscreteState q) const;
      
      /*! \brief Clear all discrete locations. */
      void clear();
      /*! \brief Returns true if every component is empty. */
      tribool empty() const;
      
      /*! \brief An iterator to the beginning of the component sets. */
      locations_iterator locations_begin();
      /*! \brief A constant iterator to the beginning of the component sets. */
      locations_const_iterator locations_begin() const;
      /*! \brief An iterator to the end of the component sets. */
      locations_iterator locations_end();
      /*! \brief A constant iterator to the end of the component sets. */
      locations_const_iterator locations_end() const;

      /*! \brief Test if the set contains a point. (Not currently implemented) */
      virtual HybridSpace space() const {
        return this->locations(); }
      /*! \brief Test if the set contains a point. (Not currently implemented) */
      virtual tribool contains(const HybridPoint<R>& pt) const {
        throw NotImplemented(__PRETTY_FUNCTION__); }
      /*! \brief Test if the set is a superset of a hybrid box. */
      virtual tribool superset(const HybridBox<R>& bx) const {
        return (*this)[bx.state()].superset(bx.set()); }
      /*! \brief Test if the set intersects a hybrid box. */
      virtual tribool intersects(const HybridBox<R>& bx) const {
        return (*this)[bx.state()].intersects(bx.set()); }
      /*! \brief Test if the set is discjoint from a hybrid box. */
      virtual tribool disjoint(const HybridBox<R>& bx) const {
        return (*this)[bx.state()].disjoint(bx.set()); }
      /*! \brief Test if the set is discjoint from a hybrid box. */
      virtual tribool subset(const HybridBox<R>& bx) const {
        throw NotImplemented(__PRETTY_FUNCTION__); }
      /*! \brief Test if the set is discjoint from a hybrid box. */
      virtual HybridBox<R> bounding_box() const {
        throw NotImplemented(__PRETTY_FUNCTION__); }

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
     private:
      std::map< DiscreteState, SetPointer > _component_sets;
    };


    template<class S> 
    std::ostream& operator<<(std::ostream& os, const HybridSet<S>& hs);

    template<class S1, class S2 > 
    tribool
    subset(const HybridSet<S1>&, const HybridSet<S2>&);

    template<class S1, class S2 > 
    tribool
    disjoint(const HybridSet<S1>&, const HybridSet<S2>&);




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


   
} // namespace Ariadne

#include "hybrid_abstract_set.inline.h"

#endif /* ARIADNE_HYBRID_ABSTRACT_SET_H */
