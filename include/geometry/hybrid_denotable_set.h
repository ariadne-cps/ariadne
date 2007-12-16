/***************************************************************************
 *            hybrid_denotable_set.h
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
 
#ifndef ARIADNE_HYBRID_DENOTABLE_SET_H
#define ARIADNE_HYBRID_DENOTABLE_SET_H

#include <set>
#include <map>

#include <string>
#include <iostream>
#include <sstream>

#include <boost/shared_ptr.hpp>

#include "base/types.h"

#include "geometry/geometrical_traits.h"
#include "geometry/set_interface.h"
#include "geometry/grid_cell_list_set.h"
#include "geometry/grid_mask_set.h"
#include "geometry/hybrid_space.h"
#include "geometry/hybrid_set_iterator.h"
#include "geometry/hybrid_basic_set.h"



namespace Ariadne {
  namespace Geometry {

    class basic_set_tag;
    class denotable_set_tag;
    class list_set_tag;



    /*! \ingroup HybridSet
     *  \brief A base class for representing subsets of hybrid state spaces.
     */
    template< class DS >
    class HybridDenotableSet
    {
     public:
      typedef denotable_set_tag set_category;
      typedef typename DS::real_type real_type;
      typedef typename DS::state_type state_type;
      typedef DiscreteState discrete_state_type;
      typedef DS continuous_state_set_type;
      typedef typename DS::basic_set_type continuous_basic_set_type;
      typedef HybridBasicSet<typename DS::basic_set_type> basic_set_type;
      typedef typename std::map<discrete_state_type,continuous_state_set_type>::iterator locations_iterator;
      typedef typename std::map<discrete_state_type,continuous_state_set_type>::const_iterator locations_const_iterator;
      typedef HybridDenotableSetIterator<DS> const_iterator;
     protected:
      /*! \brief Construct a set with no locations. */
      HybridDenotableSet();
      /*! \brief Construct from a dictionary of location identifyers and dimensions. */
      HybridDenotableSet(const HybridSpace& locations);
      /*! \brief Copy constructor. */
      HybridDenotableSet(const HybridDenotableSet<DS>& hs);
      /*! \brief Copy assignment operator. */
      HybridDenotableSet<DS>& operator=(const HybridDenotableSet<DS>& hs);
      /*! \brief Conversion constructor from another hybrid set. */
      template<class S> explicit HybridDenotableSet(const HybridDenotableSet<S>& hs);
     public:
      /*! \brief Virtual destructor. */
      virtual ~HybridDenotableSet();
      
      /*! \brief Create a new location with dimension \a d. */
      DS& new_location(discrete_state_type q, dimension_type d);
      /*! \brief Create a new location based on the set \a s. */
      DS& new_location(discrete_state_type q, const DS& s);
      /*! \brief Create a new location based on the parameter \a t. */
      template<class T> DS& new_location(discrete_state_type q, const T& t);
      /*! \brief Create a new location based on the parameters \a t1 and \a t2. */
      template<class T1, class T2> DS& new_location(discrete_state_type q, const T1& t1, const T2& t2);

      /*! \brief The discrete locations of the set. */
      HybridSpace locations() const;
      /*! \brief The space the set lies in. */
      HybridSpace space() const;
      /*! \brief The number of discrete locations or components comprising the set. */
      size_type number_of_locations() const;
      /*! \brief Check if the hybrid set has a component for discrete location \a q. */
      bool has_location(discrete_state_type q) const;
      /*! \brief A reference to the state set corresponding to discrete location \a q. */
      DS& operator[](discrete_state_type q);
      /*! \brief The state set corresponding to discrete location \a q. */
      const DS& operator[](discrete_state_type q) const;
      
      /*! \brief Clear all discrete locations. */
      void clear();
      /*! \brief Returns true if every component is empty. */
      tribool empty() const;
      /*! \brief Returns true if every component is bounded. */
      tribool bounded() const;

      /*! \brief The total number of component sets. */
      size_type size() const;
      
      /*! \brief Adjoin the set \a s to location \a q. */
      template<class S1> void adjoin(discrete_state_type q, const S1& s);
      /*! \brief Adjoin the set \a s to location \a q. */
      template<class S1> void adjoin(const HybridBasicSet<S1>& s);
      /*! \brief Adjoin another hybrid set. */
      template<class S1> void adjoin(const HybridDenotableSet<S1>& hs);
      /*! \brief Restrict to hybrid set. */
      template<class S1> void restrict(const HybridDenotableSet<S1>& hs);
      /*! \brief Remove the set \a hs. */
      template<class S1> void remove(const HybridDenotableSet<S1>& hs);
      
      /*! \brief An iterator to the beginning of the component sets. */
      locations_iterator locations_begin();
      /*! \brief A constant iterator to the beginning of the component sets. */
      locations_const_iterator locations_begin() const;
      /*! \brief An iterator to the end of the component sets. */
      locations_iterator locations_end();
      /*! \brief A constant iterator to the end of the component sets. */
      locations_const_iterator locations_end() const;

      /*! \brief A constant iterator to the beginning of the component sets. */
      const_iterator begin() const;
      /*! \brief A constant iterator to the end of the component sets. */
      const_iterator end() const;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
     protected:
      

     private:
      std::map< discrete_state_type, continuous_state_set_type > _component_sets;
    };


    template<class S> 
    std::ostream& operator<<(std::ostream& os, const HybridDenotableSet<S>& hs);

    template<class DS1, class DS2 > 
    tribool
    subset(const HybridDenotableSet<DS1>&, const HybridDenotableSet<DS2>&);

    template<class DS1, class DS2 > 
    tribool
    disjoint(const HybridDenotableSet<DS1>&, const HybridDenotableSet<DS2>&);



    /*! \ingroup HybridSet
     *  \brief A hybrid set comprising of a ListSet for every component.
     */
    template< class BS >
    class HybridListSet
      : public HybridDenotableSet< ListSet<BS> >
    {
      typedef typename BS::real_type R;
     public:
      typedef list_set_tag set_category;
      HybridListSet()
        : HybridDenotableSet< ListSet<BS> >() { }
      HybridListSet(const HybridSpace& hspc)
        : HybridDenotableSet< ListSet<BS> >(hspc) { }
      HybridListSet(const HybridListSet<BS>& hls)
        : HybridDenotableSet< ListSet<BS> >(hls) { }
      template<class HBS> void adjoin_over_approximation(const HBS& hbs) {
        this->adjoin(hbs.discrete_state(),over_approximation(hbs)); }
      virtual std::ostream& write(std::ostream&) const;
    };



    /*! \ingroup HybridSet
     *  \brief A hybrid set comprising of a single basic set for in a discrete mode.
     */
    template<class R> 
    class HybridGrid 
    {
     public:
      typedef DiscreteState discrete_state_type;
      typedef typename std::map< discrete_state_type, Grid<R> >::const_iterator locations_const_iterator;

      /*! \brief */
      HybridGrid() : _grids() { }
      /*! \brief */
      template<class GS> HybridGrid(const GS& gs) : _grids() { 
        for(typename GS::locations_const_iterator iter=gs.locations_begin(); iter!=gs.locations_end(); ++iter) {
          this->_grids.insert(std::make_pair(iter->first,iter->second.grid())); }
      }
      /*! \brief */
      void new_location(DiscreteState q, const Grid<R>& g) {
        this->_grids.insert(std::make_pair(q,g)); }
      /*! \brief */
      const Grid<R>& operator[](DiscreteState q) const {
        return this->_grids.find(q)->second; }

      HybridSpace locations() const;
      locations_const_iterator locations_begin() const { return this->_grids.begin(); }
      locations_const_iterator locations_end() const { return this->_grids.end(); }
     private:
      std::map< discrete_state_type, Grid<R> > _grids;
    };



    /*! \ingroup HybridSet
     *  \brief A hybrid set comprising of a grid cell set for in a discrete mode.
     */
    template<class R> 
    class HybridGridCell 
      : public HybridBasicSet< GridCell<R> >
    {
      typedef DiscreteState discrete_state_type;
     public:
      /*! \brief */
      HybridGridCell(const HybridBasicSet< GridCell<R> >& hbs) : HybridBasicSet< GridCell<R> >(hbs) { }
      HybridGridCell<R> operator=(const HybridBasicSet< GridCell<R> >& hbs) {
        if(this!=&hbs) { static_cast<HybridBasicSet< GridCell<R> >&>(*this)=hbs; } }
    };


    /*! \ingroup HybridSet
     *  \brief A hybrid set comprising of a GridMaskSet for every component.
     */
    template< class R >
    class HybridGridMaskSet
      : public HybridDenotableSet< GridMaskSet<R> >
    {
      typedef DiscreteState discrete_state_type;
     public:
      HybridGridMaskSet()
        : HybridDenotableSet< GridMaskSet<R> >() { }
      HybridGridMaskSet(const HybridGridMaskSet<R>& hgms)
        : HybridDenotableSet< GridMaskSet<R> >(hgms) { }
      HybridGrid<R> grid() const { return HybridGrid<R>(*this); }
      template<class HBS> void adjoin_outer_approximation(const HBS& hbs) {
        this->operator[](hbs.discrete_state()).adjoin_outer_approximation(hbs.continuous_state_set()); }
      virtual std::ostream& write(std::ostream& os) const;
    };




    /*! \ingroup HybridSet
     *  \brief A hybrid set comprising of a GridCellListSet for every component.
     */
    template< class R >
    class HybridGridCellListSet
      : public HybridDenotableSet< GridCellListSet<R> >
    {
     public:
      HybridGridCellListSet()
        : HybridDenotableSet< GridCellListSet<R> >() { }
      HybridGridCellListSet(const HybridGrid<R>& hg);
      HybridGridCellListSet(const HybridGridCellListSet<R>& hgcls)
        : HybridDenotableSet< GridCellListSet<R> >(hgcls) { }
      HybridGridCellListSet(const HybridGridMaskSet<R>& hgms)
        : HybridDenotableSet< GridCellListSet<R> >(hgms) { }
      HybridGrid<R> grid() const { return HybridGrid<R>(*this); }
      template<class HBS> void adjoin_outer_approximation(const HBS& hbs) {
        this->adjoin_outer_approximation(hbs,typename HBS::set_category()); }
      template<class HBS> void adjoin_outer_approximation(const HBS& hbs,basic_set_tag) {
        this->operator[](hbs.discrete_state()).adjoin_outer_approximation(hbs.continuous_state_set()); }
      template<class HDS> void adjoin_outer_approximation(const HDS& hds,denotable_set_tag) {
        for(typename HDS::const_iterator iter=hds.begin(); iter!=hds.end(); ++iter) { this->adjoin_outer_approximation(*iter); } }
      void unique_sort();
    };

    template<class BS> 
    HybridGridCellListSet<typename BS::real_type>
    outer_approximation(const HybridListSet<BS>& hls, const HybridGrid<typename BS::real_type>& hg);

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
   
   
  }
}


#include "hybrid_denotable_set.inline.h"

#endif /* ARIADNE_HYBRID_DENOTABLE_SET_H */
