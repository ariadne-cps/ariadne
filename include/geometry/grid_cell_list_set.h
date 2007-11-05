/***************************************************************************
 *            grid_cell_list_set.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file grid_cell_list_set.h
 *  \brief A list of grid cells
 */

#ifndef ARIADNE_GRID_CELL_LIST_SET_H
#define ARIADNE_GRID_CELL_LIST_SET_H

#include <iosfwd>

#include <boost/iterator/iterator_adaptor.hpp>

#include "../base/array.h"
#include "../base/tribool.h"
#include "../base/pointer.h"

#include "../combinatoric/lattice_set.h"

#include "../geometry/declarations.h"
#include "../geometry/set_interface.h"
#include "../geometry/grid.h"
#include "../geometry/grid_cell.h"

namespace Ariadne {
  namespace Geometry {
 
    template<class Base, class Value> class GridSetIterator;

    /*! \brief A denotable set on a grid, defined as a list of cells.
     *
     *  The data for a %GridCellListSet provided by a Grid and a LatticeCellList.
     *
     *  A cell list set is useful for storing small sparse sets, such as an 
     *  over-approximation to a basic set.
     *
     *  \ingroup DenotableSet
     *  \ingroup Grid
     */
    template<class R>
    class GridCellListSet 
      : public SetInterface<R>
    {
      friend class GridMaskSet<R>;
     public:
      typedef GridSetIterator< Combinatoric::LatticeCellListSet::const_iterator, GridCell<R> > iterator;
      typedef GridSetIterator< Combinatoric::LatticeCellListSet::const_iterator, GridCell<R> > const_iterator;

      /*! \brief A tag describing the type of set. */
      typedef denotable_set_tag set_category;
      /*! \brief The type of basic set used in the list. */
      typedef GridCell<R> basic_set_type;
      /*! \brief The type of denotable real number defining the vertices and cells of the grid. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the set. */
      typedef Point<R> state_type;
      /*! \brief The type of basic set contained by the denotable set. */
      typedef GridCell<R> value_type;
      
      /*!\brief Construct an empty set based on a Grid. */
      explicit GridCellListSet(const Grid<R>& g);

      /*!\brief Construct from a grid and a lattice cell list set. */
      explicit GridCellListSet(const Grid<R>& g, const Combinatoric::LatticeCellListSet& lcls);

      /*!\brief Convert from a GridMaskSet. */
      GridCellListSet(const GridMaskSet<R>& gms);

      /*!\brief Copy constructor. */
      GridCellListSet(const GridCellListSet<R>& gcls);

      /*!\brief Copy assignment. */
      GridCellListSet<R>& operator=(const GridCellListSet<R>& gcls);

      /*!\brief Convert to a ListSet of Rectangles. */
      operator ListSet< Rectangle<R> >() const;

      /*! \brief The underlying grid. */
      const Grid<R>& grid() const;



      //@{
      //! \name SetInterface methods
      /*! \brief Make a dynamically-allocated copy of the set. */
      virtual GridCellListSet<R>* clone() const;

      /*! \brief The space dimension of the set. */
      virtual dimension_type dimension() const;

      /*!\brief Checks if a denotable set includes a point. */
      virtual tribool contains(const Point<R>& p) const;

      /*! \brief Tests for superset of a Rectangle. */ 
      virtual tribool superset(const Rectangle<R>& r) const;

      /*! \brief Tests for intersection with a Rectangle. */
      virtual tribool intersects(const Rectangle<R>& r) const;

      /*! \brief Tests for disjointness with a Rectangle. */
      virtual tribool disjoint(const Rectangle<R>& r) const;

      /*! \brief Tests for subset of a Rectangle. */
      virtual tribool subset(const Rectangle<R>& r) const;

      /*! \brief The rectangle bounding the region of the mask. */
      virtual Rectangle<R> bounding_box() const;
      //@}


      /*! \brief True if the set is empty. */
      tribool empty() const;

      /*! \brief True if the set is bounded. */
      tribool bounded() const;

      /*! \brief The numeber of cells in the list. */
      size_type size() const;

      /*! \brief Empties the set. */
      void clear();

      /*! \brief The lattice set. */
      const Combinatoric::LatticeCellListSet& lattice_set() const;

      /*!\brief The @a i th cell in the list. */
      GridCell<R> operator[] (const size_type i) const;

      /*! \brief A constant iterator to the beginning of the list. */
      const_iterator begin() const;

      /*! \brief A constant iterator to the end of the list. */
      const_iterator end() const;

      /*! \brief Sorts the cells lexicographically, removing duplicates. */
      void unique_sort(); 

      /*!\brief Append a GridCell to the list. */
      void adjoin(const GridCell<R>& c);
      /*!\brief Append a GridBlock to the list. */
      void adjoin(const GridBlock<R>& bl);
      /*!\brief Append a GridCellListSet to the list. */
      void adjoin(const GridCellListSet<R>& cls);

      /*! \brief Adjoins an over-approximation of the rectangle \a r. */
      void adjoin_over_approximation(const Rectangle<R>& r);

      /*! \brief Adjoins an outer-approximation of the set \a s. */
      template<class S>
      void adjoin_outer_approximation(const S& s);

      /*! \brief Adjoins an inner-approximation of the set \a s. */
      template<class S>
      void adjoin_inner_approximation(const S& s);

      /*! \brief Restricts to the outer-approximation of the set \a s. */
      void restrict_outer_approximation(const SetInterface<R>& s);

      /*! \brief Restricts to the inner-approximation of the set \a s. */
      void restrict_inner_approximation(const SetInterface<R>& s);

      /*! \brief Write a summary to an output stream. */
      std::ostream& summarize(std::ostream&) const;
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     private: 
      static void _instantiate();
     private:
      Base::shared_ptr< Grid<R> > _grid_ptr;
      Combinatoric::LatticeCellListSet _lattice_set;
    };

    template<class R> tribool subset(const GridCellListSet<R>&, const GridBlock<R>&);
    template<class R> tribool subset(const GridCellListSet<R>&, const GridCellListSet<R>&);
 
    template<class R> std::ostream& operator<<(std::ostream& os, const GridCellListSet<R>& gcls);
    
    
  }
}


#include "grid_cell_list_set.inline.h"

#endif /* ARIADNE_GRID_CELL_LIST_SET_H */
