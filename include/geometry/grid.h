/***************************************************************************
 *            grid.h
 *
 *  Copyright  2005,6  Alberto Casagrande, Pieter Collins
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
 
/*! \file grid.h
 *  \brief Grids in coordinate space.
 */

#ifndef ARIADNE_GRID_H
#define ARIADNE_GRID_H

#include <iosfwd>
#include <boost/shared_ptr.hpp>

#include "../declarations.h"
#include "../base/array.h"
#include "../numeric/arithmetic.h"
#include "../linear_algebra/vector.h"
#include "../combinatoric/lattice_set.h"

namespace Ariadne {
  namespace Geometry {
    template<class R> class Grid;
    template<class R> class FiniteGrid;

    template<class R> class GridCell;
    template<class R> class GridBlock;
    
    template<class R> class Point;
    template<class R> class Rectangle;

    /*! \brief A grid of rectangles in Euclidean space with an irregular 'core' and 
     *  regular cells outside this core. */
    template<class R>
    class Grid
    { 
     public:
      /*! \brief The type of real number defining the vertices and cells of the grid. */
      typedef R real_type;
      
     public:
      /*! \brief \deprecated */
      index_type subdivision_interval(dimension_type d, const real_type& x) const { return subdivision_lower_index(d,x); }
      /*! \brief \deprecated */
      bool encloses(const Rectangle<R>& r) const { return true; }
      
     public:
      /*! \brief Destructor. */
      ~Grid();
      
      /*! \brief Construct a regular grid in dimension \a n with subdivision lengths \a l in every direction. */
      explicit Grid(const dimension_type& n, const R& l);
      
      /*! \brief Construct from a centre point and a vector of offsets. */
      explicit Grid(const Point<R>& pt, const LinearAlgebra::Vector<R>& v);
      
      /*! \brief Construct from an array of subdivision lengths \a v. */
      explicit Grid(const LinearAlgebra::Vector<R>& v);

      /*! \brief Construct from a list of subdivision coordinates in each 
       * dimension, the lower and upper lengths outsize this array, and the position of the central point of the grid. */
      explicit Grid(const array< array<R> >& sc, const array<R>& sl, const array<uint>& cp);

      /*! \brief Construct from raw pointers. */
       explicit Grid(const dimension_type& d, const size_type* nsc, const R** sc, const R* sl, const size_type* cp);

      /*! \brief Construct from a list of rectangles giving the grid points. */
      explicit Grid(const ListSet< Rectangle<R> >& ls);

      /*! \brief Join two grid, assuming they have the same behaviour for large indices. */
      Grid(const Grid<R>& g1, Grid<R>& g2);

      /*! \brief The underlying dimension of the grid. */
      dimension_type dimension() const;

      /*! \brief The coordinate of the \a n th subdivision point in dimension \a d. */
      real_type subdivision_coordinate(dimension_type d, index_type n) const;

      /*! \brief The index in dimension \a d of the subdivision 
       * point \a x. Throws an exception if \a x is not a subdivision point. */
      index_type subdivision_index(dimension_type d, const real_type& x) const;
      
      /*! \brief The index of the subdivision point below \a x. */
      index_type subdivision_lower_index(dimension_type d, const real_type& x) const;

      /*! \brief The index of the subdivision point above \a x. */
      index_type subdivision_upper_index(dimension_type d, const real_type& x) const;

      /*! \brief Tests equality of two grids. */
      bool operator==(const Grid<R>& g) const; 
      /*! \brief Tests inequality of two grids. */
      bool operator!=(const Grid<R>& g) const;

      /*! The block of irregular cells. */
      Combinatoric::LatticeBlock lattice_block() const;

      /*! The index of the cell countaining the point \a pt. */
      IndexArray index(const Point<R>& pt) const;
      /*! The index of vertex to the lower-left of \a r. */
      IndexArray lower_index(const Rectangle<R>& r) const;
      /*! The index of vertex to the upper-right of \a r. */
      IndexArray upper_index(const Rectangle<R>& r) const;

      /*! The vertex of the grid at index position \a index. */
      Point<R> point(const IndexArray& index) const;

      /*! The block of the grid bounded by index positions given by \a indices. */
      Rectangle<R> rectangle(const Combinatoric::LatticeBlock& indices) const;

      /*!\brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*!\brief Read from an input stream. */
      std::istream& read(std::istream& is);
     private:
      /*! Create cached data used to speed up computations. */
      void create();
     private:
      array< array<R> > _subdivision_coordinates;
      array<R> _subdivision_lengths;
      array<uint> _centre_positions;
      array<R*> _centre_pointers;
      array<int> _lower_indices;
      array<int> _upper_indices;
      array<R> _lower_bounds;
      array<R> _upper_bounds;
      array<R> _lower_origin;
      array<R> _upper_origin;
    };

    
    /*!\ingroup Grid
     * \brief A finite grid, suitable for defining a GridMaskSet.
     * \deprecated Use a Grid and LatticeBlock in GridMaskSet constructor instead. 
     */
    template<class R> class FiniteGrid {
     public:
      /*! \brief The type of real number defining the vertices and cells of the grid. */
      typedef R real_type;

      /*! \brief Destructor. Frees memory allocated to hold the grid if necessary. */
      ~FiniteGrid();
      
      /*! \brief A finite grid generated by restricting a grid by specifying 
       * a lattice block. */
      FiniteGrid(const Grid<R>& g, const Combinatoric::LatticeBlock& b);
      
      /*! \brief A finite grid generated by subdividing a rectangle. */
      FiniteGrid(const Rectangle<R>& g, const size_type& s);
        
      /*! \brief A finite grid generated by restricting a using a spacial 
       * bounding box. */
      FiniteGrid(const Grid<R>& g, const Rectangle<R>& bb);

      /*! \brief Construct from a list of rectangles giving the grid points. */
      explicit FiniteGrid(const ListSet< Rectangle<R> >& ls);

      
      /*! \brief The underlying grid. */
      const Grid<R>& grid() const;     
      
      /*! \brief The lattice bounds of the finite grid. */
      const Combinatoric::LatticeBlock& lattice_block() const;
     
      /*! \brief The extent of the finite grid in state space. */
      Rectangle<R> extent() const;

      /*! \brief The dimension of the grid. */
      dimension_type dimension() const;

      /*! \brief The coordinate of the \a n th subdivision point in 
       * dimension \a d. */
      real_type subdivision_coordinate(dimension_type d, index_type n) const;
  
      /*! \brief The index of interval in dimension \a d index 
       * containing \a x. */
      index_type subdivision_interval(dimension_type d, const real_type& x) const;
      
      /*! \brief The index of the subdivision point below \a x. */
      index_type subdivision_lower_index(dimension_type d, const real_type& x) const;

      /*! \brief The index of the subdivision point above \a x. */
      index_type subdivision_upper_index(dimension_type d, const real_type& x) const;
      
      /*!\brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
      // true if the class 'owns' the _grid_ptr, which means that the class is
      // responsible for deleting the Grid
      bool _own_ptr;
      // A pointer to an inifinite grid
      const Grid<R>* _grid_ptr;
      // The integer indices of the bounds of the cells covered by the grid.
      Combinatoric::LatticeBlock _lattice_block;
    };
   
    
    
    template<class R> 
    std::istream& operator>>(std::istream& is, Grid<R>& g);

    template<class R> 
    std::ostream& operator<<(std::ostream& os, const Grid<R>& g);

    template<class R>
    std::ostream& operator<<(std::ostream& os, const FiniteGrid<R>& fg);

  }
}

#include "grid.inline.h"

#endif /* ARIADNE_GRID_H */
