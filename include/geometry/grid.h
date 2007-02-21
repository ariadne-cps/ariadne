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

#ifndef _ARIADNE_GRID_H
#define _ARIADNE_GRID_H

#include <iosfwd>

#include "../declarations.h"
#include "../base/array.h"
#include "../numeric/arithmetic.h"
#include "../linear_algebra/vector.h"
#include "../combinatoric/lattice_set.h"

namespace Ariadne {
  namespace Geometry {
    template<class R> class Grid;
    template<class R> class RegularGrid;
    template<class R> class IrregularGrid;
    
    template<class R> class FiniteGrid;

    template<class R> class GridCell;
    template<class R> class GridBlock;
    
    
    /*! \ingroup Grid
     *  \brief %Base type for defining a grid.
     *
     *  We use inheritence and abstract functions here partly for ease of 
     *  development and partly since the grid coordinates only play a role when
     *  converting to rectangles and should occur with linear complexity in the 
     *  space dimension.
     */
    template<class R>
    class Grid {
    public:
      /*! \brief The type of real number defining the vertices and cells of the grid. */
      typedef R real_type;
     
      /*! \brief Destructor. */
      virtual ~Grid();

      /*! \brief Cloning operator. */
      virtual Grid<R>* clone() const = 0;

      /*! \brief Return the grid type. */
      virtual grid_type type() const = 0;
      
      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const = 0;

      /*! \brief The coordinate of the \a n th subdivision point in 
       * dimension \a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const = 0;

      /*! \brief The interval \f$[p_n,p_{n+1})\f$ in dimension \a d index 
       * containing \a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const = 0;

      /*! \brief Tests whether the grid contains the given rectangle within
       * its extent. */
      virtual bool encloses(const Rectangle<R>& r) const = 0;

      /*! \brief The index in dimension \a d of the subdivision 
       * point \a x. Throws an exception if \a x is not a subdivision point. */
      index_type subdivision_index(dimension_type d, const real_type& x) const;
      
      /*! \brief The index of the subdivision point below \a x. */
      index_type subdivision_lower_index(dimension_type d, const real_type& x) const;

      /*! \brief The index of the subdivision point above \a x. */
      index_type subdivision_upper_index(dimension_type d, const real_type& x) const;

      /*! \brief Tests equality of two grids. */
      virtual bool operator==(const Grid<R>& g) const; 
      /*! \brief Tests inequality of two grids. */
      virtual bool operator!=(const Grid<R>& g) const;

      /*! The index of the cell countaining the point \a s. */
      IndexArray index(const Point<R>& s) const;
      /*! The index of vertex to the lower-left of \a r. */
      IndexArray lower_index(const Rectangle<R>& r) const;
      /*! The index of vertex to the upper-right of \a r. */
      IndexArray upper_index(const Rectangle<R>& r) const;

      /*!\brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const = 0;
      /*!\brief Read from an input stream. */
      virtual std::istream& read(std::istream& is) = 0;
   };


   
    /*! \brief A finite, nonuniform grid of rectangles in Euclidean space. 
     *  \ingroup Grid
     */
    template<class R>
    class IrregularGrid : public Grid<R> {
     public:
      /*! \brief The type of real number defining the vertices and cells of the grid. */
      typedef R real_type;

      /*! \brief Cloning operator. */
      virtual IrregularGrid<R>* clone() const;
    
      /*! \brief Destructor. */
      virtual ~IrregularGrid();

      /*! \brief Return the grid type. */
      grid_type type() const;
      
      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const;
    
      /*! \brief The coordinate of the \a n th subdivision point in 
       * dimension \a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const;
  
      /*! \brief The index of interval in dimension \a d index 
       * containing \a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const;

      /*! \brief Tests whether the grid contains the given lattice rectangle 
       * within its extent. */
      virtual bool encloses(const Rectangle<R>& r) const;

      /*! \brief Construct from a bounding box and an equal number of 
       * subdivisions in each coordinate. */
      explicit IrregularGrid(const Rectangle<R>& r, size_type n);

      /*! \brief Construct from a bounding box and an array giving the number 
       * of subdivisions in each coordinate. */
      explicit IrregularGrid(const Rectangle<R>& r, SizeArray sz);

      /*! \brief Construct from a list of subdivision coordinates in each 
       * dimension. */
      explicit IrregularGrid(const array< std::vector<R> >& sp);

      /*! \brief Construct from a list of rectangles giving the grid points. */
      explicit IrregularGrid(const ListSet<R,Rectangle>& ls);

      /*! \brief Join two irregular grids. */
      IrregularGrid(const IrregularGrid& g1,IrregularGrid& g2);

      /*! \brief Tests equality of two irregular grids. */
      bool operator==(const Grid<R>& g) const; 
      
      /*! \brief Tests equality of two irregular grids. */
      bool operator!=(const Grid<R>& g) const;
            
      /*! \brief The lowest valid vertex index. */
      IndexArray lower() const;
      /*! \brief The highest valid vertex index. */
      IndexArray upper() const;
      /*! \brief The block of valid lattice cells. */
      Combinatoric::LatticeBlock lattice_block() const;
      /*! \brief The number of subdivision intervals in each dimension. */
      SizeArray sizes() const;
      /*! \brief The total number of cells. */
      size_type capacity() const;
      /*! \brief The number of subdivision intervals in dimension \a d. */
      size_type size(dimension_type i) const;

      /*! \brief The rectangle bounding the grid. */
      GridBlock<R> extent() const;
      
      /*! \brief Find the rule to translate elements from a grid to a 
       * refinement. */
      static array< std::vector<index_type> > index_translation(const IrregularGrid<R>& from, const IrregularGrid<R>& to);

      /*!\brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
      /*!\brief Read from an input stream. */
      virtual std::istream& read(std::istream& is);
     private:
      void create();
     private:
      array< std::vector<R> > _subdivision_coordinates;
    };



    /*! \brief An infinite, uniform grid of rectangles in Euclidean space.
     *  \ingroup Grid
     */
    template<class R> class RegularGrid : public Grid<R> {
     public:
      /*! \brief The type of real number defining the vertices and cells of the grid. */
      typedef R real_type;

      /*! \brief Cloning operator. */
      RegularGrid<R>* clone() const;
      
      /*! \brief Return the grid type. */
      grid_type type() const;
      
      /*! \brief Construct a regular grid in dimension \a n with subdivision lengths \a l in every direction. */
      explicit RegularGrid(const dimension_type& n, const R& l);

      /*! \brief Construct from an array of subdivision lengths \a sl. */
      explicit RegularGrid(const array<R>& sl);

      /*! \brief Construct from an array of subdivision lengths \a v. */
      explicit RegularGrid(const LinearAlgebra::Vector<R>& v);

      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const;

      /*! \brief The coordinate of the \a n th subdivision point in 
       * dimension \a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const;

      /*! \brief The length of the subdivision in the \a d th coordinate. */
      real_type subdivision_length(dimension_type d) const;

      /*! \brief The index of interval in dimension \a d index containing \a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const;

      /*! \brief True if \a r is contained within the region of definition of the grid. */
      virtual bool encloses(const Rectangle<R>& r) const;      
      
      /*! \brief Tests equality of two regular grids. */
      bool operator==(const Grid<R>& g) const; 
      
      /*! \brief Tests inequality of two regular grids. */
      bool operator!=(const Grid<R>& g) const; 

      /*!\brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
      /*!\brief Read from an input stream. */
      virtual std::istream& read(std::istream& is);
     private:
      array<R> _subdivision_lengths;
    };

    
    
    /*!\ingroup Grid
     * \brief A finite grid, suitable for defining a GridMaskSet.
     * \deprecated Use a Grid and LatticeBlock in GridMaskSet constructor instead. 
     */
    template<class R> class FiniteGrid {
     public:
      /*! \brief The type of real number defining the vertices and cells of the grid. */
      typedef R real_type;

      /*! \brief A finite grid generated by restricting a grid by specifying 
       * a lattice block. */
      FiniteGrid(const Grid<R>& g, const Combinatoric::LatticeBlock& b);
      
      /*! \brief A finite grid generated by subdividing a rectangle. */
      FiniteGrid(const Rectangle<R>& g, const size_type& s);
        
      /*! \brief A finite grid generated by restricting a using a spacial 
       * bounding box. */
      FiniteGrid(const Grid<R>& g, const Rectangle<R>& bb);

      
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
      const Grid<R>* _grid_ptr;
      const grid_type _grid_type; 
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

#endif /* _GRID_H */
