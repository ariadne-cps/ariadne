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
#include "../combinatoric/lattice_set.h"

namespace Ariadne {
  namespace Geometry {
    template<typename R> class Grid;
    template<typename R> class FiniteGrid;
    template<typename R> class RegularGrid;
    template<typename R> class IrregularGrid;
    
    template<typename R> class GridCell;
    template<typename R> class GridRectangle;

    template<typename R> std::ostream& operator<<(std::ostream&, const Grid<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Grid<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const FiniteGrid<R>&);
   
    /*! \ingroup Grid
     *  \brief %Base type for defining a grid.
     *
     *  We use inheritence and abstract functions here partly for ease of 
     *  development and partly since the grid coordinates only play a role when
     *  converting to rectangles and should occur with linear complexity in the 
     *  space dimension.
     */
    template<typename R>
    class Grid {
    public:
      /*! \brief The type of real number defining the vertices and cells of the grid. */
      typedef R real_type;
     
      /*! \brief Destructor. */
      virtual ~Grid() { }

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
       * its bounds. */
      virtual bool bounds_enclose(const Rectangle<R>& r) const = 0;

      /*! \brief The index in dimension \a d of the subdivision 
       * point \a x. Throws an exception if \a x is not a subdivision point. */
      index_type subdivision_index(dimension_type d, const real_type& x) const {
        index_type n=subdivision_interval(d,x);
        if(subdivision_coordinate(d,n) == x) { return n; }
        throw std::domain_error("Value is not a grid coordinate");
      }

      /*! \brief The index of the subdivision point below \a x. */
      index_type subdivision_lower_index(dimension_type d, const real_type& x) const {
        return subdivision_interval(d,x);
      }

      /*! \brief The index of the subdivision point above \a x. */
      index_type subdivision_upper_index(dimension_type d, const real_type& x) const {
        index_type n=subdivision_interval(d,x);
        return subdivision_coordinate(d,n) == x ? n : n+1;
      }

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
    template<typename R>
    class IrregularGrid : public Grid<R> {
     public:
      /*! \brief The type of real number defining the vertices and cells of the grid. */
      typedef R real_type;

      /*! \brief Cloning operator. */
      virtual IrregularGrid<R>* clone() const;
    
      /*! \brief Destructor. */
      virtual ~IrregularGrid();

      /*! \brief Return the grid type. */
      grid_type type() const {return IRREGULAR;}
      
      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const { 
        return this->_subdivision_coordinates.size();
      }
    
      /*! \brief The coordinate of the \a n th subdivision point in 
       * dimension \a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const {
        assert(d<this->dimension());
        if(!(0<=n && uint(n)<_subdivision_coordinates[d].size())) {
          std::cerr << "d=" << d << ", n=" << n << ", size=" << _subdivision_coordinates[d].size() << std::endl;
          throw std::runtime_error("index does not lie in range of finite grid");
        }
        return _subdivision_coordinates[d][n];
      }
  
      /*! \brief The index of interval in dimension \a d index 
       * containing \a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const {
        typename std::vector<R>::const_iterator pos;
        assert(d<this->dimension());
        if(x<_subdivision_coordinates[d].front() || x>_subdivision_coordinates[d].back()) {
          std::cerr << "d.front()=" << _subdivision_coordinates[d].front()
                    <<  " d.back()=" << _subdivision_coordinates[d].back()
                    <<  " x="<< x <<std::endl<<std::flush;
          throw std::runtime_error("point does not lie in extent of finite grid");
        }
        pos = std::upper_bound(_subdivision_coordinates[d].begin(),
                               _subdivision_coordinates[d].end(), x);
        return (pos - _subdivision_coordinates[d].begin()) - 1;
      }

      /*! \brief Tests whether the grid contains the given lattice rectangle 
       * within its bounds. */
      virtual bool bounds_enclose(const Rectangle<R>& r) const {
        return subset(r,Rectangle<R>(this->bounds())); }

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
      IndexArray lower() const { return IndexArray(dimension(),0); }
      /*! \brief The highest valid vertex index. */
      IndexArray upper() const { 
        IndexArray result(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          result[i]=_subdivision_coordinates[i].size()-1;
        }
        return result;
      }
      /*! \brief The block of valid lattice cells. */
      Combinatoric::LatticeRectangle block() const { return Combinatoric::LatticeRectangle(lower(),upper()); }
      /*! \brief The number of subdivision intervals in each dimension. */
      SizeArray sizes() const { return block().sizes(); }
      /*! \brief The total number of cells. */
      size_type capacity() const { return block().size(); }
      /*! \brief The number of subdivision intervals in dimension \a d. */
      size_type size(dimension_type i) const { return _subdivision_coordinates[i].size()-1; }

      /*! \brief The rectangle bounding the grid. */
      GridRectangle<R> bounds() const;
      
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
    template<typename R> class RegularGrid : public Grid<R> {
     public:
      /*! \brief The type of real number defining the vertices and cells of the grid. */
      typedef R real_type;

      /*! \brief Cloning operator. */
      RegularGrid<R>* clone() const;
      
      /*! \brief Return the grid type. */
      grid_type type() const {return REGULAR;}
      
      /*! \brief Construct from an array of subdivision lengths \a sl. */
      explicit RegularGrid(const array<R>& sl)
        : _subdivision_lengths(sl) { }

      /*! \brief Construct from an array of subdivision lengths \a sl. */
      RegularGrid(const dimension_type& n, const R& l)
        : _subdivision_lengths(n,l) { }

      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const { 
        return _subdivision_lengths.size(); }

      /*! \brief The coordinate of the \a n th subdivision point in 
       * dimension \a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const { 
        return _subdivision_lengths[d] * n; }

      /*! \brief The length of the subdivision in the \a d th coordinate. */
      real_type subdivision_length(dimension_type d) const { 
        return _subdivision_lengths[d]; }

      /*! \brief The index of interval in dimension \a d index containing \a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const {
        index_type result = Numeric::quotient(x,_subdivision_lengths[d]);
        return result;
      }

      /*! \brief True if \a r is contained within the region of definition of the grid. */
      virtual bool bounds_enclose(const Rectangle<R>& r) const { return true; }
     
      /*! \brief True if \a r is contained within the region of definition of the grid. */
      virtual bool bounds_superset(const Rectangle<R>& r) const { return true; }
     
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
     * \brief A finite grid, suitable for defining a GridMaskSet. (Deprecated)
     * \deprecated Use a Grid and LatticeRectangle in GridMaskSet constructor instead. 
     */
    template<typename R> class FiniteGrid {
     public:
      /*! \brief The type of real number defining the vertices and cells of the grid. */
      typedef R real_type;

      /*! \brief A finite grid generated by subdividing a rectangle. */
      FiniteGrid(const Rectangle<R>& g, const size_type& s);
        
      /*! \brief A finite grid generated by restricting a grid by specifying 
       * lattice bounds. */
      FiniteGrid(const Grid<R>& g, const Combinatoric::LatticeRectangle& b) 
        : _grid_ptr(&g), _grid_type(g.type()), _bounds(b)
      { }
      
      /*! \brief A finite grid generated by restricting a using a spacial 
       * bounding box. */
      FiniteGrid(const Grid<R>& g, const Rectangle<R>& bb) 
        : _grid_ptr(&g), _grid_type(g.type()), 
          _bounds(over_approximation(bb,g).lattice_set())
      { }
      
      /*! \brief Return the grid type. */
      const grid_type &type() const {return _grid_type;}
      
      /*! \brief The underlying grid. */
      const Grid<R>& grid() const { return *_grid_ptr; }
      
      /*! \brief The lattice bounds of the finite grid. */
      const Combinatoric::LatticeRectangle& bounds() const { return _bounds; }
     
      /*! \brief The bounding box of the finite grid. */
      Rectangle<R> bounding_box() const;

      /*! \brief The dimension of the grid. */
      dimension_type dimension() const;

      /*! \brief The coordinate of the \a n th subdivision point in 
       * dimension \a d. */
      real_type subdivision_coordinate(dimension_type d, index_type n) const {
        return _grid_ptr->subdivision_coordinate(d,n);
      }
  
      /*! \brief The index of interval in dimension \a d index 
       * containing \a x. */
      index_type subdivision_interval(dimension_type d, const real_type& x) const {
        return _grid_ptr->subdivision_interval(d,x);
      }
      
      /*! \brief The index of the subdivision point below \a x. */
      index_type subdivision_lower_index(dimension_type d, const real_type& x) const {
        return _grid_ptr->subdivision_lower_index(d,x);
      }

      /*! \brief The index of the subdivision point above \a x. */
      index_type subdivision_upper_index(dimension_type d, const real_type& x) const {
        return _grid_ptr->subdivision_upper_index(d,x);
      }
      
     private:
      const Grid<R>* _grid_ptr;
      const grid_type _grid_type; 
      Combinatoric::LatticeRectangle _bounds;
    };
   
  }
}

#endif /* _GRID_H */
