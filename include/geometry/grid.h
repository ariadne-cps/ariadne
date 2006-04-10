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
#include <iostream>
#include <iterator>

#include "../declarations.h"

#include "../base/array.h"
#include "../numeric/interval.h"
#include "../base/binary_word.h"


#include "../geometry/lattice_set.h"

namespace Ariadne {
  namespace Geometry {
    template<typename R> class Grid;
    template<typename R> class FiniteGrid;
    template<typename R> class UniformGrid;
    
    template<typename R> class GridCell;
    template<typename R> class GridRectangle;

    template<typename R> std::ostream& operator<<(std::ostream&, const Grid<R>&);
    
    /*!\brief Base type for defining a grid.
     * We use inheritence and abstract functions here partly for ease of development
     * and partly since the grid coordinates only play a role when converting to rectangles
     * and should occur with linear complexity in the space dimension.
     */
    template<typename R>
    class Grid {
    public:
      typedef R real_type;
     
      /*! \brief Destructor. */
      virtual ~Grid() { }

      /*! \brief Cloning operator. */
      virtual Grid<R>* clone() const = 0;
      
      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const = 0;
      /*! \brief The coordinate of the \a n th subdivision point in dimension \a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const = 0;
      /*! \brief The interval \f$[p_n,p_{n+1})\f$ in dimension \a d index containing \a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const = 0;

      /*! \brief Tests whether the grid contains the given lattice rectangle within its bounds. */
      virtual bool contains(const LatticeRectangle& r) const = 0;

      /*! \brief The index in dimension \a d of the subdivision point \a x. Throws an exception if \a x is not a subdivision point. */
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

      bool operator==(const Grid<R>& g) const { return this==&g; }
      bool operator!=(const Grid<R>& g) const { return !(*this==g); }

      IndexArray index(const Point<R>& s) const;
      IndexArray lower_index(const Rectangle<R>& r) const;
      IndexArray upper_index(const Rectangle<R>& r) const;

      virtual std::ostream& write(std::ostream& os) const = 0;
   };


    /*!\brief A finite, nonuniform grid of rectangles in Euclidean space. */
    template<typename R>
    class FiniteGrid : public Grid<R> {
      typedef R real_type;
     public:
      /*! \brief Cloning operator. */
      virtual FiniteGrid<R>* clone() const;
    
      /*! \brief Destructor. */
      virtual ~FiniteGrid();

      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const { 
        return this->_subdivision_coordinates.size();
      }
    
      /*! \brief The coordinate of the \a n th subdivision point in dimension \a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const {
        return _subdivision_coordinates[d][n];
      }
  
      /*! \brief The index of interval in dimension \a d index containing \a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const {
        typename std::vector<R>::const_iterator pos;
        if(x<_subdivision_coordinates[d].front() || x>_subdivision_coordinates[d].back()) {
          throw std::runtime_error("point does not lie in extent of finite grid");
        }
        pos = std::upper_bound(_subdivision_coordinates[d].begin(),
                               _subdivision_coordinates[d].end(), x);
        return (pos - _subdivision_coordinates[d].begin()) - 1;
      }

      /*! \brief Tests whether the grid contains the given lattice rectangle within its bounds. */
      virtual bool contains(const LatticeRectangle& r) const {
        return subset(r,this->bounds()); }

      /*! \brief Construct from a bounding box and an equal number of subdivisions in each coordinate. */
      explicit FiniteGrid(const Rectangle<R>& r, size_type n);

      /*! \brief Construct from a bounding box and an array giving the number of subdivisions in each coordinate. */
      explicit FiniteGrid(const Rectangle<R>& r, SizeArray sz);

      /*! \brief Construct from a list of subdivision coordinates in each dimension. */
      explicit FiniteGrid(const array< std::vector<R> >& sp);

      /*! \brief Construct from a list of rectangles giving the grid points. */
      explicit FiniteGrid(const ListSet<R,Rectangle>& ls);

      /*! \brief Join two finite grids. */
      FiniteGrid(const FiniteGrid& g1, FiniteGrid& g2);

      bool operator==(const FiniteGrid<R>& g) const { 
        return this->_subdivision_coordinates==g._subdivision_coordinates; }
            
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
      LatticeRectangle bounds() const { return LatticeRectangle(lower(),upper()); }
      /*! \brief The number of subdivision intervals in each dimension. */
      SizeArray sizes() const { return bounds().sizes(); }
      /*! \brief The total number of cells. */
      size_type capacity() const { return bounds().size(); }
      /*! \brief The number of subdivision intervals in dimension \a d. */
      size_type size(dimension_type i) const { return _subdivision_coordinates[i].size()-1; }

      /*! \brief The rectangle bounding the grid. */
      GridRectangle<R> bounding_box() const;
      
      /*! \brief Find the rule to translate elements from a grid to a refinement. */
      static array< std::vector<index_type> > index_translation(const FiniteGrid<R>& from, const FiniteGrid<R>& to);

      virtual std::ostream& write(std::ostream& os) const;
     private:
      void create();
     private:
      array< std::vector<R> > _subdivision_coordinates;
    };



    /*!\brief An infinite, uniform grid of rectangles in Euclidean space.
     */
    template<typename R> class InfiniteGrid : public Grid<R> {
      typedef R real_type;
     public:
      /*! \brief Cloning operator. */
      InfiniteGrid<R>* clone() const;
      
      /*! \brief Construct from an array of subdivision lengths \a sl.
       */
      explicit InfiniteGrid(const array<R>& sl)
        : _subdivision_lengths(sl) { }

      /*! \brief Construct from an array of subdivision lengths \a sl.
       */
      InfiniteGrid(const dimension_type& n, const R& l)
        : _subdivision_lengths(n,l) { }

      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const { 
        return _subdivision_lengths.size(); }

      /*! \brief The coordinate of the \a n th subdivision point in dimension \a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const { 
        return _subdivision_lengths[d] * n; }

      /*! \brief The length of the subdivision in the \a d th coordinate. */
      real_type subdivision_length(dimension_type d) const { 
        return _subdivision_lengths[d]; }

      /*! \brief The index of interval in dimension \a d index containing \a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const {
        index_type result = quotient(x,_subdivision_lengths[d]);
        return result;
      }

      virtual bool contains(const LatticeRectangle& r) const { return true; }

      virtual std::ostream& write(std::ostream& os) const;
     private:
      array<R> _subdivision_lengths;
    };

  }
}

#endif /* _GRID_H */
