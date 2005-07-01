/***************************************************************************
 *            grid_rectangle.h
 *
 *  18 January 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file grid_rectangle.h
 *  \brief Rectangles defined on a grid.
 */

#ifndef _GRID_RECTANGLE_H
#define _GRID_RECTANGLE_H

#include "binary_word.h"

#include "rectangle.h"
#include "state.h"
#include "grid.h"

namespace Ariadne {
  namespace Geometry {
    template<typename R> class Rectangle;
    template<typename R> class State;
    template<typename R> class Grid;
    template<typename R> class GridCell;
    template<typename R> class GridRectangle;
    template<typename R> class GridMaskSet;
    template<typename R> class GridCellListSet;
    template<typename R> class GridRectangleListSet;

    /*! \brief A cell in a grid.
     */
    template<typename R>
    class GridCell {
      friend class GridRectangle<R>;
      friend class GridMaskSet<R>;
      friend class GridCellListSet<R>;
     public:
      typedef R real_type;
      typedef State<R> state_type;

      /*!\brief Construct from a grid and an integer array. */
      GridCell(const Grid<R>& g, const IndexArray& pos);

      /*!\brief The grid containing the cell. */
      const Grid<R>& grid() const { return _grid; }

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const { return _position.size(); }

      /*!\brief The position of the cell in the grid. */
      const IndexArray& position() const { return _position; }

      /*!\brief Convert to an ordinary rectangle. */
      operator Rectangle<R>() const;
     private:
      const Grid<R>& _grid;
      IndexArray _position;
    };


    /*! \brief A rectangle in a grid.
     */
    template<typename R>
    class GridRectangle {
      friend class GridCell<R>;
      friend class GridMaskSet<R>;
      friend class GridRectangleListSet<R>;
     public:
      typedef R real_type;
      typedef State<R> state_type;

      /*!\brief Construct from a grid and two integer arrays giving the corners. */
      GridRectangle(const Grid<R>& g, const IndexArray& l, const IndexArray& u);
      /*!\brief Construct from a grid and an ordinary rectangle. */
      GridRectangle(const Grid<R>& g, const Rectangle<R>& r);
      /*!\brief Construct from a GridCell. */
      GridRectangle(const GridCell<R>& c);

      /*!\brief The grid containing the rectangle. */
      const Grid<R>& grid() const { return _grid; }
      /*!\brief The dimension of the rectangle. */
      dimension_type dimension() const { return _lower.size(); }

      /*!\brief The position of the rectangle in the grid. */
      IntegerRectangle position() const {
        State<index_type> lst(_lower.begin(),_lower.end());
        State<index_type> ust(_upper.begin(),_upper.end());
        return IntegerRectangle(lst,ust);
      }

      /*!\brief Convert to an ordinary rectangle. */
      operator Rectangle<R>() const;
     private:
      const Grid<R>& _grid;
      array<index_type> _lower;
      array<index_type> _upper;
    };



    /*! \brief A rectangle defined on a partition tree.
     */
    template<typename R>
    class PartitionTreeCell {
     public:
      typedef R real_type;
      typedef State<R> state_type;

      /*!\brief Construct from a grid and a binary word.
       */
      PartitionTreeCell(const PartitionGrid<R>& g,
                        const BinaryWord& w)
        : _grid(g), _word(w)
      { }

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const { return _grid.dimension(); }

      /*!\brief Convert to an ordinary rectangle. */
      operator Rectangle<R>() const;
     private:
      const PartitionGrid<R>& _grid;
      BinaryWord _word;
    };



    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g, const IndexArray& l, const IndexArray& u)
      : _grid(g), _lower(l), _upper(u)
    {
      assert(g.dimension()==l.size());
      assert(l.size()==u.size());
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g, const Rectangle<R>& r)
      : _grid(g), _lower(g.dimension()), _upper(g.dimension())
    {
      assert(g.dimension()==r.dimension());
      for(dimension_type d=0; d!=dimension(); ++d) {
        /* TODO: Catch and rethrow exceptions */
        _lower[d]=g.subdivision_index(d,r.lower(d));
        _upper[d]=g.subdivision_index(d,r.upper(d));
      }
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const GridCell<R>& c)
      : _grid(c._grid), _lower(c._position), _upper(c._position)
    {
      for(dimension_type i=0; i!=_upper.size(); ++i) {
        _upper[i]+=1;
      }
    }


    template<typename R>
    GridCell<R>::GridCell(const Grid<R>& g, const IndexArray& pos)
      : _grid(g), _position(pos)
    {
      assert(_position.size()==_grid.dimension());
    }



    template<typename R>
    inline
    GridCell<R>::operator Rectangle<R>() const {
      Rectangle<R> result(dimension());

      for(dimension_type i=0; i!=dimension(); ++i) {
        result.set_lower(i, _grid.subdivision_coordinate(i,_position[i]));
        result.set_upper(i, _grid.subdivision_coordinate(i,_position[i]+1));
      }

      return result;
    }


    template<typename R>
    inline
    GridRectangle<R>::operator Rectangle<R>() const {
      Rectangle<R> result(dimension());

      for(size_type i=0; i!=dimension(); ++i) {
        result.set_lower(i, _grid.subdivision_coordinate(i,_lower[i]));
        result.set_upper(i, _grid.subdivision_coordinate(i,_upper[i]));
      }

      return result;
    }



    template<typename R>
    inline
    PartitionTreeCell<R>::operator Rectangle<R>() const {
      Rectangle<R> res(_grid.bounding_box());
      sequence<unsigned short>::const_iterator coord_iter=_grid.subdivision_coordinates().begin();
      BinaryWord::const_iterator word_iter=_word.begin();

      while(word_iter!=_word.end()) {
        size_type i=(*coord_iter);
        R centre = ( res.lower(i) + res.upper(i) ) / 2;
        if( (*word_iter)==0 ) {
          res.set_upper(i,centre); }
        else {
          res.set_lower(i,centre);
        }
      }

      return res;
    }

    template<typename R>
    inline
    std::ostream&
    operator<<(std::ostream& os, const GridRectangle<R>& r) {
      return os << Rectangle<R>(r);
    }

    template<typename R>
    inline
    std::ostream&
    operator<<(std::ostream& os, const GridCell<R>& c) {
      return os << Rectangle<R>(c);
    }


  }
}

#endif /* _GRID_RECTANGLE_H */
