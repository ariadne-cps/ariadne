/***************************************************************************
 *            grid.h
 *
 *  18 January 2005
 *  Copyright  2004,2005  Alberto Casagrande, Pieter Collins
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

/*! \file grid.h
 *  \brief Cuboidal grids.
 */

#ifndef _GRID_H
#define _GRID_H

#include <algorithm>
#include <vector>
#include <iostream>

#include "array.h"
#include "sequence.h"
#include "utility.h"

#include "grid_operations.h"

#include "list_set.h"

namespace Ariadne {
  int quotient(double x, double y) {
    assert(y!=0.0);
    int q = int(x/y+0.5);
    if(x < q*y) { q=q-1; }
    assert(q*y <= x && x < (q+1)*y);
    return q;
  }

  int quotient(Rational x, Rational y) {
    assert(y!=Rational(0));
    Rational d = x/y;
    Integer qz = d.get_num() / d.get_den();
    int q = int(qz.get_si());
    assert(q*y <= x && x < (q+1)*y);
    return q;
  }

  int quotient(Dyadic x, Dyadic y) {
    assert(y!=Dyadic(0));
    Dyadic d = x/y + Dyadic(0.5);
    Dyadic qd = floor(d);
    int q = int(qd.get_si());
    if(x < q*y) { q=q-1; }
    assert(q*y <= x && x < (q+1)*y);
    return q;
  }

  namespace Geometry {
    template<typename R> class Rectangle;
    template<typename R, template<typename> class BS> class ListSet;

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

      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const = 0;
      /*! \brief The coordinate of the @a n th subdivision point in dimension @a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const = 0;
      /*! \brief The interval @math [p_n,p_{n+1}] in dimension @a d index containing @a x. */
     virtual index_type subdivision_interval(dimension_type d, const real_type& x) const = 0;

      /*! \brief The index in dimension @a d of the subdivision point @a x. Throws an exception if @a x is not a subdivision point. */
      index_type subdivision_index(dimension_type d, const real_type& x) const {
        index_type n=subdivision_interval(d,x);
        if(subdivision_coordinate(d,n) == x) { return n; }
        throw std::domain_error("Value is not a grid coordinate");
      }

      /*! \brief The index of the subdivision point below x. */
      index_type subdivision_lower_index(dimension_type d, const real_type& x) const {
        return subdivision_interval(d,x);
      }

      /*! \brief The index of the subdivision point above x. */
      index_type subdivision_upper_index(dimension_type d, const real_type& x) const {
        index_type n=subdivision_interval(d,x);
        return subdivision_coordinate(d,n) == x ? n : n+1;
      }

      bool operator==(const Grid<R>& g) const { return this==&g; }

      IndexArray lower_index(const Rectangle<R>& r) const {
        IndexArray res(r.dimension());
        for(size_type i=0; i!=res.size(); ++i) {
          res[i]=subdivision_lower_index(i,r.lower(i));
        }
        return res;
      }

      IndexArray upper_index(const Rectangle<R>& r) const {
        IndexArray res(r.dimension());
        for(size_type i=0; i!=res.size(); ++i) {
          res[i]=subdivision_upper_index(i,r.upper(i));
        }
        return res;
      }

      virtual std::ostream& write(std::ostream& os) const = 0;
   };


    /*!\brief A finite, nonuniform grid of rectangles in Euclidean space.
     */
    template<typename R>
    class FiniteGrid : public Grid<R> {
      typedef R real_type;
     public:
      /*! \brief Destructor. */
      virtual ~FiniteGrid() { }

      /*! \brief Construct from a list of subdivision coordinates in each dimension. */
      explicit FiniteGrid(const array< std::vector<R> >& sp);

      /*! \brief Construct from a list of rectangles giving the grid points. */
      explicit FiniteGrid(const ListSet<R,Rectangle>& ls);

      /*! \brief Join two finite grids. */
      FiniteGrid(const FiniteGrid& g1, FiniteGrid& g2);

      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const { return _subdivision_coordinates.size(); }
      /*! \brief The total number of cells. */
//      size_type size() const { return _strides[dimension()]; }
      /*! \brief The number of subdivision intervals in dimension @a d. */
      size_type size(dimension_type d) const { return _subdivision_coordinates[d].size(); }

//      const array< std::vector<R> >& subdivision_coordinates() const { return _subdivision_coordinates; }
//      const std::vector<R>& subdivision_coordinates(dimension_type d) const { return _subdivision_coordinates[d]; }
      /*! \brief The coordinate of the @a n th subdivision point in dimension @a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const {
        return _subdivision_coordinates[d][n];
      }
      /*! \brief The index of interval in dimension @a d index containing @a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const {
        typename std::vector<R>::const_iterator pos;
        pos = std::upper_bound(_subdivision_coordinates[d].begin(),
                                    _subdivision_coordinates[d].end(), x);
        return (pos - _subdivision_coordinates[d].begin()) - 1;
      }

      /*! \brief Find the rule to translate elements from a grid to a refinement. */
      static array< std::vector<index_type> > index_translation(const FiniteGrid<R>& from, const FiniteGrid<R>& to);

      virtual std::ostream& write(std::ostream& os) const {
        return os << *this;
      }
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
      /*! \brief Construct from an array of subdivision lengths \a sl.
       */
      InfiniteGrid(const array<R>& sl)
        : _subdivision_lengths(sl) { }

      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const { return _subdivision_lengths.size(); }

//      const array<R>& subdivision_lengths() const { return _subdivision_lengths; }
      /*! \brief The length of the subdivision in the \a d th coordinate. */
      real_type subdivision_length(dimension_type d) const { return _subdivision_lengths[d]; }
      /*! \brief The coordinate of the @a n th subdivision point in dimension @a d. */
      virtual real_type subdivision_coordinate(dimension_type d, size_type n) const { return _subdivision_lengths[d] * n; }
      /*! \brief The index of interval in dimension @a d index containing @a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const {
        index_type result = quotient(x,_subdivision_lengths[d]);
        return result;
      }
     private:
      array<R> _subdivision_lengths;
    };

    template<typename R>
    FiniteGrid<R>::FiniteGrid(const array< std::vector<R> >& sp)
      : _subdivision_coordinates(sp)
    {
      create();
    }

    template<typename R>
    FiniteGrid<R>::FiniteGrid(const ListSet<R,Rectangle>& ls)
      : _subdivision_coordinates(ls.dimension())
    {
      for(typename ListSet<R,Rectangle>::const_iterator riter=ls.begin(); riter!=ls.end(); ++riter) {
        for(dimension_type n=0; n!=ls.dimension(); ++n) {
          _subdivision_coordinates[n].push_back(riter->lower(n));
          _subdivision_coordinates[n].push_back(riter->upper(n));
        }
      }
      create();
    }

    template<typename R>
    FiniteGrid<R>::FiniteGrid(const FiniteGrid<R>& g1, FiniteGrid<R>& g2)
      : _subdivision_coordinates(g1.dimension())
    {
      for(dimension_type d=0; d!=dimension(); ++d) {
        std::vector<R>& sc(_subdivision_coordinates[d]);
        const std::vector<R>& sc1(g1._subdivision_coordinates[d]);
        const std::vector<R>& sc2(g2._subdivision_coordinates[d]);
        sc.resize(sc1.size()+sc2.size());
        std::merge(sc1.begin(),sc1.end(),sc2.begin(),sc2.end(),sc.begin());
      }
      create();
    }

    template<typename R>
    void
    FiniteGrid<R>::create()
    {
       for(dimension_type i=0; i!=dimension(); ++i) {
        std::vector<R>& pos=_subdivision_coordinates[i];
        std::sort(pos.begin(),pos.end());
        typename std::vector<R>::iterator newend=std::unique(pos.begin(),pos.end());
        pos.resize(std::distance(pos.begin(),newend));
      }
    }

    template<typename R>
    array< std::vector<index_type> >
    FiniteGrid<R>::index_translation(const FiniteGrid<R>& from, const FiniteGrid<R>& to)
    {
      assert(from.dimension()==to.dimension());
      array< std::vector<index_type> > result(from.dimension());
      for(dimension_type d=0; d!=from.dimension(); ++d) {
        for(size_type n=0; n!=from.size(d); ++n) {
          index_type i=to.subdivision_index(d,from.subdivision_coordinate(d,n));
          result[d].push_back(i);
        }
      }
      return result;
    }


   /*!\brief A partition tree grid of rectangles in Euclidean space.
     */
    template<typename R> class PartitionGrid {
      typedef R real_type;
    public:
      /*! \brief Construct from a bounding box \a bb and a sequence of subdivision coordinates \a sc. */
      PartitionGrid(const Rectangle<R>& bb, const sequence<dimension_type>& sc)
        : _bounding_box(bb), _subdivision_coordinates(sc) { }

      /*! \brief The underlying dimension of the grid. */
      dimension_type dimension() const { return _bounding_box.dimension(); }

      /*! \brief The outer bounding box of the grid. */
      const Rectangle<R>& bounding_box() const { return _bounding_box; }
      const sequence<dimension_type>& subdivision_coordinates() const { return _subdivision_coordinates; }
      /*! \brief The coordinate in which the \a n th subdivision is performed. */
      dimension_type subdivision_coordinate(size_type n) const { return _subdivision_coordinates[n]; }
     private:
      Rectangle<R> _bounding_box;
      sequence<dimension_type> _subdivision_coordinates;
    };




    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const Grid<R>& g)
    {
      return g.write(os);
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const FiniteGrid<R>& g)
    {
      os << "[ ";
      for(dimension_type n=0; n!=g.dimension(); ++n) {
        if(n!=0) { os << " x "; }
        os << "[";
        for(size_type i=0; i!=g.size(n); ++i) {
          if(i!=0) { os << ", "; }
          os << g.subdivision_coordinate(n,i);
        }
        os << "]";
      }
      os << " ]";
      return os;
    }

  }
}

#endif /* _GRID_H */
