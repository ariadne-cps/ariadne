/***************************************************************************
 *            grid.tpl
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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

#include "../base/arithmetic.h"

#include "../geometry/grid.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"

namespace Ariadne {
  namespace Geometry {

    template<typename R>
    FiniteGrid<R>*
    FiniteGrid<R>::clone() const
    {
      std::cerr << "WARNING: Cloning FiniteGrid<R> causes memory leak" << std::endl;
      return new FiniteGrid<R>(*this);
    }

    template<typename R>
    FiniteGrid<R>::~FiniteGrid()
    {
    }

    template<typename R>
    FiniteGrid<R>::FiniteGrid(const Rectangle<R>& r, size_type n)
      : _subdivision_coordinates(r.dimension())
    {
      for(dimension_type i=0; i!=r.dimension(); ++i) {
        R lower(r.lower_bound(i));
        R upper(r.upper_bound(i));
        R step=(upper-lower)/n;
        _subdivision_coordinates[i].push_back(lower);
        for(size_type j=1; j!=n; ++j) {
          _subdivision_coordinates[i].push_back(lower+j*step);
        }
        _subdivision_coordinates[i].push_back(upper);
      }
      create();
    }

    template<typename R>
    FiniteGrid<R>::FiniteGrid(const Rectangle<R>& r, SizeArray sz)
      : _subdivision_coordinates(r.dimension())
    {
      for(dimension_type i=0; i!=r.dimension(); ++i) {
        R lower(r.lower_bound(i));
        R upper(r.upper_bound(i));
        R step=(upper-lower)/sz[i];
        _subdivision_coordinates[i].push_back(lower);
        for(size_type j=1; j!=sz[i]; ++j) {
          _subdivision_coordinates[i].push_back(lower+j*step);
        }
        _subdivision_coordinates[i].push_back(upper);
      }
      create();
    }

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
          _subdivision_coordinates[n].push_back(riter->lower_bound(n));
          _subdivision_coordinates[n].push_back(riter->upper_bound(n));
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
      _strides=compute_strides(sizes());
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


    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const Grid<R>& g)
    {
      return g.write(os);
    }

    template<typename R>
    std::ostream&
    FiniteGrid<R>::write(std::ostream& os) const
    {
      return os << "FiniteGrid(" << this->_subdivision_coordinates << ")";
    }

    template<typename R>
    std::ostream&
    InfiniteGrid<R>::write(std::ostream& os) const
    {
        return os << "InfiniteGrid";
    }

  }
}
