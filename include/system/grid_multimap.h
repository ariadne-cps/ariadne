/***************************************************************************
 *            grid_multimap.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
/*! \file grid_multimap.h
 *  \brief Denotable multivalued maps on grids.
 */

#ifndef ARIADNE_GRID_MULTIMAP_H
#define ARIADNE_GRID_MULTIMAP_H

#include <iosfwd>

#include <boost/iterator/iterator_adaptor.hpp>

#include "../base/tribool.h"
#include "../base/array.h"

#include "../combinatoric/lattice_map.h"

#include "../geometry/exceptions.h"
#include "../geometry/rectangle_expression.h"
#include "../geometry/grid.h"
#include "../geometry/grid_set.h"

#include "exceptions.h"

namespace Ariadne {
  namespace System {

    template<class R> class GridMultiMapIterator { };
    template<class R> class GridMultiMap;

    /*! \brief A multivalued map on grid cells.
     *  \ingroup Grid
     */
    template<class R>
    class GridMultiMap
    {
     public:
      /*! \brief The type of denotable real number defining the vertices and cells of the grid. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the set. */
      typedef Geometry::Point<R> state_type;
      /*! \brief The type of denotable point contained by the set. */
      typedef Geometry::GridCell<R> basic_set_type;
      
      typedef GridMultiMapIterator<R> const_iterator;

      /*!\brief Construct from the result grid and the argument grid. */
      GridMultiMap(const Geometry::Grid<R>& ag, const Geometry::Grid<R>& rg)
        : _argument_grid_ref(ag), _result_grid_ref(rg), _lattice_map(ag.dimension(),rg.dimension()) { }

      /*!\brief Construct from the result grid, the argument grid, and a lattice denotable multivalued map. */
      GridMultiMap(const Geometry::Grid<R>& ag, const Geometry::Grid<R>& rg, const Combinatoric::LatticeMultiMap& lmm)
        : _argument_grid_ref(ag), _result_grid_ref(rg), _lattice_map(lmm) 
      {
        ARIADNE_CHECK_RESULT_DIMENSION(lmm,ag,"GridMultiMap::GridMultiMap(Grid ag, Grid rg, LatticeMultiMap lmm)");
        ARIADNE_CHECK_RESULT_DIMENSION(lmm,rg,"GridMultiMap::GridMultiMap(Grid ag, Grid rg, LatticeMultiMap lmm)");
      }

       /*! \brief Adjoin set \a igc to the image of set \a agc. */
      void adjoin_to_image(const Geometry::GridCell<R>& agc, const Geometry::GridCell<R>& igc) {
        ARIADNE_CHECK_ARGUMENT_DIMENSION(*this,agc,"void GridMultiMap::adjoin_to_image(GridCell agc, GridCell igc)");
        ARIADNE_CHECK_RESULT_DIMENSION(*this,igc,"void GridMultiMap::adjoin_to_image(GridCell agc, GridCell igc)");
        this->_lattice_map.adjoin_to_image(agc.lattice_set(),igc.lattice_set());
      }

      /*! \brief Adjoin set \a igcls to the image of set \a agc. */
      void adjoin_to_image(const Geometry::GridCell<R>& agc, const Geometry::GridCellListSet<R>& igcls) {
        ARIADNE_CHECK_ARGUMENT_DIMENSION(*this,agc,"void GridMultiMap::adjoin_to_image(GridCell agc, GridCell igclc)");
        ARIADNE_CHECK_RESULT_DIMENSION(*this,igcls,"void GridMultiMap::adjoin_to_image(GridCell agc, GridCell igcls)");
        this->_lattice_map.adjoin_to_image(agc.lattice_set(),igcls.lattice_set());
      }
        
      /*! \brief  The map applied to a cell. */
      Geometry::GridCellListSet<R> image(const Geometry::GridCell<R>& gc) const {
        ARIADNE_CHECK_ARGUMENT_DIMENSION(*this,gc,"GridCellListSet GridMultiMap::image(GridCell agc)");
        return Geometry::GridCellListSet<R>(this->_result_grid_ref,this->_lattice_map.image(gc.lattice_set()));
      }

      /*! \brief  The map applied to a grid mask set. */
     Geometry::GridCellListSet<R> image(const Geometry::GridMaskSet<R>& gms) const {
       ARIADNE_CHECK_ARGUMENT_DIMENSION(*this,gms,"GridCellListSet GridMultiMap::image(GridMaskSet gms)");
        return Geometry::GridCellListSet<R>(this->_result_grid_ref,this->_lattice_map.image(gms.lattice_set()));
      }

      /*! \brief  The map applied to a cell. */
      Geometry::GridCellListSet<R> operator() (const Geometry::GridCell<R>& gc) const {
        ARIADNE_CHECK_ARGUMENT_DIMENSION(*this,gc,"GridCellListSet GridMultiMap::operator()(GridCell gc)");
        return Geometry::GridCellListSet<R>(this->_result_grid_ref,this->_lattice_map(gc.lattice_set()));
      }
        
      /*! \brief  The map applied to a lattice rectangle. */
      Geometry::GridCellListSet<R> operator() (const Geometry::GridBlock<R>& gb) const {
        ARIADNE_CHECK_ARGUMENT_DIMENSION(*this,gb,"GridCellListSet GridMultiMap::operator()(GridBlock gb)");
        return Geometry::GridCellListSet<R>(this->_result_grid_ref,this->_lattice_map(gb.lattice_set()));
      }
        
      /*! \brief  The map applied to a cell list set. */
      Geometry::GridCellListSet<R> operator() (const Geometry::GridCellListSet<R>& gcls) const {
        ARIADNE_CHECK_ARGUMENT_DIMENSION(*this,gcls,"GridCellListSet GridMultiMap::operator()(GridCellListSet gcls)");
        return Geometry::GridCellListSet<R>(this->_result_grid_ref,this->_lattice_map(gcls.lattice_set()));
      }
                
      /*! \brief  The map applied to a mask set. */
      Geometry::GridCellListSet<R> operator() (const Geometry::GridMaskSet<R>& gms) const {
        ARIADNE_CHECK_ARGUMENT_DIMENSION(*this,gms,"GridCellListSet GridMultiMap::operator()(GridMaskSet gms)");
        return Geometry::GridCellListSet<R>(this->_result_grid_ref,this->_lattice_map(gms.lattice_set()));
      };

      /*! \brief  The set of cells which map into \a lms. */
      Geometry::GridCellListSet<R> strong_preimage (const Geometry::GridMaskSet<R>& gms) const {
        ARIADNE_CHECK_RESULT_DIMENSION(*this,gms,"GridCellListSet GridMultiMap::strong_preimage(GridMaskSet gms)");
        return Geometry::GridCellListSet<R>(this->_argument_grid_ref,this->_lattice_map.strong_preimage(gms.lattice_set()));
      };

      /*! \brief  The set of cells which map over a cell in \a lms. */
      Geometry::GridCellListSet<R> weak_preimage (const Geometry::GridMaskSet<R>& gms) const {
        ARIADNE_CHECK_RESULT_DIMENSION(*this,gms,"GridCellListSet GridMultiMap::weak_preimage(GridMaskSet gms)");
        return Geometry::GridCellListSet<R>(this->_argument_grid_ref,this->_lattice_map.weak_preimage(gms.lattice_set()));
      };

      /*! \brief  The dimension of the argument. */
      dimension_type argument_dimension() const {
        return this->_lattice_map.argument_dimension();
      }
      
      /*! \brief The dimension of the result. */
      dimension_type result_dimension() const {
        return this->_lattice_map.result_dimension();
      }
      
      /*! brief  The inverse of the map. */
      GridMultiMap<R> inverse() const {
        return GridMultiMap<R>(this->_result_grid_ref,this->_argument_grid_ref,this->_lattice_map.inverse());
      }

      
      /*! brief  A constant iterator to the first non-default image. */
      typename GridMultiMap<R>::const_iterator begin() const;
      
      /*! brief  A constant iterator to the last non-default image. */
      typename GridMultiMap<R>::const_iterator end() const;
      
      /*! \brief The name of the system. */
      std::string name() const { return "GridMultiMap"; }
      /*! \brief The name of the system. */
      std::ostream& write(std::ostream&) const;

#ifdef DOXYGEN
      /*! \brief Write to an output stream. */
      friend std::ostream& operator<<(std::ostream&, const GridMultiMap<R>&);
#endif
     private:
      const Geometry::Grid<R>& _argument_grid_ref;
      const Geometry::Grid<R>&  _result_grid_ref;
      Combinatoric::LatticeMultiMap _lattice_map;
    };
    
    template<class R>
    std::ostream&
    GridMultiMap<R>::write(std::ostream& os) const
    {
      return os << "GridMultiMap( " << this->_argument_grid_ref
                << ", result_grid=" << this->_result_grid_ref
                << ", lattice_map=" << this->_lattice_map << " )";
    }
  
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const GridMultiMap<R>& gmm)
    {
      return gmm.write(os);
    }

  }
}


#endif /* ARIADNE_GRID_MULTIMAP_H */
