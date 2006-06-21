/***************************************************************************
 *            lattice_map.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it,  Pieter.Collins@cwi.nl
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
 
/*! \file lattice_map.h
 *  \brief Maps on a lattice, implemented by listing the image of each cell.
 */

#ifndef _ARIADNE_LATTICE_MAP_H
#define _ARIADNE_LATTICE_MAP_H

#include <map>

#include "../geometry/lattice_set.h"


namespace Ariadne {
  namespace System {

    /*! \brief A combinatorial map on a lattice. */
    class LatticeMap {
     public:
     
      explicit LatticeMap(const dimension_type arg_dim, const dimension_type res_dim) 
        : _argument_dimension(arg_dim), _result_dimension(res_dim) { }
        
      /*! \brief Adjoin set \a img to the image of set \a lc. */
      void adjoin_to_image(const Geometry::LatticeCell& lc, const Geometry::LatticeCell& img);
      /*! \brief Adjoin set \a img to the image of set \a lc. */
      void adjoin_to_image(const Geometry::LatticeCell& lc, const Geometry::LatticeCellListSet& img);
        
      /*! \brief  The map applied to a cell. */
      Geometry::LatticeCellListSet apply(const Geometry::LatticeCell& lc) const;

      /*! \brief  The map applied to a cell. */
      Geometry::LatticeCellListSet operator() (const Geometry::LatticeCell& lc) const;
        
      /*! \brief  The map applied to a lattice rectangle. */
      Geometry::LatticeCellListSet operator() (const Geometry::LatticeRectangle& lr) const;
        
      /*! \brief  The map applied to a cell list set. */
      Geometry::LatticeCellListSet operator() (const Geometry::LatticeCellListSet& lcsl) const;
        
      /*! \brief  The map applied to a cell list set. */
      Geometry::LatticeCellListSet operator() (const Geometry::LatticeRectangleListSet& lrls) const;
        
      /*! \brief  The map applied to a cell list set. */
      Geometry::LatticeCellListSet operator() (const Geometry::LatticeMaskSet& lms) const;

      /*! \brief  The dimension of the argument. */
      dimension_type argument_dimension() const {
        return _argument_dimension;
      }
      
      /*! \brief The dimension of the result. */
      dimension_type result_dimension() const {
        return _result_dimension;
      }
      
      LatticeMap inverse() const;
      
      std::string name() const { return "LatticeMap"; }
     private:
      friend std::ostream& operator<<(std::ostream&, const LatticeMap&);
     private:
      dimension_type _argument_dimension;
      dimension_type _result_dimension;
      mutable std::map<Geometry::LatticeCell,Geometry::LatticeCellListSet> _map;
    };
      
    std::ostream& operator<<(std::ostream&, const LatticeMap&);
    
  }
}


#endif /* _ARIADNE_LATTICE_MAP_H */
