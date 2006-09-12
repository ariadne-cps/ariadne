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

#include "../combinatoric/lattice_set.h"


namespace Ariadne {
  namespace Combinatoric {

    /*! \brief A combinatorial map on a lattice. */
    class LatticeMultiMap {
     public:
     
      /*! \brief Construct a map with argumbent dimension \a arg_dim and result dimension \a res_dim,
       *  which maps all cells to the empty set. */
      explicit LatticeMultiMap(const dimension_type arg_dim, const dimension_type res_dim) 
        : _argument_dimension(arg_dim), _result_dimension(res_dim) { }
        
      /*! \brief Adjoin set \a img to the image of set \a lc. */
      void adjoin_to_image(const LatticeCell& lc, const LatticeCell& img);
      /*! \brief Adjoin set \a img to the image of set \a lc. */
      void adjoin_to_image(const LatticeCell& lc, const LatticeCellListSet& img);
        
      /*! \brief  The map applied to a cell. */
      LatticeCellListSet apply(const LatticeCell& lc) const;

      /*! \brief  The map applied to a cell. */
      LatticeCellListSet operator() (const LatticeCell& lc) const;
        
      /*! \brief  The map applied to a lattice rectangle. */
      LatticeCellListSet operator() (const LatticeRectangle& lr) const;
        
      /*! \brief  The map applied to a cell list set. */
      LatticeCellListSet operator() (const LatticeCellListSet& lcsl) const;
        
      /*! \brief  The map applied to a cell list set. */
      LatticeCellListSet operator() (const LatticeRectangleListSet& lrls) const;
        
      /*! \brief  The map applied to a cell list set. */
      LatticeCellListSet operator() (const LatticeMaskSet& lms) const;

      /*! \brief  The dimension of the argument. */
      dimension_type argument_dimension() const {
        return _argument_dimension;
      }
      
      /*! \brief The dimension of the result. */
      dimension_type result_dimension() const {
        return _result_dimension;
      }
      
      /*! brief  The inverse of the map. */
      LatticeMultiMap inverse() const;
      
      /*! \brief The name of the system. */
      std::string name() const { return "LatticeMultiMap"; }
     private:
      friend std::ostream& operator<<(std::ostream&, const LatticeMultiMap&);
     private:
      dimension_type _argument_dimension;
      dimension_type _result_dimension;
      mutable std::map<LatticeCell,LatticeCellListSet> _map;
    };
      
    std::ostream& operator<<(std::ostream&, const LatticeMultiMap&);
    
  }
}


#endif /* _ARIADNE_LATTICE_MAP_H */
