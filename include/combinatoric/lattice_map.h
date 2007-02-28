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

#ifndef ARIADNE_LATTICE_MAP_H
#define ARIADNE_LATTICE_MAP_H

#include <map>

#include "../combinatoric/lattice_set.h"


namespace Ariadne {
  namespace Combinatoric {

    /*! \brief A combinatorial map on a lattice. */
    class LatticeMultiMap {
     public:
      typedef std::map<LatticeCell,LatticeCellListSet>::const_iterator const_iterator;
     public:
     
      /*! \brief Construct a map with argument dimension \a arg_dim and result dimension \a res_dim,
       *  which maps all cells to the empty set. */
      explicit LatticeMultiMap(const dimension_type arg_dim, const dimension_type res_dim) 
        : _argument_dimension(arg_dim), _result_dimension(res_dim) { }
        
      /*! \brief Adjoin set \a img to the image of set \a lc. */
      void adjoin_to_image(const LatticeCell& lc, const LatticeCell& img);
      /*! \brief Adjoin set \a img to the image of set \a lc. */
      void adjoin_to_image(const LatticeCell& lc, const LatticeCellListSet& img);
        
      /*! \brief  The map applied to a cell. \deprecated */
      LatticeCellListSet apply(const LatticeCell& lc) const;

      /*! \brief  The map applied to a cell. */
      LatticeCellListSet image(const LatticeCell& lc) const;

      /*! \brief  The map applied to a lattice mask set. */
      LatticeCellListSet image(const LatticeMaskSet& lc) const;

      /*! \brief  The map applied to a cell. */
      LatticeCellListSet operator() (const LatticeCell& lc) const;
        
      /*! \brief  The map applied to a lattice rectangle. */
      LatticeCellListSet operator() (const LatticeBlock& lr) const;
        
      /*! \brief  The map applied to a cell list set. */
      LatticeCellListSet operator() (const LatticeCellListSet& lcsl) const;
                
      /*! \brief  The map applied to a mask set. */
      LatticeCellListSet operator() (const LatticeMaskSet& lms) const;

      /*! \brief  The set of cells which map into \a lms. */
      LatticeCellListSet strong_preimage (const LatticeMaskSet& lms) const;

      /*! \brief  The set of cells which map over a cell in \a lms. */
      LatticeCellListSet weak_preimage (const LatticeMaskSet& lms) const;

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
      
      /*! brief  A constant iterator to the first non-default image. */
      LatticeMultiMap::const_iterator begin() const;
      
      /*! brief  A constant iterator to the last non-default image. */
      LatticeMultiMap::const_iterator end() const;
      
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
    
    /*! \brief A discrete control system on a lattice. */
    class LatticeSystem {
     public:
      typedef std::map<LatticeCell,LatticeCellListSet>::const_iterator const_iterator;
     public:
     
      /*! \brief Construct a system with space dimension \a sd and input dimension \a id,
       *  which maps all cells to the empty set. */
      explicit LatticeSystem(const dimension_type sd, const dimension_type id) 
        : _space_dimension(sd), _input_dimension(id) { }
        
      /*! \brief . */
      void set_control_values(const LatticeCell& lc, const LatticeCellListSet& img);
      /*! \brief . */
      void set_noise_values(const LatticeCell& lc, const LatticeCellListSet& img);
      /*! \brief . */
      void set_noise_values(const LatticeCell& lc, const LatticeCell& lc, const LatticeCellListSet& img);
        
      /*! \brief  The set of cells which can be controlled into lms regardless of the noise. */
      LatticeCellListSet preimage (const LatticeMaskSet& lms) const;

      /*! \brief  The dimension of the argument. */
      dimension_type space_dimension() const {
        return _space_dimension;
      }
      
      /*! \brief The dimension of the result. */
      dimension_type input_dimension() const {
        return _input_dimension;
      }
      
      /*! brief  The inverse of the map. */
      LatticeMultiMap inverse() const;
      
      /*! brief  A constant iterator to the first non-default image. */
      LatticeMultiMap::const_iterator begin() const;
      
      /*! brief  A constant iterator to the last non-default image. */
      LatticeMultiMap::const_iterator end() const;
      
      /*! \brief The name of the system. */
      std::string name() const { return "LatticeMultiMap"; }
     private:
      friend std::ostream& operator<<(std::ostream&, const LatticeMultiMap&);
     private:
      dimension_type _space_dimension;
      dimension_type _input_dimension;
      mutable std::map<LatticeCell,LatticeCellListSet> _control_map;
      mutable std::map<LatticeCell,LatticeCellListSet> _noise_map;
    };
      
    std::ostream& operator<<(std::ostream&, const LatticeSystem&);
    
  }
}


#endif /* ARIADNE_LATTICE_MAP_H */
