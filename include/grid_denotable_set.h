/***************************************************************************
 *            grid_denotable_set.h
 *
 *  10 January 2005
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
 
/*! \file grid_denotable_set.h
 *  \brief Denotable sets.
 */

#ifndef _GRID_DENOTABLE_SET_H
#define _GRID_DENOTABLE_SET_H

#include "interval.h"
#include "binary_word.h"
#include "rectangle.h"
#include "state.h"
#include "grid.h"
#include "grid_rectangle.h"

namespace Ariadne {	
  namespace Geometry {	
    template<typename R> class Rectangle;
    template<typename R> class State;
    
    /*! \brief A denotable set on a finite grid, defined using a mask.
     */
    template<typename R>
    class FiniteGridMaskDenotableSet {
     public:
      typedef R real_type;
      typedef State<R> state_type;
   
      /*!\brief Construct from a grid and an integer array.
       */
      FiniteGridMaskDenotableSet(const FiniteGrid<R>& g)
	: _grid(g), _mask(g.size()) { } 
      
      void increment(array<size_type>& index) {
	while(true) {
	  ++index[i];
	  if(index[i] != _grid.size(i)) {
	    return; 
	  }
	  if(i == (_grid.dimension()-1)) {
	    return;
	  }
	  index[i]=0;
	  ++i;
	}
      }
     private:
      const FiniteGrid<R>& _grid;
      std::vector<bool> _mask;
    };

    std::vector<bool>
    operator&(const std::vector<bool>& v1, const std::vector<bool>& v2) 
    {
      std::vector<bool> result(v1.size());

      assert(v1.size()==v2.size());
      typedef std::vector<bool>::const_iterator const_iterator;
      typedef std::vector<bool>::iterator iterator;
      
      iterator res_iter=result.begin();
      const_iterator v1_iter=v1.begin();
      const_iterator v2_iter=v2.begin();
      iterator res_end=result.end();
      
      while(res_iter!=res_end) {
	(*res_iter) = ( (*v1_iter) & (*v2_iter) );
	++res_iter;
	++v1_iter;
	++v2_iter;
      }
    }

    std::vector<bool> 
    operator|(const std::vector<bool>& v1, const std::vector<bool>& v2) 
    {
      std::vector<bool> result(v1.size());

      assert(v1.size()==v2.size());
      typedef std::vector<bool>::const_iterator const_iterator;
      typedef std::vector<bool>::iterator iterator;

      iterator res_iter=result.begin();
      const_iterator v1_iter=v1.begin();
      const_iterator v2_iter=v2.begin();
      iterator res_end=result.end();
      
      while(res_iter!=res_end) {
	(*res_iter) = ( (*v1_iter) & (*v2_iter) );
	++res_iter;
	++v1_iter;
	++v2_iter;
      }
    }


  }
}

#endif /* _GRID_DENOTABLE_SET_H */
