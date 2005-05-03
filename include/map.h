/***************************************************************************
 *            map.h
 *
 *  Wed Feb  2 18:33:10 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
#ifndef _MAP_H
#define _MAP_H

#include <denotable_set.h>
#include <linear_algebra.h>
#include <rectangle.h>

namespace Ariadne {	
namespace Map{
	
enum MapKind {
	LINEAR,
	AFFINE,
	GENERAL
};

enum MapResultKind {
	SINGLE_VALUE,
	MULTI_VALUE
};

}
}

#include "affine_map.h"

namespace Ariadne {	
namespace Map{
	
template <typename M>
class Map {
	
	public:
		typedef M MapT;
		typedef typename Map::DenotableSet DenotableSet;
		typedef typename DenotableSet::BasicSet BasicSet;
		typedef typename BasicSet::State State;
		typedef typename State::Real Real;
	
		Map(const MapT &T):_map(T) {}
		
		inline BasicSet operator() (const BasicSet &A) const { return _map(A); }
		
		inline DenotableSet operator() (const DenotableSet &A) const{ 
			
			DenotableSet trans_ds(A.dim());
			
			for (size_t i=0; i< A.size(); i++) {
				trans_ds.inplace_union(_map(A[i]));
			}
			
			return trans_ds;
		}
		
		inline Map<MapT> &operator=(const Map<MapT>  &A) {
			
			this->_map=A._map;
		}
		
		inline size_t dim() const {
			return (this->_map).dim();
		}
		
		inline bool invertible() const {return this->_map.invertible();}
		
	private:
	
		MapT _map;
};


/* WARNING!!!! Is it the same of an inclusion map? Here I can set
threshold in inclusion not, but the map seem similar */ 
template <typename MAP>
class ThresholdMap {
	
	public:
		typedef MAP Map;
		typedef typename Map::DenotableSet DenotableSet;
		typedef typename DenotableSet::BasicSet BasicSet;
		typedef typename BasicSet::State State;
		typedef typename State::Real Real;
	
		ThresholdMap(const Map &T, const BasicSet &threshold):_map(T), 
					_threshold(threshold) {}
		
		ThresholdMap(const Map &T, Real threshold = 0):_map(T), 
					_threshold(T.dim(), threshold) {}
						
		inline BasicSet operator() (const BasicSet &A) const { 
			return _map(A)+this->_threshold; 
		}
		
		inline void set_threshold(const Real &threshold) { 
			BasicSet new_th((this._map).dim(), threshold);
			
			this->_threshold=new_th;
		}
		
		inline DenotableSet operator() (const DenotableSet &A) const{ 
			
			DenotableSet trans_ds(A.dim());
			
			for (size_t i=0; i< A.size(); i++) {
				trans_ds.inplace_union(_map(A[i])+this->threshold);
			}
			
			return trans_ds;
		}
		
		inline ThresholdMap<MAP> &operator=(const ThresholdMap<MAP>  &A) {
			
			this->_map=A._map;
			this->_threshold=A._threshold;

		}
		
		inline size_t dim() const {
			return (this->_map).dim();
		}
		
	private:
	
		Map _map;
	
		BasicSet _threshold;
};

}

}

#endif /* _MAP_H */
