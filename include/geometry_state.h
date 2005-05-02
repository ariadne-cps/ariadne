/***************************************************************************
 *            geometry_state.h
 *
 *  Sun Jan 23 18:00:21 2005
 *  Copyright  2005  Alberto Casagrande
 *  Email casagrande@dimi.uniud.it
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
 
#ifndef _GEOMETRY_STATE_H
#define _GEOMETRY_STATE_H

#include <vector>
#include <iostream>
#include <stdexcept>

#include "utility.h"
#include "linear_algebra.h"

namespace Ariadne {
namespace Geometry {

template <typename R = AriadneRationalType> class AriadneState;

template <typename R> 
class AriadneState {
	public:	
		typedef R Real;
	
	private:
		/*! \brief The vector defining the state */
		 boost::numeric::ublas::vector<Real> _vector;
	
	public:
		AriadneState() : _vector(0) {
		}
		
		AriadneState(size_t dim) : _vector(dim) {
			if (dim<0) 
				throw std::domain_error("States should have at least dimension 0.");
		}
		
		AriadneState(size_t dim, const Real default_value) : 
			_vector(dim) {
			if (dim<0) 
				throw std::domain_error("States should have at least dimension 0.");
			

			for (size_t i=0; i<dim; i++) 
				this->_vector[i]=default_value;
		}

		AriadneState(const AriadneState& original){
			_vector=original._vector;
		}
		
		inline Real &operator[](size_t index) {
			
			if (((this->_vector).size() <= index)||(index<0))
				throw std::out_of_range("Out of the vector's range.");
			
			return  (this->_vector[index]);
		}
	
		inline Real operator[](size_t index) const{
			
			if (((this->_vector).size() <= index)||(index<0))
				throw std::out_of_range("Out of the vector's range.");
			
			return  (this->_vector[index]);
		}

		
		/*! \brief Checks equivalence between two states. */
		inline bool operator==(const AriadneState<Real> &A) const {
			
		    /* Return false if states have different dimensions */

		    if (this->dim()!=A.dim()) { 
			return false; 
		    }
			// throw std::domain_error("The parameters have different space dimentions");
			
			/* for each dimension i */
			for (size_t i=0; i<this->dim(); i++) {
				
				if (this->_vector[i]!=A._vector[i])  return false;
			}
			
			return true;
		}
		
		/*! \brief Checks equivalence between two states. */
		inline bool operator!=(const AriadneState<Real> &A) const {
		    return !( *this == A );
		}

		/*! \brief Returns the state's dimention. */
		inline size_t dim() const {
			return (this->_vector).size();
		}
		
		inline AriadneState<Real> &operator=(const AriadneState<Real> &A){
			
			this->_vector=A._vector;
			
		  	return *this;	
		}
		
		template <typename RType>
		friend std::ostream& operator<<(std::ostream &os, const AriadneState<RType> &state);

		template <typename RType>
                friend std::istream& operator>>(std::istream &is, AriadneState<RType> &state);

};

template <typename R> 
std::ostream& operator<< (std::ostream &os, const AriadneState<R> &state)
{
	os << "[";
	if(state.dim() > 0) {
	    os << state._vector[0] ;
			
	    for (size_t i=1; i<state.dim(); i++) {			
		os << ", " << state._vector[i];
	    }
	}
	os << "]" ;

	return os;
}

template <typename R> 
std::istream& operator>> (std::istream &is, AriadneState<R> &state)
{
    static size_t last_size;
    
    std::vector<R> v;
    v.reserve(last_size);
    Ariadne::operator>> (is,v);
    // is >> v;
    last_size = v.size();

    state._vector.resize(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
	state._vector[i]=v[i];
    }
    return is;
}

}
}

#endif /* _GEOMETRY_STATE_H */
