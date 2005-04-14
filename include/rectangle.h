/***************************************************************************
 *            rectangle.h
 *
 *  Sun Jan 23 15:31:54 2005
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
 
#ifndef _RECTANGLE_H
#define _RECTANGLE_H

#include<list>

#include <geometry_state.h>

namespace Ariadne {	
namespace Geometry {

template <typename S = AriadneState<AriadneRationalType> > 
class AriadneRectangle;

namespace IO_Operators{	

template <typename S>
std::ostream& operator<<(std::ostream &os, const AriadneRectangle<S> &r) {
			
	os << "[ " << (r._lower_corner) << " , " ;

	os << (r._upper_corner) << " ]" ;
			
	return os;
}

}

/*! \brief Tests disjointness */
template <typename S>
bool disjoint(const AriadneRectangle<S> &A, const AriadneRectangle<S> &B) {
	
	if (A.dim()!=B.dim()) 
		throw std::domain_error("The two parameters have different space dimentions");
			

	for (size_t i=0; i< A.dim(); i++) {
		if ((A._upper_corner[i]<B._lower_corner[i])|| 
			(B._upper_corner[i]<A._lower_corner[i])) return true;
	}
		
	return false;
}

/*! \brief Tests intersection of interiors */
template <typename S>
bool intersects_interior(const AriadneRectangle<S> &A, 
		const AriadneRectangle<S> &B) {
			
	if (A.dim()!=B.dim()) 
		throw std::domain_error("The two parameters have different space dimentions");

	if (A.empty()||B.empty()) return false;
	
	for (size_t i=0; i< A.dim(); i++) {
		if ((A._upper_corner[i]<=B._lower_corner[i])|| 
			(B._upper_corner[i]<=A._lower_corner[i])) return false;
	}
		
	return true;
}



/*! \brief Tests inclusion of interiors */
template <typename S>
bool subset_of_interior(const AriadneRectangle<S> &A, 
				const AriadneRectangle<S> &B) {
	
	if (A.dim()!=B.dim()) 
		throw std::domain_error("The two parameters have different space dimentions");
			
	if (A.empty()||B.empty()) return false;
	
	for (size_t i=0; i< A.dim(); i++) {
		if ((A._upper_corner[i] >= B._upper_corner[i])|| 
			(B._lower_corner[i] >= A._lower_corner[i])) return false;
	}
		
	return true;				
}

/*! \brief Tests intersection of interiors */
template <typename S>
bool interiors_intersect(const AriadneRectangle<S> &A, 
		const AriadneRectangle<S> &B) {
			
	if (A.dim()!=B.dim()) 
		throw std::domain_error("The two parameters have different space dimentions");

	if (A.empty()||B.empty()) return false;
	
	for (size_t i=0; i< A.dim(); i++) {
		if ((A._upper_corner[i]<=B._lower_corner[i])|| 
			(B._upper_corner[i]<=A._lower_corner[i])) return false;
	}
		
	return true;
}



/*! \brief Tests inclusion of interiors */
template <typename S>
bool interior_subset_of_interior(const AriadneRectangle<S> &A, 
				const AriadneRectangle<S> &B) {
	
	if (A.dim()!=B.dim()) 
		throw std::domain_error("The two parameters have different space dimentions");
			
	if (A.empty()||B.empty()) return false;
	
	for (size_t i=0; i< A.dim(); i++) {
		if ((A._upper_corner[i] > B._upper_corner[i])|| 
			(B._lower_corner[i] > A._lower_corner[i])) return false;
	}
		
	return true;				
}


//*! \brief Tests inclusion */
/*
template <typename S >
bool subset_of_open_cover(const AriadneRectangle<S> &A, 
				const std::list< AriadneRectangle<S> > &list) {

}
*/

/*! \brief Makes intersection of interior */
template <typename S>
AriadneRectangle<S> closure_of_intersection_of_interior(
		const AriadneRectangle<S> &A, const AriadneRectangle<S> &B){

	if (A.dim()!=B.dim()) 
		throw std::domain_error("The two parameters have different space dimentions");
	
	AriadneRectangle<S> C(A.dim());

	if (A.empty()||B.empty()) return C;
	
	for (size_t i=0; i< A.dim(); i++) {
	
		if (A._upper_corner[i] > B._upper_corner[i]) {
			C._upper_corner[i]=B._upper_corner[i];
		} else {
			C._upper_corner[i]=A._upper_corner[i];
		}
		
		if (B._lower_corner[i] > A._lower_corner[i]) {
			C._lower_corner[i]=B._lower_corner[i];
		} else {
			C._lower_corner[i]=A._lower_corner[i];
		}
	}
	
	return C;
			
}

template <typename S>
class AriadneRectangle {
	
		/*! \brief Tests disjointness */
		friend bool disjoint <> (const AriadneRectangle<S> &A, 
				const AriadneRectangle<S> &B);

		/*! \brief Tests intersection of interiors */
		friend bool intersects_interior <> (const AriadneRectangle<S> &A, 
				const AriadneRectangle<S> &B);

		/*! \brief Tests inclusion of interiors */
		friend bool subset_of_interior <> (const AriadneRectangle<S> &A, 
				const AriadneRectangle<S> &B);

		/*! \brief Tests intersection of interiors */
		friend bool interiors_intersect <> (const AriadneRectangle<S> &A, 
				const AriadneRectangle<S> &B);

		/*! \brief Tests inclusion of interiors */
		friend bool interior_subset_of_interior <> (const AriadneRectangle<S> &A, 
				const AriadneRectangle<S> &B);
		
		//*! \brief Tests inclusion */
		//friend bool subset_of_open_cover <> (const AriadneRectangle<S> &A, 
		//		const std::list< AriadneRectangle<S> > &list);

		/*! \brief Makes intersection of interiors */
		friend AriadneRectangle<S> closure_of_intersection_of_interior <> (
				const AriadneRectangle<S> &A, const AriadneRectangle<S> &B);
	
	public:
		typedef typename S::Real Real;	
		typedef S State;
	
	private:
		/*! \brief Rectangle's upper corner */
		State _upper_corner;
		
		/*! \brief Rectangle's lower corner */
		State _lower_corner;

		/*! \brief The emptyness flag */
		bool _empty;
	
	public:
		AriadneRectangle(size_t dim = 0): _upper_corner(dim), 
						  _lower_corner(dim), _empty(true) {}
			
		AriadneRectangle(const AriadneRectangle<State> &original): 
			_upper_corner(original._upper_corner),
			_lower_corner(original._lower_corner),_empty(original._empty) {}
			
		AriadneRectangle(const State &u_corner, const State &l_corner):
			_upper_corner(u_corner.dim()), _lower_corner(u_corner.dim()) {
								
			if (u_corner.dim()!=l_corner.dim()) 
				throw std::domain_error("The parameters have different space dimentions");
			
			/* for each dimension i */
			for (size_t i=0; i<this->dim(); i++) {
				
				if (l_corner[i] > u_corner[i]) {
					this->_empty=true;

					return;
				}
				
				this->_upper_corner[i]=u_corner[i];
				this->_lower_corner[i]=l_corner[i];
			}

			this->_empty=false;
									
		};
	
		/*! \brief Returns rectangle's space dimensions */
		inline size_t dim() const {
			return (this->_upper_corner).dim();
		}
		
		/*! \brief Returns $true$ if the rectangle is empty */
		inline bool empty() const {
			return (this->_empty);
		}
		
		/*! \brief Tests if a state is included into a rectangle. */
		inline bool contains(const State& state) const {
			
			if (state.dim()!=this->dim()) 
				throw std::domain_error("This object and parameter have different space dimentions");
			
			if (this->empty()) return false;
			
			/* for each dimension i */
			for (size_t i=0; i<this->dim(); i++) {
				
				/* if the i dim of the state is greater than the one of the 
				 * rectangle's upper corner then state is not contained into 
				 * this object */
				if (state[i] > this->_upper_corner[i]) return false;
					
				/* if the i dim of the state is smaller than the one of the 
				 * rectangle's lower corner then state is not contained into 
				 * this object */
				if (state[i] < this->_lower_corner[i]) return false;
			}
			
			return true;
			
		}
		
		/*! \brief Tests if a state is included into the interior a rectangle. */
		inline bool interior_contains(const State& state) const {
			
			if (state.dim()!=this->dim()) 
				throw std::domain_error("This object and parameter have different space dimentions");
			
			if (this->empty()) return false;
			
			/* for each dimension i */
			for (size_t i=0; i<this->dim(); i++) {
				
				/* if the i dim of the state is greater or equal than the one 
				 * of the rectangle's upper corner then state is not contained 
				 * in the interior of this object */
				if (state[i] >= this->_upper_corner[i]) return false;
					
				/* if the i dim of the state is smaller or equal than the one 
				 * of the rectangle's lower corner then state is not contained
				 * in the interior of this object */
				if (state[i] <= this->_lower_corner[i]) return false;
			}
			
			return true;
			
		}

		inline State upper_corner() const {
			return this->_upper_corner;
		}

		inline State lower_corner() const {
			return this->_lower_corner;
		}

		inline AriadneRectangle<State> find_quadrant(size_t q) const {
			
			size_t j;
			
			AriadneRectangle<State> quadrant(this->dim());
			
			for (j=0; j< this->dim(); j++) {
				
				if (q%2) {
					quadrant._lower_corner[j]=(this->_upper_corner[j]+
							this->_lower_corner[j])/2;
					quadrant._upper_corner[j]=this->_upper_corner[j];
				} else {
					quadrant._upper_corner[j]=(this->_upper_corner[j]+
							this->_lower_corner[j])/2;
					quadrant._lower_corner[j]=this->_lower_corner[j];
				}
				q=q/2;
			}

			return quadrant;
		}

		inline AriadneRectangle<State> &expand_by(const Real &delta) {
			
			for (size_t j=0; j< this->dim(); j++) {
				
				this->_upper_corner[j]+=delta;
				this->_lower_corner[j]-=delta;
			}
			
			return *this;
		}
		
		inline AriadneRectangle<State> &expand_by(const State &delta) {
			
			for (size_t j=0; j< this->dim(); j++) {
				
				this->_upper_corner[j]+=delta[j];
				this->_lower_corner[j]-=delta[j];
			}
			
			return *this;
		}
		
		inline bool operator!=(const AriadneRectangle<State> &A) const {
			
			if (this->dim()!= A.dim()) return true;
		
			if ((A.empty()&&!this->empty())||(!A.empty()&&this->empty())) return false;
			
			for (size_t j=0; j< this->dim(); j++) {
				
				if (this->_upper_corner[j]!=A._upper_corner[j]) return true;
				if (this->_lower_corner[j]!=A._lower_corner[j]) return true;
			}
			
			return false;
		}
		
		template <typename STATE>
		friend std::ostream& IO_Operators::operator<<(std::ostream &os, 
				const AriadneRectangle<STATE> &r);

};

}
}

#endif /* _RECTANGLE_H */
