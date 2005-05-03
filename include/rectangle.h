/***************************************************************************
 *            rectangle.h
 *
 *  Mon 2 May 2005
 *  Copyright 2005  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include "ariadne.h"
#include "utility.h"
#include "geometry_state.h"

namespace Ariadne {	

namespace Geometry {

template <typename S = State<Rational> > 
class Rectangle;


/*! \brief Tests disjointness */
template <typename S>
bool disjoint(const Rectangle<S> &A, const Rectangle<S> &B) {
	
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
bool intersects_interior(const Rectangle<S> &A, 
		const Rectangle<S> &B) {
			
	if (A.dim()!=B.dim()) 
		throw std::domain_error("The two parameters have different space dimensions");

	if (A.empty()||B.empty()) return false;
	
	for (size_t i=0; i< A.dim(); i++) {
		if ((A._upper_corner[i]<=B._lower_corner[i])|| 
			(B._upper_corner[i]<=A._lower_corner[i])) return false;
	}
		
	return true;
}


/*! \brief Tests inclusion of interiors */
template <typename S>
bool subset_of_interior(const Rectangle<S> &A, 
				const Rectangle<S> &B) {
	
	if (A.dim()!=B.dim()) 
		throw std::domain_error("The two parameters have different space dimentions");
			
	if (A.empty()||B.empty()) return false;
	
	for (size_t i=0; i< A.dim(); i++) {
		if ((A._upper_corner[i] >= B._upper_corner[i])|| 
			(B._lower_corner[i] >= A._lower_corner[i])) return false;
	}
		
	return true;				
}



/*! \brief Tests inclusion of interiors */
template <typename S>
bool is_subset_of(const Rectangle<S> &A, 
		  const Rectangle<S> &B) {
	
	if (A.dim()!=B.dim()) 
		throw std::domain_error("The two parameters have different space dimentions");
			
	if (A.empty()||B.empty()) return false;
	
	for (size_t i=0; i< A.dim(); i++) {
		if ((A._upper_corner[i] > B._upper_corner[i])|| 
			(B._lower_corner[i] > A._lower_corner[i])) return false;
	}
		
	return true;				
}


//*! \brief Tests inclusion in an open cover.  */
template <typename S >
bool subset_of_open_cover(const Rectangle<S> &A, 
			  const std::list< Rectangle<S> > &list);


/*! \brief Makes intersection of interior */
template <typename S>
Rectangle<S> closure_of_intersection_of_interior(
		const Rectangle<S> &A, const Rectangle<S> &B){

	if (A.dim()!=B.dim()) 
		throw std::domain_error("The two parameters have different space dimentions");
	
	Rectangle<S> C(A.dim());

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

/*! \brief A cuboid of arbitrary dimension.
 *  COMMENT: I think we should disallow empty rectangles. 
 *  This should make the algorithms faster, since we don't have to check emptyness every time.
 */
template <typename S>
class Rectangle {
	
		/*! \brief Tests disjointness */
		friend bool disjoint <> (const Rectangle<S> &A, 
				const Rectangle<S> &B);

		/*! \brief Tests intersection of interiors */
		friend bool intersects_interior <> (const Rectangle<S> &A, 
				const Rectangle<S> &B);

                /*! \brief Tests if \a A is a subset of the interior of \a B. */
		friend bool subset_of_interior <> (const Rectangle<S> &A, 
				const Rectangle<S> &B);

		/*! \brief Tests inclusion */
		friend bool is_subset_of <> (const Rectangle<S> &A, 
				const Rectangle<S> &B);
		
		//*! \brief Tests inclusion in an open cover. */
		friend bool subset_of_open_cover <> (const Rectangle<S> &A, 
						     const std::list< Rectangle<S> > &list);

		/*! \brief Makes intersection of interiors */
		friend Rectangle<S> closure_of_intersection_of_interior <> (
				const Rectangle<S> &A, const Rectangle<S> &B);
	
	public:
                /*! \brief The type of denotable real number used for the corners. */
                typedef typename S::Real Real;	
                /*! \brief The type of denotable state contained by the rectangle. */
		typedef S State;
	private:
                typedef Ariadne::interval<Real> Interval;
	private:
		/* Rectangle's lower corner */
		State _lower_corner;

		/* Rectangle's upper corner */
		State _upper_corner;
		
		/* The emptyness flag */
		bool _empty;
	
	public:
                /*! \brief Default constructor construcs an empty rectangle of dimension 0. */
		Rectangle(size_t dim = 0)
		    : _lower_corner(dim),  _upper_corner(dim), _empty(true) {}
    
                /*! \brief Copy constructor. */
		Rectangle(const Rectangle<State> &original) 
		    : _lower_corner(original._lower_corner),
		      _upper_corner(original._upper_corner),
		      _empty(original._empty) {}
			
                /*! \brief Construct from an array of intervals. */
		Rectangle(size_t dim, const Interval* intervals) 
		    : _lower_corner(dim), _upper_corner(dim), _empty(true) 
                {
		    for(size_t i=0; i!=dim; ++i) {
			if(intervals[i].lower() >= intervals[i].upper()) {
			    this->_empty=true;
			}
		    }
		    this->_empty=false;

		    for(size_t i=0; i!=dim; ++i) {
			this->_lower_corner[i] = intervals[i].lower();
			this->_upper_corner[i] = intervals[i].upper();
		    }    
		}

                /*! \brief Construct from lower and upper corners. */
                Rectangle(const State &l_corner, const State &u_corner):
		    _lower_corner(l_corner.dim()), _upper_corner(l_corner.dim()), _empty(true) 
                {
		    /* Test to see if corners have same dimensions */
		    if (l_corner.dim()!=u_corner.dim()) {
				throw std::domain_error("The parameters have different space dimensions");
		    }
			
		    /* Test for emptyness */
		    for (size_t i=0; i<this->dim(); i++) {
			if (l_corner[i] > u_corner[i]) {
			    this->_empty=true;
			    return;
			}
		    }
		    this->_empty=false;
				
		    /* Set coordinates */
		    for (size_t i=0; i<this->dim(); i++) {
			this->_lower_corner[i]=l_corner[i];
			this->_upper_corner[i]=u_corner[i];
		    }
		    
		}
	
		/*! \brief Returns rectangle's space dimensions */
		inline size_t dim() const {
			return (this->_lower_corner).dim();
		}
		
		/*! \brief Returns \c true if the rectangle is empty */
		inline bool empty() const {
			return (this->_empty);
		}
		
                /*! \brief The lower corner. */
		inline State lower_corner() const {
			return this->_lower_corner;
		}

                /*! \brief The upper corner. */
		inline State upper_corner() const {
			return this->_upper_corner;
		}

                /*! \brief Returns the projection onto the \a n th coordinate. */
                inline interval<Real> interval(size_t n) const {
		    return Interval(this->_lower_corner[n],this->_upper_corner[n]);
		}

                /*! \brief Returns the projection onto the \a n th coordinate. */
                inline Interval operator[] (size_t n) const {
		    return Interval(this->_lower_corner[n],this->_upper_corner[n]);
		}

                /*! \brief Returns the lower bound of the \a n th coordinate */
                inline Real lower(size_t n) const {
		    return this->_lower_corner[n];
		}

                /*! \brief Returns the upper bound of the \a n th coordinate */
                inline Real upper(size_t n) const {
		    return this->_upper_corner[n];
		}

		/*! \brief Tests if \a state is included into a rectangle. */
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
		
		/*! \brief Tests if \a state is included into the interior a rectangle. */
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

		inline Rectangle<State> find_quadrant(size_t q) const {
			
			size_t j;
			
			Rectangle<State> quadrant(this->dim());
			
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

		inline Rectangle<State> &expand_by(const Real &delta) {
			
			for (size_t j=0; j< this->dim(); j++) {
				
				this->_upper_corner[j]+=delta;
				this->_lower_corner[j]-=delta;
			}
			
			return *this;
		}
		
		inline Rectangle<State> &expand_by(const State &delta) {
			
			for (size_t j=0; j< this->dim(); j++) {
				
				this->_upper_corner[j]+=delta[j];
				this->_lower_corner[j]-=delta[j];
			}
			
			return *this;
		}
		
		inline bool operator!=(const Rectangle<State> &A) const {
			
			if (this->dim()!= A.dim()) return true;
		
			if ((A.empty()&&!this->empty())||(!A.empty()&&this->empty())) return false;
			
			for (size_t j=0; j< this->dim(); j++) {
				
				if (this->_upper_corner[j]!=A._upper_corner[j]) return true;
				if (this->_lower_corner[j]!=A._lower_corner[j]) return true;
			}
			
			return false;
		}
		
		friend std::ostream& 
		operator<< <> (std::ostream &os, 
			       const Rectangle<S> &r);

		friend std::istream& 
		operator>> <> (std::istream &is, 
			       Rectangle<S> &r);

};

template <typename S>
std::ostream& 
operator<<(std::ostream &os, const Rectangle<S> &r) 
{

    /*
    os << "{ lower=" << (r._lower_corner) << ", " ;
    os << "upper=" << (r._upper_corner) << " }" ;
    */

    os << "[ ";
    if(r.dim() > 0) {
	os << r[0];
	for(size_t i=1; i!=r.dim(); ++i) {
	    os << ", " << r[i];
	}
    }
    os << " ]";

    return os;
}

template <typename S>
std::istream& 
operator>>(std::istream &is, Rectangle<S> &r) 
{
    typedef typename Rectangle<S>::Real Real;
    typedef typename Ariadne::interval<Real> Interval;
    
    char c;
    is >> c;
    is.putback(c);
    if(c=='[') {
	/* Representation as list of intervals */ 
	std::vector< Interval > v;
	Ariadne::operator>>(is,v);
	r=Rectangle<S>(v.size(),&v[0]);
    }
    else {
	/* representation as lower and upper corners */
	/* FIXME */
	throw invalid_input("Not implemented");
    }
    return is;
}

}
}

#endif /* _RECTANGLE_H */
