/***************************************************************************
 *            point.h
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

/*! \file point.h
 *  \brief A point in Euclidean space.
 */

#ifndef _ARIADNE_POINT_H
#define _ARIADNE_POINT_H

#include <vector>
#include <iostream>
#include <stdexcept>

#include "../declarations.h"

#include "../utility/stlio.h"
#include "../linear_algebra/vector.h"


namespace Ariadne {
  namespace Geometry {
    
    template<typename R> typename Point<R>::difference_type operator-(const Point<R>&, const Point<R>&);
    template<typename R> Point<R> operator+(const Point<R>&, const typename Point<R>::difference_type&);
    template<typename R> Point<R> operator-(const Point<R>&, const typename Point<R>::difference_type&);
    template<typename R> Point<R> project_on_dimensions(const Point<R>& A,
		    			const Base::array<bool> &dims);
    template<typename R> Point<R> project_on_dimensions(const Point<R>& A,
		    			const size_t &x, const size_t &y);
    template<typename R> Point<R> project_on_dimensions(const Point<R>& A,
		    			const size_t &x, const size_t &y,
					const size_t &z);
    template<typename R> std::ostream& operator<<(std::ostream&, const Point<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Point<R>&);

    /*! \brief A point in Euclidean space. */
    template <typename R>
    class Point {
     public:
      typedef Ariadne::LinearAlgebra::vector<R> vector_type;

      typedef R real_type;
      typedef real_type value_type;
      typedef vector_type difference_type;
     private:
      /*! \brief The vector defining the state */
      vector_type _vector;

     public:
      Point() : _vector(0) { }

      Point(size_t dim) : _vector(dim) {
        for(size_t i=0; i!=dimension(); ++i) {
          _vector[i]=real_type(0);
        }
      }

      Point(size_t dim, const real_type default_value) : _vector(dim) {
        for(size_t i=0; i!=dimension(); ++i) {
          _vector[i]=default_value;
        }
      }

      template<class ForwardIterator>
      Point(ForwardIterator b, ForwardIterator e) : _vector(std::distance(b,e))
      {
        for(size_t i=0; i!=dimension(); ++i) {
          _vector[i]=*b;
          ++b;
        }

      }

      explicit Point(const vector_type& position) : _vector(position) { }

      Point(const Point& original) : _vector(original._vector) { }

      inline real_type& operator[] (size_t index) {
        if ((this->_vector).size() <= index) { 
	 /* Since index has type  size_t there is no need to check whever
	      index < 0 or not */ 
          throw std::out_of_range("Out of the vector's range.");
        }
        return  (this->_vector[index]);
      }

      inline const real_type& operator[](size_t index) const {
        if ((this->_vector).size() <= index) { 
	   /* Since index has type  size_t there is no need to check whever
	      index < 0 or not */ 	
            throw std::out_of_range("Out of the vector's range.");
        }
        return  (this->_vector[index]);
      }


      /*! \brief Checks equivalence between two states. */
      inline bool operator==(const Point<real_type>& A) const {
        /* Return false if states have different dimensions */
        if (this->dimension()!=A.dimension()) { return false; }

        /* for each dimension i */
        for (size_t i=0; i<this->dimension(); i++) {
          if (this->_vector[i]!=A._vector[i]) { return false; }
        }

        return true;
      }

      /*! \brief Returns the position vector of the point. */
      inline const vector_type& position_vector() const {
        return this->_vector; 
      }
      
      /*! \brief Checks equivalence between two states. */
      inline bool operator!=(const Point<real_type>& A) const {
        return !( *this == A );
      }

      /*! \brief The dimension of the Euclidean space the state lies in. */
      inline size_t dimension() const {
        return (this->_vector).size();
      }

      inline real_type get(size_t index) const {
        if ((this->_vector).size() <= index) { 
	   /* Since index has type  size_t there is no need to check whever
	      index < 0 or not */ 
            throw std::out_of_range("Out of the vector's range.");
        }
        return  (this->_vector[index]);
      }

      inline void set(size_t index, const real_type& r) {
        if ((this->_vector).size() <= index) {
	    /* Since index has type  size_t there is no need to check whever
	      index < 0 or not */ 
            throw std::out_of_range("Out of the vector's range.");
        }
        this->_vector[index]=r;
      }

      inline Point<real_type>& operator=(const Point<real_type>& A) {
        this->_vector=A._vector;
        return *this;
      }

       /*! \brief Difference of two states. */
      friend vector_type operator- <> (const Point<real_type>& s1,
                                  const Point<real_type>& s2);

       /*! \brief Difference of two states. */
      friend Point<real_type> operator+ <> (const Point<real_type>& s,
                                       const vector_type& v);

       /*! \brief Difference of two states. */
      friend Point<real_type> operator- <> (const Point<real_type>& s,
                                       const vector_type& v);

       /*! \brief Project a point. */
      friend Point<real_type> project_on_dimensions<>(const Point<real_type>& A,
		    			const Base::array<bool> &dims);
       /*! \brief Project a point. */
      friend Point<real_type> project_on_dimensions<>(const Point<real_type>& A,
		                        const size_t &x, const size_t &y); 
       /*! \brief Project a point. */
      friend Point<real_type> project_on_dimensions<>(const Point<real_type>& A,
		    	                const size_t &x, const size_t &y, 
				        const size_t &z);

      friend std::ostream& operator<< <>(std::ostream& os, const Point<real_type>& state);

      friend std::istream& operator>> <> (std::istream& is, Point<real_type>& state);
    };


    template <typename R>
    typename Point<R>::difference_type 
    operator-(const Point<R>& s1, const Point<R>& s2)
    {
      return s1._vector - s2._vector;
    }

    template <typename R>
    Point<R> 
    operator+(const Point<R>& s, const typename Point<R>::difference_type& v)
    {
      return Point<R>(s._vector + v);
    }

    template <typename R>
    Point<R> 
    operator-(const Point<R>& s, const typename Point<R>::difference_type& v)
    {
      return Point<R>(s._vector - v);
    }

    template <typename R>
    std::ostream& operator<<(std::ostream& os, const Point<R>& state)
    {
      os << "(";
      if(state.dimension() > 0) {
        os << state[0] ;
        for (size_t i=1; i<state.dimension(); i++) {
          os << ", " << state[i];
        }
      }
      os << ")" ;

      return os;
    }

    template <typename R>
    std::istream& operator>>(std::istream& is, Point<R>& state)
    {
      static size_t last_size;

      std::vector<R> v;
      v.reserve(last_size);
      is >> v;
      last_size = v.size();

      state._vector.resize(v.size());
      for(size_t i=0; i!=v.size(); ++i) {
        state._vector[i]=v[i];
      }
      return is;
    }

    template <typename R>
    inline Point<R> project_on_dimensions(const Point<R> &A, 
		    const Base::array<bool> &dims) {
	    
      if (A.dimension()!=dims.size())
         throw "project_on_dimensions(const Point & ,...): the two parameters have different dimension";
 
      size_t new_dim=0;

      for (size_t i=0; i< A.dimension(); i++)
        if (dims[i])
          new_dim++;

      Point<R> new_point(new_dim);
      
      size_t new_i=0;
      for (size_t i=0; i<dims.size(); i++) {
        if (dims[i]) {
           new_point.set(new_i,A[i]);
	   new_i++;
	}
      }

      return new_point;
    }

    template <typename R>
    inline Point<R> project_on_dimensions(const Point<R> &A, 
		    const size_t &x, const size_t&y, const size_t&z) {
	    
      if ((A.dimension()<=x)||(A.dimension()<=y)||(A.dimension()<=z))
         throw "project_on_dimensions(const Point& ,...): one of the projection dimensions is greater than the Point dimension";
      
      Point<R> new_point(3);
      
      new_point.set(0,A[x]);
      new_point.set(1,A[y]);
      new_point.set(2,A[z]);

      return new_point;
    }

    template <typename R>
    inline Point<R> project_on_dimensions(const Point<R> &A, 
		    const size_t &x, const size_t&y) {
	    
      if ((A.dimension()<=x)||(A.dimension()<=y))
         throw "project_on_dimensions(const Point& ,...): one of the projection dimensions is greater than the Point dimension";
      
      Point<R> new_point(2);
      
      new_point.set(0,A[x]);
      new_point.set(1,A[y]);

      return new_point;
    }

  }
}

#endif /* _ARIADNE_POINT_H */
