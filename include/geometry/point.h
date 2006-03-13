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

#include "../base/utility.h"

#include "../linear_algebra/vector.h"

#include "../geometry/geometry_declarations.h"

namespace Ariadne {
  namespace Geometry {
    
    template<typename R> typename Point<R>::difference_type operator-(const Point<R>&, const Point<R>&);
    template<typename R> Point<R> operator+(const Point<R>&, const typename Point<R>::difference_type&);
    template<typename R> Point<R> operator-(const Point<R>&, const typename Point<R>::difference_type&);
    template<typename R> std::ostream& operator<<(std::ostream&, const Point<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Point<R>&);

    /*! \brief A point in Euclidean space. */
    template <typename R>
    class Point {
     public:
      typedef Ariadne::LinearAlgebra::vector<R> Vector;

      typedef R Real;
      typedef Real value_type;
      typedef Vector difference_type;
     private:
      /*! \brief The vector defining the state */
      Vector _vector;

     public:
      Point() : _vector(0) { }

      Point(size_t dim) : _vector(dim) {
        for(size_t i=0; i!=dimension(); ++i) {
          _vector[i]=Real(0);
        }
      }

      Point(size_t dim, const Real default_value) : _vector(dim) {
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

      explicit Point(const Vector& position) : _vector(position) { }

      Point(const Point& original) : _vector(original._vector) { }

      inline Real& operator[] (size_t index) {
        if ((this->_vector).size() <= index) { 
	 /* Since index has type  size_t there is no need to check whever
	      index < 0 or not */ 
          throw std::out_of_range("Out of the vector's range.");
        }
        return  (this->_vector[index]);
      }

      inline const Real& operator[](size_t index) const {
        if ((this->_vector).size() <= index) { 
	   /* Since index has type  size_t there is no need to check whever
	      index < 0 or not */ 	
            throw std::out_of_range("Out of the vector's range.");
        }
        return  (this->_vector[index]);
      }


      /*! \brief Checks equivalence between two states. */
      inline bool operator==(const Point<Real>& A) const {
        /* Return false if states have different dimensions */
        if (this->dimension()!=A.dimension()) { return false; }

        /* for each dimension i */
        for (size_t i=0; i<this->dimension(); i++) {
          if (this->_vector[i]!=A._vector[i]) { return false; }
        }

        return true;
      }

      /*! \brief Returns the position vector of the point. */
      inline const Vector& position_vector() const {
        return this->_vector; 
      }
      
      /*! \brief Checks equivalence between two states. */
      inline bool operator!=(const Point<Real>& A) const {
        return !( *this == A );
      }

      /*! \brief The dimension of the Euclidean space the state lies in. */
      inline size_t dimension() const {
        return (this->_vector).size();
      }

      inline Real get(size_t index) const {
        if ((this->_vector).size() <= index) { 
	   /* Since index has type  size_t there is no need to check whever
	      index < 0 or not */ 
            throw std::out_of_range("Out of the vector's range.");
        }
        return  (this->_vector[index]);
      }

      inline void set(size_t index, const Real& r) {
        if ((this->_vector).size() <= index) {
	    /* Since index has type  size_t there is no need to check whever
	      index < 0 or not */ 
            throw std::out_of_range("Out of the vector's range.");
        }
        this->_vector[index]=r;
      }

      inline Point<Real>& operator=(const Point<Real>& A) {
        this->_vector=A._vector;
        return *this;
      }

       /*! \brief Difference of two states. */
      friend Vector operator- <> (const Point<Real>& s1,
                                  const Point<Real>& s2);

       /*! \brief Difference of two states. */
      friend Point<Real> operator+ <> (const Point<Real>& s,
                                       const Vector& v);

       /*! \brief Difference of two states. */
      friend Point<Real> operator- <> (const Point<Real>& s,
                                       const Vector& v);

      friend std::ostream& operator<< <>(std::ostream& os, const Point<Real>& state);

      friend std::istream& operator>> <> (std::istream& is, Point<Real>& state);
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

  }
}

#endif /* _ARIADNE_POINT_H */
