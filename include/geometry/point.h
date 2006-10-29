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

#include <iosfwd>
#include <stdexcept>

#include "../declarations.h"
#include "../base/exceptions.h"
#include "../base/array.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace Geometry {
   
    template<class> class Rectangle;
    template<class> class Parallelotope;
    template<class> class Zonotope;
    template<class> class Sphere;
    template<class> class Ellipsoid;

    template<template<class> class BS1, template<class> class BS2> 
    inline bool is_a(){ return false; }

    template<> inline bool is_a<Point, Point>(){ return true; }

    /* Forward declaration of friends. */
    template<class R> std::ostream& operator<<(std::ostream&, const Point<R>&);
    template<class R> std::istream& operator>>(std::istream&, Point<R>&);

    /*!\brief A point in Euclidean space. 
     *
     * A point is defined by its coordinates or type #real_type, which are accessed by 
     * operator[](dimension_type). 
     *
     * We consider Euclidean space as an affine space, so it is allowed to
     * add a vector (of type #vector_type or LinearAlgebra::Vector < R >) to a point, or subtract two points to obtain a vector,
     * but not add two points.
     *
     * Points can be constructed from string literals of the form
     * "(p0,p1,...)". Note the use of round brackets in the
     * literal expression.
     *
     * Points are ordered by the lexicographic order.
     *
     * Operations on points which cannot be computed exactly in the arithmetic
     * used, give results of type FuzzyPoint< R >, which is a Point< Numeric::Interval < R > >.
     * Currently, %Ariadne does not support other space types, such as annuli
     * or general differential manifolds.
     */
    template<class R>
    class Point {
      typedef typename Numeric::traits<R>::arithmetic_type F;
     public:
      /*!\brief The type of denotable real number giving the point's values. */
      typedef R real_type;
      /*!\brief The type of denotable vector in the space the point lies in. */
      typedef LinearAlgebra::Vector<R> vector_type;

      /*! \brief Default constructor. */
      Point() : _vector(0) { }

      /*! \brief The origin in dimension \a dim. */
      explicit Point(dimension_type d) : _vector(d) {
        for(size_type i=0; i!=dimension(); ++i) {
          _vector[i]=real_type(0);
        }
      }

      /*! \brief Construct a point a strided array of data. */
      template<class Rl>
      Point(dimension_type d, const Rl* data, size_type inc=1)
        : _vector(d,data,inc) { }
      
      /*! \brief Construct a point from a range of values. */
      template<class ForwardIterator>
      Point(ForwardIterator b, ForwardIterator e) : _vector(std::distance(b,e))
      {
        for(size_type i=0; i!=dimension(); ++i) {
          _vector[i]=*b;
          ++b;
        }
      }

      /*! \brief Construct a point from a position vector. */
      explicit Point(const LinearAlgebra::Vector<R>& position) : _vector(position) { }

      /*! \brief Construct a point from a string literal. */
      explicit Point(const std::string& s);

      /*! \brief Convert from a point which possibly lies in a different space. */
      template<class R2>
      Point(const Point<R2>& original) : _vector(original.position_vector()) { }
      
      /*! \brief Copy constructor. */
      Point(const Point<R>& original) : _vector(original._vector) { }

      /*! \brief Assignment operator. */
      Point<R>& operator=(const Point<R>& original) {
        if(this!=&original) { this->_vector=original._vector; }
        return *this; 
      }
      
 
      /*! \brief Checks equivalence between two states. */
      bool operator==(const Point<R>& A) const {
        /* Return false if states have different dimensions */
        if (this->dimension()!=A.dimension()) { return false; }
        for (size_type i=0; i<this->dimension(); i++) {
          if (this->_vector[i]!=A._vector[i]) { return false; }
        }
        return true;
      }
      
      /*! \brief Inequality operator. */
      bool operator!=(const Point<real_type>& A) const {
        return !( *this == A );
      }

      /*! \brief The dimension of the Euclidean space the state lies in. */
      dimension_type dimension() const {
        return this->_vector.size();
      }

      /*! \brief Subcripting operator. */
      real_type& operator[] (dimension_type index) {
        check_index(*this,index,"Point<R>::operator[]");
        return  (this->_vector[index]);
      }

      /*! \brief Subcripting operator. */
      const real_type& operator[](dimension_type index) const {
        check_index(*this,index,"Point<R>::operator[]");
        return  (this->_vector[index]);
      }

      /*! \brief The position vector of the point. */
      const LinearAlgebra::Vector<R>& position_vector() const {
        return this->_vector; 
      }
      
     
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);

      friend std::ostream& operator<< <>(std::ostream& os, const Point<real_type>& state);
      friend std::istream& operator>> <> (std::istream& is, Point<real_type>& state);
     private:
      LinearAlgebra::Vector<R> _vector;
    };

    template<class R> inline
    Point<R> approximation(const Point< Interval<R> >& ipt) 
    {
      Point<R> result(ipt.dimension());
      for(dimension_type i=0; i!=ipt.dimension(); ++i) {
        result[i]=ipt[i].centre();
      }
      return result;
    }
    
    template<class R>
    inline
    Point<typename Numeric::traits<R>::arithmetic_type>
    minkowski_sum(const Point<R>& pt1, const Point<R>& pt2) {
      check_dimension(pt1,pt2,"minkowski_sum(Point<R>,Point<R>)");
      return Point<typename Numeric::traits<R>::arithmetic_type>(pt1.position_vector()+pt2.position_vector());
    }
    
    template<class R>
    inline
    Point<typename Numeric::traits<R>::arithmetic_type>
    minkowski_difference(const Point<R>& pt1, const Point<R>& pt2) {
      check_dimension(pt1,pt2,"minkowski_difference(Point<R>,Point<R>)");
      return Point<typename Numeric::traits<R>::arithmetic_type>(pt1.position_vector()-pt2.position_vector());
    }
    
    template<class R1,class R2>
    inline
    LinearAlgebra::Vector<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator-(const Point<R1> pt1, const Point<R2>& pt2) 
    {
      check_dimension(pt1,pt2,"operator-(Point<R>,Point<R>)");
      return pt1.position_vector()-pt2.position_vector();
    }
    
    template<class R1,class R2>
    inline
    Point<typename Numeric::traits<R1,R2>::arithmetic_type> 
    operator+(const Point<R1>& pt, const LinearAlgebra::Vector<R2>& v)
    {
      check_dimension(pt,v,"operator+(Point<R>,Vector<R>)");
      return Point<typename Numeric::traits<R1,R2>::arithmetic_type>(pt.position_vector() + v);
    }


    template<class R1,class R2>
    inline
    Point<typename Numeric::traits<R1,R2>::arithmetic_type> 
    operator-(const Point<R1>& pt, const LinearAlgebra::Vector<R2>& v)
    {
      check_dimension_size(pt,v,"operator-(Point<R>,Vector<R>)");
      return Point<typename Numeric::traits<R1,R2>::arithmetic_type>(pt.position_vector() - v);
    }

    template<class R>
    inline
    Point<R> 
    add_approx(const Point<R>& pt, const LinearAlgebra::Vector<R>& v)
    {
      check_dimension_size(pt,v,"add_approx(Point<R>,Vector<R>)");
      return Point<R>(add_approx(pt.position_vector(),v));
    }

    template<class R>
    inline
    Point<R> 
    sub_approx(const Point<R>& pt, const LinearAlgebra::Vector<R>& v)
    {
      check_dimension_size(pt,v,"sub_approx(Point<R>,Vector<R>)");
      return Point<R>(sub_approx(pt.position_vector(),v));
    }

    
    template<class R>
    inline
    Point<R>
    approximate_value(const Point< Interval<R> >& pt) 
    {
      Point<R> result(pt.dimension());
      for(dimension_type i=0; i!=result.dimension(); ++i) {
        result[i]=approximate_value(pt[i]);
      }
      return result;
    }
    
    template<class R>
    inline 
    Point<R> 
    project_on_dimensions(const Point<R> &A, const Base::array<bool>& dims) 
    {
      if (A.dimension()!=dims.size()) {
        throw std::runtime_error("project_on_dimensions(const Point & ,...): the two parameters have different dimension");
      }
      size_type new_dim=0;

      for (size_type i=0; i< A.dimension(); i++) {
        if (dims[i]) {
          new_dim++;
        }
      }
      Point<R> new_point(new_dim);
      
      size_type new_i=0;
      for (size_t i=0; i<dims.size(); i++) {
        if (dims[i]) {
           new_point.set(new_i,A[i]);
           new_i++;
        } 
      }

      return new_point;
    }

    template<class R>
    inline 
    Point<R> 
    project_on_dimensions(const Point<R> &A, 
                          const size_type& x, const size_type& y, const size_type& z) 
    {

      if ((A.dimension()<=x)||(A.dimension()<=y)||(A.dimension()<=z)) {
        throw std::runtime_error("project_on_dimensions(const Point & ,...): the two parameters have different dimension");
      }
      
      Point<R> new_point(3);
      
      new_point.set(0,A[x]);
      new_point.set(1,A[y]);
      new_point.set(2,A[z]);

      return new_point;
    }

    template<class R>
    inline 
    Point<R> 
    project_on_dimensions(const Point<R> &A, 
                          const size_type &x, const size_type&y) 
    {
      if ((A.dimension()<=x)||(A.dimension()<=y)) {
         throw "project_on_dimensions(const Point& ,...): one of the projection dimensions is greater than the Point dimension";
      }
      Point<R> new_point(2);
      
      new_point.set(0,A[x]);
      new_point.set(1,A[y]);

      return new_point;
    }
    
    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const Point<R>& pt) {
      return pt.write(os);
    }
    
    template<class R> inline
    std::istream& operator>>(std::istream& is, Point<R>& pt) {
      return pt.read(is);
    }


  }
}

#endif /* _ARIADNE_POINT_H */
