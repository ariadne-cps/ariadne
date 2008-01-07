/***************************************************************************
 *            point.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#ifndef ARIADNE_POINT_H
#define ARIADNE_POINT_H

#include <iosfwd>
#include <stdexcept>

#include "base/array.h"
#include "linear_algebra/declarations.h"
#include "linear_algebra/vector.h"

#include "exceptions.h"

namespace Ariadne {
  namespace Geometry {
   
    template<class> class Box;
    template<class> class Sphere;
    template<class> class Ellipsoid;

    /*!\brief A point in Euclidean space.
     *
     * Points may be instantiated with both numeric and interval coefficients. 
     * An interval point can be automatically converted to and from a rectangle with numerical coefficients.
     * However, an interval point canonically represents an approximation to a single point, whereas a rectangle represents a set of points.
     */
    template<class X>
    class Point {
      typedef typename Numeric::traits<X>::number_type R;
      typedef typename Numeric::traits<X>::arithmetic_type F;
      typedef typename Numeric::traits<X>::interval_type I;
     public:
      /*!\brief The type of denotable real number giving the point's values. */
      typedef R real_type;
      /*!\brief The type of denotable vector in the space the point lies in. */
      typedef LinearAlgebra::Vector<X> vector_type;

      /*! \brief Default constructor. */
      Point();

      /*! \brief The origin in dimension \a dim. */
      explicit Point(dimension_type d);

      /*! \brief Construct a point from an array of data. */
      template<class XX>
      explicit Point(const array<XX>& ary);
      
      /*! \brief Construct a point a strided array of data. */
      template<class XX>
      explicit Point(dimension_type d, const XX* data, size_type inc=1);
      
      /*! \brief Construct a point from a range of values. */
      template<class ForwardIterator>
      explicit Point(ForwardIterator b, ForwardIterator e);

      /*! \brief Construct a point from a position vector. */
      explicit Point(const LinearAlgebra::Vector<X>& position);

      /*! \brief Construct a point from a string literal. */
      explicit Point(const std::string& s);

      /*! \brief Convert from a point which possibly lies in a different space. */
      template<class XX> Point(const Point<XX>& original);
      
      /*! \brief Copy constructor. */
      Point(const Point<X>& original);

      /*! \brief Assignment operator. */
      Point<X>& operator=(const Point<X>& original);
      
 
      /*! \brief Checks equivalence between two states. */
      bool operator==(const Point<X>& A) const;

      /*! \brief Inequality operator. */
      bool operator!=(const Point<X>& A) const;

      /*! \brief The array containing the point's values. */
      const array<X>& data() const;

      /*! \brief The dimension of the Euclidean space the state lies in. */
      dimension_type dimension() const;

      /*! \brief Change the dimension of the space the point lies in. */
      void resize(dimension_type d);
      
      /*! \brief Subcripting operator (unchecked). */
      X& operator[](dimension_type index);

      /*! \brief Subcripting operator (unchecked). */
      const X& operator[](dimension_type index) const;

      /*! \brief Subcripting operator (checked). */
      X& at(dimension_type index);

      /*! \brief Subcripting operator (checked). */
      const X& at(dimension_type index) const;

      /*! \brief The position vector of the point. */
      const LinearAlgebra::Vector<X>& position_vector() const;
      
     
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);

     private:
      LinearAlgebra::Vector<X> _vector;
    };

    
    template<class R>
    Point<R> approximation(const Point<R>& pt);
    
    template<class R>
    Point<R> approximation(const Point< Numeric::Interval<R> >& ipt);

    
    template<class R> 
    Point<R> midpoint(const Point< Numeric::Interval<R> >& ipt);
    
    template<class R> 
    R radius(const Point< Numeric::Interval<R> >& ipt);

    template<class R> 
    bool encloses(const Point< Numeric::Interval<R> >& ipt, const Point<R>& pt);
    
    template<class R> 
    bool refines(const Point< Numeric::Interval<R> >& ipt1, const Point< Numeric::Interval<R> >& ipt2);


    
    template<class X>
    Point<typename Numeric::traits<X>::arithmetic_type>
    minkowski_sum(const Point<X>& pt1, const Point<X>& pt2);
    
    template<class X>
    Point<typename Numeric::traits<X>::arithmetic_type>
    minkowski_difference(const Point<X>& pt1, const Point<X>& pt2);
    
    template<class X1,class X2>
    LinearAlgebra::Vector<typename Numeric::traits<X1,X2>::arithmetic_type>
    operator-(const Point<X1>& pt1, const Point<X2>& pt2);
    
    template<class X1,class X2>
    Point<typename Numeric::traits<X1,X2>::arithmetic_type> 
    operator+(const Point<X1>& pt, const LinearAlgebra::Vector<X2>& v);


    template<class X1,class X2>
    Point<typename Numeric::traits<X1,X2>::arithmetic_type> 
    operator-(const Point<X1>& pt, const LinearAlgebra::Vector<X2>& v);

    template<class X>
    Point<X> add_approx(const Point<X>& pt, const LinearAlgebra::Vector<X>& v);

    template<class X>
    Point<X> sub_approx(const Point<X>& pt, const LinearAlgebra::Vector<X>& v);

    
  //template<class X>
  //Point<X> 
  //project_on_dimensions(const Point<X> &A, const Base::array<bool>& dims);

    template<class X>
    Point<X> 
    project_on_dimensions(const Point<X> &A, 
                          const size_type& x, const size_type& y, const size_type& z);

    template<class X>
    Point<X>  
    project_on_dimensions(const Point<X> &A, 
                          const size_type &x, const size_type&y);
    

    template<class X>  
    std::ostream& operator<<(std::ostream& os, const Point<X>& pt);
    
    template<class X> 
    std::istream& operator>>(std::istream& is, Point<X>& pt);


  }
}

#include "point.inline.h"

#endif /* ARIADNE_POINT_H */
