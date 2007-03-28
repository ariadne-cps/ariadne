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

#include "../declarations.h"
#include "../exceptions.h"
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

    template<class R>
    class Point {
      typedef typename Numeric::traits<R>::arithmetic_type F;
     public:
      /*!\brief The type of denotable real number giving the point's values. */
      typedef R real_type;
      /*!\brief The type of denotable vector in the space the point lies in. */
      typedef LinearAlgebra::Vector<R> vector_type;

      /*! \brief Default constructor. */
      Point();

      /*! \brief The origin in dimension \a dim. */
      explicit Point(dimension_type d);

      /*! \brief Construct a point a strided array of data. */
      template<class Rl>
      Point(dimension_type d, const Rl* data, size_type inc=1);
      
      /*! \brief Construct a point from a range of values. */
      template<class ForwardIterator>
      Point(ForwardIterator b, ForwardIterator e);

      /*! \brief Construct a point from a position vector. */
      explicit Point(const LinearAlgebra::Vector<R>& position);

      /*! \brief Construct a point from a string literal. */
      explicit Point(const std::string& s);

      /*! \brief Convert from a point which possibly lies in a different space. */
      template<class R2> Point(const Point<R2>& original);
      
      /*! \brief Copy constructor. */
      Point(const Point<R>& original);

      /*! \brief Assignment operator. */
      Point<R>& operator=(const Point<R>& original);
      
 
      /*! \brief Checks equivalence between two states. */
      bool operator==(const Point<R>& A) const;

      /*! \brief Inequality operator. */
      bool operator!=(const Point<R>& A) const;

      /*! \brief The array containing the point's values. */
      const array<R>& data() const;

      /*! \brief The dimension of the Euclidean space the state lies in. */
      dimension_type dimension() const;

      /*! \brief Subcripting operator. */
      R& operator[](dimension_type index);

      /*! \brief Subcripting operator. */
      const R& operator[](dimension_type index) const;

      /*! \brief The position vector of the point. */
      const LinearAlgebra::Vector<R>& position_vector() const;
      
     
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);

     private:
      LinearAlgebra::Vector<R> _vector;
    };

    template<class R> 
    Point<R> approximation(const Point<R>& pt);
    
    template<class R> 
    Point<R> approximation(const Point< Interval<R> >& ipt);
    
    template<class R> 
    bool contains_value(const Point< Interval<R> >& ipt, const Point<R>& pt);
    
    template<class R> 
    R error_bound(const Point< Interval<R> >& ipt);
    

    
    template<class R>
    Point<typename Numeric::traits<R>::arithmetic_type>
    minkowski_sum(const Point<R>& pt1, const Point<R>& pt2);
    
    template<class R>
    Point<typename Numeric::traits<R>::arithmetic_type>
    minkowski_difference(const Point<R>& pt1, const Point<R>& pt2);
    
    template<class R1,class R2>
    LinearAlgebra::Vector<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator-(const Point<R1> pt1, const Point<R2>& pt2);
    
    template<class R1,class R2>
        Point<typename Numeric::traits<R1,R2>::arithmetic_type> 
    operator+(const Point<R1>& pt, const LinearAlgebra::Vector<R2>& v);


    template<class R1,class R2>
    Point<typename Numeric::traits<R1,R2>::arithmetic_type> 
    operator-(const Point<R1>& pt, const LinearAlgebra::Vector<R2>& v);

    template<class R>
    Point<R> add_approx(const Point<R>& pt, const LinearAlgebra::Vector<R>& v);

    template<class R>
    Point<R> sub_approx(const Point<R>& pt, const LinearAlgebra::Vector<R>& v);

    
    template<class R>
    Point<R>
    approximate_value(const Point< Interval<R> >& pt);
    
  //template<class R>
  //Point<R> 
  //project_on_dimensions(const Point<R> &A, const Base::array<bool>& dims);

    template<class R>
    Point<R> 
    project_on_dimensions(const Point<R> &A, 
                          const size_type& x, const size_type& y, const size_type& z);

    template<class R>
    Point<R>  
    project_on_dimensions(const Point<R> &A, 
                          const size_type &x, const size_type&y);
    

    template<class R>  
    std::ostream& operator<<(std::ostream& os, const Point<R>& pt);
    
    template<class R> 
    std::istream& operator>>(std::istream& is, Point<R>& pt);


  }
}

#include "point.inline.h"

#endif /* ARIADNE_POINT_H */
