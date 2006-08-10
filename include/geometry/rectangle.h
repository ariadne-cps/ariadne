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
 
/*! \file rectangle.h
 *  \brief Rectangles and cuboids.
 */

#ifndef _ARIADNE_RECTANGLE_H
#define _ARIADNE_RECTANGLE_H

#include <iosfwd>

#include "../declarations.h"

#include "../numeric/arithmetic.h"
#include "../numeric/interval.h"

#include "../linear_algebra/interval_vector.h"

#include "../geometry/point.h"

namespace Ariadne {
  namespace Geometry {

    template<> 
    inline bool is_a<Rectangle,Rectangle>() { return true; }
    template<> 
    inline bool is_a<Rectangle,Parallelotope>() { return true; }
    template<> 
    inline bool is_a<Rectangle,Zonotope>() { return true; }
    template<> 
    inline bool is_a<Rectangle,Polyhedron>() { return true; }

    /*! \brief A cuboid of arbitrary dimension.
     */
    template <typename R>
    class Rectangle {
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the rectangle. */
      typedef Point<R> state_type;
      /*! \brief An iterator through the vertices of the rectangle. */
      typedef typename array< Point<R> >::const_iterator const_vertex_iterator;
    
     public:
      /*! \brief Construct a rectangle with coordinates zero of dimension \a n. */
      explicit Rectangle(size_type n=0)
        : _bounds(n)
      { 
      }
      
      /*! \brief Construct a rectangle from a range of interval values. */
      template<class ForwardIterator>
      Rectangle(ForwardIterator b, ForwardIterator e)
        : _bounds(std::distance(b,e))
      {
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          this->_bounds[i]=*b;
          ++b;
        }
      }
      
      /*! \brief Construct from an array of intervals. */
      explicit Rectangle(const array< Interval<R> >& a)
        : _bounds(a.size())
      {
        for(dimension_type i=0; i!=a.size(); ++i) {
          this->_bounds[i]=a[i];
        }
      }
      
      /*! \brief Construct from a std::vector of intervals. */
      explicit Rectangle(const std::vector< Interval<R> >& v)
        : _bounds(v.size())
      {
        for(dimension_type i=0; i!=v.size(); ++i) {
          this->_bounds[i]=v[i];
        }
      }

      /*! \brief Construct a degenerate rectangle from a single point. */
      explicit Rectangle(const Point<R>& p) 
        : _bounds(p.dimension())
      {
        for(dimension_type i=0; i!=p.dimension(); ++i) {
          this->_bounds[i]=p[i];
        }
      }
      
      /*! \brief Construct from two corners. */
      explicit Rectangle(const Point<R>& p1, const Point<R>& p2) 
        : _bounds(p1.dimension())
      {
        if (p1.dimension()!=p2.dimension()) {
          throw std::domain_error("The parameters have different space dimensions");
        }
        for (size_type i=0; i!=this->dimension(); ++i) {
          this->_bounds[i]=Interval<R>(std::min(p1[i],p2[i]),std::max(p1[i],p2[i]));
        }
      }
      
      /*! \brief Construct from a string literal. */
      explicit Rectangle(const std::string& s);
      
      /*! \brief Construct from an interval vector. */
      explicit Rectangle(const LinearAlgebra::IntervalVector<R>& iv)
        : _bounds(iv)
      {
      }

      /*! \brief Copy constructor. */
      Rectangle(const Rectangle<R>& original)
        : _bounds(original._bounds)
      { }
    
      /*! \brief Copy assignment operator. */
      Rectangle<R>& operator=(const Rectangle<R>& A) {
        if(this != &A) {
          this->_bounds = A._bounds;
        }
        return *this;
      }
      
      /*! \brief The equality operator */
      bool operator==(const Rectangle<R>& A) const
      {
        if (A.empty() && this->empty()) { return true; }
        if (A.empty() || this->empty()) { return false; }
        if(this->dimension()!=A.dimension()) { return false; }
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          if (this->lower_bound(i)!=A.lower_bound(i)) { return false; }
          if (this->upper_bound(i)!=A.upper_bound(i)) { return false; }
        }
        return true;
      }
      
      /*! \brief The inequality operator */
      bool operator!=(const Rectangle<R>& A) const {
        return !(*this == A);
      }

      /*! \brief The dimension of the Euclidean space the rectangle lies in. */
      size_type dimension() const {
        return this->_bounds.size();
      }
      
      /*! \brief Makes the rectangle empty. */
      void clear() {
        if(this->_bounds.size()!=0) {
          this->_bounds[0]=Interval<R>();
        }
      }
      
      /*! \brief True if the rectangle is empty. */
      bool empty() const {
        if(this->dimension()==0) {
          return true;
        }
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          if(this->interval(i).empty()) {
            return true;
          }
        }
        return false;
      }
      
      /*! \brief True if the rectangle has empty interior. */
      bool empty_interior() const {
        if(this->dimension()==0) {
          return true;
        }
        for(size_type i=0; i!=this->dimension(); ++i) {
          if(this->lower_bound(i) >= this->upper_bound(i)) {
            return true;
          }
        }
        return false;
      }
      
      /*! \brief A rectangle containing the given rectangle; returns a copy. */
      Rectangle bounding_box() const {
        return *this;
      }
        
      /*! \brief The lower corner. */
      Point<R> lower_corner() const {
        Point<R> result(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          result[i]=this->lower_bound(i);
        }
        return result;
      }
      
      /*! \brief The upper corner. */
      Point<R> upper_corner() const {
        Point<R> result(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          result[i]=this->upper_bound(i);
        }
        return result;
      }
      
      /*! \brief The centre. */
      Point<R> centre() const {
        return Point<R>(this->_bounds.centre());
      }
      
      /*! \brief The radius in the sup norm. */
      R radius() const {
        R diameter=0;
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          diameter=max(diameter,this->_bounds[i].length());
        }
        return diameter/2;
      }
      
      /*! \brief Returns the projection onto the \a i th coordinate. */
      IntervalReference<R> operator[] (dimension_type i) {
        return IntervalReference<R>(this->_bounds[i]);
      }
      //Interval<R>& operator[] (size_type i) {
      //  return this->_bounds[i];
      //}
      
      /*! \brief Returns the projection onto the \a i th coordinate. */
      Interval<R> operator[] (dimension_type i) const {
        return this->_bounds[i];
      }
      
      /*! \brief Returns the projection onto the \a i th coordinate. */
      Interval<R> interval(size_type i) const {
        return this->_bounds[i];
      }
      
      /*! \brief Returns the lower bound of the \a i th coordinate */
      const R& lower_bound(size_type i) const {
        return this->_bounds[i].lower();
      }
      
      /*! \brief Returns the upper bound of the \a i th coordinate */
      const R& upper_bound(size_type i) const {
        return this->_bounds[i].upper();
      }
      
      /*! \brief Sets the \a i th interval. */
      void set_interval(size_type i, Interval<R> x) {
        this->_bounds[i]=x;
      }
      
      /*! \brief Sets the lower bound of the \a i th coordinate to \a r. */
      void set_lower_bound(size_type i, const R& l) {
        this->_bounds[i]=Interval<R>(l,this->_bounds[i].upper());
      }
      
      /*! \brief Sets the upper bound of the \a i th coordinate to \a u. */
      void set_upper_bound(size_type i, const R& u) {
        this->_bounds[i]=Interval<R>(this->_bounds[i].lower(),u);
      }
      
      /*! \brief Convert to a rational polyhedron. */
      operator Polyhedron<Rational> () const;
            
      /*! \brief The set of position vectors of the rectangle. */
      const LinearAlgebra::IntervalVector<R>& position_vectors() const {
        return this->_bounds;
      }
      
      /*! \brief Expand the Rectangle by \a delta in each direction. */
      Rectangle<R>& expand_by(const R& delta);
      
      /*! \brief Compute a quadrant of the Rectangle determined by \a q.
       *
       *  \a q is a binary word such that the ith bit of q is 0 if the lower half
       *  of the rectangle in the ith coordinate is used, and 1 if the upper
       *  half is used.
       */
      Rectangle<R> quadrant(const BinaryWord& q) const;
      /*! \brief Subdivide into smaller pieces. */
      ListSet<R,Geometry::Rectangle> subdivide() const;
      
      /*! The vertices of the rectangle. */
      std::vector< Point<R> > vertices() const;
      /*! \brief The \a i th vertex. */
      Point<R> vertex(size_type i) const;
      
      
      /*! \brief Tests if \a point is included into a rectangle. */
      bool contains(const Point<R>& pt) const;
      /*! \brief Tests if \a point is included into the interior a rectangle. */
      bool interior_contains(const Point<R>& pt) const;

      /*! \brief Set equality operator. */
      static bool equal(const Rectangle<R>& r1, const Rectangle<R>& r2);
     
      /*! \brief Tests disjointness with \a r. */
      static bool disjoint(const Rectangle<R>& r1, const Rectangle<R>& r2);
      /*! \brief Tests overlap with \a r. */
      static bool interiors_intersect(const Rectangle<R>& r1, const Rectangle<R>& r2);
      /*! \brief Tests inclusion in the interior of \a set. */
      static bool inner_subset(const Rectangle<R>& r1, const Rectangle<R>& r2);
      /*! \brief Tests if the rectangle is a subset of another rectangle. */
      static bool subset(const Rectangle<R>& r1, const Rectangle<R>& r2);
      
      /*! \brief The intersection of \a A and \a B. */
      static Rectangle<R> intersection(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The closure of the intersection of the interiors of \a A and \a B. */
      static Rectangle<R> regular_intersection(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The smallest rectangle containing \a A and \a B. */
      static Rectangle<R> rectangular_hull(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The componentwise sum of rectangles \a A and \a B. */
      static Rectangle<R> minkowski_sum(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The componentwise difference of rectangles \a A and \a B. */
      static Rectangle<R> minkowski_difference(const Rectangle<R>& A, const Rectangle<R>& B); 
      
      std::ostream& write(std::ostream& os) const;
      std::istream& read(std::istream& is);
     private:
      LinearAlgebra::IntervalVector<R> _bounds;
    };
  

    
    template <typename R>
    inline
    bool 
    Rectangle<R>::contains(const Point<R>& p) const 
    {
      const Rectangle<R>& self=*this;
      if (p.dimension()!=self.dimension()) {
        throw std::domain_error("The point has a different dimension to the rectangle.");
      }
      for (size_type i=0; i!=self.dimension(); ++i) {
        if(!self[i].contains(p[i])) {
          return false;
        }
      }
      return true;
    }

    template <typename R>
    inline
    bool 
    Rectangle<R>::interior_contains(const Point<R>& p) const 
    {
      const Rectangle<R>& self=*this;
      if (p.dimension()!=self.dimension()) {
        throw std::domain_error("The point has a different dimension to the rectangle.");
      }
      for (size_type i=0; i!=self.dimension(); ++i) {
        if(!self[i].interior_contains(p[i])) {
          return false;
        }
      }
      return true;
    }
      
    template <typename R>
    inline
    bool 
    Rectangle<R>::equal(const Rectangle<R>& r1, const Rectangle<R>& r2)
    {
      if(r1.dimension()!=r2.dimension()) {
        return false;
      }
      for(size_type i=0; i!=r1.dimension(); ++i) {
        if(!Interval<R>::equal(r1[i],r2[i])) {
          return false;
        }
      }
      return true;
    }
      
    template <typename R>
    inline
    bool 
    Rectangle<R>::disjoint(const Rectangle<R>& r1, const Rectangle<R>& r2)
    {
      if(r1.dimension()!=r2.dimension()) {
        throw std::domain_error("disjoint(Rectangle,Rectangle): The two rectangles have different dimensions");
      }
      for(size_type i=0; i!=r1.dimension(); ++i) {
        if(Interval<R>::disjoint(r1[i],r2[i])) {
          return true;
        }
      }
      return false;
    }
   
    template <typename R>
    inline
    bool 
    Rectangle<R>::interiors_intersect(const Rectangle<R>& r1, const Rectangle<R>& r2)
    {
      if (r1.dimension()!=r2.dimension()) {
        throw std::domain_error("interiors_intersect(Rectangle,Rectangle): The two rectangles have different dimensions");
      }
      for(size_type i=0; i!=r1.dimension(); ++i) {
        if(!Interval<R>::interiors_intersect(r1[i],r2[i])) {
          return false;
        }
      }
      return true;
    }
    
    template <typename R>
    inline
    bool 
    Rectangle<R>::inner_subset(const Rectangle<R>& r1, const Rectangle<R>& r2)
    {
      if (r1.dimension()!=r2.dimension()) {
        throw std::domain_error("inner_subset(Rectangle,Rectangle): The two rectangles have different dimensions");
      }
      for (size_type i=0; i!=r1.dimension(); ++i) {
        if(!Interval<R>::inner_subset(r1[i],r2[i])) {
          return false;
        }
      }
      return true;
    }
    
    template <typename R>
    inline
    bool 
    Rectangle<R>::subset(const Rectangle<R>& r1, const Rectangle<R>& r2)
    {
      if (r1.dimension()!=r2.dimension()) {
        throw std::domain_error("subset(Rectangle,Rectangle): The two rectangles have different dimensions");
      }
      for (size_type i=0; i< r1.dimension(); ++i) {
        if(!Interval<R>::subset(r1[i],r2[i])) {
          return false;
        }
      }
      return true;
    }
    
    
    template <typename R>
    inline
    Rectangle<R> 
    Rectangle<R>::intersection(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("Cannot intersect with a rectangle of a different dimension.");
      }
      for(size_type i=0; i != C.dimension(); ++i) {
        C[i]=Interval<R>::intersection(A[i],B[i]);
      }
      return C;
    }

    template <typename R>
    inline
    Rectangle<R> 
    Rectangle<R>::regular_intersection(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("Cannot intersect with a rectangle of a different dimension.");
      }
      for(size_type i=0; i != C.dimension(); ++i) {
        C[i]=Interval<R>::regular_intersection(A[i],B[i]);
      }
      return C;
    }

    template <typename R>
    inline
    Rectangle<R>
    Rectangle<R>::rectangular_hull(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("The two parameters have different space dimensions");
      }      
      for(size_type i=0; i != C.dimension(); ++i) {
        C[i]=Interval<R>::hull(A[i],B[i]);
      }
      return C;
    }

    template<typename R> 
    Rectangle<R> 
    Rectangle<R>::minkowski_sum(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error(
              "minkowski_sum: the two rectangles have different dimension.");
      }
      for(dimension_type i=0; i!=C.dimension(); ++i) {
        C[i]=A[i]+B[i];
      }
      return C;
    }

    template<typename R> 
    Rectangle<R> 
    Rectangle<R>::minkowski_difference(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error(
              "minkowski_difference: the two rectangles have different dimension.");
      }
      for(dimension_type i=0; i!=C.dimension(); ++i) {
          C[i]=A[i]-B[i];
      }
      return C;
    }
    

    
    /*! \brief Tests equality of \a A and \a B. */
    template <typename R>
    inline
    bool 
    equal(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      return Rectangle<R>::equal(A,B);
    }

    /*! \brief Tests disjointness of \a A and \a B. */
    template <typename R>
    inline
    bool 
    disjoint(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      return Rectangle<R>::disjoint(A,B);
    }

   
    /*! \brief Tests intersection of the interior of \a A with the interior of \a B. */
    template <typename R>
    inline
    bool 
    interiors_intersect(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      return Rectangle<R>::interiors_intersect(A,B);
    }
    
    
    /*! \brief Tests if \a A is a subset of the interior of \a B. */
    template <typename R>
    inline
    bool 
    inner_subset(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      return Rectangle<R>::inner_subset(A,B);
    }
    
    
    /*! \brief Tests if \a A is a subset of \a B. */
    template <typename R>
    inline
    bool 
    subset(const Rectangle<R>& A, const Rectangle<R>& B) 
    {
      return Rectangle<R>::subset(A,B);
    }
    
    
    
    /*! \brief Tests inclusion of \a A in \a B. */    
    template <typename R> 
    inline
    bool 
    subset(const Rectangle<R>& A, ListSet<R,Geometry::Rectangle>& B);
    
    /*! \brief Tests inclusion in the interior of \a set. */
    template <typename R> 
    inline
    bool 
    inner_subset(const Rectangle<R>& A, 
                 const ListSet<R,Geometry::Rectangle>& B);
    
    /*! \brief Tests inclusion in an open cover. */
    template <typename R> 
    inline
    bool 
    subset_of_open_cover(const Rectangle<R>& A, 
                         const ListSet<R,Geometry::Rectangle>& U);




    
    /*! \brief The intersection of \a A and \a B. */
    template <typename R>
    inline
    Rectangle<R> 
    intersection(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      return Rectangle<R>::intersection(A,B);
    }
    
    /*! \brief The closure of the intersection of the interiors of \a A and \a B. */
    template <typename R>
    inline
    Rectangle<R>
    regular_intersection(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      return Rectangle<R>::regular_intersection(A,B);
    }
    
    /*! \brief The smallest rectangle containing \a A and \a B. */
    template <typename R>
    inline
    Rectangle<R>
    rectangular_hull(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      return Rectangle<R>::rectangular_hull(A,B);
    }
    
    /*! \brief The componentwise sum of rectangles \a A and \a B. */
    template <typename R>
    inline
    Rectangle<R> 
    minkowski_sum(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      return Rectangle<R>::minkowski_sum(A,B);
    }
    
    /*! \brief The componentwise difference of rectangles \a A and \a B. */
    template<typename R> 
    inline
    Rectangle<R> 
    minkowski_difference(const Rectangle<R>& A, const Rectangle<R>& B)
    { 
      return Rectangle<R>::minkowski_difference(A,B);
    }
    

    
    
    
    /*! \brief The difference between two rectangles. */
    template<typename R>
    inline
    LinearAlgebra::IntervalVector<R> 
    operator-(const Geometry::Rectangle<R>& r1, 
              const Geometry::Rectangle<R>& r2)
    {
      assert(r1.dimension()==r2.dimension());
      
      return r1.position_vectors()-r2.position_vectors();
    }

    /*! \brief Adds a vector to a rectangle. */
    template<typename R>
    inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector<R>& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      assert(r.dimension()==v.size());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]+v(i));
      }
      return result;
    }

    /*! \brief Adds an interval vector to a rectangle. */
    template<typename R>
    inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::IntervalVector<R>& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      assert(r.dimension()==v.size());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]+v(i));
      }
      return result;
    }

    /*! \brief Subtracts a vector from a rectangle. */
    template<typename R>
    inline
    Geometry::Rectangle<R> 
    operator-(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector<R>& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      assert(r.dimension()==v.size());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]-v(i));
      }
      return result;
    }

    /*! \brief Subtracts an interval vector from a rectangle. */
    template<typename R>
    inline
    Geometry::Rectangle<R> 
    operator-(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::IntervalVector<R>& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      assert(r.dimension()==v.size());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]-v(i));
      }
      return result;
    }

    /*! \brief Scale a rectangle. */
    template<typename R>
    inline
    Geometry::Rectangle<R> 
    scale(const Geometry::Rectangle<R>& r, const R& scale_factor) {

      Geometry::Rectangle<R> result(r.dimension());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,scale_factor*r[i]);
      }

      return result;
    }
    
    template<typename R> inline 
    std::ostream& operator<<(std::ostream& os, const Rectangle<R>& r) {
      return r.write(os);
    }
    
    template<typename R> inline
    std::istream& operator>>(std::istream& is, Rectangle<R>& r) {
      return r.read(is);
    }

    
  }
}

#endif /* _ARIADNE_RECTANGLE_H */
