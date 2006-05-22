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
#include "../geometry/point.h"

namespace Ariadne {
  namespace Geometry {

    /* Forward declaration of friends. */
    template<typename R> Rectangle<R> rectangular_hull(const Rectangle<R>&, const Rectangle<R>&);
    template<typename R> Rectangle<R> regular_intersection(const Rectangle<R>&, const Rectangle<R>&);
    template<typename R> Rectangle<R> intersection(const Rectangle<R>&, const Rectangle<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const Rectangle<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Rectangle<R>&);

    /* Value holder to allow expressions of the form r[i]=[a,b]. */
    template <typename R>
    class RectangleInterval {
     public:
      RectangleInterval(Rectangle<R>& r, const size_type& n) : _r(r), _n(n) { }
      void operator=(const Interval<R>& i) { _r.set_interval(_n,i); }
      operator Interval<R> () const { return _r[_n]; }
      R lower() const { return _r.lower_bound(_n); }
      R upper() const { return _r.upper_bound(_n); }
     private:
      Rectangle<R>& _r; const size_type& _n;
    };
    
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
      /*! \brief Construct an empty rectangle of dimension \a n. */
      explicit Rectangle(size_type n=0)
        : _lower_corner(n),  _upper_corner(n) 
      { 
        if(n>0) {
          set_lower_bound(0,R(1));
          set_upper_bound(0,R(0));
        }
      }
      
      /*! \brief Construct a rectangle from a range of interval values. */
      template<class ForwardIterator>
      Rectangle(ForwardIterator b, ForwardIterator e)
        : _lower_corner(std::distance(b,e)), _upper_corner(std::distance(b,e))
      {
        for(size_type i=0; i!=dimension(); ++i) {
          _lower_corner[i]=b->lower();
          _upper_corner[i]=b->upper();
          ++b;
        }
      }
      
      /*! \brief Construct from an array of intervals. */
      explicit Rectangle(const array< Interval<real_type> >& a);
      
      /*! \brief Construct from a std::vector of intervals. */
      explicit Rectangle(const std::vector< Interval<real_type> >& v);
      
      /*! \brief Construct from two corners. */
      explicit Rectangle(const state_type& p1, const state_type& p2);
      
      /*! \brief Construct from a string literal. */
      explicit Rectangle(const std::string& s);
      
      /*! \brief Construct from an interval vector. */
      explicit Rectangle(const LinearAlgebra::IntervalVector<R>& iv);
      
      /*! \brief Copy constructor. */
      Rectangle(const Rectangle<R>& original)
        : _lower_corner(original._lower_corner),
          _upper_corner(original._upper_corner)
      { }
    
      /*! \brief The \a i th vertex. */
      state_type vertex(size_type i) const 
      {
        size_type dim=this->_lower_corner.dimension();
	state_type output(dim); 
	      
        if (i >= (size_type)(1<<dim))
           throw std::domain_error("Rectangle::vertex(i): Wrong vertex index.");
	     
	for (size_type j=0; j<dim; j++) {
          if (i%2)
	    output[j]=this->_lower_corner[j];
	  else
	    output[j]=this->_upper_corner[j];
	  i=i/2;
        }

	return output;
      }
      
      /*! \brief Convert to a paralletope. */
      operator Parallelotope<R> () const;
      
      /*! \brief Convert to a zonotope. */
      operator Zonotope<R> () const;
      
      /*! \brief Convert to a polyhedron. */
      operator Polyhedron<R> () const;
      
      /*! \brief Copy assignment operator. */
      Rectangle<R>& operator=(const Rectangle<R>& original) {
        if(this != &original) {
          this->_lower_corner = original._lower_corner;
          this->_upper_corner = original._upper_corner;
        }
        return *this;
      }
      
      /*! \brief The equality operator */
      bool operator==(const Rectangle<real_type>& A) const
      {
        if (A.empty() && this->empty()) { return true; }
        if (A.empty() || this->empty()) { return false; }
        if (this->dimension() != A.dimension()) { return false ; }
        for (size_type j=0; j != this->dimension(); ++j) {
          if (this->_lower_corner[j] != A._lower_corner[j]) { return false; }
          if (this->_upper_corner[j] != A._upper_corner[j]) { return false; }
        }
        return true;
      }
      
      /*! \brief The inequality operator */
      bool operator!=(const Rectangle<real_type>& A) const {
        return !(*this == A);
      }

      /*! \brief The dimension of the Euclidean space the rectangle lies in. */
      size_type dimension() const {
        return (this->_lower_corner).dimension();
      }
      
      /*! \brief True if the rectangle is empty. */
      bool empty() const {
        if(this->dimension()==0) {
          return true;
        }
        for(size_type i=0; i!=this->dimension(); ++i) {
          if(this->lower_bound(i) > this->upper_bound(i)) {
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
      state_type lower_corner() const {
        return this->_lower_corner;
      }
      
      /*! \brief The upper corner. */
      state_type upper_corner() const {
        return this->_upper_corner;
      }
      
      /*! \brief The centre. */
      state_type centre() const {
        return state_type((this->_lower_corner.position_vector()+this->_upper_corner.position_vector())/2);
      }
      
      /*! \brief The radius in the sup norm. */
      real_type radius() const {
        real_type diameter=0;
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          diameter=max(diameter,this->interval(i).length());
        }
        return diameter/2;
      }
      
      /*! \brief Returns the projection onto the \a n th coordinate. */
      Interval<real_type> operator[] (size_type n) const {
        return Interval<real_type>(this->_lower_corner[n],this->_upper_corner[n]);
      }
      
      /*! \brief Returns the projection onto the \a n th coordinate. */
      Interval<real_type> interval(size_type n) const {
        return Interval<real_type>(this->_lower_corner[n],this->_upper_corner[n]);
      }
      
      /*! \brief Returns the lower bound of the \a n th coordinate */
      const real_type& lower_bound(size_type n) const {
        return this->_lower_corner[n];
      }
      
      /*! \brief Returns the upper bound of the \a n th coordinate */
      const real_type& upper_bound(size_type n) const {
        return this->_upper_corner[n];
      }
      
      /*! \brief Returns a smart reference projection onto the \a n th coordinate. */
      RectangleInterval<real_type> operator[] (size_type n) {
        return RectangleInterval<real_type>(*this,n);
      }
      
      /*! \brief Sets the \a n th interval. */
      void set_interval(size_type n, Interval<real_type> i) {
         this->_lower_corner[n]=i.lower();
         this->_upper_corner[n]=i.upper();
      }
      
      /*! \brief Sets the lower bound of the \a n th coordinate to \a r. */
      void set_lower_bound(size_type n, const real_type& r) {
        this->_lower_corner[n] = r;
      }
      
      /*! \brief Sets the upper bound of the \a n th coordinate to \a r. */
      void set_upper_bound(size_type n, const real_type& r) {
        this->_upper_corner[n] = r;
      }
      
      
      /*! \brief The set of position Vectors of the rectangle. */
      LinearAlgebra::IntervalVector<R> position_vectors() const;
      
      /*! \brief Expand the Rectangle by \a delta in each direction. */
      Rectangle<R>& expand_by(const real_type& delta);
      

      
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
      array<state_type> vertices() const;
      
      /*! \brief Tests if \a point is included into a rectangle. */
      bool contains(const state_type& p) const;
     
      /*! \brief Tests if the rectangle contains a \a rectangle. */
      bool contains(const Rectangle<R>& rect) const; 
      
      /*! \brief Tests if \a point is included into the interior a rectangle. */
      bool interior_contains(const state_type& p) const;

      /*! \brief Tests inclusion in the interior of \a set. */
      bool inner_subset(const ListSet<R,Geometry::Rectangle>& set) const;
      /*! \brief Tests inclusion in \a set. */
      bool subset(const ListSet<R,Geometry::Rectangle>& set) const;
      /*! \brief Tests inclusion in an open cover. */
      bool subset_of_open_cover(const ListSet<R,Geometry::Rectangle>& cover) const;
     private:
      friend Rectangle<R> rectangular_hull<>(const Rectangle<R>&, const Rectangle<R>&);
      friend Rectangle<R> regular_intersection<>(const Rectangle<R>&, const Rectangle<R>&);
      friend Rectangle<R> intersection<>(const Rectangle<R>&, const Rectangle<R>&);
      friend std::ostream& operator<< <> (std::ostream& os, const Rectangle<R>& r);
      friend std::istream& operator>> <> (std::istream& is, Rectangle<R>& r);
     private:
      /* Lower corner */
      state_type _lower_corner;
      /* Upper corner */
      state_type _upper_corner;
    };
  

    template <typename R>
    inline
    bool 
    Rectangle<R>::contains(const state_type& p) const 
    {
      if (this->empty()) {
        return false;
      }
      if (p.dimension()!=this->dimension()) {
        throw std::domain_error("This object and parameter have different space dimensions");
      }
      for (size_type i=0; i<this->dimension(); i++) {
        if (p[i] < this->_lower_corner[i]) { return false; }
        if (p[i] > this->_upper_corner[i]) { return false; }
      }
      return true;
    }

    template <typename R>
    inline
    bool 
    Rectangle<R>::contains(const Rectangle& r) const 
    {
      if (this->empty()) {
        return false;
      }
      if (r.dimension()!=this->dimension()) {
        throw std::domain_error("This object and parameter have different space dimensions");
      }
      for (size_type i=0; i<this->dimension(); i++) {
        if (r._lower_corner[i] < this->_lower_corner[i]) { return false; }
        if (r._upper_corner[i] > this->_upper_corner[i]) { return false; }
      }
      return true;
    }

    template <typename R>
    inline
    bool 
    Rectangle<R>::interior_contains(const state_type& p) const 
    {
      if (this->empty()) {
        return false;
      }
      if (p.dimension()!=this->dimension()) {
        throw std::domain_error("This object and parameter have different space dimensions");
      }
      for (size_type i=0; i<this->dimension(); i++) {
        if (p[i] >= this->_upper_corner[i]) { return false; }
        if (p[i] <= this->_lower_corner[i]) { return false; }
      }
      return true;
    }
      
    template <typename R>
    inline
    bool 
    disjoint(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      if (A.empty() || B.empty()) {
        return true;
      }
      if(A.dimension()!=B.dimension()) {
        throw std::domain_error("The two parameters have different space dimensions");
      }
      for(size_type i=0; i< A.dimension(); i++) {
        if ((A.upper_bound(i)<B.lower_bound(i))|| 
            (B.upper_bound(i)<A.lower_bound(i))) 
        {
          return true;
        }
      }
      return false;
    }
   
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline
    bool 
    interiors_intersect(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      if (A.empty() || B.empty()) {
        return false;
      }
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("The two parameters have different space dimensions");
      }
      for(size_type i=0; i< A.dimension(); i++) {
        if ((A.upper_bound(i)<=B.lower_bound(i))|| 
            (B.upper_bound(i)<=A.lower_bound(i))) 
        {
          return false;
        }
      }
      return true;
    }
    
    
    template <typename R>
    inline
    bool 
    inner_subset(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      if(A.empty()) {
        return true;
      }
      if(B.empty()) {
        return false;
      }
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("The two parameters have different space dimensions");
      }
      for (size_type i=0; i< A.dimension(); i++) {
        if((A.lower_bound(i) >= B.lower_bound(i))||
            (B.upper_bound(i) >= A.upper_bound(i))) 
        {
          return false;
        }
      }
      return true;
    }

    

    /*! \brief Tests inclusion */
    template <typename R>
    inline
    bool 
    subset(const Rectangle<R>& A, const Rectangle<R>& B) 
    {
      if(A.empty()) {
        return true;
      }
      if(B.empty()) {
        return false;
      }
      if(A.dimension()!=B.dimension()) {
        throw std::domain_error("The two parameters have different space dimensions");
      }
      for(size_type i=0; i< A.dimension(); i++) {
        if((A.lower_bound(i) > B.lower_bound(i))||
            (B.upper_bound(i) > A.upper_bound(i))) 
        {
          return false;
        }
      }
      return true;
    }
    
    template <typename R>
    inline
    bool 
    inner_subset(const Rectangle<R>& A, const ListSet<R,Rectangle>& B)
    {
      return A.inner_subset(B);
    }
    
    template <typename R>
    bool 
    subset(const Rectangle<R>& A, const ListSet<R,Rectangle>& B)
    { 
      return A.subset(B);
    }
    
    /*! \brief The intersections of \a A and \a B. */
    template <typename R>
    inline
    Rectangle<R> 
    intersection(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("The two parameters have different space dimensions");
      }

      Rectangle<R> C(A.dimension());

      for(size_type i=0; i != C.dimension(); ++i) {
        C._lower_corner[i] = max(A._lower_corner[i],B._lower_corner[i]);
        C._upper_corner[i] = min(A._upper_corner[i],B._upper_corner[i]);
        if(C._lower_corner[i] > C._upper_corner[i]) {
          C._lower_corner[0]=1;
          C._upper_corner[0]=0;
          return C;
        }
      }
      return C;
    }

    /*! \brief The closure of the intersection of the interiors of \a A and \a B. */
    template <typename R>
    inline
    Rectangle<R>
    regular_intersection(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("The two parameters have different space dimensions");
      }

      Rectangle<R> C(A.dimension());

      if (A.empty() || B.empty()) {
        return C;
      }

      for(size_type i=0; i != C.dimension(); ++i) {
        C._lower_corner[i] = std::max(A._lower_corner[i],B._lower_corner[i]);
        C._upper_corner[i] = std::min(A._upper_corner[i],B._upper_corner[i]);
        if(C._lower_corner[i] >= C._upper_corner[i]) {
          C._lower_corner[0]=1;
          C._upper_corner[0]=0;
          return C;
        }
      }
      return C;
    }

    /*! \brief The smallest rectangle containing \a A and \a B. */
    template <typename R>
    inline
    Rectangle<R>
    rectangular_hull(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("The two parameters have different space dimensions");
      }
      if(A.empty()) { 
        return B;
      }
      if(B.empty()) {
        return A;
      }
      
      Rectangle<R> C(A.dimension());
      for(size_type i=0; i != C.dimension(); ++i) {
        C._lower_corner[i] = std::min(A._lower_corner[i],B._lower_corner[i]);
        C._upper_corner[i] = std::max(A._upper_corner[i],B._upper_corner[i]);
      }
      return C;
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

  }
}

#endif /* _ARIADNE_RECTANGLE_H */
