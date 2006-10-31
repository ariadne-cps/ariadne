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

#include "../base/array.h"
#include "../base/iterator.h"
#include "../base/tribool.h"
#include "../exceptions.h"

#include "../numeric/arithmetic.h"
#include "../numeric/function.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"

#include "../geometry/point.h"
#include "../geometry/rectangle_expression.h"

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

    template<class R> class RectangleVerticesIterator;
    
    /*! \ingroup BasicSet
     *  \brief A cuboid of arbitrary dimension.
     * 
     * The most important geomtric object in %Ariadne, the %Rectangle class 
     * describes a rectangle, cuboid or hypercuboid in Euclidean space.
     * Rectangles are closed under most major geometric operations, and geometric
     * predicates and operations are simple to compute.
     * Rectangular sets require little data to describe, and are used as cells in grid-based and lattice based sets,
     * making them ubiquitous in storage-representations.
     * Further, rectangles are easily and exactly convertible to the other major polyhedral
     * set representations, namely, Zonotope,  Polytope and Polyhedron.
     * 
     * Rectangles are decribed by the lower and upper bounds in each of their
     * dimensions, as accessed by the lower_bound(dimension_type) const and upper_bound(dimension_type) const methods. 
     * 
     * Rectangles are by default ordered by the lexicographic order on lower-left corner.
     * If these are equal, then points are ordered by the lexicographic order 
     * on the upper-right corner.
     *
     * A rectangle is described by the literal "[a0,b0]x[a1,b1]x...".
     * 
     * Rectangles are always bounded, but may be empty. A zero-dimensional rectangle is considered to be non-empty.
     *
     * \b Storage: A %Rectangle of dimension \a d is specified by \a 2d real 
     * number giving the lower and upper bounds in each coordinate.
     * The lower bound in the \a i th coordinate is the \a 2i th element, and the 
     * upper bound is the \a 2i+1 st element.
     */
    template<class R>
    class Rectangle 
      : public RectangleExpression< Rectangle<R> >
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the rectangle. */
      typedef Point<R> state_type;
      /*! \brief An iterator to the vertices of the rectangle. */
      typedef RectangleVerticesIterator<R> vertices_iterator;
      typedef RectangleVerticesIterator<R> vertices_const_iterator;
     public:
      //@{
      //! \name Constructors
      /*! \brief Construct an empty rectangle with dimension \a d. */
      explicit Rectangle(size_type d=0)
        : _bounds(2*d)
      { 
        if(d!=0) { this->_bounds[0]=1; this->_bounds[1]=0; }
      }
      
      /*! \brief Construct a rectangle from a range of interval values. */
      template<class ForwardIterator>
      Rectangle(ForwardIterator b, ForwardIterator e)
        : _bounds(2*std::distance(b,e))
      {
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          this->set_lower_bound(i,b->lower());
          this->set_upper_bound(i,b->upper());
          ++b;
        }
      }
      
      /*! \brief Construct from an array of intervals. */
      explicit Rectangle(const array< Interval<R> >& a)
        : _bounds(2*a.size())
      {
        for(dimension_type i=0; i!=a.size(); ++i) {
          this->set_lower_bound(i,a[i].lower());
          this->set_upper_bound(i,a[i].upper());
        }
      }
      
      /*! \brief Construct from a std::vector of intervals. */
      explicit Rectangle(const std::vector< Interval<R> >& v)
        : _bounds(2*v.size())
      {
        for(dimension_type i=0; i!=v.size(); ++i) {
          this->set_lower_bound(i,v[i].lower());
          this->set_upper_bound(i,v[i].upper());
        }
      }

      /*! \brief Construct a degenerate rectangle from a single point. */
      explicit Rectangle(const Point<R>& pt)
        : _bounds(2*pt.dimension())
      {
        for(dimension_type i=0; i!=pt.dimension(); ++i) {
          this->set_lower_bound(i,pt[i]);
          this->set_upper_bound(i,pt[i]);
        }
      }
      
      /*! \brief Construct a rectangle from an interval point. */
      explicit Rectangle(const Point<I>& pt)
        : _bounds(2*pt.dimension())
      {
        for(dimension_type i=0; i!=pt.dimension(); ++i) {
          this->set_lower_bound(i,pt[i].lower());
          this->set_upper_bound(i,pt[i].upper());
        }
      }
      
      /*! \brief Construct from two corners. */
      explicit Rectangle(const Point<R>& pt1, const Point<R>& pt2) 
        : _bounds(2*pt1.dimension())
      {
        check_equal_dimensions(pt1,pt2,__PRETTY_FUNCTION__);
        for (size_type i=0; i!=this->dimension(); ++i) {
          this->set_lower_bound(i,Numeric::min_exact(pt1[i],pt2[i]));
          this->set_upper_bound(i,Numeric::max_exact(pt1[i],pt2[i]));
        }
      }
      
      /*! \brief Construct from a string literal. */
      explicit Rectangle(const std::string& s);
      
      /*! \brief Construct from an interval vector. */
      explicit Rectangle(const LinearAlgebra::Vector< Interval<R> >& iv)
        : _bounds(2*iv.size())
      {
        for (size_type i=0; i!=this->dimension(); ++i) {
          this->set_lower_bound(i,iv[i].lower());
          this->set_upper_bound(i,iv[i].upper());
        }
      }

      /*! \brief Convert from a rectangle expression. */
      template<class E>
      Rectangle(const RectangleExpression<E>& original)
        : _bounds(2*original().dimension())
      {         
        const E& expression=original();
        for (size_type i=0; i!=this->dimension(); ++i) {
          this->_bounds[2*i]=expression.lower_bound(i);
          this->_bounds[2*i+1]=expression.upper_bound(i);
        }
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

      /*! \brief Assign from a rectangle expression. */
      template<class E>
      Rectangle<R>& operator=(const RectangleExpression<E>& original)
      {         
        const E& expression=original();
        this->_bounds.resize(2*expression.dimension());
        for (size_type i=0; i!=this->dimension(); ++i) {
          this->set_lower_bound(i,expression.lower_bound(i));
          this->set_upper_bound(i,expression.upper_bound(i));
        }
        return *this;
      }
    
      //@}
      
      
      //@{
      //! \name Conversion operators
      /*! \brief Convert to an interval point. */
      operator Point< Interval<R> >() const {
        return Point< Interval<R> >(this->dimension(),reinterpret_cast<const Interval<R>*>(this->_bounds.begin()));
      }
      //@}
      
      
      //{@
      //! \name Comparison operators
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
      //@}
      

      //@{
      //! \name Data access
      /*! \brief Returns the projection onto the \a i th coordinate. */
      Interval<R>& operator[] (dimension_type i) {
        check_coordinate(*this,i,__PRETTY_FUNCTION__);
        return reinterpret_cast<Interval<R>&>(this->_bounds[2*i]);
      }
      //IntervalReference<R> operator[] (dimension_type i) {
      //  return IntervalReference<R>(this->_bounds[i]);
      //}
      
      /*! \brief The lower bound of the \a i th coordinate */
      const R& lower_bound(dimension_type i) const {
        check_coordinate(*this,i,__PRETTY_FUNCTION__);
        return this->_bounds[2*i];
      }
      
      /*! \brief A reference to the lower bound of the \a i th coordinate */
      R& lower_bound(dimension_type i) {
        check_coordinate(*this,i,__PRETTY_FUNCTION__);
        return this->_bounds[2*i];
      }
      
      /*! \brief The upper bound of the \a i th coordinate */
      const R& upper_bound(dimension_type i) const {
        check_coordinate(*this,i,__PRETTY_FUNCTION__);
        return this->_bounds[2*i+1];
      }
      
      /*! \brief A reference to the upper bound of the \a i th coordinate */
      R& upper_bound(dimension_type i) {
        check_coordinate(*this,i,__PRETTY_FUNCTION__);
        return this->_bounds[2*i+1];
      }
      
      /*! \brief The projection onto the \a i th coordinate. */
      const Interval<R>& operator[] (dimension_type i) const {
        check_coordinate(*this,i,__PRETTY_FUNCTION__);
        return reinterpret_cast<const Interval<R>&>(this->_bounds[2*i]);
      }
      
      /*! \brief The interval of values in the \a i th coordinate. */
      const Interval<R>& interval(dimension_type i) const {
        check_coordinate(*this,i,__PRETTY_FUNCTION__);
        return reinterpret_cast<const Interval<R>&>(this->_bounds[2*i]);
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
      
      /*! \brief The set of position vectors of the rectangle. */
      LinearAlgebra::Vector< Interval<R> > position_vectors() const {
        LinearAlgebra::Vector< Interval<R> > result(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          result[i]=this->interval(i);
        }
        return result;
      }
      //@}
      
      
      //@{ 
      //! \name Modifying operations
      /*! \brief Makes the rectangle empty. */
      void clear() {
        if(this->_bounds.size()!=0) {
          this->_bounds[0]=1;
          this->_bounds[1]=0;
        }
      }
      
      /*! \brief Sets the \a i th interval. */
      void set_interval(dimension_type i, Interval<R> x) {
        check_coordinate(*this,i,__PRETTY_FUNCTION__);
        this->set_lower_bound(i,x.lower());
        this->set_upper_bound(i,x.upper());
      }
      
      /*! \brief Sets the lower bound of the \a i th coordinate to \a r. */
      void set_lower_bound(dimension_type i, const R& l) {
        check_coordinate(*this,i,__PRETTY_FUNCTION__);
        this->_bounds[2*i]=l;
      }
      
      /*! \brief Sets the upper bound of the \a i th coordinate to \a u. */
      void set_upper_bound(dimension_type i, const R& u) {
        check_coordinate(*this,i,__PRETTY_FUNCTION__);
        this->_bounds[2*i+1]=u;
      }

      /*! \brief Expand the Rectangle by \a delta in each direction. */
      Rectangle<R>& expand_by(const R& delta);
      //@}
      
      
      //@{
      //! \name Rectangle geometric operations
      /*! \brief The dimension of the Euclidean space the rectangle lies in. */
      size_type dimension() const {
        return this->_bounds.size()/2;
      }
      
      /*! \brief True if the rectangle is empty. A zero-dimensional rectangle is considered empty. */
      tribool empty() const {
        tribool result=false;
        if(this->dimension()==0) {
          return true;
        }
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          if(this->lower_bound(i) > this->upper_bound(i)) {
            return true;
          }
          if(this->lower_bound(i)== this->upper_bound(i)) {
            result=indeterminate;
          }
        }
        return result;
      }
      
      /*! \brief The centre. */
      Point<R> centre() const {
        Point<R> result(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          result[i]=Numeric::div_approx(Numeric::add_approx(this->lower_bound(i),this->upper_bound(i)),R(2));
        }
        return result;
      }
      
      /*! \brief The radius in the sup norm. */
      R radius() const {
        R diameter=0;
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          diameter=Numeric::max_up(diameter,Numeric::sub_up(this->upper_bound(i),this->lower_bound(i)));
        }
        return div_up(diameter,R(2));
      }
      
      /*! An approximation to the volume. */
      R volume() const {
        R result=1;
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          result=mul_approx(result,sub_approx(this->upper_bound(i),this->lower_bound(i)));
        }
        return result;
      }
      
      /*! \brief Compute a quadrant of the Rectangle determined by \a q.
       *  \a q is a binary word such that the ith bit of q is 0 if the lower half
       *  of the rectangle in the ith coordinate is used, and 1 if the upper
       *  half is used.
       */
      Rectangle<R> quadrant(const Combinatoric::BinaryWord& q) const;
      /*! \brief Subdivide into smaller pieces. */
      ListSet<R,Geometry::Rectangle> subdivide() const;
      
      /*! \brief The vertices of the rectangle. */
      PointList<R> vertices() const;

      /*! \brief An iterator to the first vertex of the rectangle. */
      vertices_const_iterator vertices_begin() const;
      /*! \brief An iterator to the end vertex of the rectangle. */
      vertices_const_iterator vertices_end() const;

      /*! \brief The number of vertices. */
      size_type number_of_vertices() const;
        
      /*! \brief The \a i th vertex. */
      Point<R> vertex(size_type i) const;
        
      /*! \brief Checks for boundedness. */
      tribool bounded() const { return true; }
      
      /*! \brief Tests if \a point is included into a rectangle. */
      tribool contains(const Point<R>& pt) const;
      
      /*! \brief A rectangle containing the given rectangle; returns a copy. */
      Rectangle bounding_box() const {
        return *this;
      }
      //@}
      
#ifdef DOXYGEN
      //@{ 
      //! \name Binary geometric predicates
      /*! \brief Set equality operator. */
      friend tribool equal(const Rectangle<R>& A, const Rectangle<R>& B) const;
     
      /*! \brief Tests disjointness with \a r. */
      friend tribool disjoint(const Rectangle<R>& A,__PRETTY_FUNCTION__
 const Rectangle<R>& B) const;
      /*! \brief Tests if the rectangle is a subset of another rectangle \a r. */
      friend tribool subset(const Rectangle<R>& A, const Rectangle<R>& B) const;
      //@{ 
      //! \name Binary geometric operations
      /*! \brief The intersection of \a A and \a B. */
      friend Rectangle<R> intersection(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The closure of the intersection of the interiors of \a A and \a B. */
      friend Rectangle<R> regular_intersection(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The smallest rectangle containing \a A and \a B. */
      friend Rectangle<R> rectangular_hull(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The componentwise sum of rectangles \a A and \a B. */
      friend Rectangle<R> minkowski_sum(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The componentwise difference of rectangles \a A and \a B. */
      friend Rectangle<R> minkowski_difference(const Rectangle<R>& A, const Rectangle<R>& B); 
      
      /*! \brief The difference between two rectangles. */
      friend LinearAlgebra::Vector< Interval<R> > operator-(const Rectangle<R>& A, 
                                                        const Rectangle& B);
      /*! \brief Adds a vector to a rectangle. */
      friend Rectangle<R> operator+(const Rectangle<R>& r, 
                                    const LinearAlgebra::Vector<R>& v);
      /*! \brief Adds an interval vector to a rectangle. */
      friend Rectangle<R> operator+(const Rectangle<R>& r, 
                                    const LinearAlgebra::Vector< Interval<R> >& v);
      /*! \brief Subtracts a vector from a rectangle. */
      friend Rectangle<R> operator-(const Rectangle<R>& r, 
                                    const LinearAlgebra::Vector<R>& v);
      /*! \brief Subtracts an interval vector from a rectangle. */
      friend Rectangle<R> operator-(const Rectangle<R>& r, 
                                    const LinearAlgebra::Vector< Interval<R> >& v);
      //@}
#endif
      
      //@{ 
      //! \name Input/output operations
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}
     private:
      array<R> _bounds;
    };
  

    template<class R>
    class Rectangle< Interval<R> > 
      : public RectangleExpression< Rectangle< Interval<R> > >
    {
      typedef Interval<R> I;
     public:
      Rectangle(dimension_type d) : _bounds(2*d) { }
      template<class E> Rectangle(const RectangleExpression<E>& e)
        : _bounds(2*e().dimension()) { this->assign(e()); }
      template<class E> Rectangle< Interval<I> > operator=(const RectangleExpression<E>& e) {
         this->_bounds.resize(e().dimension()); this->assign(e); }
      dimension_type dimension() const { return _bounds.size()/2; }
      const I& lower_bound(const dimension_type& i) const { return _bounds[2*i]; }
      const I& upper_bound(const dimension_type& i) const { return _bounds[2*i+1]; }
      void set_lower_bound(const dimension_type& i, const I& x) { _bounds[2*i]=x; }
      void set_upper_bound(const dimension_type& i, const I& x) { _bounds[2*i+1]=x; }
     private:
      template<class RE> void assign(const RE& re) { 
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          this->_bounds[2*i]=re.lower_bound(i); this->_bounds[2*i+1]=re.upper_bound(i);
        }
      }
     private:
      array<I> _bounds;
    };
    
    
    template<class R> inline
    Rectangle<R> over_approximation(const Rectangle< Numeric::Interval<R> >& ir) {
      Rectangle<R> result(ir.dimension());
      for(dimension_type i=0; i!=result.dimension(); ++i) {
        result.set_lower_bound(i,ir.lower_bound(i).lower());
        result.set_upper_bound(i,ir.upper_bound(i).upper());
      }
      return result;
    }
    
    template<class R> inline
    Rectangle<R> under_approximation(const Rectangle< Numeric::Interval<R> >& ir) {
      Rectangle<R> result(ir.dimension());
      for(dimension_type i=0; i!=result.dimension(); ++i) {
        result.set_lower_bound(i,ir.lower_bound(i).upper());
        result.set_upper_bound(i,ir.upper_bound(i).lower());
      }
      return result;
    }
    
    template<class R>
    inline
    tribool 
    Rectangle<R>::contains(const Point<R>& p) const 
    {
      tribool result=true;
      check_equal_dimensions(*this,p,__PRETTY_FUNCTION__);
      const Rectangle<R>& self=*this;
      for (size_type i=0; i!=self.dimension(); ++i) {
        if(self.lower_bound(i)>p[i] || p[i]>self.upper_bound(i)) {
          return false;
        }
        if(self.lower_bound(i)==p[i] || p[i]==self.upper_bound(i)) { 
          result=indeterminate;
        }
      }
      return result;
    }

      
    template<class R>
    inline
    tribool 
    equal(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(size_type i=0; i!=A.dimension(); ++i) {
        if(A.lower_bound(i)!=B.lower_bound(i) || A.upper_bound(i)!=B.upper_bound(i)) {
          return false;
        }
      }
      return indeterminate;
    }
      
    template<class R>
    inline
    tribool 
    disjoint(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      tribool result=false;
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(size_type i=0; i!=A.dimension(); ++i) {
        if(A.lower_bound(i)>B.upper_bound(i) || A.upper_bound(i)<B.lower_bound(i)) {
          return true;
        }
        if(A.lower_bound(i)==B.upper_bound(i) || A.upper_bound(i)==B.lower_bound(i)) {
          result=indeterminate;
        }
      }
      return result;
    }
       

    template<class R>
    inline
    tribool 
    subset(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      tribool result=true;
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for (size_type i=0; i!=A.dimension(); ++i) {
        if(A.lower_bound(i)<B.lower_bound(i) || A.upper_bound(i)>B.upper_bound(i)) {
          return false;
        }
        if(A.lower_bound(i)==B.lower_bound(i) || A.upper_bound(i)==B.upper_bound(i)) {
          result=indeterminate;
        }
      }
      return result;
    }
    
    
    template<class R>
    inline
    Rectangle<R> 
    closed_intersection(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(size_type i=0; i != C.dimension(); ++i) {
        C[i]=Numeric::intersection(A[i],B[i]);
      }
      return C;
    }

    template<class R>
    inline
    Rectangle<R> 
    open_intersection(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(size_type i=0; i != C.dimension(); ++i) {
        C[i]=Numeric::intersection(A[i],B[i]);
        if(C[i].lower()>=C[i].upper()) {
          C[i]=Interval<R>();
        }
      }
      return C;
    }

    template<class R>
    inline
    Rectangle<R>
    rectangular_hull(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(size_type i=0; i != C.dimension(); ++i) {
        C[i]=Numeric::hull(A[i],B[i]);
      }
      return C;
    }

    template<class R1, class R2> 
    Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> 
    minkowski_sum(const Rectangle<R1>& A, const Rectangle<R2>& B)
    {
      Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> C(A.dimension());
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(dimension_type i=0; i!=C.dimension(); ++i) {
        C.set_lower_bound(i,A.lower_bound(i)+B.lower_bound(i));
        C.set_lower_bound(i,A.upper_bound(i)+B.upper_bound(i));
      }
      return C;
    }

    template<class R1, class R2> 
    Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> 
    minkowski_difference(const Rectangle<R1>& A, const Rectangle<R2>& B)
    {
      Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> C(A.dimension());
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(dimension_type i=0; i!=C.dimension(); ++i) {
        C.set_lower_bound(i,A.lower_bound(i)-B.lower_bound(i));
        C.set_upper_bound(i,A.upper_bound(i)-B.upper_bound(i));
      }
      return C;
    }
    
    
    template<class R> 
    inline
    tribool 
    subset(const Rectangle<R>& A, ListSet<R,Geometry::Rectangle>& B);
        

    
    
    
    template<class R>
    inline
    LinearAlgebra::Vector< Interval<R> > 
    operator-(const Geometry::Rectangle<R>& r1, 
              const Geometry::Rectangle<R>& r2)
    {
      check_equal_dimensions(r1,r2,__PRETTY_FUNCTION__);
       
      return r1.position_vectors()-r2.position_vectors();
    }

    template<class R, class E>
    inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const boost::numeric::ublas::vector_expression<E>& v)
    {
      const E& ev=v();
      check_dimension(r,ev.size(),__PRETTY_FUNCTION__);
      LinearAlgebra::Vector< Interval<R> > iv=ev;
      return r+iv; 
    }
      
    template<class R>
    inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector<R>& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      check_dimension(r,v.size(),__PRETTY_FUNCTION__);
       
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]+v(i));
      }
      return result;
    }

    template<class R>
    inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector< Interval<R> >& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      check_dimension(r,v.size(),__PRETTY_FUNCTION__);
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]+v(i));
      }
      return result;
    }

    template<class R>
    inline
    Geometry::Rectangle<R> 
    operator-(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector<R>& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      check_dimension(r,v.size(),__PRETTY_FUNCTION__);
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]-v(i));
      }
      return result;
    }

    template<class R>
    inline
    Geometry::Rectangle<R> 
    operator-(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector< Interval<R> >& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      check_dimension(r,v.size(),__PRETTY_FUNCTION__);
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]-v(i));
      }
      return result;
    }

    template<class R>
    inline
    Geometry::Rectangle<R> 
    scale(const Geometry::Rectangle<R>& r, const R& scale_factor) 
    {
      Geometry::Rectangle<R> result(r.dimension());
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,scale_factor*r[i]);
      }
      return result;
    }
    
    template<class R>
    class RectangleVerticesIterator 
      : public boost::iterator_facade<RectangleVerticesIterator<R>,
                                      Point<R>,
                                      boost::forward_traversal_tag,
                                      Point<R> const&,
                                      Point<R> const*
                                     >
    {
     public:
      RectangleVerticesIterator(const Rectangle<R>& r, const bool end);
      bool equal(const RectangleVerticesIterator<R>& other) const;
      const Point<R>& dereference() const;
      void increment();
     private:
      const Rectangle<R>* _r; long unsigned int _i; bool _parity; Point<R> _pt;
    };
    
    template<class R> inline
    RectangleVerticesIterator<R>::RectangleVerticesIterator(const Rectangle<R>& r, const bool end)
      : _r(&r), _i(end==true ? (1<<(r.dimension()-1))*3 : 0), _parity(0), _pt(r.lower_corner()) { }
      
    template<class R> inline
    bool RectangleVerticesIterator<R>::equal(const RectangleVerticesIterator<R>& other) const {
      return this->_i==other._i && this->_r==other._r; }
      
    template<class R> inline
    const Point<R>& RectangleVerticesIterator<R>::dereference() const { 
      return this->_pt; }
      
    template<class R> inline
    void RectangleVerticesIterator<R>::increment() { 
      uint j=0; uint m=1; if(this->_parity) { while(!(m&(this->_i))) { ++j; m*=2u; } ++j; m*=2u; }
      this->_parity=!this->_parity;
      if(j==this->_r->dimension()) { this->_i+=m; return; }
      if(m&(this->_i)) { this->_pt[j]=this->_r->lower_bound(j); this->_i-=m; }
      else { this->_pt[j]=this->_r->upper_bound(j); this->_i+=m; }
    }
    
    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const Rectangle<R>& r) {
      return r.write(os);
    }
    
    template<class R> inline
    std::istream& operator>>(std::istream& is, Rectangle<R>& r) {
      return r.read(is);
    }

    
  }
}

#endif /* _ARIADNE_RECTANGLE_H */
