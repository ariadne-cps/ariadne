/***************************************************************************
 *            box.h
 *
 *  Copyright 2005-7  Alberto Casagrande, Pieter Collins
 *
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
 
/*! \file box.h
 *  \brief Boxes in Euclidean space.
 */

#ifndef ARIADNE_GEOMETRY_BOX_H
#define ARIADNE_GEOMETRY_BOX_H

#include <iosfwd>
#include <boost/concept_check.hpp>

#include "base/array.h"
#include "base/iterator.h"
#include "base/tribool.h"

#include "numeric/interval.h"

#include "linear_algebra/vector.h"

#include "combinatoric/declarations.h"

#include "geometry/exceptions.h"
#include "geometry/point.h"
#include "geometry/basic_set_concept.h"
#include "geometry/box_expression.h"

namespace Ariadne {
  

    class basic_set_tag;
    class EuclideanSpace;
    template<class R> class BoxVerticesIterator;
    template<class R> class PointList;
    template<class BS> class ListSet;
    
    /*! \brief A box of arbitrary dimension.
     * 
     * The most important geomtric object in %Ariadne, the %Box class 
     * describes a box, cuboid or hypercuboid in Euclidean space.
     * Boxes are closed under most major geometric operations, and geometric
     * predicates and operations are simple to compute.
     * Boxes sets require little data to describe, and are used as cells in grid-based and lattice based sets,
     * making them ubiquitous in storage-representations.
     * Further, boxes are easily and exactly convertible to the other major polyhedral
     * set representations, namely, Box, Zonotope,  Polytope and Polyhedron.
     * 
     * Boxes are decribed by the lower and upper bounds in each of their
     * dimensions, as accessed by the lower_bound(dimension_type) const and upper_bound(dimension_type) const methods. 
     * 
     * Boxes are by default ordered by the lexicographic order on lower-left corner.
     * If these are equal, then points are ordered by the lexicographic order 
     * on the upper-right corner.
     *
     * A box is described by the literal "[a0,b0]x[a1,b1]x...".
     * 
     * Boxes are always bounded, but may be empty. A zero-dimensional box is considered to be non-empty.
     *
     * \b Storage: A %Box of dimension \a d is specified by \a 2d real 
     * number giving the lower and upper bounds in each coordinate.
     * The lower bound in the \a i th coordinate is the \a 2i th element, and the 
     * upper bound is the \a 2i+1 st element.
     */
    template<class R>
    class Box 
      : public BoxExpression< Box<R> >
    {
      typedef Box<R> Self;
      BOOST_CLASS_REQUIRE(Self,Ariadne,BasicSetConcept);
      typedef typename traits<R>::interval_type I;
      typedef typename traits<R>::arithmetic_type A;
     private:
      array<R> _data;
     public:
      /*! \brief A tag describing the type of set. */
      typedef basic_set_tag set_category;
      /*! \brief The type used to describe the space the set lies in. */
      typedef EuclideanSpace space_type;
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type used for the corners. */
      typedef R value_type;
      /*! \brief The type of denotable point contained by the box. */
      typedef Point<R> state_type;
      /*! \brief An iterator to the vertices of the box. */
      typedef BoxVerticesIterator<R> vertices_const_iterator;
     public:
      //@{
      //! \name Constructors
      /*! \brief Construct an empty box with dimension \a d. */
      explicit Box(size_type d=0);
      
      /*! \brief Construct a box from a range of interval values. */
      template<class ForwardIterator> explicit Box(ForwardIterator b, ForwardIterator e);
      
      /*! \brief Construct from a one-dimensional C-style array of the form <code>{ l1,u2, l2,u2, ... , ln,un }</code>. */
      template<class RR> explicit Box(const dimension_type& d, const RR* ptr);
      
      /*! \brief Construct from a two-dimensional C-style array of the form <code>{ {l1,u2}, {l2,u2}, ... , {ln,un} }</code>. */
      template<class RR> explicit Box(const dimension_type& d, const RR ary[][2]);
      
      /*! \brief Construct from an array of intervals. */
      explicit Box(const array< Interval<R> >& a);
      
      /*! \brief Construct from a std::vector of intervals. */
      explicit Box(const std::vector< Interval<R> >& v);

      /*! \brief Construct a degenerate box from a single point. */
      explicit Box(const Point<R>& pt);
      
      /*! \brief Construct a box from an interval point. */
      explicit Box(const Point< Interval<R> >& pt);;
      
      /*! \brief Construct from two corners. */
      explicit Box(const Point<R>& pt1, const Point<R>& pt2);
      
      /*! \brief Construct from a string literal. */
      explicit Box(const std::string& s);
      
      /*! \brief Construct from an interval vector. */
      explicit Box(const Vector< Interval<R> >& iv);
      
      /*! \brief Convert from a box expression. */
      template<class E> Box(const BoxExpression<E>& r);
    
    
      /*! \brief Copy constructor. */
      Box(const Box<R>& bx);
    
      /*! \brief Copy assignment operator. */
      Box<R>& operator=(const Box<R>& bx);

      /*! \brief Assign from a box expression. */
      template<class E> Box<R>& operator=(const BoxExpression<E>& r);
    

      /*! \brief Make an empty box. */
      static Box<R> empty_box(dimension_type d);
      /*! \brief Make a unit box. */
      static Box<R> unit_box(dimension_type d);
      /*! \brief Make the positive orthant in \f$\mathbb{R}^d\f$. */
      static Box<R> positive_orthant(dimension_type d);
      /*! \brief Make entire Euclidean space. */
      static Box<R> entire_space(dimension_type d);

      //@}
      
      
      //@{
      //! \name Conversion operators
      /*! \brief Convert to an interval point. */
      operator Point< Interval<R> >() const;
      //@}
      
      
      //{@
      //! \name Comparison operators
      /*! \brief The equality operator */
      bool operator==(const Box<R>& bx) const;
      
      /*! \brief The inequality operator */
      bool operator!=(const Box<R>& bx) const;
      //@}
      

      //@{
      //! \name Data access
      /*! \brief Returns a reference to the array of data. */
      array<R>& data();
     
      /*! \brief Returns a constant reference to the array of data. */
      const array<R>& data() const;
     
      /*! \brief The lower bound of the \a i th coordinate */
      const R& lower_bound(dimension_type i) const;
      
      /*! \brief A reference to the lower bound of the \a i th coordinate */
      R& lower_bound(dimension_type i);
      
      /*! \brief The upper bound of the \a i th coordinate */
      const R& upper_bound(dimension_type i) const;
      
      /*! \brief A reference to the upper bound of the \a i th coordinate */
      R& upper_bound(dimension_type i);
      
      /*! \brief Returns the projection onto the \a i th coordinate (unchecked). */
      Interval<R>& operator[] (dimension_type i);
     
      /*! \brief The projection onto the \a i th coordinate (unchecked). */
      const Interval<R>& operator[] (dimension_type i) const;
      
      /*! \brief The interval of values in the \a i th coordinate. */
      const Interval<R>& interval(dimension_type i) const;
      
      /*! \brief The lower corner. */
      Point<R> lower_corner() const;
      
      /*! \brief The upper corner. */
      Point<R> upper_corner() const;
      
      /*! \brief The set of position vectors of the box. */
      Vector< Interval<R> > position_vectors() const;
      //@}
      
      
      //@{ 
      //! \name Modifying operations
      /*! \brief Makes the box empty. */
      void clear();
      
      /*! \brief Sets the \a i th interval. */
      void set_interval(dimension_type i, Interval<R> x);
      
      /*! \brief Sets the lower bound of the \a i th coordinate to \a r. */
      void set_lower_bound(dimension_type i, const R& l);
      
      /*! \brief Sets the upper bound of the \a i th coordinate to \a u. */
      void set_upper_bound(dimension_type i, const R& u);

      /*! \brief Return a copy of the %Box expanded by \a delta in each direction. */
      Box<R> neighbourhood(const R& delta) const;
      //@}
      
      
      //@{
      //! \name Box geometric operations
      /*! \brief The dimension of the Euclidean space the box lies in. */
      dimension_type dimension() const;
      
      /*! \brief The centre. */
      Point<R> centre() const;
      
      /*! \brief The radius in the sup norm. */
      R radius() const;
      
      /*! \brief An approximation to the volume. */
      R volume() const;
      
      /*! \brief Determines whether the box is empty. */
      tribool empty() const;
      
      /*! \brief Determines whether the box is bounded. */
      tribool bounded() const;
      
      /*! \brief A bounding box for the box; returns the identity. */
      const Box<R>& bounding_box() const;
      
      /*! \brief Determines whether the box contains a point. */
      tribool contains(const Point<R>& pt) const;
      tribool contains(const Point<I>& pt) const;
      
      /*! \brief Determines whether the box is disjoint from another box. */
      tribool disjoint(const Box<R>& bx) const;
      
      /*! \brief Determines whether the box intersects another box. */
      tribool intersects(const Box<R>& bx) const;
      
      /*! \brief Determines whether the box is a subset of another box. */
      tribool subset(const Box<R>& bx) const;
      
      /*! \brief Determines whether the box is a superset of another box. */
      tribool superset(const Box<R>& bx) const;
      
      /*! \brief Compute a quadrant of the Box determined by \a q.
       *  \a q is a binary word such that the ith bit of q is 0 if the lower half
       *  of the box in the ith coordinate is used, and 1 if the upper
       *  half is used.
       */
      Box<R> quadrant(const BinaryWord& q) const;
      
      /*! \brief The number of vertices. */
      size_type number_of_vertices() const;
      /*! \brief The \a i th vertex. */
      Point<R> vertex(size_type i) const;
      /*! \brief A list of all vertices. */
      PointList<R> vertices() const;

      /*! \brief An iterator to the first vertex of the box. */
      vertices_const_iterator vertices_begin() const;
      /*! \brief An iterator to the end vertex of the box. */
      vertices_const_iterator vertices_end() const;

      /*! \brief Split into two smaller boxes along the longest edge. */
      ListSet< Box<R> > split() const;
      /*! \brief Split into boxes along all longest edge. */
      ListSet< Box<R> > subdivide() const;
      //@}
      
#ifdef DOXYGEN
      //@{ 
      //! \name Binary geometric predicates
      /*! \brief Tests if box \a bx contains point \a pt. */
      friend tribool contains(const Box<R>& bx, const Point<R>& pt);

      /*! \brief Tests disjointness of \a bx1 and \a bx2. */
      friend tribool disjoint(const Box<R>& bx1, const Box<R>& bx2);
      /*! \brief Tests if \a bx1 and \a bx2 intersect. */
      friend tribool intersect(const Box<R>& bx1, const Box<R>& bx2);
      /*! \brief Tests if box \a bx1 is a subset of another box \a bx2. */
      friend tribool subset(const Box<R>& bx1, const Box<R>& bx2);
      /*! \brief Tests if box \a bx1 is a superset of another box \a bx2. */
      friend tribool superset(const Box<R>& bx1, const Box<R>& bx2);
      //@}

      //@{ 
      //! \name Binary geometric operations
      /*! \brief The intersection of \a bx1 and \a bx2. */
      friend Box<R> closed_intersection(const Box<R>& bx1, const Box<R>& bx2); 
      /*! \brief The closure of the intersection of the interiors of \a bx1 and \a bx2. */
      friend Box<R> open_intersection(const Box<R>& bx1, const Box<R>& bx2); 

      /*! \brief The smallest box containing \a bx1 and \a bx2. */
      friend Box<R> rectangular_hull(const Box<R>& bx1, const Box<R>& bx2); 
      //@{ 
      //! \name Input/output operations
      /*! \brief Write to an output stream. */
      friend std::ostream& operator<<(std::ostream& os, const Box<R>&);
      /*! \brief Read from an input stream. */
      friend std::istream& operator>>(std::istream& is, const Box<R>&);
      //@}
#endif
    };
    
    
    
    template<class R>
    class BoxVerticesIterator 
      : public boost::iterator_facade<BoxVerticesIterator<R>,
                                      Point<R>,
                                      boost::forward_traversal_tag,
                                      Point<R> const&,
                                      Point<R> const*
                                     >
    {
     public:
      BoxVerticesIterator(const Box<R>& bx, const bool end);
      bool equal(const BoxVerticesIterator<R>& other) const;
      const Point<R>& dereference() const;
      void increment();
     private:
      const Box<R>* _bx; 
      long unsigned int _i; 
      bool _parity; 
      Point<R> _pt;
    };
    
    
    
    
    
    template<class R, class X> tribool contains (const Box<R>& bx, const Point<X>& pt);
    template<class R> tribool disjoint(const Box<R>& bx1, const Box<R>& bx2);
    template<class R> tribool intersect(const Box<R>& bx1, const Box<R>& bx2);
    template<class R> tribool subset(const Box<R>& bx1, const Box<R>& bx2);
    template<class R> tribool superset(const Box<R>& bx1, const Box<R>& bx2);
    template<class R> Box<R> bounding_box(const Box<R>& bx);
    template<class R> ListSet< Box<R> > split(const Box<R>& bx);
    template<class R> ListSet< Box<R> > subdivide(const Box<R>& bx);
    template<class R> Box<R> closed_intersection(const Box<R>& bx1, const Box<R>& bx2);
    template<class R> Box<R> open_intersection(const Box<R>& bx1, const Box<R>& bx2);
    template<class R> Box<R> rectangular_hull(const Box<R>& bx1, const Box<R>& bx2);
      
    template<class R> Box<R> operator+(const Box<R>& bx, const Vector< Interval<R> >& iv);
 
    template<class R> std::istream& operator>>(std::istream& is, Box<R>& bx);
    template<class R> std::ostream& operator<<(std::ostream& os, const Box<R>& bx);

    

} // namespace Ariadne

#include "box.inline.h"

#endif /* ARIADNE_RECTANGLE_H */
