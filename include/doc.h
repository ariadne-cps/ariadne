/***************************************************************************
 *            doc.h
 *
 *  Wed 16 Sept 2004
 *  Copyright  2004  Pieter Collins
 *  Pieter.Collins@cwi.nl
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

#ifndef _DOC_H
#define _DOC_H

/*! \file doc.h
 * \brief Miscellaneous documentation pages
 */
	
/*! \page geometric Geometric Representation
 *
 * \section Introduction
 * 
 * One of the main concerns of Ariadne is the computation of points and sets in Euclidean space.
 * Unfortunately, the set of points in Euclidean space is uncountable (it has continuum cardinality) and so there
 * is no way of modelling all points in Euclidean space using binary words. The set of subsets of Euclidean space
 * is also uncountable, and has an even greater cardinality!
 *
 * To represent an uncountable set, we can use approximations.
 * We take a set of elements which we model exactly, such as @a double reals or @a rational reals, which we call <em>denotable</em> elements.
 * This denotable set may be finite or countably infinite, the former being denoted by words of a fixed length, and the latter by arbitrary length words.
 * An <em>approximation</em> to an element is then a pair consisting of a denotable element and an <em>error</em>.
 * Ideally, we wish to be able to approximate to arbitrary precision;
 * this means that the set of denotable elements must be <em>dense</em> in the set of all elements.
 * We can then represent arbitrary elements be a <em>convergent sequence</em> of approximations.
 * In terms of a concrete representation on a computer, we can think of a such a sequence as a neverending <em>stream</em> of data.
 *
 * The material in this section is heavily influenced by the book @ref weihrauch2000 "Computational Analysis" by Klaus Weihrauch.
 *
 * \section real Representation of real numbers
 * 
 * Real numbers are traditionally described by floating-point types @a float and @a double.
 * The set of denotable elements is <em>finite</em>.
 * Arithmetic operations can only be performed approximately to maximal precision determined by the data type.
 * However, these types have the advantage of requiring a known amount of memory, which means they can be statically
 * allocated, and having hardware-supported arithmetical approximations. This makes them especially suitable in
 * situations where execution speed and memory usage are more important than computational accuracy.
 *
 * Alternatively, real numbers may be described by rational types such as rational and dyadic, or even by more
 * complicated types, such as @a algebraic. We have a <em>countable</em> set of denotable elements which are dense
 * in the set of all real numbers. Real-valued functions can be approximated to arbitrary precision, which may be 
 * by the user. Arithmetic operations can usually be performed exactly. However, these types require arbitrarily
 * large amounts of memory which is typically dynamically allocated, and hardware support for arithmetic does not
 * exist, so are typically less efficient than the floating-point types.
 *
 * COMMENT FOR ALBERTO
 * 
 * Sometime when designing class hierarchies, we are going to have to use templates. 
 * In particular, if we want to be able to use either (say) double or rational
 * as the underlying representation of real numbers, almost everything will be
 * templated by this type. The question is where do these types occur in the
 * class hierarchy? I would suggest that since the real_type is an implementation
 * issue and not a abstraction issue, that we define pure abstract interfaces without any
 * templates, but denotable instances use templates.
 *
 * This also raises the issue of whether to have an abstract class representing real numbers. 
 * The advantage of this is that is will allow our abstract state types to return their coordinates.
 * The disadvantage is that for every class representing a (denotable) real number, we will need to 
 * introduce it into the class hierarchy. I think this is a bad idea -- lots of wasted time with no reward.
 *
 * Better, is to have a system-wide default  denotable_real_type  (e.g. rational) which is returned by
 * abstract functions such as volume(). Of course, even then the volume is only an approximation in some cases.
 * Even better then is to have a volume() function return a RealApproximation<rational> or interval<rational> or 
 * interval<denotable_real_type> object which approximates the volume().
 * 
 *
 * \section state Representation of points in Euclidean space
 *
 * Points in Euclidean space can be represented by the templated class @a DenotablePoint<T>, where @a T is a numeric type
 * such as double or rational. From this, we can define a PointApproximation class, which gives an approximation of a 
 * state, plus an error bound in terms of the sup norm or Euclidean norm.
 * 
 * @code
 * // Abstract point class representing possibly 
 * // uncomputable points in Euclidean space!
 * class Point 
 * {
 *   // This class is so abstract, we can't even define it's coordinates,
 *   // since they depend on the approximation!!
 *   virtual ~Point() = 0;
 *   virtual Point* clone() const = 0;
 * };
 * 
 * // Concrete point class for efficiency.
 * template<class R>
 * class DenotablePoint
 * {
 *   typedef R real_type;
 *
 *   DenotablePoint(const DenotablePoint& p);
 *   DenotablePoint& operator=(const DenotablePoint& p);
 *
 *   real_type operator[] (uint n) const;
 * ...
 * };
 * 
 * // A point in Euclidean space represented by real numbers of type R.
 * template<class R>
 * class EuclideanPoint : public Point, private DenotablePoint<R> { 
 *   typedef R real_type;
 *   
 *   EuclideanPoint(const DenotablePoint<R>& p) : DenotablePoint<R>(p) { }
 *   virtual EuclideanPoint* clone() const { return new EuclideanPoint(*this); }
 *   real_type operator[] (uint n) const { return DenotablePoint<R>(n); };'
 *  ...
 * };
 *
 * template<class P, class R=P::real_type>
 * class PointApproximation { 
 *   typedef R point_type;
 *   typedef R real_type;
 *
 *   point_type& approximation() const;
 *   real_type error() const;
 * };
 * @endcode
 *
 * \section set Representation of (closed) sets.
 *
 * Basic sets.
 *
 * @code
 * // Abstract base for Set class.
 * // There's not much we can do with such a general Set, 
 * // since we don't have type information, and even 
 * // membership of a general set may be undefined.
 * class Set 
 * {
 *   virtual ~Set() = 0;
 * };
 * 
 * 
 * // An abstract base class for BasicSets, to allow mixing of 
 * // of a BasicRectangle, or as a template parameter for DenotableSet.
 * class BasicSet
 *   : public Set
 * {
 *   // A copy of the BasicSet.
 *   virtual BasicSet* clone() const = 0; 
 *   
 *   // No assignment operator, since we only work through pointers.
 *   // BasicSet& operator=(const BasicSet&);
 *
 *   // Reference or pointer? Since we need to make this a general interface, and 
 *   // hence independent of representation, we can't assume that we can return a 
 *   // reference to an existing object, and hence the result must be allocated 
 *   // on the heap and requires deletion.
 *   virtual Point* centre() const = 0; 
 *   virtual rational min_radius() const = 0; //alternatively, return an object of default_real_type
 *   virtual rational max_radius() const = 0; 
 *   virtual rational volume() const = 0; 
 *
 *   // Set up dynamic dispatching to show errors if is_subset and is_disjoint are not defined for actual types.
 *   virtual bool is_subset(const BasicSet& b) const throw undefined_operation { throw undefined_operation("is_subset"); }
 *   virtual bool is_disjoint(const BasicSet& b) const throw undefined_operation { throw undefined_operation("is_disjoint"); } = 0;
 *   
 *   // Note that there are no union or intersection operations defined!
 *   // This is because the intersection of two BasicSets need not be a BasicSet.
 *   // There are hence only very few operations supported by all BasicSets which
 *   // we can use in a library, unless we use unsafe "fat interfaces".
 * };
 * 
 * // A DenotableSet<BS> implements a set as a union of basic sets type BS.
 * // Maybe we want a further abstraction level here 
 * template<class BS>
 * class DenotableSet
 *   : public Set
 * {
 *   typedef BS::real_type real_type;
 *   typedef BS::point_type point_type;
 *   typedef BS basic_set_type;
 *   typedef R real_type;
 *
 *   // Must return a new pointer here as there is no 
 *   // guarentee that the basic sets are stored 
 *   // in a form in which they can be referenced.
 *   basic_set_type* operator[] (uint n) const = 0;
 *
 *   // A DenotableSet<BS> can determine whether is contains a denotable point.
 *   bool contains(const point_type& p) const = 0;
 * };
 * 
 * // A rectangle in Euclidean space, which can be used as an implementation
 * // of a BasicRectangle, or as a template parameter for DenotableSet.
 * templace<class R>
 * class Rectangle
 * {
 *   typedef R real_type;
 *   typedef DenotablePoint<R> point_type;
 *
 *   // Simple operations
 *   real_type lower(uint n) const;
 *   real_type upper(uint n) const;
 *   Interval<real_type> interval(uint n) const;
 *   bool contains(const point_type& p) const;
 *   
 *   // I suggest we don't define binary operations as members,
 *   // we also need binary operations between different types of
 *   // concrete sets.
 *   // is_subset(const Rectangle& r) const; // NOT DEFINED
 *   // is_disjoint(const Rectangle& r) const; // NOT DEFINED
 *   // intersect(const Rectangle& r) const; // NOT DEFINED
 *
 *   // friend declarations for efficiency
 *   friend Rectangle intersection(const Rectangle& r1, const Rectangle& r2); 
 *   friend bool is_subset(const Rectangle& r, const Ellipsoid& e);
 *   friend bool are_disjoint(const Rectangle& r1, const Rectangle& r2); 
 * };
 * 
 * // An ellipsoid in Euclidean space
 * templace<class R>
 * class Ellispsoid {
 *  ...
 * };
 *
 * // Four is_subset and is_disjoint operations for just two classes!
 * bool is_subset(const Rectangle& r, const Rectangle& r);
 * bool is_subset(const Rectangle& r, const Ellipsoid& e);
 * bool is_subset(const Ellipsoid& e, const Rectangle& r);
 * bool is_subset(const Ellipsoid& e, const Ellipsoid& e);
 * 
 * bool is_disjoint(const Rectangle& r, const Rectangle& r);
 * bool is_disjoint(const Rectangle& r, const Ellipsoid& e);
 * bool is_disjoint(const Ellipsoid& e, const Rectangle& r);
 * bool is_disjoint(const Ellipsoid& e, const Ellipsoid& e);
 * 
 * // The intersection of two rectangles is a rectangle, but intersections
 * // involving ellipsoids cannot be defined explicitly at the level
 * // of concrete basic sets.
 * Rectangle intersection(const Rectangle& r1, const Rectangle& r2); 
 * 
 *
 * templace<class R>
 * class BasicRectangle
 *   : public BasicSet, private Rectangle<R>
 * {
 *   typedef R real_type;
 *   typedef DenotablePoint<R> point_type;
 *  
 *   // Construct from a concrete Rectangle<R>
 *   BasicRectangle(const Rectangle<R>& r) : Rectangle(r) { };
 *
 *   // Copy constructor and assignment
 *   BasicRectangle(const BasicRectangle<R>& r) : Rectangle(r) { };
 *   BasicRectangle<R> operator=(const const BasicRectangle<R>& r);
 *
 *   // Implement virtual functions
 *   virtual BasicRectangle<R>* clone() const { return new BasicRectangle(*this); } 
 *   virtual Point* centre() const { return Rectangle::centre(); }
 *   virtual rational min_radius() const { return Rectangle::min_radius(); }
 *   virtual rational max_radius() const { return Rectangle::max_radius(); } 
 *   virtual rational volume() const { return Rectangle::volume(); }; 
 *
 *   // Implement double-dispatched operations.
 *   virtual bool is_subset(const BasicSet& b) const throw undefined_operation { throw undefined_operation("is_subset"); } 
 *   virtual bool is_subset(const BasicRectangle& b) const { return is_subset(*this, b); } // Dispatch to  is_subset(Rectangle&, Rectangle&);
 *   virtual bool is_subset(const BasicEllipsoid& b) const { return is_subset(*this, b); } // Dispatch to  is_subset(Rectangle&, Ellipsoid&);

 *   // Implement Rectangle-specific functions
 *   Interval<real_type> interval(uint n) const { return Rectangle<R>::upper(n); };
 *   real_type lower(uint n) const { return Rectangle<R>::lower(n); }
 *   real_type upper(uint n) const { return Rectangle<R>::upper(n); }
 * };
 * 
 * template<class S>
 * class SetApproximation { 
 *   typedef S::real_type real_type;
 *   typedef S::point_type point_type;
 *   typedef S set_type;
 *   set_type& approximation() const;
 *   real_type error() const;
 * };
 * @endcode
 */

#endif /* _DOC_H */
