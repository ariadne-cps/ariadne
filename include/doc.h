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

#ifndef _ARIADNE_DOC_H
#define _ARIADNE_DOC_H

/*! \file doc.h
 * \brief Miscellaneous documentation pages
 */

/*! \page real Real Number Types
 *
 * \section Introduction
 *
 * Ariadne supports three different real number types, \a double, \a dyadic and \a rational.
 * The \a double type is a \emph finite-precision type, and the \a dyadic and \a rational types are \emph arbitrary-precision types.
 *
 * \section Finite-Precision Real Numbers
 *
 * Real numbers are traditionally described by floating-point types \a float and \a double.
 * The set of denotable elements is <em>finite</em>.
 * Arithmetic operations can only be performed approximately to maximal precision determined by the data type.
 * However, these types have the advantage of requiring a known amount of memory, which means they can be statically
 * allocated, and having hardware-supported arithmetical approximations. This makes them especially suitable in
 * situations where execution speed and memory usage are more important than computational accuracy.
 *
 * A finite-precision real number type \c FPReal must satisfy the following requirements
 *
 * \code
 * concept FPReal {
 *  public:
 *   FPReal(int);
 *   FPReal(const FPReal &);
 *   FPReal& operator=(const FPReal&);
 * };
 *
 * // Arithmetical operators.
 * FPReal operator-(FPReal);
 * FPReal operator+(FPReal,FPReal);
 * FPReal operator-(FPReal,FPReal);
 * FPReal operator*(FPReal,FPReal);
 * FPReal operator/(FPReal,FPReal);
 * \endcode
 *
 * We also define the following functions.
 *
 * \code
 * // Arithmetical functions.
 * FPReal neg(FPReal);
 * FPReal add(FPReal,FPReal);
 * FPReal sub(FPReal,FPReal);
 * FPReal mul(FPReal,FPReal);
 * FPReal recip(FPReal);
 * FPReal div(FPReal,FPReal);
 *
 * // Conversions.
 * FPReal floor(FPReal);
 * FPReal ceil(FPReal);
 *
 * // Comparisons.
 * FPReal min(FPReal, FPReal)
 * FPReal max(FPReal, FPReal)
 *
 * // Absolute value.
 * FPReal abs(FPReal)
 *
 * // Algebraic functions.
 * FPReal pow(FPReal,int)
 * FPReal square(FPReal)
 * FPReal sqrt(FPReal)
 *
 * // Transcendental functions.
 * FPReal exp(FPReal);
 * FPReal log(FPReal);
 *
 * FPReal sin(FPReal);
 * FPReal cos(FPReal);
 * FPReal tan(FPReal);
 *
 * FPReal asin(FPReal);
 * FPReal acos(FPReal);
 * FPReal atan(FPReal);
 *
 * FPReal sinh(FPReal);
 * FPReal cosh(FPReal);
 * FPReal tanh(FPReal);
 *
 * FPReal asinh(FPReal);
 * FPReal acosh(FPReal);
 * FPReal atanh(FPReal);
 * \endcode
 *
 * The arithmetical operators should follow IEEE standards.
 * The arithmetical functions are synonyms for the arithmetical operators.
 * The precision of the algebraic and transcendental functions is not guarenteed,
 * except that the functions must give values in the correct range.
 *
 * In order to support interval arithmetic, we also define the following functions.
 *
 * \code
 * FPReal f_approx(FPReal); // Synonym for FPReal f(FPReal)
 * FPReal f_down(FPReal);
 * FPReal f_up(FPReal)
 * \endcode
 * An alternative syntax is
 * \code
 * enum ApproximationKind { EXACT, APPROX, LOWER, UPPER };
 *
 * FPReal f(FPReal,ApproximationKind);
 * \endcode
 * where \p ApproximationKind describes whether the function return value is exact, and approximation, a lower bound or an upper bound.
 *
 * Further, we define arithmetical operations on Interval<FPReal>.
 *
 * The upper and lower bounds satisfy the mathematical postcondition \f$ f(x) \in [\textrm{f\_down}(x),\textrm{f\_up}(x)] \f$
 * and the implementation postcondition <tt>f_down(x) <= f_approx(x) <= f_up(x)</tt>.
 *
 * \note
 * No other conditions are required of \c f, \c f_down, and \c f_up.
 * Hence a valid (but useless) implementation of \f$ \sin(x) \f$ is
 * \code
 * FPReal sin_approx(FPReal x) { return 0.0; }
 * FPReal sin_down(FPReal x) { return -1.0; }
 * FPReal sin_up(FPReal x) { return 1.0; }
 * \endcode
 *
 * \internal
 * The names are chosen for compatibility with \c \<cmath\> and Boost.
 * The only name change is that \c multiplicative_inverse from Boost has become \c recip.
 *
 * Should we only support explicit approximation functions, like \c exp_approx for finite-precision types?
 * While it is in some ways appealing to make the approximation explicit,
 * I think it's probably better to allow <tt> double exp(double)</tt> in Ariadne, for two reasons.
 * The first reason is that this is common useage, and people might well include <cmath> anyway.
 * The second is operator overloading; we should insist that \c operator+ and add give the same answers,
 * and by disallowing add, we also need to disallow all operator overloading on double-precision types.
 *
 * With \c floor and \c ceil, it would be useful to have these return an \c int.
 * Unfortunately, this is not the behaviour in \c \<cmath\>, and it's probably not a good idea to change this.
 * A solution for \c double is to use the conversion to \c int, but this doesn't work for \c dyadic and \c rational types.
 * Solutions are
 * \li provide functions <tt>int ifloor(Real)</tt> and <tt>int iceil(Real)</tt>
 * \li wrap the GnuMP classes to provide conversion operator to int (rounding to zero)
 * \li provide a function <tt>int convert_to_int(Real)</tt> (ugly)
 *
 * Should we use the \c f_down, \c f_up syntax, or the \c f(Real,ApproximationKind) syntax?
 * Boost uses \c f_down and \c f_up.
 * I currently have a slight preference for \c f_down and \c f_up,
 * but this might change when we fix the names for approximations of sets and functions.
 *
 *
 * \section arbitraryprecision Arbitrary-precision types.
 *
 * Alternatively, real numbers may be described by rational types such as \c rational and \c dyadic, or even by more
 * complicated types, such as \c algebraic. We have a <em>countable</em> set of denotable elements which are dense
 * in the set of all real numbers. Real-valued functions can be approximated to arbitrary precision, which may be
 * by the user. Arithmetic operations can usually be performed exactly. However, these types require arbitrarily
 * large amounts of memory which is typically dynamically allocated, and hardware support for arithmetic does not
 * exist, so are typically less efficient than the floating-point types.
 *
 * \code
 * concept APReal {
 *  public:
 *   APReal(int);
 *   APReal(double);
 *   APReal(const APReal &);
 *   APReal& operator=(const APReal&);
 * };
 *
 * // Arithmetical operators.
 * APReal operator-(APReal);
 * APReal operator+(APReal,APReal);
 * APReal operator-(APReal,APReal);
 * APReal operator*(APReal,APReal);
 * rational operator/(rational,rational);
 * dyadic operator/(dyadic,int i); // Behaviour undefined unless i is a power of 2.
 *
 * // Arithmetical functions.
 * APReal neg(APReal);
 * APReal add(APReal,APReal);
 * APReal sub(APReal,APReal);
 * APReal mul(APReal,APReal);
 *
 * rational recip(rational);
 * rational div(rational,rational);
 * dyadic div(dyadic,int); // Behaviour undefined unless i is a power of 2.
 *
 * // Conversions.
 * APReal floor(APReal);
 * APReal ceil(APReal);
 *
 * // Comparisons.
 * APReal min(APReal, APReal)
 * APReal max(APReal, APReal)
 *
 * // Absolute value.
 * APReal abs(APReal)
 *
 * // Algebraic functions.
 * APReal pow(APReal,int)
 * APReal square(APReal)
 *
 * // Approximate algebraic functions.
 * APReal div_approx(APReal x, APReal y, APReal d, APReal& e); // d gives the error in y.
 * APReal recip_approx(APReal x, APReal d, APReal& e);
 * APReal sqrt_approx(APReal x, APReal d, APReal& e);
 *
 * // Approximate transcendental functions.
 * APReal exp_approx(APReal x, APReal d, APReal& e);
 * APReal log_approx(APReal x, APReal d, APReal& e);
 *
 * APReal sin_approx(APReal x, APReal d, APReal& e);
 *  ...
 * \endcode
 *
 * Unlike a finite-precision type, arbitrary-precision types are intended for use when precise error specifications are required.
 * Hence, whenever a function returns an approximation, the suffix \c _approx is \em always added to the function name.
 * This ensures that the user is always aware of the use of approximations.
 *
 * The \c dyadic type does not in general support exact division.
 * The exception is that division by a power of 2 yields an exact answer.
 * Generic code intended for use with all arbitrary precision types should \em not use division, except by a power of two.
 *
 * \b Warning: The compiler may not catch errors resulting from inexact division involving \c dyadic numbers!!
 *
 * The semantics of \c f_approx is as follows:
 *
 * The parameter \a d is the maximum error in the argument \a x.
 * If \c d==0, then \a x is exact.
 * The parameter \a e is an input/output parameter giving the error in the output.
 * The function \c f_approx tries to find an output with a maximum error of \a e, if possible.
 * After the function call, \a e is overwritten with the best guarenteed error bound.
 *
 * If \c e==0, then \c f_approx is guarenteed to converge to the exact answer as \a d tends to 0.
 * If \c d==0 and \c e==0, then either an approximation is computed to some arbitrary error, or an exception is thrown.
 * This choice is implementation-dependent.
 *
 * \section interval Interval Arithmetic
 *
 * Interval arithmetic is commonly used in finite-precision computations to obtain guarenteed bounds for the result.
 * However, it can also be used in arbtitrary-precision computations to compute the image of basic sets.
 * The semantics is as follows:
 *
 * \subsection finite_precision_interval Finite-precision interval arithmetic
 *
 * Finite precision functions are of the form
 * \code
 * Interval<FPReal> f_approx(Interval<FPReal>);
 * \endcode
 * and satisfy the mathematical postcondition \f$ \forall x\in I,\ f(x)\in\textrm{f\_approx}(I) \f$,
 * and the implementation postcondition <tt>I.contains(x)</tt> implies <tt>f_approx(I).contains(f_approx(x))</tt>.
 *
 * \subsection arbitrary_precision_interval Arbitrary-precision interval arithmetic
 *
 * Exact arbitrary-precision interval functions are given by
 * \code
 * Interval<APReal> f(Interval<APReal>);
 * \endcode
 * and approximate arbitrary-precision interval functions by
 * \code
 * Interval<APReal> f_approx(Interval<APReal>);
 * \endcode
 * where \c f(I) is the exact image of the interval \c I, and \c f_approx satisfies the mathematical postcondition
 * \f$ \textrm{f\_approx}(I) \supset f(I)\f$, and the convergence criterion that the length of \c f_approx(I) approaches 0
 * as the length of \c I approaches 0.
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
 * We take a set of elements which we model exactly, such as \a double reals or \a rational reals, which we call <em>denotable</em> elements.
 * This denotable set may be finite or countably infinite, the former being denoted by words of a fixed length, and the latter by arbitrary length words.
 * An <em>approximation</em> to an element is then a pair consisting of a denotable element and an <em>error</em>.
 * Ideally, we wish to be able to approximate to arbitrary precision;
 * this means that the set of denotable elements must be <em>dense</em> in the set of all elements.
 * We can then represent arbitrary elements be a <em>convergent sequence</em> of approximations.
 * In terms of a concrete representation on a computer, we can think of a such a sequence as a neverending <em>stream</em> of data.
 *
 * The material in this section is heavily influenced by the book \ref weihrauch2000 "Computational Analysis" by Klaus Weihrauch.
 *
 * \section state Representation of points in Euclidean space
 *
 * Points in Euclidean space can be represented by the templated class \a DenotablePoint<T>, where \a T is a numeric type
 * such as double or rational. From this, we can define a PointApproximation class, which gives an approximation of a
 * state, plus an error bound in terms of the sup norm or Euclidean norm.
 *
 * \code
 * template<class R>
 * class State
 * {
 *   typedef R Real;
 *
 *   State(const State &);
 *   State& operator=(const State &);
 *
 *   Real operator[] (size_type) const;
 * };
 * \endcode
 *
 * \section basicset Basic sets.
 *
 * A \c BasicSet provides a building block for representing more complicated sets.
 * Mathematically, a \c BasicSet type represents elements of a countable base of a topological space.
 * More precisely, a \c BasicSet represents the \em closure of a basic set for the topology.
 *
 * \code
 * // The basic set concept.
 * concept BasicSet
 * {
 *   typename Real;
 *   typename State;
 *
 *   BasicSet(const BasicSet &);
 *   BasicSet & operator=(const BasicSet &);
 *
 *   bool operator==(const BasicSet &) const;
 *   bool operator!=(const BasicSet &) const;
 *
 *   size_type dimension() const;
 *   bool empty() const;
 *   bool contains(const State &) const;
 *
 *   //Optional
 *   State centre() const;
 *   Real min_radius() const;
 *   Real max_radius() const;
 *   Real volume() const;
 *
 *   //Optional
 *   Rectangle<Real> bounding_box() const;
 *   Sphere<Real> bounding_sphere() const;
 * };

 *  // Optional, since the intersection of two basic sets need not be a basic set.
 *  BasicSet regular_intersection(const BasicSet &, const BasicSet &);
 *
 *  bool interiors_intersect(const BasicSet &, const BasicSet &);
 *  bool disjoint(const BasicSet &, const BasicSet &);
 *  bool inner_subset(const BasicSet &, const BasicSet &);
 *
 *  // Optional, but highly recommended.
 *  bool subset(const BasicSet &, const BasicSet &);
 * \endcode
 *
 * Classes fulfilling the \c BasicSet concept are \c Rectangle (or \c Cuboid), \c Simplex, \c Parallelopiped, \c Polytope and \c Ellipsoid.
 * Actually, these are templates, parameterised by the real number type \c Real.
 *
 * \subsection rectangle Rectangles
 * A \c Rectangle describes a cuboid in Euclidean space.
 * \code
  * templace<class R>
 * class Rectangle
 * {
 *   typedef R Real;
 *   typedef State<R> State;
 *
 *   Rectangle(const State &, const State &);
 *
 *   Rectangle(const Rectangle &);
 *   Rectangle & operator=(const Rectangle &);
 *
 *   bool operator==(const Rectangle &) const;
 *   bool operator!=(const Rectangle &) const;
 *
 *   size_type dimension() const;
 *   bool empty() const;
 *   bool contains(const State &) const;
 *
 *   // Simple operations
 *   Real lower(size_type) const;
 *   Real upper(size_type) const;
 *   Interval<Real> interval(size_type) const;
 *
 *   State lower_corner() const;
 *   State upper_corner() const;
 *   Interval<Real> operator[] (size_type) const;
 * };
 * 
 * Rectangle regular_intersection(const Rectangle&, const Rectangle&);
 *
 * bool interiors_intersect(const Rectangle&, const Rectangle&);
 * bool disjoint(const Rectangle&, const Rectangle&);
 * bool subset(const Rectangle&, const Rectangle&);
 * bool inner_subset(const Rectangle&, const Rectangle&);
 * \endcode
 *
 * \subsection ellipsoids Ellipsoids
 *
 * \code
 * // An ellipsoid in Euclidean space
 * templace<class R>
 * class Ellispsoid {
 *  ...
 * };
 *
 * \\ No regular_intersection function.
 *
 * bool interiors_intersect(const Ellipsoid &, const Ellipsoid &);
 * bool disjoint(const Ellipsoid &, const Ellipsoid &);
 * bool inner_subset(const Ellipsoid &, const Ellipsoid &);
 * bool subset(const Ellipsoid &, const Ellipsoid &);
 * \endcode
 *
 * Additionally, we have mixed comparison operators.
 * \code
 * bool interiors_intersect(const Rectangle &, const Ellipsoid &);
 * bool interiors_intersect(const Ellipsoid &, const Rectangle &);
 * bool disjoint(const Rectangle &, const Ellipsoid &);
 * bool disjoint(const Ellipsoid &, const Rectangle &);
 * bool inner_subset(const Rectangle &, const Ellipsoid &);
 * bool inner_subset(const Ellipsoid &, const Rectangle &);
 * bool subset(const Rectangle &, const Ellipsoid &);
 * bool subset(const Ellipsoid &, const Rectangle &);
 * \endcode
 *
 *
 * \section denotable_set Denotable Sets
 *
 * A DenotableSet implements a set as a union of basic sets type \c DenotableSet::BasicSet.
 * \code
 * concept DenotableSet
 * {
 *   typename Real;
 *   typename State;
 *   typename BasicSet;
 *
 *   typename const_iterator; // Must satisfy the requirements of a ForwardIterator.
 *
 *   // No default constructor required.
 *
 *   DenotableSet(const DenotableSet &);
 *   DenotableSet & operator=(const DenotableSet &);
 *
 *   // No equality operator required.
 *
 *   // Set-theoretic operations
 *   size_type dimension() const;
 *   bool empty() const;
 *   bool contains(const State &) const;
 *
 *   void adjoin(const BasicSet &);
 *   void adjoin(const DenotableSet &);
 *
 *   // List operations
 *   const_iterator begin() const;
 *   const_iterator end() const;
 *
 *   size_type size() const; // Only required if the iterator is a RandomAccessIterator.
 *   BasicSet operator[] (size_type) const; // Only required if the iterator is a RandomAccessIterator.
 *
 *   void push_back(const BasicSet &); // Only used if the DenotableSet is an ordered list.
 *   void insert(const BasicSet &); // Only used if the DenotableSet is an unordered or sorted list.
 *
 *   // Miscellaneous operations
 *   Rectangle<Real> bounding_box() const; // Optional.
 * };
 *
 * DenotableSet join(const DenotableSet &, const DenotableSet &);
 *
 * bool inner_subset(const BasicSet &, const DenotableSet &);
 * bool subset(const BasicSet &, const DenotableSet &); // Optional, but highly recommended.
 *
 * bool interiors_intersect(const DenotableSet &, const DenotableSet &);
 * bool disjoint(const DenotableSet &, const DenotableSet &);
 * bool inner_subset(const DenotableSet &, const DenotableSet &);
 * bool subset(const BasicSet &, const DenotableSet &); // Optional, but highly recommended.
 *
 * DenotableSet regular_intersection(const BasicS &, const DenotableSet &); // Optional.
 * \endcode
 *
 * We can define approximations as follows.
 * \code
 * template<class P, class R=P::real_type>
 * class StateApproximation {
 *   typedef R point_type;
 *   typedef R real_type;
 *
 *   point_type& approximation() const;
 *   real_type error() const;
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
 * \endcode
 */

#endif /* _ARIADNE_DOC_H */
