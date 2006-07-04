/***************************************************************************
 *            documentation.h
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

#ifndef _ARIADNE_DOCUMENTATION_H
#define _ARIADNE_DOCUMENTATION_H

/*! \file documentation.h
 * \brief Miscellaneous documentation pages
 */

/*! \mainpage
 *
 * \section Introduction
 * 
 * Ariadne is a C++ package for set-based analysis of dynamical and control systems, including reachability analysis and verification.
 * 
 * \section Requirements
 * 
 * To compile Ariadne, you will need a C++ compiler (we recommend g++ version 4.0.2 or higher, which can be downloaded from <a href="http://gcc.gnu.org/">http://gcc.gnu.org</a>).
 *
 * You will also need the following libraries:
 *   - The GNU Multiple-Precision Library (version 4.1.2 or higher)  <a href="http://www.swox.com/gmp/">http://www.swox.com/gmp/</a>).
 *   - MPFR Library (version 2.0.1 or higher) <a href="http://www.mpfr.org/">http://www.mpfr.org/</a>.
 *   - The Boost C++ Libraries (version 1.33.1) <a href="http://www.boost.org/">http://www.boost.org/</a>.
 *   - The Parma Polyhedra Library (version 0.8) <a href="http://www.cs.unipr.it/ppl/">http://www.cs.unipr.it/ppl/</a>.
 *
 * For the Python interface, you will also need:
 *   - Python (version 2.4) <a href="http://www.python.org/">http://www.python.org/</a>.
 *
 * To make this documentation, you have used:
 *   - Doxygen (version 1.4.6 or higher is recommended) 
 *        <a href="http://www.stack.nl/~dimitri/doxygen/">http://www.stack.nl/~dimitri/doxygen/</a>
 *
 * \section Installation
 *
 * From the ariadne/trunk/ directory of the main tree, type\n
 * <tt>> make install</tt>\n
 * to make and install the library. 
 *
 * The default installation directory for the library is $HOME/lib/, and for the Python interface is $HOME/python/.
 * These defaults can be changed by editing the PREFIX variable in the file config.mk, or, for more control, the 
 * variables LIBPREFIX and PYTHONPREFIX.
 */



 
 
/*! \page real Real Number Types
 *
 * \section Introduction
 *
 * Ariadne supports three different real number types, \a double, \a dyadic and \a rational.
 * The \a double type is a \em finite-precision type, and the \a dyadic and \a rational types are \em arbitrary-precision types.
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
 * Should we only support explicit approximation functions, like \c exp_approx for finite-precision types?
 * While it is in some ways appealing to make the approximation explicit,
 * I think it's probably better to allow <tt> double exp(double)</tt> in Ariadne, for two reasons.
 * The first reason is that this is common useage, and people might well include \c <cmath> anyway.
 * The second is operator overloading; we should insist that \c operator+ and add give the same answers,
 * and by disallowing add, we also need to disallow all operator overloading on double-precision types.
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
 * dyadic div(dyadic,int i); // Behaviour undefined unless i is a power of 2.
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
 * APReal div_approx(APReal x, APReal y, APReal e); 
 * APReal recip_approx(APReal x, APReal e);
 * APReal sqrt_approx(APReal x, APReal e);
 *
 * // Approximate transcendental functions.
 * APReal exp_approx(APReal x, APReal e);
 * APReal log_approx(APReal x, APReal e);
 *
 * APReal sin_approx(APReal x, APReal e);
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
 * The parameter \a e is an input parameter giving the desired error in the output.
 *
 * \internal Alternatively, since we always need an error bound, we don't need to add the suffix \c _approx 
 * to function names for exact arithmetic.
 *
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
 * Interval<FPReal> f(Interval<FPReal>);
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
 * Interval<APReal> f(Interval<APReal>); // (Possible) syntactic sugar for f_approx.
 * \endcode
 * where \c f(I) is the exact image of the interval \c I, and \c f_approx satisfies the mathematical postcondition
 * \f$ \textrm{f\_approx}(I) \supset f(I)\f$, and the convergence criterion that the length of \c f_approx(I) approaches 0
 * as the length of \c I approaches 0.
 *
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
 * The material in this section is heavily influenced by the book "Computational Analysis" by Klaus Weihrauch.
 *
 * \section state Representation of points in Euclidean space
 *
 * Points in Euclidean space can be represented by the templated class \a DenotablePoint<T>, where \a T is a numeric type
 * such as double or rational. From this, we can define a PointApproximation class, which gives an approximation of a
 * state, plus an error bound in terms of the sup norm or Euclidean norm.
 *
 * \code
 * concept State
 * {
 *   typename real_type;
 *
 *   State(const State &);
 *   State& operator=(const State &);
 *
 *   real_type operator[] (size_type) const;
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
 *   typename real_type; // The type of denotable real number used for the representation.
 *   typename state_type; // The type of denotable point the set contains.
 *
 *   BasicSet(const BasicSet &); // Copy constructor.
 *   BasicSet & operator=(const BasicSet &); // Assignment operator.
 *
 *   bool operator==(const BasicSet &) const; // Equality operator.
 *   bool operator!=(const BasicSet &) const; // Inequality operator.
 *
 *   dimension_type dimension() const; // The dimension of the set.
 *   bool empty() const; // Returns true if the set is empty.
 *   bool empty_interior() const; // Returns true if the interior of the set is empty.
 *   bool contains(const State &) const; // Returns true if the set contains a given state.
 *   bool interior_contains(const State &) const; // Returns true if the interior of the set contains a given state.
 *
 *   state_type centre() const; // A point in the set (typically, the "centre" point, if this makes sense).
 *   State radius() const; // The maximum distance from the centre to another point in the set in an appropriate metric. 
 *   Real volume() const; // The volume of the set. (Optional).
 *
 *   Rectangle<real_type> bounding_box() const; // A rectangle containing the set.
 *   Sphere<real_type> bounding_sphere() const; // A sphere containing the set. (Optional)
 * };
 *
 *  // Optional, since the intersection of two basic sets need not be a basic set.
 *  BasicSet regular_intersection(const BasicSet &, const BasicSet &); // The closure of the intersection of the interiors of the two sets.
 *
 *  bool interiors_intersect(const BasicSet &, const BasicSet &); // Returns true if the interiors of the two sets intersect.
 *  bool disjoint(const BasicSet &, const BasicSet &); // Returns true if the two sets are disjoint.
 *  bool inner_subset(const BasicSet &, const BasicSet &); // Returns true if the first argument is a subset of the interior of the second.
 *  bool subset(const BasicSet &, const BasicSet &); // Returns true if the first argument is a subset of the second.
 * \endcode
 *
 * Classes fulfilling the \c BasicSet concept are \c Rectangle (or \c Cuboid), \c Simplex, \c Parallelotope, \c Zonotope, \c Polytope and \c Ellipsoid.
 * Actually, these are templates, parameterised by the real number type \c real_type.
 *
 * \subsection rectangle Rectangles
 * A \c Rectangle describes a cuboid in Euclidean space.
 * \code
  * templace<class R>
 * class Rectangle
 * {
 *   typedef R real_type;
 *   typedef Point<R> state_type;
 *
 *   Rectangle(const state_type &, const state_type &);
 *   Rectangle(const IntervalVector<real_type> &);
 *
 *   Rectangle(const Rectangle &);
 *   Rectangle & operator=(const Rectangle &);
 *
 *   bool operator==(const Rectangle &) const;
 *   bool operator!=(const Rectangle &) const;
 *
 *   size_type dimension() const;
 *   bool empty() const;
 *   bool empty_interior() const;
 *   bool contains(const State &) const;
 *   bool interior_contains(const State &) const;
 *
 *   state_type centre() const;
 *   real_type radius() const;
 *   Rectangle<real_type> bounding_box() const;
 *
 *   // Simple operations
 *   real_type lower_bound(dimension_type) const;
 *   real_type upper_bound(dimension_type) const;
 *   Interval<real_type> operator[] (dimension_type) const;
 *
 *   state_type lower_corner() const;
 *   state_type upper_corner() const;
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
 * A DenotableSet implements a set as a union of basic sets type \c DenotableSet::basic_set_type.
 * \code
 * concept DenotableSet
 * {
 *   typename real_type;
 *   typename state_type;
 *   typename basic_set_type;
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
 *   dimension_type dimension() const;
 *   bool empty() const;
 *   bool contains(const state_type &) const;
 *
 *   Rectangle<real_type> bounding_box() const; // Optional.
 *
 *   void adjoin(const basic_set_type &);
 *   void adjoin(const DenotableSet &);
 *
 *   // List operations
 *   const_iterator begin() const;
 *   const_iterator end() const;
 *
 *   size_type size() const; // Only required if the iterator is a RandomAccessIterator.
 *   basic_set_type operator[] (size_type) const; // Only required if the iterator is a RandomAccessIterator.
 *
 *   void push_back(const basic_set_type &); // Only used if the DenotableSet is an ordered list. (Optional)
 *   basic_set_type pop_back(); // Only used if the DenotableSet is an ordered list. (Optional)
 *
 *   void insert(const basic_set_type &); // Only used if the DenotableSet is an unordered or sorted list. (Optional)
 *   void remove(const basic_set_type &); // Only used if the DenotableSet is an unordered or sorted list. (Optional)
 * };
 *
 * DenotableSet join(const DenotableSet &, const DenotableSet &);
 * DenotableSet regular_intersection(const BasicSet &, const DenotableSet &); // Optional.
 *
 * bool subset(const BasicSet &, const DenotableSet &); // Optional, but highly recommended.
 *
 * bool interiors_intersect(const DenotableSet &, const DenotableSet &); // Optional, but highly recommended.
 * bool disjoint(const DenotableSet &, const DenotableSet &);
 * bool subset(const DenotableSet &, const DenotableSet &); // Optional, but highly recommended.
 *
 * \endcode
 */

/*! \page evaluation Function Evaluation
 * 
 * Ariadne is primarily a module for set-based computationan.
 * For this reason, functions are best defined by their actions on sets.
 * However, we sometimes also want to compute function values on points, and to evaluate real-valued functions.
 * For this reason, we also allow definition of functions on points.
 * 
 * We distinguish between computations on \em fixed-precision and \em arbitrary-precision types, 
 * between computations on \em points and \em sets, and between \em exact and \em approximate computations.
 *
 * The basic computation on sets is to compute \em over-approximations to the image of <em>basic sets</em>
 * From these computations, arbitrarily-accurate computations can be performed.
 * The basic computation on points can be \em exact and \em approximate, as appropriate.
 *
 * \section set_functions  Computations on sets.
 *
 * A valid arbitrary-precision computation on sets is defined as follows
 * If \f$f\f$ is a mathematical function, and \c BS is a basic set type,
 * then a valid representation \c f of \f$f\f$ is a function which, 
 * for any basic set \f$A\f$, returns a set \f$B\f$ such that \f$f(A)\subset B\f$,
 * and such that whenever \f$A_n\f$ is a decreasing sequence of sets with \f$\bigcap_{n=1}^{\infty} A_n=\{x\}\f$,
 * then \f$\bigcap_{n=1}^{\infty} B_n=\{y\}\f$, where \f$y=f(x)\f$.
 *
 * A valid fixed-precision computation on sets is an over-approximation. 
 * In other words, a valid representation \c f of \f$f\f$ is a function which,
 * for any basic set \f$A\f$, returns a set \f$B\f$ such that \f$f(A)\subset B\f$.
 * No guarentees on the accuracy are required.
 * Note that it does not make sense to consider a sequence \f$A_n\f$ converging to a point for fixed-precision types. 
 *
 * \section point_functions Computations on points.
 *
 * A arbitrary-precision computation on points may be \em exact or \em approximate.
 * Examples of exact operations are polynomial functions on a ring (e.g. dyadic numbers)
 * and rational functions on a field (e.g. rational numbers). If \f$f\f$ is a mathematical function, 
 * then \c f(x) computes \f$f(x)\f$ exactly if possible, and is undefined (compile-time or run-time error) otherwise.
 * \c f(x,e) computes \f$f(x)\f$ with an error of at most \c e.
 * 
 * Note that even if \c f is exact, it is impossible to compute 
 * the action of \f$f\f$ on a set just from the action of \c f on points,
 * unless a modulus of continuity for \f$f\f$ is known.
 *
 * All fixed-precision computations on points are approximate. Further, the accuracy of the approximation is often unknown,
 * or hard to compute. Since the aim of the fixed-precision computation is rapid computation, no guarentees are given about 
 * the error. We use the syntax \c f(x) to compute an approximation to \f$f(x)\f$ with no controls on the accuracy.
 *
 *
 * \section real_functions Computations on real numbers.
 *
 * Most continuous functions used in science and engineering are built up from elementary real-valued functions,
 * either arithmetical or defined using Taylor series expansions. 
 * For this reason, Ariadne provides extended operations for computation on real-valued functions.
 * 
 * Arbitrary-precision computations may be exact or approximate. 
 * The function \c f(x) computes \f$f(x)\f$ exactly, if possible, 
 * and gives an error (at compile-time or run-time) if the result cannot be computed exactly.
 * The function \c f_approx(x,e) computes \f$f(x)\f$ with an error of at most \f$e\f$.
 *
 * Although it is not, in general, possible to perform evaluation of functions on sets from their definitions on points,
 * in many cases such a computation can be extracted. 
 * In particular, we can construct interval computations from pointwise computations in the following cases:
 * <ul>
 *    <li>The function is Lipschitz with a known Lipschitz constant e.g. \f$\sin,\ \cos\f$.</li>
 *    <li>The function is monotone e.g. \f$\exp\f$.</li>
 *    <li>The function is piecewise-monotone with known branches e.g. arithmetic.</li>
 * </ul>
 *
 *
 * Fixed-precision computations are inherantly approximate, 
 * though in some cases they may happen to give exact result. 
 * The function \c f_approx(x) computes \f$f(x)\f$ approximately, with no control on the error bound.
 * The function \c f_down(x) computes a lower-approximation to \f$f(x)\f$, 
 * and \c f_up(x) computes an upper-approximation.
 * For consistency with existing practise, we also allow \c fnc(x) as a valid alternative to \c fnc_approx(x).
 *
 * Note that in many cases, including arithmetic and simple functions, it is possible to compute an interval \f$J\f$ containing \f$f(I)\f$ 
 * using \c f_down and \c f_upp. This allows an implementation of the standard set-based function \c f(I).
 * 
 * \section function_syntax Syntax for continuous functions and function objects.
 *
 * Computable functions are represented as function objects, i.e. objects in a
 * class with an overloaded operator().
 * Elementary and arithmetical functions are represented as ordinary functions
 * with the same syntax.
 *
 * \code
 * // Syntax for built-in real-valued functions.
 * FPReal f(FPReal x); // A synonym for f_approx(x)
 * FPReal f_approx(FPReal x); // Compute an approximation to f(x).
 * FPReal f_down(FPReal x); // Compute a lower-approximation to f(x). 
 * FPReal f_up(FPReal x); // Compute an upper-approximation to f(x). 
 *
 * APReal f(APReal x); // Compute f(x) exactly, if possible.
 * APReal f(APReal x, APReal e); // A synonym for f_approx(x,e).
 * APReal f_approx(APReal x, APReal e); // Compute f(x) with an error of at most e.
 *
 * // Syntax for user-defined single-valued continuous functions with fixed-precision arithmetic.
 * Point<FPReal> f(Point<FPReal> A); // A synonym for f.approx(x)
 * Point<FPReal> f.approx(Point<FPReal A); // Compute an approximation to f(x).
 * Point<FPReal> f.approximate(Point<FPReal> A); // Compute an approximation to f(x).
 *
 * BasicSet<FPReal> f(BasicSet<FPReal> A); // Compute an over-approximation to f(A).
 *
 * // Syntax for user-defined single-valued continuous functions with arbitrary-precision arithmetic.
 * Point<APReal> f(Point<APReal> x); // Compute f(x) exactly, if possible.
 * Point<APReal> f(Point<APReal> x, APReal e); // Compute f(x) with an error of at most e.
 * Point<APReal> f.approx(Point<APReal> x, APReal e); // A more descriptive syntax for f(x,e)
 * Point<APReal> f.approximate(Point<APReal> x, APReal e); // A more descriptive syntax for f(x,e)
 *
 * BasicSet<APReal> f(BasicSet<APReal> A); // Compute a convergent over-approximation to f(A).
 * BasicSet<APReal> f(BasicSet<APReal> A, OverApproximation); // Compute a convergent over-approximation to f(A).
 * BasicSet<APReal> f.approx(BasicSet<APReal> A); // A more descriptive syntax for f(A).
 * BasicSet<APReal> f.approximate(BasicSet<APReal> A); // A more descriptive syntax for f(A).
 * BasicSet<APReal> f.approximation(BasicSet<APReal> A); // A more descriptive syntax for f(A).
 * BasicSet<APReal> f.over(BasicSet<APReal> A); // A more descriptive syntax for f(A).
 * BasicSet<APReal> f.over_approximate(BasicSet<APReal> A); // A more descriptive syntax for f(A).
 *
 * // Definition for fixed-precision continuous functions.
 * template<class FPReal>
 * class ContinuousFunction {
 *   // Compute an approximation to f(p) with no guarentees on error.
 *   Point<FPReal> operator() (Point<FPReal> p); 
 *   Point<FPReal> approximate() (Point<FPReal> p); 
 *   Point<FPReal> approx() (Point<FPReal> p); 
 *
 *   // Compute an over-approximation to f(A) with no guarentees on error.
 *   BasicSet<FPReal> operator() (BasicSet<FPReal> A); 
 * };
 *
 *
 * // Definition for arbitrary-precision continuous functions.
 * template<class APReal>
 * class ContinuousFunction {
 *   // Compute f(p) exactly, but only if possible.
 *   Point<APReal> operator() (Point<APReal> p); 
 *
 *   // Compute an approximation to f(p) with error at most e.
 *   Point<APReal> operator() (Point<APReal> p, APReal e); 
 *   Point<APReal> approximate() (Point<APReal> p, APReal e); 
 *
 *   // Compute an over-approximation to f(A). The result is guarenteed to 
 *   // converge to a one-point set as the argument converges to a one-point set.
 *   BasicSet<APReal> operator() (BasicSet<APReal> A); 
 * };
 *
 * \endcode
 *
 * \internal 
 *    Use \c f.approx(p) or f.approximate(p) for function objects?
 *    Use \c f_approx(p) or f_approximate(p) for functions?     
 * 
 */
 
/*! \page integration Integration methods
 *
 * \section taylor Taylor methods 
 * All integration methods are based on Taylor expansion of solutions curves.
 * 
 * \f[ \begin{array}{rl} \displaystyle
 *      \frac{dx}{dt} &= f(x(t)) \\
 *  \frac{d^2x}{dt^2} &= Df(x(t))f(x(t)) \\
 *  \frac{d^3x}{dt^3} &= D^2f(x(t))f(x(t))f(x(t))+Df(x(t))Df(x(t))f(x(t)) \\
 *  \frac{d^4x}{dt^4} &= D^3f(x(t))f(x(t))f(x(t))f(x(t)) + 4D^2f(x(t))Df(x(t))f(x(t))f(x(t)) + Df(x(t))Df(x(t))Df(x(t))f(x(t)) \\
 *  \frac{d^5x}{dt^5} &= D^4f(x(t))f(x(t))f(x(t))f(x(t))f(x(t)) + 7D^3f(x(t))Df(x(t))f(x(t))f(x(t))f(x(t)) \\
 *                    &\ \qquad + 4D^2f(x(t))D^2f(x(t))f(x(t))f(x(t))f(x(t)) + 11D^2f(x(t))Df(x(t))Df(x(t))f(x(t))f(x(t)) \\
 *                    &\ \qquad + Df(x(t))Df(x(t))Df(x(t))Df(x(t))f(x(t)) \end{array} \f]
 *
 *
 * \subsection euler Euler method
 * \f[ x_1 = x_0 + hf(x_0) \approx x(h) \f]
 * One-step error
 * \f[ ||x_1-x(h)|| = h ||f(x_0)-f(\xi)|| = \frac{h^2}{2} || Df(\xi)f(\xi) || \f]
 *
 * \subsection second_order_taylor 2nd Order Taylor Method
 * \f[ x_1 = x_0 + hf(x_0) + \frac{h^2}{2} Df(x_0)f(x_0) \approx x(h) \f]
 * One-step error
 * \f[ ||x_1-x(h)|| = \frac{h^2}{2} || Df(x_0)f(x_0)-Df(\xi)f(\xi)|| = \frac{h^3}{6} || D^2f(\xi)f(\xi)f(\xi) + Df(\xi)Df(\xi)f(\xi) || \f]
 *
 * \subsection second_order_rk 2nd Order Runge-Kutta Method
 * 
 * \f[ x_1 = x_0 + \frac{h}{2}\left( f(x_0)+f(x_0+hf(x_0))\right) \approx x_0 + hf(x_0) + \frac{h^2}{2} Df(x_0)f(x_0) 
 *      + \frac{h^3}{4} D^2f(\xi)f(\xi)f(\xi) \f]
 */

/*! \page zonotope Zonotopic reduction methods
 *
 * Throughout this sections, we use the supremum norm on \f$R^n\f$, and the correspoinding operator norm on \f$\mathbb{R}^{m\times n}\f$.
 * 
 * Given a zonotope \f$ Z=\{ c+Ae \mid ||e||\leq 1 \}\subset \mathbb{R}^n\f$, where \f$A\in \mathbb{R}^{n\times p}\f$, 
 * we wish to compute a zonotope \f$Z' = \{ c + A' e' \mid ||e'||\leq 1\}\f$ with fewer generators 
 * i.e. \f$A'\in \mathbb{R}^{n\times p'}\f$ with \f$p'<p\f$.
 * The general reduction method is to choose \f$ A'\f$ such that \f$ A = A' B\f$ with \f$||B||\leq 1\f$.
 * The key to zonotopic reduction is to choose a method with good properties.
 *
 * A simple criterion to note is that if \f$\sum_{j=1}^{p} |b_{ij}|<1\f$ for some \f$i\f$, then we can improve the approximation by taking
 * \f$D=\mathrm{diag}(d_{i})\f$ with \f$d_{i}=\sum_{j=1}^{p} |b_{ij}|\f$, 
 * and \f$B'= D^{-1}B\f$ which has \f$b'_{ij}=b_{ij}/\sum_{k=1}^{p} |b_{ik}|\f$.
 * 
 * It is clear that if the rows of \f$B\f$ are close to a set of mutually orthogonal coordinate vectors, then the approximation is good, 
 * since the image of \f$B\f$ is close to the unit ball. This suggests the following algorithm:
 * 
 *
 * \section interval_zonotope Interval zonotopic reduction
 *
 * An <em>interval zonotope</em> is a set of the form \f$ \{ y = c + A e \mid c\in R,\ A\in\mathcal{A} \text{ and } ||e||\leq1 \} 
 * = R + \mathcal{A} B\f$.
 * To reduce an interval zonotope, we first write \f$ R = \{ c + Be\mid ||e||\leq 1 \}\f$ and combine this in \f$\mathcal{A}\f$.
 * To reduce \f$ \mathcal{A} \f$, write \f$ \mathcal{A} = A \mathcal{B} \mathcal{A} \f$ where \f$ A\in\mathcal{A}\f$ and \f$A\mathcal{B}\ni I\f$.
 * Then take \f$ \mathcal{C} = \mathcal{B} \mathcal{A} \f$ and \f$ || \mathcal{C} || 
 *   = \sup_{i]1}^{n} \sum_{j=1}^{p} \max |\mathcal{C}_{ij}| \f$, 
 * where \f$ \max |\mathcal{C}_{ij}| = \max\{ |x| \mid x\in \mathcal{C}_{ij}\f$.
 * 
 */
 
/*! \page references References
 *
 * \section computable_analysis_references Computable Analysis
 *
 * Ker-I Ko <em>Complexity Theory of Real Functions</em>, Birkh\"aser, 1991, ISBN 3-7643-3586-6.
 * 
 * Klaus Weihrauch, <em>Computable Analysis</em>, Springer, 2000.
 * 
 * \section interval_references Interval Arithmetic
 * Ramon E. Moore, <em>Methods and applications of interval analysis</em>,
 *   SIAM Studies in Applied Mathematics, 2.
 * Society for Industrial and Applied Mathematics (SIAM), Philadelphia, Pa., 1979. xi+190 pp. ISBN 0-89871-161-4 
 *
 * Baker R. Kearfott, "Interval computations: introduction, uses, and resources",
 *   <em>Euromath Bull.</em> <b>2</b> (1996), no. 1, 95--112. 
 *
 * Marcel Gavriliu, "Towards more efficient interval analysis: corner forms and a remainder Newton method", 
 *   Ph.D. Thesis, California Institute of Technology, 2005. <br>
 *
 * R. Krawczyk, "A class of interval-Newton-operators",
 *   <em>Computing</em> <b>37</b> (1986), no. 2, 179--183.
 *
 * Arnold Neumaier, <em>Interval methods for systems of equations</em>,
 *   Encyclopedia of Mathematics and its Applications, 37.
 *   Cambridge University Press, Cambridge, 1990. xvi+255 pp. ISBN 0-521-33196-X 
 *
 * A. Neumaier, "The wrapping effect, ellipsoid arithmetic, stability and confidence regions",
 *   <em>Computing Supplementum</em> <b>9</b> (1993), 175-190.
 * 
 * \section zonotope_references Zonotopes
 *
 * K. Fukuda, "From the zonotope construction to the Minkowski addition of convex polytopes", Preprint, 2003. 
 *
 * Leonidas J. Guibas, An Nguyen and Li Zhang, "Zonotopes as bounding volumes", Preprint.
 *
 * Ari Ingimundarson, Jose Manuel Bravo, Vicenc Puig and Teodoro Alama, 
 *   "Robust Fault Diagnosis using Parallelotope-based Set-membership Consistency Tests",
 *   In <em>Proceedings of CDC-ECC 2005</em>.
 *
 * \section integration_references Integration
 *
 * Rudolf J. Lohner, "Enclosing the solutions of ordinary initial and boundary value problems",
 *   <em>Computer Arithmetic</em>, 255--286, Teubner, Stuttgart, 1987. 
 *
 * Rudolf J. Lohner, "Computation of guaranteed enclosures for the solutions of ordinary initial and boundary value problems",
 *   <em>Computational ordinary differential equations (London, 1989)</em>,  425--435.
 *
 * N. S. Nedialkov, K. R. Jackson and G. F. Corliss,
 *   "Validated solutions of initial value problems for ordinary differential equations",
 *   <em>Appl. Math. Comput.</em> <b>105</b> (1999), no. 1, 21--68.
 *
 * Piotr Zgliczynski, "C^1 Lohner algorithm", <i>Found. Comput. Math.</i> <b>2</b> (2002), no. 4, 429--465.
 *
 * \section reachability_references Reachability Analysis and Control
 *
 * Antoine Girard, Colas Le Guernic and Oded Maler, 
 * "Efficient Computation of Reachable Sets of Linear Time-Invariant Systems with Inputs",
 * 
 * Antoine Girard, "Reachability of Uncertain Linear Systems using Zonotopes," 
 *   in <em>Proceedings of HSCC 2005</em>, LNCS 3414, pp 291--305, 2005.
 * 
 * Alex Kurzhanskiy and Pravin Varaiya, "Ellipsoidal Techniques for Reachability Analysis of Discrete-Time Linear Systems",
 *
 * Alexander Kurzhanski and Pravin Varaiya, "On ellipsoidal techniques for reachability analysis", 
 *   <em>Optim. Methods Softw.</em> <b>17</b> (2002), no. 2, 207--237
 *
 * S. V. Rakovic and D. Q. Mayne, "Set Robust Control Invariance for Linear Discrete Time Systems", 
 *   in <em>Proceedings of CDC-ECC 2005</em>.
 *
 * F. Lydoire and P. Poignet, "Nonlinear Model Predictive Control via Interval Analysis", 
 *   in <em>Proceedings of CDC-ECC 2005</em>.
 *
 * Dietmar Szolnoki, "Set oriented methods for computing reachable sets and control sets",
 *   <em>Discrete Contin. Dyn. Syst. Ser. B</em> <b>3</b> (2003), no. 3, 361--382.
 *
 * \section spacial_data_structure_references Spacial Data Structures
 *
 * Hanan Samet, <em>The Design and Analysis of Spacial Data Structures</em>, Addison-Wesley, 1990, ISBN 0-201-50255-0.
 *
 * \section algebraic_topology_references Algebraic Topology
 *
 * Tomasz Kaczynski, Konstantin Mischaikow, Marian Mrozek, <em>Computational Homology</em>, Springer-Verlag, 2004, ISBN 0-387-40853-3.
 *
 * Afra J. Zomorodian, <em>Topology for Computing</em>, Cambridge University Press, 2005, ISBN 0-521-83666-2.
 *
 * \section automatic_differentiation_references Automatic Differentiation
 *
 * Andreas Giewank, <em>Evaluating Derivatives</em>, SIAM, 2000, ISBN 0-89871-451-6.
 */
 
#endif /* _ARIADNE_DOCUMENTATION_H */
