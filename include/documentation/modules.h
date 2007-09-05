/***************************************************************************
 *            modules.h
 *
 *  Copyright  2006  Pieter Collins
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

#ifndef _ARIADNE_MODULES_H
#define _ARIADNE_MODULES_H

/*! \file modules.h
 * \brief Documentation for %Ariadne modules.
 */

/*!\addtogroup Base Base Module
 * \brief Fundamental classes and operations, mostly implemented as wrappers 
 * around other libraries.
 *
 * \defgroup Storage Storage
 * \ingroup Base
 * \brief Array classes for data storage
 * \addtogroup Storage
 * \class array<bool>
 * 
 *
 * \defgroup Traversal Traversal
 * \ingroup Base
 * \brief Iterator classes for data traversal
 * 
 * \defgroup Numeric Numeric
 * \ingroup Base
 * \brief Numerical types and intervals.
 *
 * In this module, we define %Ariadne's base numerical types, and interval 
 * arithmetic. The main classes are the various <em>real number types</em>
 * which must either conform to the ExactArithmetic interface or the 
 * ApproximateArithmetic interface. This use of standard operations allows
 * new classes to be easily dropped in as wrappers. Currently, the Float64 and
 * FloatMP classes support approximate arithmetic, algebraic and transcendental
 * functions, and the Dyadic and Rational classes support exact arithmetic, but
 * no algebraic or transcendental operations. Further, Dyadic does not support
 * division, but does support the median() operation.
 *
 * Interval arithmetic and functions are provided through the Interval<Real>
 * template. Note that Interval<Real> is not part of the Geometry module, so
 * does not support the same geometric operations.
 *
 * Additionally, rounded operations are provided using the round_down
 * and round_up rounding types. Approximate arithmetic without control of the
 * rounding mode can be specified by the approx rounding type, and exact 
 * arithmetic (the default) can be specified explicitly using the exact 
 * rounding type.
 *
 * Note that to prevent accidental calls to approximate functions, all 
 * non-exact operations must be explicitly qualified. This means that the
 * operation operator*(Float64,Float64) generates an error at compile-time,  
 * even though Float64 is based on the C++ double type, which has an 
 * (approximate) multiplication operator. Use 
 * mul(Float64,Float64,RoundApproximate) instead.
 * 
 * We use the following abbreviations for template parameters 
 * for classes and concepts in this module:
 *  - \b Q: Rational
 *  - \b F: Float 
 *  - \b I: Interval<Float>
 *  - \b A: Rational, Interval<Float>
 *  - \b R: Rational, Float
 *  - \b X: Rational, Float, Interval<Float>
 *
 * Here, Float refers to any floating-point type, such as Float64 or FloatMP.
 * 
 * \internal The above template names are currently not used consistently
 *
 * See the page \ref realnumbertypes for more information.
 *
 * \defgroup LinearAlgebra Linear Algebra
 * \ingroup Base
 * \brief Vector, matrix and tensor classes, solution of linear equations,
 * LU and QR factorisation and singular value decomposition.
 * 
 * Basic types are represented by the classes Vector, Matrix
 * and Tensor. Specialist types include PermutationMatrix and DiagonalMatrix.
 * Operations such as vector sums and matrix-vector products are
 * accessed using operator overloading. Expression templates and slices are 
 * used to improve efficiency.
 * 
 * More complicated matrix operations, such as factorisation, are represented
 * by sub-classes which contain the factorised data, from which the factors 
 * can easily be extracted.
 *
 * Currently sparse and structured (e.g. triangular, symmetric, banded) 
 * matrices are not directly supported, except to the extent needed for the
 * matrix decomposition classes, or for derivative tensors.
 *
 * All classes support exact arithmetic and interval arithmetic.
 *
 * We use the following abbreviations for template parameters for classes 
 * in this module:
 *  - \b SCA: Scalar
 *  - \b VEC: Vector
 *  - \b MX: Matrix
 *  - \b TEN: Tensor
 * 
 * \defgroup LinearProgramming Linear Programming
 * \ingroup Base
 * \brief Classes and functions for solving linear programming problems (including feasibility problems).
 *
 * The LinearProgram class supports construction and solution of linear 
 * programming problems. Tests for feasibility only are also possible. The 
 * computations are performed by a BLAS/LAPACK style lpslv() routing.
 *
 */

/*!\defgroup Combinatoric Combinatoric Module
 * \brief Combinatoric sets and maps, based either on integer lattices or
 * partition trees.
 *
 * 
 * \defgroup ArrayOperations Array Operations
 * \ingroup Combinatoric
 * \brief Operations on integer arrays
 *
 * \defgroup BinaryTree Binary Words and Trees
 * \ingroup Combinatoric
 * \brief Classes for binary words and binary trees.
 *
 * 
 * \defgroup Lattice Lattice Sets and Functions
 * \ingroup Combinatoric
 * \brief Classes for integer lattices.
 *
 * An <em>integer lattice</em> is a subset of \f$\mathbb{Z}^n\f$. Integer 
 * lattices provide a representation for sets in \f$\mathbb{R}^n\f$ based on
 * coordinate-aligned grids.
 *
 * As well as sets of cells, %Ariadne also supports maps between lattices.
 * The most useful form takes a single grid cell to a set of grid cells, and
 * so is a type of multivalued map. Hence, no LatticeMap class is given.
 *
 * In a future version, methods will be provided to export integer lattice sets
 * and maps to the CHomP family of computational homology programs.
 *
 * 
 * \defgroup SubdivisionTree Subdivision Tree Sets and Functions
 * \ingroup Combinatoric
 * \brief Classes for binary subdivision trees.
 *
 * A <em>subdivision tree</em> is a binary tree describing a subdivision of 
 * a rectangle (<em>unit cell</em>) into smaller cells. The subdivision is 
 * described by a SubdivisionSequence giving the coordinate of the ith 
 * subdivision, and a BinaryTree giving the subdivisions.
 *
 * e.g. The cell described by the SubdivisionSequence '[0,0,1,...]' and the 
 * BinaryWord '[0,1,1]' is given by subdividing the unit interval twice in the 
 * \f$x_0\f$ direction, first taking the lower half, and then the upper half, 
 * and by subdividing once in the \f$x_1\f$ direction, taking the upper half.
 * This gives the dyadic rectangle '[0.25,0.5]x[0.5,1.0]'.
 */
 
/*!\defgroup Function Function Module
 * \brief Functions on Euclidean space.
 *
 * In this module are the classes used to represent concrete functions in Euclidean space.
 * These functions can be used as the building blocks for geometric objects and dynamic systems.
 *
 * \defgroup FunctionTypes Function Types
 * \ingroup Function
 * \brief General and specialised function classes.
 *
 * \defgroup FunctionModel Function Models
 * \ingroup Function
 * \brief Approximate models of general (differentiable) functions.
 */

/*!\defgroup Geometry Geometry Module
 * \brief Geometric calculus module.
 *
 * In this module are the classes used to represent sets and approximations to 
 * sets. Sets fall into two classes, <em>basic sets</em> which are typically
 * simple sets which form a base for the topology, and <em>denotable sets</em>
 * which are unions of basic sets.
 * 
 * The fundamental geometric operations are \c contains(Set,State), 
 * \a disjoint(Set,Set) and subset(Set,Set). All these return \a tribool values,
 * where the indeterminate value is used to indicate a result which is either
 * not robust (e.g. the boundaries of two sets intersect but their interiors do
 * not) or which cannot be computed to the given precision.
 *
 * %Ariadne uses <em>fuzzy basic sets</em>, which are sets of sets defined using
 * interval data, to store intermediate results if these cannot be computed 
 * exactly. These sets can be converted to ordinary basic sets by over- or
 * under-approximation. Alternatively, the fundamental binary predicates can
 * be computed directed for the fuzzy set types.
 * Denotable sets are never fuzzy.
 *
 * We use the following abbreviations for template parameters for classes in this module:
 *  - \b S: SetInterface
 *  - \b BS: BasicSet
 *  - \b DS: DenotableSet
 *  - \b LS: ListSet
 *  - \b GS: GridSet
 *  - \b PTS: PartitionTreeSet
 *
 * See \ref geometric for more information.
 * 
 * \defgroup BasicSet Basic Sets
 * \ingroup Geometry
 * \brief Basic set classes.
 *
 * Basic sets are so-called because they form a base for the topology of the
 * space. Typically, basic sets are (a subclass of) convex polytopes or 
 * ellipsoids. Basic sets must support the fundamental geometric predicates,
 * both within the class and with the Rectangle class. Additionally, basic sets
 * may support the optional geometric operations, but only if the class of set
 * under consideration is closed under that operation. The result may be 
 * exactly computable if it involves no arithmetic (e.g. intersection of two
 * rectangles) or may need to be represented by a fuzzy set (e.g. Minkowski sum
 * of two rectangles).
 *
 *
 * \defgroup DenotableSet Denotable Sets
 * \ingroup Geometry
 * \brief Denotable set classes.
 *
 * Denotable sets are unions of basic sets of a particular kind. In addition to
 * the fundamental geometric predicates, basic sets must also support iteration
 * through their elements, and union (join) with basic sets and denotable sets
 * of the same kind. 
 * 
 *
 * \defgroup List List Sets
 * \ingroup DenotableSet
 * \brief Denotable sets based on arbitrary lists of elements.
 *
 *
 * \defgroup Grid Grid Sets 
 * \ingroup DenotableSet
 * \brief Sets based on grids.
 *
 *
 * \defgroup PartitionTree Partition Tree Sets
 * \ingroup DenotableSet
 * \brief Sets based on partition trees.
 *
 * \defgroup SetInterface Set interfaces
 * \ingroup Geometry
 * \brief Interfaces for various classes of sets.
 *
 * \defgroup ExactSet Exact sets
 * \ingroup Geometry
 * \brief Sets which can be defined exactly.
 *
 * Any set can be uniquely specified by the results geometric predicates involving
 * rectangles. Exact sets are defined by an interface in which these predicates are
 * defined.
 *
 * \defgroup HybridSet Hybrid sets
 * \ingroup Geometry
 * \brief Sets in spaces with many components.
 *
 * The state space of a hybrid system consists of many disconnected components,
 * one for each discrete mode. Hybrid sets consist of a union of many different
 * sets, one for each component of the state space.
 */

/*!\defgroup System System Module
 * \brief Abstract base classes for system interface, and some commonly used systems.
 *
 * Systems fall into three basic classes, namely <em>discrete-time 
 * systems</em>, <em>continuous-time systems</em>, and <em>hybrid-time 
 * systems</em>, which are a combination of dicrete-time and hybrid-time 
 * systems. There is a further subdivision into
 * <em>dynamical systems</em>, which are deterministic, <em>multivalued
 * systems</em>, which are nondeterministic, and
 * <em>control systems</em> in which control inputs and noise may be present.
 * 
 * A system is defined by an algorithm to approximate the image of a basic set,
 * and possibly to approximate the derivative over a basic set.
 * It is <em>not</em> the purpose of the classes in the system module to
 * compute the image of denotable sets; this is the purpose of the Evaluation 
 * module.
 *
 * We use the following abbreviations for template parameters for classes in this module:
 *  - \b FN: %Function
 *  - \b MP: %Map
 *  - \b VF: %VectorField
 *
 * See the page on \ref function for more information.
 *
 * \defgroup DiscreteTime Discrete-Time Systems
 * \ingroup System
 * \brief Discrete-time systems.
 *
 * \defgroup ContinuousTime Continuous-Time Systems
 * \ingroup System 
 * \brief Continuous-time systems.
 *
 * \defgroup HybridTime Hybrid-Time Systems
 * \ingroup System 
 * \brief Hybrid-time systems.
 */

/*! \defgroup Evaluation Evaluation Module
 *  \brief Functions and methods for computing the evolution of a system and solving equations.
 *
 * 
 *  \defgroup Solve Solve
 *  \ingroup Evaluation
 *  \brief Functions for solving systems of equations.
 *
 *  \defgroup Apply Apply
 *  \ingroup Evaluation
 *  \brief Functions for iterating forward discrete-time systems.
 *
 *  \defgroup Integrate Integrate
 *  \ingroup Evaluation
 *  \brief Classes for integrating continuous-time systems.
 *
 *  \defgroup Evolve Evolve
 *  \ingroup Evaluation
 *  \brief Classes for computing the evolution of hybrid-time systems.
 *
 */

/*! \defgroup Input Input Module
 *
 *  \defgroup Modelica Modelica Input
 *  \ingroup Input 
 *  \brief Modelica input.
 */

/*! \defgroup Output Output Module
 *
 *  \defgroup Postscript Postscript Output
 *  \ingroup Output 
 *  \brief Encapsulated Postscript output.
 *
 *  \defgroup CHomP Homology Output
 *  \ingroup Output 
 *  \brief Homology Output.
 */


 
#endif /* _ARIADNE_MODULES_H */
