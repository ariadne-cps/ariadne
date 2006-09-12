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

/*! \addtogroup Base Base
 *
 *  \defgroup LinearAlgebra Linear Algebra
 *  \ingroup Base
 *  \brief Vector, Matrix and Tensor classes, LU and QR factorisation, singular value decomposition, and linear programming.
 *
 */

/*! \defgroup Numeric Numeric
 *  \ingroup Base
 *  \brief Numerical types and intervals.
 */

/*! \defgroup Combinatoric Combinatoric
 *  \brief Combinatoric sets and maps.
 *
 *  \defgroup Lattice Lattice Sets and Functions
 *  \ingroup Combinatoric
 *  \brief Classes for integer lattices.
 *
 *  \defgroup SubdivisionTree Subdivision Tree Sets and Functions
 *  \ingroup Combinatoric
 *  \brief Classes for binary subdivision trees.
 *
 */
 
/*! \defgroup Geometry Geometry
 *  \brief Geometric calculus module.
 *
 *
 *  \defgroup BasicSet Basic Sets
 *  \ingroup Geometry
 *  \brief Basic set classes.
 *
 *
 *  \defgroup DenotableSet Denotable Sets
 *  \ingroup Geometry
 *  \brief Denotable set classes.
 *
 *
 *  \defgroup List List Sets 
 *  \ingroup DenotableSet
 *  \brief Denotable sets based on arbitrary lists of elements.
 *
 *
 *  \defgroup Grid Grid Sets 
 *  \ingroup DenotableSet
 *  \brief Sets based on grids.
 *
 *
 *  \defgroup PartitionTree Partition Tree Sets
 *  \ingroup DenotableSet
 *  \brief Sets based on partition trees.
 */

/*! \defgroup System System
 *  \brief Abstract base classes for system interface, and some commonly used systems.
 *
 *  \defgroup DiscreteTime Discrete-Time Systems
 *  \ingroup System
 *  \brief Discrete-time systems.
 *
 *  \defgroup ContinuousTime Continuous-Time Systems
 *  \ingroup System 
 *  \brief Continuous-time systems.
 *
 *  \defgroup HybridTime Hybrid-Time Systems
 *  \ingroup System 
 *  \brief Hybrid-time systems.
 */

/*! \defgroup Evaluation Evaluation
 *  \brief Functions and methods for computing the evolution of a system and solving equations.
 *
 *  \defgroup Apply Apply
 *  \ingroup Evaluation
 *  \brief Functions for iterating forward discrete-time systems.
 *
 *  \defgroup Integrate Integrate
 *  \ingroup Evaluation
 *  \brief Classes for integrating discrete-time systems.
 *
 *  \defgroup Solve Solve
 *  \ingroup Evaluation
 *  \brief Functions for solving systems of equations.
 */

/*! \defgroup Output Output
 *
 *  \defgroup Postscript Postscript Output
 *  \ingroup Output 
 *  \brief Encapsulated Postscript output.
 */


 
#endif /* _ARIADNE_MODULES_H */
