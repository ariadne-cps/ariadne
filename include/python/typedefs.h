/***************************************************************************
 *            typedefs.h
 *
 *  06 Feb 2006
 *  Copyright  2005  Pieter Collins
 *  pieter.collins@cwi.nl
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
 
/*! \file typedefs.h
 *  \brief Typedefs for python modules.
 */

#ifndef _ARIADNE_PYTHON_TYPEDEFS_H
#define _ARIADNE_PYTHON_TYPEDEFS_H

#include "../linear_algebra/linear_algebra_declarations.h"
#include "../geometry/geometry_declarations.h"
#include "../evaluation/evaluation_declarations.h"

#include "../python/real_typedef.h"

typedef Ariadne::Dyadic Dyadic;
typedef Ariadne::Rational Rational;

typedef Ariadne::Interval<Real> RInterval;

typedef Ariadne::LinearAlgebra::vector<Rational> QVector;
typedef Ariadne::LinearAlgebra::matrix<Rational> QMatrix;

typedef Ariadne::LinearAlgebra::vector<Real> RVector;
typedef Ariadne::LinearAlgebra::matrix<Real> RMatrix;
typedef Ariadne::LinearAlgebra::interval_vector<Real> RIntervalVector;
typedef Ariadne::LinearAlgebra::interval_matrix<Real> RIntervalMatrix;

typedef Ariadne::Geometry::Point<Real> RPoint;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;
typedef Ariadne::Geometry::Parallelotope<Real> RParallelotope;
typedef Ariadne::Geometry::Zonotope<Real> RZonotope;
typedef Ariadne::Geometry::Simplex<Real> RSimplex;
typedef Ariadne::Geometry::Polyhedron<Real> RPolyhedron;

typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Rectangle> RRectangleListSet;
typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Parallelotope> RParallelotopeListSet;

typedef Ariadne::Geometry::Grid<Real> RGridBase;
typedef Ariadne::Geometry::FiniteGrid<Real> RFiniteGrid;
typedef Ariadne::Geometry::GridCell<Real> RGridCell;
typedef Ariadne::Geometry::GridRectangle<Real> RGridRectangle;
typedef Ariadne::Geometry::GridCellListSet<Real> RGridCellListSet;
typedef Ariadne::Geometry::GridRectangleListSet<Real> RGridRectangleListSet;
typedef Ariadne::Geometry::GridMaskSet<Real> RGridMaskSet;

typedef Ariadne::Geometry::PartitionScheme<Real> RPartitionScheme;
typedef Ariadne::Geometry::PartitionTree<Real> RPartitionTree;
typedef Ariadne::Geometry::PartitionTreeCell<Real> RPartitionTreeCell;
typedef Ariadne::Geometry::PartitionTreeSet<Real> RPartitionTreeSet;

typedef Ariadne::Evaluation::Map<Real> RMap;
typedef Ariadne::Evaluation::AffineMap<Real> RAffineMap;
typedef Ariadne::Evaluation::VectorField<Real> RVectorField;
typedef Ariadne::Evaluation::AffineVectorField<Real> RAffineVectorField;

typedef Ariadne::Evaluation::HenonMap<Real> RHenonMap;
typedef Ariadne::Evaluation::LorenzSystem<Real> RLorenzSystem;

#endif /* _ARIADNE_PYTHON_TYPEDEFS_H */
