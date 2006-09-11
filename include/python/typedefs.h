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

#include "../real_typedef.h"

#include "../numeric/float64.h"
#include "../numeric/mpfloat.h"
#include "../numeric/dyadic.h"
#include "../numeric/rational.h"

#include "../declarations.h"

namespace Ariadne {
  
typedef Numeric::Float64 Float64;
typedef Numeric::MPFloat MPFloat;
typedef Numeric::Dyadic Dyadic;
typedef Numeric::Rational Rational;

typedef Interval<Real> RInterval;
typedef Interval<Field> FInterval;

typedef LinearAlgebra::Vector<Real> RVector;
typedef LinearAlgebra::Matrix<Real> RMatrix;
typedef LinearAlgebra::IntervalVector<Real> RIntervalVector;
typedef LinearAlgebra::IntervalMatrix<Real> RIntervalMatrix;

typedef LinearAlgebra::Vector<Field> FVector;
typedef LinearAlgebra::Matrix<Field> FMatrix;
typedef LinearAlgebra::LinearProgram<Field> FLinearProgram;

typedef Geometry::Point<Real> RPoint;
typedef Geometry::PointList<Real> RPointList;
typedef Geometry::Rectangle<Real> RRectangle;
typedef Geometry::Parallelotope<Real> RParallelotope;
typedef Geometry::Zonotope<Real> RZonotope;
typedef Geometry::Simplex<Real> RSimplex;
typedef Geometry::Polyhedron<Real> RPolyhedron;
typedef Geometry::Polytope<Real> RPolytope;

typedef Geometry::ListSet<Real,Geometry::Rectangle> RRectangleListSet;
typedef Geometry::ListSet<Real,Geometry::Parallelotope> RParallelotopeListSet;
typedef Geometry::ListSet<Real,Geometry::Zonotope> RZonotopeListSet;
typedef Geometry::ListSet<Real,Geometry::Polyhedron> RPolyhedronListSet;

typedef Geometry::LatticeCell LatticeCell;
typedef Geometry::LatticeRectangle LatticeRectangle;
typedef Geometry::LatticeMaskSet LatticeMaskSet;
typedef Geometry::LatticeCellListSet LatticeCellListSet;
typedef Geometry::LatticeRectangleListSet LatticeRectangleListSet;

typedef Geometry::grid_type grid_type;
typedef Geometry::Grid<Real> RGrid;
typedef Geometry::RegularGrid<Real> RRegularGrid;
typedef Geometry::IrregularGrid<Real> RIrregularGrid;
typedef Geometry::FiniteGrid<Real> RFiniteGrid;
typedef Geometry::GridCell<Real> RGridCell;
typedef Geometry::GridRectangle<Real> RGridRectangle;
typedef Geometry::GridCellListSet<Real> RGridCellListSet;
typedef Geometry::GridRectangleListSet<Real> RGridRectangleListSet;
typedef Geometry::GridMaskSet<Real> RGridMaskSet;

typedef Geometry::PartitionScheme<Real> RPartitionScheme;
typedef Geometry::PartitionTree<Real> RPartitionTree;
typedef Geometry::PartitionTreeCell<Real> RPartitionTreeCell;
typedef Geometry::PartitionTreeSet<Real> RPartitionTreeSet;

typedef System::LatticeMultiMap LatticeMultiMap;

typedef System::Map<Real> RMapBase;
typedef System::AffineMap<Real> RAffineMap;
typedef System::AffineMultiMap<Real,Geometry::Rectangle> RAffineMultiMapRect;
typedef System::AffineMultiMap<Real,Geometry::Parallelotope> RAffineMultiMapPltp;
typedef System::AffineMultiMap<Real,Geometry::Zonotope> RAffineMultiMapZntp;
typedef System::Monomial<Real> RMonomial;
typedef System::Polynomial<Real> RPolynomial;
typedef System::PolynomialMap<Real> RPolynomialMap;
typedef System::PolynomialMatrix<Real> RPolynomialMatrix;
typedef System::VectorField<Real> RVectorFieldBase;
typedef System::AffineVectorField<Real> RAffineVectorField;

typedef System::HenonMap<Real> RHenonMap;
typedef System::LorenzSystem<Real> RLorenzSystem;

}

#endif /* _ARIADNE_PYTHON_TYPEDEFS_H */
