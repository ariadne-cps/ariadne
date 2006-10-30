/***************************************************************************
 *            declarations.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file declarations.h
 *  \brief Forward declarations of classes.
 */

#ifndef _ARIADNE_DECLARATIONS_H
#define _ARIADNE_DECLARATIONS_H

#include <iosfwd>

#define NO_CBLAS


namespace Ariadne { namespace Base {
template<class T, size_t N=0> class array;
template<class T> class sequence;
}}

namespace Ariadne { namespace Numeric {
class Integer;
class Float64;
class MPFloat;
class Rational;
template<class R> class Interval;
template<class T1, class T2=T1> class traits;
}}

namespace Ariadne { namespace LinearAlgebra {
class MultiIndex;
template<class R> class Vector;
template<class R> class VectorSlice;
template<class R> class Matrix;
template<class R> class MatrixSlice;
template<class R> class Tensor;
template<class R> class LinearProgram;
}}

namespace Ariadne { namespace Combinatoric {
class BinaryWord;
class BinaryTree;  

class LatticeCell;
class LatticeBlock;
class LatticeMaskSet;
class LatticeCellListSet;
class LatticeBlockListSet;
class LatticeTransformation;
class LatticeMultiMap;
class LatticeSystem;
}}

namespace Ariadne { namespace Geometry {
template<class R> class Point;
template<class R> class PointList;

template<class R> class Rectangle;
template<class R> class Parallelotope;
template<class R> class Zonotope;
template<class R> class Simplex;
template<class R> class Polyhedron;
template<class R> class Polytope;
template<class R> class Sphere;
template<class R> class Ellipsoid;

template<class R, template<class> class BS> class ListSet;

typedef enum {
 REGULAR,
 IRREGULAR
} grid_type;

template<class R> class Grid;
template<class R> class RegularGrid;
template<class R> class IrregularGrid;
template<class R> class FiniteGrid;
template<class R> class GridCell;
template<class R> class GridBlock;
template<class R> class GridMaskSet;
template<class R> class GridCellListSet;
template<class R> class GridBlockListSet;

template<class R> class PartitionScheme;
template<class R> class PartitionTree;
template<class R> class PartitionTreeCell;
template<class R> class PartitionTreeSet;

template<class R> class Set;

template<class R> class HybridGridMaskSet;
}}


namespace Ariadne { namespace System {
template<class R> class Map;
template<class R> class AffineMap;
template<class R, template<class> class BS > class AffineMultiMap;
template<class R> class Monomial;
template<class R> class Polynomial;
template<class R> class PolynomialMap;
template<class R> class PolynomialMatrix;
template<class R> class VectorField;
template<class R> class AffineVectorField;

template<class R> class DiscreteLocation;
template<class R> class DiscreteTransition;
template<class R> class HybridAutomaton;
  
}}


namespace Ariadne { namespace Evaluation {
template<class R> class Solver;
template<class R> class Applicator;
template<class R> class Integrator;
template<class R> class HybridEvolver;
}}

namespace Ariadne { namespace Base {
/*! \brief An unsigned integral type used to represent a coordinate in state space. */
typedef unsigned short dimension_type;
/*! \brief An unsigned integral type used to represent the size of a list. */
typedef size_t size_type;
/*! \brief An signed integral type used to represent the position in a list with positive and negative indices. */
typedef int index_type;
/*! \brief The type of a machine byte. */
typedef unsigned char byte_type; 
/*! \brief The type used for a unique identifyer or key. */
typedef size_type id_type;
/*! \brief The type used for a unique identifyer or key. */
typedef Numeric::Rational time_type;

/*! \brief An array of boolean values. */
typedef array<bool> BooleanArray;
/*! \brief An array of unsigned integer values. */
typedef array<size_type> SizeArray;
/*! \brief An array of integer values. */
typedef array<index_type> IndexArray;
/*! \brief An array of integer values. */
typedef LinearAlgebra::MultiIndex multi_index_type;
}}

namespace Ariadne {
using namespace Base;
using namespace Numeric;
}

#endif /* _ARIADNE_DECLARATIONS_H */
