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
template<typename T, size_t N=0> class array;
template<typename T> class sequence;
}}

namespace Ariadne { namespace Numeric {
template<typename R> class Interval;
template<typename T> class numerical_traits;
}}

namespace Ariadne { namespace LinearAlgebra {
class MultiIndex;
template<typename R> class Vector;
template<typename R> class Matrix;
template<typename R> class Tensor;
template<typename R> class LinearProgram;
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
template<typename R> class Point;
template<typename R> class PointList;

template<typename R> class Rectangle;
template<typename R> class Parallelotope;
template<typename R> class Zonotope;
template<typename R> class Simplex;
template<typename R> class Polyhedron;
template<typename R> class Polytope;
template<typename R> class Sphere;
template<typename R> class Ellipsoid;

template<typename R, template<typename> class BS> class ListSet;

typedef enum {
 REGULAR,
 IRREGULAR
} grid_type;

template<typename R> class Grid;
template<typename R> class RegularGrid;
template<typename R> class IrregularGrid;
template<typename R> class FiniteGrid;
template<typename R> class GridCell;
template<typename R> class GridBlock;
template<typename R> class GridMaskSet;
template<typename R> class GridCellListSet;
template<typename R> class GridBlockListSet;

template<typename R> class PartitionScheme;
template<typename R> class PartitionTree;
template<typename R> class PartitionTreeCell;
template<typename R> class PartitionTreeSet;
}}


namespace Ariadne { namespace System {
template <typename R> class Map;
template <typename R> class AffineMap;
template <typename R, template<typename> class BS > class AffineMultiMap;
template <typename R> class Monomial;
template <typename R> class Polynomial;
template <typename R> class PolynomialMap;
template <typename R> class PolynomialMatrix;
template <typename R> class VectorField;
template <typename R> class AffineVectorField;

template <typename R> class HenonMap;
template <typename R> class LorenzSystem;
}}


namespace Ariadne { namespace Evaluation {
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
