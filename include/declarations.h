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
 
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace Ariadne { namespace LinearAlgebra {
using boost::numeric::ublas::vector;
using boost::numeric::ublas::matrix;

template<typename R> class interval_vector;
template<typename R> class interval_matrix;
}}

namespace Ariadne { namespace Geometry {
template<typename R> class Point;

template<typename R> class Rectangle;
template<typename R> class Parallelotope;
template<typename R> class Zonotope;
template<typename R> class Simplex;
template<typename R> class Polyhedron;

template<typename R, template<typename> class BS> class ListSet;

template<typename R> class Grid;
template<typename R> class FiniteGrid;
template<typename R> class GridCell;
template<typename R> class GridRectangle;
template<typename R> class GridMaskSet;
template<typename R> class GridCellListSet;
template<typename R> class GridRectangleListSet;

template<typename R> class PartitionScheme;
template<typename R> class PartitionTree;
template<typename R> class PartitionTreeCell;
template<typename R> class PartitionTreeSet;
}}


namespace Ariadne { namespace Evaluation {
template <typename R> class Map;
template <typename R> class AffineMap;
template <typename R> class VectorField;
template <typename R> class AffineVectorField;

template <typename R> class HenonMap;
template <typename R> class LorenzSystem;
}}

#endif /* _ARIADNE_DECLARATIONS_H */
