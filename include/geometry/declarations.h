/***************************************************************************
 *            geometry/declarations.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file geometry/declarations.h
 *  \brief Forward declarations of classes in the Geometry module.
 */

#ifndef ARIADNE_GEOMETRY_DECLARATIONS_H
#define ARIADNE_GEOMETRY_DECLARATIONS_H

#include <iosfwd>

namespace Ariadne { 
  
    
    class ExactTag;
    class UniformErrorTag;
    class IntervalTag;

    class EuclideanSpace;
    template<class X> class Point;
    template<class R> class PointList;
    template<class R> class Box;
    template<class R> class BoxListSet;
    template<class E> class BoxExpression;

    template<class R> class Zonotope;
    template<class R> class Simplex;
    template<class R> class Halfspace;
    template<class R> class Polyhedron;
    template<class R> class Polytope;
    template<class R> class Sphere;
    template<class R> class Ellipsoid;

    template<class R> class TaylorSet;
    
    template<class BS> class ListSet;

    template<class R> class Grid;
    template<class R> class FiniteGrid;
    template<class R> class GridCell;
    template<class R> class GridBlock;
    template<class R> class GridCellListSet;
    template<class R> class GridMaskSet;

    template<class R> class IrregularGrid;
    template<class R> class IrregularGridMaskSet;

    template<class R> class PartitionScheme;
    template<class R> class PartitionTree;
    template<class R> class PartitionTreeCell;
    template<class R> class PartitionTreeSet;

    template<class R> class BasicSetInterface;
    template<class R> class DenotableSetInterface;
    template<class R> class SetInterface;

    template<class R> class SetReference;

    template<class R> class CurveInterface;
    template<class R> class ConstraintInterface;
  
    template<class R> class ImageSet;
    template<class R> class ConstraintSet;

    template<class T> class OrbitInterface;
    template<class T, class ES, class RS=ES> class Orbit;

    class HybridSpace;
    template<class S> class HybridSet;
    template<class R> class HybridGrid;
    template<class R> class HybridGridCell;
    template<class R> class HybridGridMaskSet;
    template<class R> class HybridGridCellListSet;
    template<class R> class HybridBox;
    template<class R> class HybridBoxListSet;
    template<class BS> class HybridBasicSet;
    template<class BS> class HybridListSet;

    template<class T, class BS> class TimedSet;

    class basic_set_tag;
    class denotable_set_tag;

    template<class BS> std::string name();


} // namespace Ariadne

#endif /* ARIADNE_GEOMETRY_DECLARATIONS_H */
