/***************************************************************************
 *            python/export_grid.cc
 *
 *  21 October 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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

#include "numerical_type.h"
#include "linear_algebra.h"
#include "grid_set.h"

#include "parallelopiped.h"

#include <boost/python.hpp>

#include "real_typedef.h"

using Ariadne::BooleanArray;
using Ariadne::IndexArray;
using Ariadne::Geometry::IntegerRectangle;

typedef Ariadne::Geometry::Point<Real> RPoint;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;
typedef Ariadne::Geometry::Parallelopiped<Real> RParallelopiped;
typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Rectangle> RRectangleListSet;

typedef Ariadne::Geometry::Grid<Real> RGrid;
typedef Ariadne::Geometry::FiniteGrid<Real> RFiniteGrid;
typedef Ariadne::Geometry::GridCell<Real> RGridCell;
typedef Ariadne::Geometry::GridCellListSet<Real> RGridCellListSet;
typedef Ariadne::Geometry::GridRectangle<Real> RGridRectangle;
typedef Ariadne::Geometry::GridRectangleListSet<Real> RGridRectangleListSet;
typedef Ariadne::Geometry::GridMaskSet<Real> RGridMaskSet;

using Ariadne::Geometry::regular_intersection;
using Ariadne::Geometry::interiors_intersect;
using Ariadne::Geometry::disjoint;
using Ariadne::Geometry::inner_subset;
using Ariadne::Geometry::subset;

using boost::python::class_;
using boost::python::init; 
using boost::python::self;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;
using boost::python::def;
using boost::python::self_ns::str;

inline RGridCellListSet over_approximation_parallelopiped(const RFiniteGrid& g, const RParallelopiped& p) {
  return Ariadne::Geometry::over_approximation(g,p);
}

void export_grid() {
  class_<RFiniteGrid>("FiniteGrid",init<RFiniteGrid>())
    .def(init<RRectangleListSet>())
    .def(init<RRectangle,uint>())
    .def("dimension", &RFiniteGrid::dimension)
    .def("subdivision_coordinate", &RFiniteGrid::subdivision_coordinate)
    .def("subdivision_interval", &RFiniteGrid::subdivision_interval)
    .def("subdivision_index", &RFiniteGrid::subdivision_index)
    .def("subdivision_lower_index", &RFiniteGrid::subdivision_lower_index)
    .def("subdivision_upper_index", &RFiniteGrid::subdivision_upper_index)
    .def(str(self))    // __str__
    ;

  class_<RGridCell>("RGridCell",init<RFiniteGrid,IndexArray>())
    .def("dimension", &RGridCell::dimension)
    .def(str(self))    // __str__
    ;
  
  class_<RGridCellListSet>("GridCellListSet",init<RFiniteGrid>())
    .def(init<RGridMaskSet>())
    .def(init<RGridCellListSet>())
    .def("dimension", &RGridCellListSet::dimension)
//    .def("adjoin", &RGridRectangleListSet::adjoin)
    .def("__len__", &RGridCellListSet::size)
    .def(str(self))    // __str__
    ;
    
  class_<RGridRectangle>("RGridCell",init<RFiniteGrid,IndexArray,IndexArray>())
    .def("dimension", &RGridCell::dimension)
    .def(str(self))    // __str__
    ;
  
  
  class_<RGridRectangleListSet>("GridRectangleListSet",init<RRectangleListSet>())
    .def(init<RGridRectangleListSet>())
    .def("dimension", &RGridRectangleListSet::dimension)
//    .def("adjoin", &RGridRectangleListSet::adjoin)
    .def("__len__", &RGridRectangleListSet::size)
    .def(str(self))    // __str__
    ;
    
  class_<RGridMaskSet>("GridMaskSet",init<RFiniteGrid,IntegerRectangle>())
    .def("dimension", &RGridMaskSet::dimension)
//    .def("adjoin", &RGridMaskSet::adjoin)
    .def("__len__", &RGridMaskSet::size)
    .def(str(self))    // __str__
    ;
    
    
  def("over_approximation",&over_approximation_parallelopiped);
}
