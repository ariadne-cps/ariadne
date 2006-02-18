/***************************************************************************
 *            python/export_postscript_output.cc
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
 
#include "geometry/rectangle.h"
#include "geometry/parallelopiped.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"

#include "utility/epsfstream.h"

#include "python/real_typedef.h"

typedef Ariadne::Geometry::Rectangle<Real> RRectangle;
typedef Ariadne::Geometry::GridRectangle<Real> RGridRectangle;
typedef Ariadne::Geometry::GridCell<Real> RGridCell;
typedef Ariadne::Geometry::Parallelopiped<Real> RParallelopiped;
typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Rectangle> RRectangleListSet;
typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Parallelopiped> RParallelopipedListSet;
typedef Ariadne::Geometry::GridMaskSet<Real> RGridMaskSet;
typedef Ariadne::Geometry::GridCellListSet<Real> RGridCellListSet;
typedef Ariadne::Geometry::GridRectangleListSet<Real> RGridRectangleListSet;
typedef Ariadne::Geometry::PartitionTreeSet<Real> RPartitionTreeSet;

using Ariadne::Postscript::epsfstream;

inline void write_rectangle(epsfstream& eps, const RRectangle& r) { eps << r; }
inline void write_grid_rectangle(epsfstream& eps, const RGridRectangle& r) { eps << RRectangle(r); }
inline void write_grid_cell(epsfstream& eps, const RGridCell& r) { eps << RRectangle(r); }
inline void write_parallelopiped(epsfstream& eps, const RParallelopiped& p) { eps << p; }
inline void write_rectangle_list_set(epsfstream& eps, const RRectangleListSet& r) { eps << r; }
inline void write_parallelopiped_list_set(epsfstream& eps, const RParallelopipedListSet& s) { eps << s; }
inline void write_grid_mask_set(epsfstream& eps, const RGridMaskSet& s) { eps << s; }
inline void write_grid_cell_list_set(epsfstream& eps, const RGridCellListSet& s) { eps << s; }
inline void write_grid_rectangle_list_set(epsfstream& eps, const RGridRectangleListSet& s) { eps << s; }
inline void write_partition_tree_set(epsfstream& eps, const RPartitionTreeSet& s) { eps << s; }
inline void epsfstream_open(epsfstream& eps, const Ariadne::Geometry::Rectangle<Real>& bbox) { eps.open(bbox); }

#include <boost/python.hpp>

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;
using boost::python::def;
using boost::python::self_ns::str;

void export_postscript_output()
{
  class_<epsfstream>("EpsPlot",init<const char*,RRectangle>())
    .def("open",&epsfstream_open)
    .def("close",&epsfstream::close)
    .def("set_pen_colour",&epsfstream::set_pen_colour)
    .def("set_fill_colour",&epsfstream::set_fill_colour)
    .def("set_line_style",&epsfstream::set_line_style)
    .def("write",&write_rectangle)
    .def("write",&write_grid_rectangle)
    .def("write",&write_grid_cell)
    .def("write",&write_parallelopiped)
    .def("write",&write_rectangle_list_set)
    .def("write",&write_parallelopiped_list_set)
    .def("write",&write_grid_mask_set)
    .def("write",&write_grid_cell_list_set)
    .def("write",&write_grid_rectangle_list_set)
    .def("write",&write_partition_tree_set)
  ;
  
}
