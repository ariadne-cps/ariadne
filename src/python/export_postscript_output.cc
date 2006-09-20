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
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"

#include "output/epsfstream.h"

#include "python/typedefs.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

typedef Postscript::epsfstream epsfstream;

inline void write_rectangle(epsfstream& eps, const RRectangle& r) { eps << r; }
inline void write_grid_rectangle(epsfstream& eps, const RGridRectangle& r) { eps << RRectangle(r); }
inline void write_grid_cell(epsfstream& eps, const RGridCell& r) { eps << RRectangle(r); }
inline void write_parallelotope(epsfstream& eps, const RParallelotope& p) { eps << p; }
inline void write_zonotope(epsfstream& eps, const RZonotope& z) { eps << z; }
inline void write_polytope(epsfstream& eps, const RPolytope& p) { eps << p; }
inline void write_rectangle_list_set(epsfstream& eps, const RRectangleListSet& r) { eps << r; }
inline void write_parallelotope_list_set(epsfstream& eps, const RParallelotopeListSet& s) { eps << s; }
inline void write_zonotope_list_set(epsfstream& eps, const RZonotopeListSet& s) { eps << s; }
inline void write_polytope_list_set(epsfstream& eps, const RPolytopeListSet& s) { eps << s; }
inline void write_grid_mask_set(epsfstream& eps, const RGridMaskSet& s) { eps << s; }
inline void write_grid_cell_list_set(epsfstream& eps, const RGridCellListSet& s) { eps << s; }
inline void write_grid_rectangle_list_set(epsfstream& eps, const RGridRectangleListSet& s) { eps << s; }
inline void write_partition_tree(epsfstream& eps, const RPartitionTree& s) { eps << s; }
inline void write_partition_tree_set(epsfstream& eps, const RPartitionTreeSet& s) { eps << s; }
inline void epsfstream_open(epsfstream& eps, const Ariadne::Geometry::Rectangle<Real>& bbox) { eps.open("Ariadne",bbox); }
inline void epsfstream_close(epsfstream& eps) { eps.close(); }

void export_postscript_output()
{
  class_<epsfstream>("EpsPlot",init<const char*,RRectangle>())
    .def(init<const char*,RRectangle, const unsigned int&, const unsigned int&>())
    .def(init<const char*,RRectangle, const unsigned int&, const unsigned int&, const char*, const char*>())
    .def("open",&epsfstream_open)
    .def("close",&epsfstream_close)
    .def("set_pen_colour",&epsfstream::set_pen_colour)
    .def("set_fill_colour",&epsfstream::set_fill_colour)
    .def("set_line_style",&epsfstream::set_line_style)
    .def("set_fill_style",&epsfstream::set_fill_style)
    .def("write",&write_rectangle)
    .def("write",&write_grid_rectangle)
    .def("write",&write_grid_cell)
    .def("write",&write_parallelotope)
    .def("write",&write_zonotope)
    .def("write",&write_parallelotope)
    .def("write",&write_polytope)
    .def("write",&write_rectangle_list_set)
    .def("write",&write_parallelotope_list_set)
    .def("write",&write_zonotope_list_set)
    .def("write",&write_polytope_list_set)
    .def("write",&write_grid_mask_set)
    .def("write",&write_grid_cell_list_set)
    .def("write",&write_grid_rectangle_list_set)
    .def("write",&write_partition_tree)
    .def("write",&write_partition_tree_set)
  ;
  
}
