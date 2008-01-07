/***************************************************************************
 *            python/export_postscript_output.cc
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *
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


#include "python/float.h"

#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/box_list_set.h"
#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"

#include "geometry/orbit.h"
#include "output/epsstream.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class S> inline void write(epsfstream& eps, const S& s) { eps << s; }
template<class R> inline void epsfstream_open(epsfstream& eps, const Ariadne::Geometry::Box<R>& bbox, int ix, int iy) { eps.open("Ariadne",bbox,ix,iy); }
template<class R> inline void epsfstream_open_with_defaults(epsfstream& eps, const Ariadne::Geometry::Box<R>& bbox) { eps.open("Ariadne",bbox); }
inline void epsfstream_close(epsfstream& eps) { eps.close(); }

void export_postscript_output()
{

  class_<PlanarProjectionMap>("PlanarProjectionMap",init<dimension_type, dimension_type, dimension_type>())
    .def(self_ns::str(self))
  ;
    
  class_<epsfstream, boost::noncopyable>("EpsPlot",init<>())
    .def("open",(void(epsfstream::*)(const char* fn,const Box<FloatPy>&))&epsfstream::open<FloatPy>)
    .def("open",(void(epsfstream::*)(const char* fn,const Box<FloatPy>&,uint,uint))&epsfstream::open<FloatPy>)
    .def("open",(void(epsfstream::*)(const char* fn,const Box<FloatPy>&,const PlanarProjectionMap&))&epsfstream::open<FloatPy>)
    .def("open",(void(epsfstream::*)(const char*,const Rectangle2d&,const PlanarProjectionMap&))&epsfstream::open)
    .def("open",&epsfstream_open_with_defaults<FloatPy>)
    .def("close",&epsfstream_close)
    .def("set_pen_colour",(void(epsfstream::*)(const char*))&epsfstream::set_pen_colour)
    .def("set_fill_colour",(void(epsfstream::*)(const char*))&epsfstream::set_fill_colour)
    .def("set_pen_colour",(void(epsfstream::*)(const Colour&))&epsfstream::set_pen_colour)
    .def("set_fill_colour",(void(epsfstream::*)(const Colour&))&epsfstream::set_fill_colour)
    .def("set_line_style",(void(epsfstream::*)(bool))&epsfstream::set_line_style)
    .def("set_fill_style",(void(epsfstream::*)(bool))&epsfstream::set_fill_style)
    .def("write",&write< Box<FloatPy> >)
    .def("write",&write< GridCell<FloatPy> >)
    .def("write",&write< GridBlock<FloatPy> >)
    .def("write",&write< GridCellListSet<FloatPy> >)
    .def("write",&write< GridMaskSet<FloatPy> >)
    .def("write",&write< PartitionTreeSet<FloatPy> >)
    .def("write",&write< FiniteGrid<FloatPy> >)
    .def("write",&write< PartitionTree<FloatPy> >)
    .def("write",&write< Zonotope<FloatPy> >)
    .def("write",&write< Polytope<FloatPy> >)
    .def("write",&write< Polyhedron<FloatPy> >)
    .def("write",&write< ListSet< Polytope<FloatPy> > >)
    .def("write",&write< ListSet< Zonotope<FloatPy> > >)
    .def("write",&write< BoxListSet<FloatPy> >)
    .def("write",&write< RectangularSet<FloatPy> >)
    ;
}
