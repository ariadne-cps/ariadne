/***************************************************************************
 *            python/export_grid_multimap.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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

#include "real_typedef.h"

#include "linear_algebra/vector.h"

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "system/grid_multimap.h"


#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Combinatoric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;

#include <boost/python.hpp>
using namespace boost::python;


template<class R>
void export_grid_multimap() 
{
  typedef Interval<R> I;

  class_< GridMultiMap<R> >("GridMultiMap",init<const Grid<R>&,const Grid<R>&>())
    .def(init<const Grid<R>&,const Grid<R>&,const LatticeMultiMap&>())
    .def("argument_dimension", &GridMultiMap<R>::argument_dimension)
    .def("result_dimension", &GridMultiMap<R>::result_dimension)
    .def("adjoin_to_image", (void(GridMultiMap<R>::*)(const GridCell<R>&,const GridCell<R>&))(&GridMultiMap<R>::adjoin_to_image))
    .def("adjoin_to_image", (void(GridMultiMap<R>::*)(const GridCell<R>&,const GridCellListSet<R>&))(&GridMultiMap<R>::adjoin_to_image))
    .def("image", (GridCellListSet<R>(GridMultiMap<R>::*)(const GridCell<R>&)const)(&GridMultiMap<R>::image))
    .def("image", (GridCellListSet<R>(GridMultiMap<R>::*)(const GridMaskSet<R>&)const)(&GridMultiMap<R>::image))
    .def("__call__", (GridCellListSet<R>(GridMultiMap<R>::*)(const GridCell<R>&)const)(&GridMultiMap<R>::operator()))
    .def("__call__", (GridCellListSet<R>(GridMultiMap<R>::*)(const GridBlock<R>&)const)(&GridMultiMap<R>::operator()))
    .def("__call__", (GridCellListSet<R>(GridMultiMap<R>::*)(const GridCellListSet<R>&)const)(&GridMultiMap<R>::operator()))
    .def("__call__", (GridCellListSet<R>(GridMultiMap<R>::*)(const GridMaskSet<R>&)const)(&GridMultiMap<R>::operator()))
    .def("strong_preimage", (GridCellListSet<R>(GridMultiMap<R>::*)(const GridMaskSet<R>&)const)(&GridMultiMap<R>::strong_preimage))
    .def("weak_preimage", (GridCellListSet<R>(GridMultiMap<R>::*)(const GridMaskSet<R>&)const)(&GridMultiMap<R>::weak_preimage))
    .def("inverse", &GridMultiMap<R>::inverse)
    .def(self_ns::str(self))    // __str__
    ;

}

template void export_grid_multimap<Real>();
