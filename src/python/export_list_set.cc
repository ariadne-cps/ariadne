/***************************************************************************
 *            python/export_list_set.cc
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


#include "linear_algebra/matrix.h"
#include "geometry/box.h"
#include "geometry/zonotope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"

#include "python/utilities.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;


template<class R>
void export_list_set() 
{
  typedef Interval<R> I;
  
  typedef Box<R> RBox;
  typedef Polytope<R> RPolytope;
  typedef Zonotope<R> RZonotope;

  typedef BoxListSet<R> RBoxListSet;
  typedef ListSet<RPolytope> RPolytopeListSet;
  typedef ListSet<RZonotope> RZonotopeListSet;

  typedef GridCellListSet<R> RGridCellListSet;
  typedef GridMaskSet<R> RGridMaskSet;
  typedef PartitionTreeSet<R> RPartitionTreeSet;
  

  
  class_<RZonotopeListSet>("ZonotopeListSet",init<int>())
    .def(init<RZonotope>())
    .def(init<RZonotopeListSet>())
    .def("dimension", &RZonotopeListSet::dimension)
    .def("pop", (RZonotope(RZonotopeListSet::*)())(&RZonotopeListSet::pop))
    .def("adjoin", (void(RZonotopeListSet::*)(const RZonotope&))(&RZonotopeListSet::adjoin))
    .def("adjoin", (void(RZonotopeListSet::*)(const RZonotopeListSet&))(&RZonotopeListSet::adjoin))
    .def("size", &RZonotopeListSet::size)
    .def("__len__", &RZonotopeListSet::size)
    .def("__getitem__", &__getitem__<RZonotopeListSet>)
    .def("__iter__", iterator<RZonotopeListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;


}

template void export_list_set<FloatPy>();
