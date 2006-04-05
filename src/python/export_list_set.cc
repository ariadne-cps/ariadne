/***************************************************************************
 *            python/export_list_set.cc
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

#include "base/numerical_type.h"

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"

#include <boost/python.hpp>

#include "python/typedefs.h"
#include "python/python_utilities.h"

using Ariadne::Geometry::regular_intersection;
using Ariadne::Geometry::interiors_intersect;
using Ariadne::Geometry::disjoint;
using Ariadne::Geometry::inner_subset;
using Ariadne::Geometry::subset;

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::iterator;
using boost::python::def;
using boost::python::self_ns::str;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;

void export_list_set() {
  typedef bool (*RectLSBinPred) (const RRectangleListSet&, const RRectangleListSet&);
  typedef RRectangleListSet (*RectLSBinFun) (const RRectangleListSet&, const RRectangleListSet&);

  def("regular_intersection", RectLSBinFun(&regular_intersection));
  def("interiors_intersect", RectLSBinPred(&interiors_intersect));
  def("disjoint", RectLSBinPred(&disjoint));
  def("inner_subset", RectLSBinPred(&inner_subset));
  def("subset", RectLSBinPred(&subset));

  class_<RRectangleListSet>("RectangleListSet",init<int>())
    .def(init<RRectangle>())
    .def(init<RRectangleListSet>())
    .def(init<RGridRectangleListSet>())
    .def(init<RGridCellListSet>())
    .def(init<RGridMaskSet>())
    .def(init<RPartitionTreeSet>())
    .def("dimension", &RRectangleListSet::dimension)
    .def("push_back", &RRectangleListSet::push_back)
    .def("size", &RRectangleListSet::size)
    .def("__len__", &RRectangleListSet::size)
    .def("__getitem__", &RRectangleListSet::get, return_value_policy<copy_const_reference>())
    .def("__setitem__", &RRectangleListSet::set)
    .def("__iter__", iterator<RRectangleListSet>())
    .def(str(self))    // __str__
  ;
  
  class_<RParallelotopeListSet>("ParallelotopeListSet",init<int>())
    .def(init<RParallelotope>())
    .def(init<RParallelotopeListSet>())
    .def("dimension", &RParallelotopeListSet::dimension)
    .def("push_back", &RParallelotopeListSet::push_back)
    .def("size", &RParallelotopeListSet::size)
    .def("__len__", &RParallelotopeListSet::size)
    .def("__getitem__", &RParallelotopeListSet::get, return_value_policy<copy_const_reference>())
    .def("__setitem__", &RParallelotopeListSet::set)
    .def("__iter__", iterator<RParallelotopeListSet>())
    .def(str(self))    // __str__
  ;

}
