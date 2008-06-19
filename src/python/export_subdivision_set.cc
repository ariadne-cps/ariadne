/***************************************************************************
 *            python/export_subdivision_set.cc
 *
 *  Copyright  2007  Pieter Collins
 *
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



#include "combinatoric/subdivision_tree_set.h"

using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

void export_subdivision_set() {

  class_<SubdivisionCell>("SubdivisionCell",init<uint,uint,uint>())
    .def("__eq__",&SubdivisionCell::operator==)
    .def("dimension", &SubdivisionCell::dimension)
    .def("depth", &SubdivisionCell::depth)
    .def("subdivisions", &SubdivisionCell::subdivisions)
    .def("box", &SubdivisionCell::box)
    .def("split", &SubdivisionCell::split)
    .def(self_ns::str(self))
    ;

  class_<SubdivisionBox>("SubdivisionBox",init<SubdivisionCell>())
    .def(self_ns::str(self))
    .def("dimension", &SubdivisionBox::dimension)
    .def("depth", &SubdivisionBox::depth)
    .def("lower_bound", &SubdivisionBox::lower_bound)
    .def("upper_bound", &SubdivisionBox::upper_bound)
    ;
 
  class_<SubdivisionCellListSet>("SubdivisionCellListSet",init<>())
    .def("dimension", &SubdivisionCellListSet::dimension)
    .def("size", &SubdivisionCellListSet::size)
    .def("pop", &SubdivisionCellListSet::pop)
    .def("adjoin", (void(SubdivisionCellListSet::*)(const SubdivisionCell&))(&SubdivisionCellListSet::adjoin))
    .def("unique_sort", &SubdivisionCellListSet::unique_sort)
    .def("__iter__", iterator<SubdivisionCellListSet>())
    .def(self_ns::str(self))
    ;

  class_<SubdivisionMaskSet>("SubdivisionMaskSet",init<uint,uint>())
    .def("dimension", &SubdivisionMaskSet::dimension)
    .def("depth", &SubdivisionMaskSet::dimension)
    .def("adjoin", (void(SubdivisionMaskSet::*)(const SubdivisionCell&))(&SubdivisionMaskSet::adjoin))
    .def("__iter__", iterator<SubdivisionMaskSet>())
    .def(self_ns::str(self))
    ;

  def("hull",&hull);

}
