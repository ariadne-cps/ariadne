/***************************************************************************
 *            python/export_array.cc
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
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "base/array.h"
#include "combinatoric/array_operations.h"
#include "base/stlio.h"

#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Base;

#include <boost/python.hpp>
using namespace boost::python;
  
void export_array() {
  class_<BooleanArray>("BooleanArray",init<uint>())
    .def(init<BooleanArray>())
    .def("__len__", &BooleanArray::size)
    .def("__getitem__", &get_item<BooleanArray>)
    .def("__setitem__", &set_item<BooleanArray>)
    .def(self_ns::str(self))    // __self_ns::str__
  ;

  class_<IndexArray>("IndexArray",init<uint>())
    .def(init<IndexArray>())
    .def("__len__", &IndexArray::size)
    .def("__getitem__", &get_item<IndexArray>)
    .def("__setitem__", &set_item<IndexArray>)
    .def(self_ns::str(self))    // __self_ns::str__
  ;

  class_<SizeArray>("SizeArray",init<uint>())
    .def(init<SizeArray>())
    .def("__len__", &SizeArray::size)
    .def("__getitem__", &get_item<SizeArray>)
    .def("__setitem__", &set_item<SizeArray>)
    .def(self_ns::str(self))    // __self_ns::str__
  ;

}
