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
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "base/array.h"
#include "base/array_operations.h"
#include "base/utility.h"
#include "base/basic_type.h"

#include <boost/python.hpp>

#include "python/real_typedef.h"
#include "python/python_utilities.h"

using Ariadne::BooleanArray;
using Ariadne::IndexArray;
using Ariadne::SizeArray;

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::self_ns::str;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;
  
void export_array() {
  class_<BooleanArray>("BooleanArray",init<uint>())
    .def(init<BooleanArray>())
    .def("__len__", &BooleanArray::size)
    .def("__getitem__", &get<BooleanArray>)
    .def("__setitem__", &set<BooleanArray>)
    .def(str(self))    // __str__
  ;

  class_<IndexArray>("IndexArray",init<uint>())
    .def(init<IndexArray>())
    .def("__len__", &IndexArray::size)
    .def("__getitem__", &get<IndexArray>)
    .def("__setitem__", &set<IndexArray>)
    .def(str(self))    // __str__
  ;

  class_<SizeArray>("SizeArray",init<uint>())
    .def(init<SizeArray>())
    .def("__len__", &SizeArray::size)
    .def("__getitem__", &get<SizeArray>)
    .def("__setitem__", &set<SizeArray>)
    .def(str(self))    // __str__
  ;
}
