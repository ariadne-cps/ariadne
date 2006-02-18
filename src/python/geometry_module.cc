/***************************************************************************
 *            python/geometry_module.cc
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

#include "base/binary_word.h"
#include "base/binary_tree.h"

#include <boost/python.hpp>

#include "python/real_typedef.h"
#include "python/python_utilities.h"


using Ariadne::BooleanArray;
using Ariadne::BinaryWord;
using Ariadne::BinaryTree;
using Ariadne::BinaryWord;
using Ariadne::SizeArray;
using Ariadne::SubdivisionSequence;

void export_point();
void export_rectangle();
void export_parallelopiped();
void export_simplex();
void export_polyhedron();
void export_list_set();
void export_grid();
void export_partition_tree();

void export_postscript_output();

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::self_ns::str;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;
  
BOOST_PYTHON_MODULE(geometry)
{
  class_<BooleanArray>("BooleanArray")
    .def(init<int>())
    .def("__len__", &BooleanArray::size)
    .def("__getitem__", &get<BooleanArray>)
    .def("__setitem__", &set<BooleanArray>)
    .def(str(self))    // __str__
  ;

  class_<BinaryTree>("BinaryTree")
    .def(init<BooleanArray>())
    .def("__len__", &BinaryTree::size)
    .def(str(self))    // __str__
  ;

  class_<BinaryWord>("BinaryWord")
    .def(init<BinaryWord>())
    .def("__len__", &BinaryWord::size)
    .def("__getitem__", &BinaryWord::get)
    .def(str(self))    // __str__
  ;

  class_<SizeArray>("SizeArray")
    .def(init<SizeArray>())
    .def("__len__", &SizeArray::size)
    .def("__getitem__", &get<SizeArray>)
    .def(str(self))    // __str__
  ;

  class_<SubdivisionSequence>("SubdivisionSequence")
    .def("__len__", &BinaryTree::size)
    .def("__getitem__", &SubdivisionSequence::get)
    .def(str(self))    // __str__
  ;

  export_point();
  export_rectangle();
  export_parallelopiped();
  export_simplex();
  export_polyhedron();
  export_list_set();
  export_grid();
  export_partition_tree();
  
  export_postscript_output();

}
