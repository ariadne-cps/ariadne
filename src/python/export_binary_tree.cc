/***************************************************************************
 *            python/export_binary_tree.cc
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

#include "base/binary_word.h"
#include "base/binary_tree.h"
#include "base/basic_type.h"

#include <boost/python.hpp>

#include "python/real_typedef.h"
#include "python/python_utilities.h"

using Ariadne::BooleanArray;
using Ariadne::Base::BinaryWord;
using Ariadne::Base::BinaryTree;

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::self_ns::str;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;
  
void export_binary_tree() {
  class_<BinaryWord>("BinaryWord",init<>())
    .def(init<BinaryWord>())
    .def("back", &BinaryWord::back)
    .def("set_back", &BinaryWord::set_back)
    .def("pop_back", &BinaryWord::pop_back)
    .def("push_back", &BinaryWord::push_back)
    .def("__len__", &BinaryWord::size)
    .def("__getitem__", &BinaryWord::get)
    .def(str(self))    // __str__
  ;

  class_<BinaryTree>("BinaryTree",init<>())
    .def(init<uint>())
    .def(init<BooleanArray>())
    .def("__len__", &BinaryTree::size)
    .def(str(self))    // __str__
  ;
}
