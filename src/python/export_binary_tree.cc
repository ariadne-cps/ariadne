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
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "combinatoric/binary_word.h"
#include "combinatoric/binary_tree.h"


#include "python/typedefs.h"
#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Combinatoric;

#include <boost/python.hpp>
using namespace boost::python;
  
void export_binary_tree() {
  class_<BinaryWord>("BinaryWord",init<>())
    .def(init<BinaryWord>())
    .def("back", &BinaryWord::back)
    .def("set_back", &BinaryWord::set_back)
    .def("pop_back", &BinaryWord::pop_back)
    .def("push_back", &BinaryWord::push_back)
    .def("__len__", &BinaryWord::size)
    .def("__getitem__", &BinaryWord::get)
    .def(self_ns::str(self))    // __self_ns::str__
  ;

  class_<BinaryTree>("BinaryTree",init<>())
    .def(init<uint>())
    .def(init<BooleanArray>())
    .def("__len__", &BinaryTree::size)
    .def(self_ns::str(self))    // __self_ns::str__
  ;
}
