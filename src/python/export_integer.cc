/***************************************************************************
 *            python/export_integer.cc
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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


#include "python/python_utilities.h"

#include "numeric/integer.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;

#include <boost/python.hpp>
using namespace boost::python;

void export_integer() {
  class_<Integer>("Integer")
    .def(init<int>())
    .def(init<Integer>())
    .def("__neg__", &neg<Integer,Integer>)
    .def("__add__", &add<Integer,Integer,Integer>)
    .def("__add__", &add<Integer,Integer,int>)
    .def("__radd__", &add<Integer,Integer,int>)
    .def("__sub__", &sub<Integer,Integer,Integer>)
    .def("__sub__", &sub<Integer,Integer,int>)
    .def("__rsub__", &rsub<Integer,Integer,int>)
    .def("__mul__", &mul<Integer,Integer,Integer>)
    .def("__mul__", &mul<Integer,Integer,int>)
    .def("__rmul__", &mul<Integer,Integer,int>)
    .def("__eq__", &eq<bool,Integer,int>)
    .def("__eq__", &eq<bool,Integer,Integer>)
    .def("__ne__", &ne<bool,Integer,int>)
    .def("__ne__", &ne<bool,Integer,Integer>)
    .def("__lt__", &lt<bool,Integer,int>)
    .def("__lt__", &lt<bool,Integer,Integer>)
    .def("__gt__", &gt<bool,Integer,int>)
    .def("__gt__", &gt<bool,Integer,Integer>)
    .def("__le__", &le<bool,Integer,int>)
    .def("__le__", &le<bool,Integer,Integer>)
    .def("__ge__", &ge<bool,Integer,int>)
    .def("__ge__", &ge<bool,Integer,Integer>)
    .def(self_ns::str(self))
  ;

  
}
