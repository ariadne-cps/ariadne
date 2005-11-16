/***************************************************************************
 *            python/state.cc
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

#include <iostream>
#include "numerical_type.h"
#include "state.h"

#include <boost/python.hpp>

#include "real_typedef.h"
#include "container_utilities.h"

typedef Ariadne::Geometry::State<Real> RState;

using boost::python::class_;
using boost::python::init; 
using boost::python::self;
using boost::python::def;
using boost::python::self_ns::str;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;

void export_state() {
  class_<RState>("State",init<int>())
    .def(init<int,Real>())
    .def(init<RState>())
    .def("dimension", &RState::dimension)
    .def("__len__", &RState::dimension)
    .def("__getitem__", &RState::get)
    .def("__setitem__", &RState::set)
    .def("__eq__", &RState::operator==)
    .def("__ne__", &RState::operator!=)
    .def(str(self))    // __str__
  ;
}
