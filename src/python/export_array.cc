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

#include "python/subscripting.h"
#include "python/read_array.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;

template<class X> inline
array<X>*
make_array(const boost::python::object& obj)
{
  array<X>* ary=new array<X>();
  read_array(*ary,obj);
  return ary;
}

template<class T>
void export_array(const char* name) {

  class_< array<T> > array_class(name,init<uint>());
  array_class.def("__init__", make_constructor(&make_array<T>));
  array_class.def(init< array<T> >());
  array_class.def("__len__", &array<T>::size);
  array_class.def("__getitem__", &__getitem__< array<T> >);
  array_class.def("__setitem__", &__setitem__< array<T>, T >);
  array_class.def(self_ns::str(self));    // __self_ns::str__

}

void export_arrays() {
  export_array<bool>("BooleanArray");
  export_array<index_type>("IndexArray");
  export_array<size_type>("SizeArray");
  export_array<Integer>("IntegerArray");
  export_array<Rational>("RationalArray");
  export_array<FloatPy>("FloatArray");
}

