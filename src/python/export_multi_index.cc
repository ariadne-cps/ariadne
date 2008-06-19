/***************************************************************************
 *            python/export_multi_index.cc
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

#include "base/array.h"
#include "function/multi_index.h"

#include "python/operators.h"
#include "python/read_array.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;

MultiIndex*
make_multi_index(const boost::python::object& obj) 
{
  array<uint> ary;
  read_array(ary,obj);
  return new MultiIndex(ary.size(),ary.begin());
}


void export_multi_index()
{
  
  class_<MultiIndex> multi_index_class("MultiIndex",no_init);
  multi_index_class.def("__init__", make_constructor(&make_multi_index));
  multi_index_class.def(init<size_type>());
  multi_index_class.def(init<size_type,smoothness_type>());
  multi_index_class.def(init<MultiIndex>());
  multi_index_class.def("degree",&MultiIndex::degree);
  multi_index_class.def("increment",&MultiIndex::operator++,return_value_policy<reference_existing_object>());
  multi_index_class.def("__getitem__",&MultiIndex::operator[],return_value_policy<copy_const_reference>());
  multi_index_class.def("__setitem__",&MultiIndex::set_index);
  multi_index_class.def(self_ns::str(self));
}

