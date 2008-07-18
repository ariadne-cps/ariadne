/***************************************************************************
 *            python/export_slice.cc
 *
 *  Copyright  2008  Pieter Collins
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

#include "linear_algebra/slice.h"

#include "python/operators.h"
#include "python/read_array.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;


void export_slice()
{
  
  class_<Slice> slice_class("Slice",init<uint,uint,uint>());
  slice_class.def("start",&Slice::start);
  slice_class.def("stop",&Slice::stop);
  slice_class.def("size",&Slice::size);
  slice_class.def("stride",&Slice::stride);

  def("slice",&slice);
  def("Range",&Ariadne::range);
}

