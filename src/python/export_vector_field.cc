/***************************************************************************
 *            python/export_vector_field.cc
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

#include "system/vector_field.h"

#include "python/typedefs.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

struct RVectorField : RVectorFieldBase, wrapper<RVectorFieldBase>
{
  dimension_type dimension() const { return this->get_override("dimension")(); }
  std::string name() const { return this->get_override("name")(); }
};

void export_vector_field() {
  class_<RVectorField, boost::noncopyable>("VectorField")
    .def("dimension", pure_virtual(&RVectorFieldBase::dimension))
    .def("name", pure_virtual(&RVectorFieldBase::name))
  ;
}
