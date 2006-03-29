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
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "evaluation/vector_field.h"

#include <boost/python.hpp>

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::bases;
using boost::python::wrapper;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;
using boost::python::pure_virtual;

#include "python/real_typedef.h"
  
using Ariadne::dimension_type;
using Ariadne::Interval;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;

typedef Ariadne::Evaluation::VectorField<Real> RVectorFieldBase;

struct RVectorField : RVectorFieldBase, wrapper<RVectorFieldBase>
{
  dimension_type dimension() const { return this->get_override("dimension")(); }
};

void export_vector_field() {
  class_<RVectorField, boost::noncopyable>("VectorField")
    .def("dimension", pure_virtual(&RVectorFieldBase::dimension))
  ;
}
