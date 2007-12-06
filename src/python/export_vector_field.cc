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

#include "python/float.h"

#include "linear_algebra/vector.h"
#include "system/vector_field.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
class VectorFieldWrapper
  : public VectorFieldInterface<R>, public wrapper< VectorFieldInterface<R> >
{
  typedef typename VectorFieldInterface<R>::F F;
 public: 
  VectorFieldInterface<R>* clone() const { return this->get_override("clone")(); }
  Vector<F> image(const Point<F>&) const { return this->get_override("image")(); }
  dimension_type dimension() const { return this->get_override("dimension")(); }
  smoothness_type smoothness() const { return this->get_override("smoothness")(); }
  std::string name() const { return this->get_override("name")(); }
};

template<class R>
void export_vector_field() 
{
  class_<VectorFieldWrapper<R>, boost::noncopyable>("VectorFieldInterface")
    .def("dimension", pure_virtual(&VectorFieldInterface<R>::dimension))
    .def("smoothness", pure_virtual(&VectorFieldInterface<R>::smoothness))
    .def("name", pure_virtual(&VectorFieldInterface<R>::name))
  ;
}

template void export_vector_field<FloatPy>();
