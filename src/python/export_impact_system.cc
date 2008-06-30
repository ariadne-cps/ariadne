/***************************************************************************
 *            python/export_impact_system.cc
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
#include "linear_algebra/matrix.h"
#include "differentiation/taylor_derivative.h"
#include "function/function_interface.h"
#include "geometry/point.h"
#include "geometry/box.h"

#include "system/impact_system.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;


template<class R>
void export_impact_system() 
{
  typedef typename traits<R>::arithmetic_type A;


  class_< ImpactSystem<R> >("ImpactSystem",
                            init<const FunctionInterface<R>&, const FunctionInterface<R>&, const FunctionInterface<R>&>())
    .def("state_space", &ImpactSystem<R>::state_space)
    .def("vector_field", &ImpactSystem<R>::vector_field,return_value_policy<copy_const_reference>())
    .def("guard_condition", &ImpactSystem<R>::guard_condition,return_value_policy<copy_const_reference>())
    .def("impact_map", &ImpactSystem<R>::impact_map,return_value_policy<copy_const_reference>())
    .def(self_ns::str(self))
  ;

}

template void export_impact_system<FloatPy>();
