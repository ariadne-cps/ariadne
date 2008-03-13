/***************************************************************************
 *            python/export_evolution_parameters.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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

#include "geometry/box.h"
#include "geometry/grid.h"
#include "evaluation/evolution_parameters.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_evolution_parameters() 
{
  class_< EvolutionParameters<R> >("EvolutionParameters",init<>())
    .def("minimum_step_size",&EvolutionParameters<R>::minimum_step_size)
    .def("maximum_step_size",&EvolutionParameters<R>::maximum_step_size)
    .def("minimum_basic_set_radius",&EvolutionParameters<R>::minimum_basic_set_radius)
    .def("maximum_basic_set_radius",&EvolutionParameters<R>::maximum_basic_set_radius)
    .def("lock_to_grid_time",&EvolutionParameters<R>::lock_to_grid_time)
    .def("lock_to_grid_steps",&EvolutionParameters<R>::lock_to_grid_steps)
    .def("grid_length",&EvolutionParameters<R>::grid_length)
    .def("argument_grid_length",&EvolutionParameters<R>::argument_grid_length)
    .def("result_grid_length",&EvolutionParameters<R>::result_grid_length)
    .def("bounding_domain_size",&EvolutionParameters<R>::bounding_domain_size)
		.def("hybrid_bounding_domain",&EvolutionParameters<R>::hybrid_bounding_domain)

    .def("grid",&EvolutionParameters<R>::grid)
    .def("finite_grid",&EvolutionParameters<R>::finite_grid)

    .def("set_minimum_step_size",(void(EvolutionParameters<R>::*)(double))&EvolutionParameters<R>::set_minimum_step_size)
    .def("set_maximum_step_size",(void(EvolutionParameters<R>::*)(double))&EvolutionParameters<R>::set_maximum_step_size)
    .def("set_minimum_enclosure_radius",(void(EvolutionParameters<R>::*)(double))&EvolutionParameters<R>::set_minimum_basic_set_radius)
    .def("set_maximum_enclosure_radius",(void(EvolutionParameters<R>::*)(double))&EvolutionParameters<R>::set_maximum_basic_set_radius)
    .def("set_lock_to_grid_time",(void(EvolutionParameters<R>::*)(double))&EvolutionParameters<R>::set_lock_to_grid_time)
    .def("set_lock_to_grid_steps",(void(EvolutionParameters<R>::*)(double))&EvolutionParameters<R>::set_lock_to_grid_steps)
    .def("set_grid_length",(void(EvolutionParameters<R>::*)(double))&EvolutionParameters<R>::set_grid_length)
    .def("set_argument_grid_length",(void(EvolutionParameters<R>::*)(double))&EvolutionParameters<R>::set_argument_grid_length)
    .def("set_result_grid_length",(void(EvolutionParameters<R>::*)(double))&EvolutionParameters<R>::set_result_grid_length)
    .def("set_bounding_domain_size",(void(EvolutionParameters<R>::*)(double))&EvolutionParameters<R>::set_bounding_domain_size)
    .def("set_verbosity",&EvolutionParameters<R>::set_verbosity)
    .def("set_hybrid_grid",&EvolutionParameters<R>::set_hybrid_grid)
    .def("set_grid",&EvolutionParameters<R>::set_grid)
		.def("set_hybrid_bounding_domain",&EvolutionParameters<R>::set_hybrid_bounding_domain)

    .def(self_ns::str(self))
  ;
 
}

template void export_evolution_parameters<FloatPy>();
