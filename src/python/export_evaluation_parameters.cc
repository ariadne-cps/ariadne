/***************************************************************************
 *            python/export_evaluation_parameters.cc
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

#include "python/python_float.h"

#include "geometry/rectangle.h"
#include "geometry/grid.h"
#include "evaluation/evaluation_parameters.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::Evaluation;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_evaluation_parameters() 
{
  class_< EvaluationParameters<R> >("EvaluationParameters",init<>())
    .def("minimum_step_size",&EvaluationParameters<R>::minimum_step_size)
    .def("maximum_step_size",&EvaluationParameters<R>::maximum_step_size)
    .def("minimum_basic_set_radius",&EvaluationParameters<R>::minimum_basic_set_radius)
    .def("maximum_basic_set_radius",&EvaluationParameters<R>::maximum_basic_set_radius)
    .def("lock_to_grid_time",&EvaluationParameters<R>::lock_to_grid_time)
    .def("grid_length",&EvaluationParameters<R>::grid_length)
    .def("argument_grid_length",&EvaluationParameters<R>::argument_grid_length)
    .def("result_grid_length",&EvaluationParameters<R>::result_grid_length)
    .def("bounding_domain_size",&EvaluationParameters<R>::bounding_domain_size)

    .def("grid",&EvaluationParameters<R>::grid)
    .def("finite_grid",&EvaluationParameters<R>::finite_grid)

    .def("set_minimum_step_size",&EvaluationParameters<R>::set_minimum_step_size)
    .def("set_maximum_step_size",&EvaluationParameters<R>::set_maximum_step_size)
    .def("set_minimum_basic_set_radius",&EvaluationParameters<R>::set_minimum_basic_set_radius)
    .def("set_maximum_basic_set_radius",&EvaluationParameters<R>::set_maximum_basic_set_radius)
    .def("set_lock_to_grid_time",&EvaluationParameters<R>::set_lock_to_grid_time)
    .def("set_grid_length",&EvaluationParameters<R>::set_grid_length)
    .def("set_argument_grid_length",&EvaluationParameters<R>::set_argument_grid_length)
    .def("set_result_grid_length",&EvaluationParameters<R>::set_result_grid_length)
    .def("set_bounding_domain_size",&EvaluationParameters<R>::set_bounding_domain_size)
  ;
 
}

template void export_evaluation_parameters<Float>();
