/***************************************************************************
 *            python/export_solve.cc
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

#include "evaluation/newton.h"
#include "python/utilities.h"

using namespace Ariadne;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_solver() 
{
  return_value_policy<copy_const_reference> copy_const_reference;
  
  class_< IntervalNewtonSolver<R> >("IntervalNewtonSolver",init<R,uint>())
    .def(init<double,uint>())
    .def("maximum_error",&SolverInterface<R>::maximum_error,copy_const_reference)
    .def("maximum_number_of_steps",&SolverInterface<R>::maximum_number_of_steps,copy_const_reference)
    .def("solve",&IntervalNewtonSolver<R>::solve)
    .def("fixed_point",&SolverInterface<R>::fixed_point)
  ;
  
}

template void export_solver<FloatPy>();
