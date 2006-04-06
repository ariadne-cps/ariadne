/***************************************************************************
 *            python/export_linear_program.cc
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

#include "base/numerical_type.h"

#include "linear_algebra/linear_program.h"

#include <boost/python.hpp>

#include "python/typedefs.h"
#include "python/python_utilities.h"

#include <utility>  //for std::pair

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::self_ns::str;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;

using boost::python::tuple;
using boost::python::extract;

inline FMatrix flinear_program_tableau(const FLinearProgram& lp) {
  return lp.tableau();
}


void export_linear_program() {

  class_<FLinearProgram>("FLinearProgram",init<FMatrix,FVector,FVector>())
    .def(init<FMatrix>())
    .def("solve",&FLinearProgram::solve)
    .def("is_satisfiable",&FLinearProgram::is_satisfiable)
    .def("optimizing_point",&FLinearProgram::optimizing_point)
    .def("optimal_value",&FLinearProgram::optimal_value)
    .def("tableau",&flinear_program_tableau)
    .def(str(self))    // __str__
  ;
}
