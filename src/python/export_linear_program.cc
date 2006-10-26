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
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/linear_program.h"

#include "python/python_utilities.h"
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_linear_program() 
{
  class_< LinearProgram<R> >(python_name<R>("LinearProgram").c_str(),
                             init< Matrix<R>, Vector<R>, Vector<R> >())
    .def(init< Matrix<R> >())
    .def("solve",&LinearProgram<R>::solve)
    .def("is_feasible",&LinearProgram<R>::is_feasible)
    .def("optimizing_point",&LinearProgram<R>::optimizing_point)
    .def("optimal_value",&LinearProgram<R>::optimal_value)
//    .def("tableau",&LinearProgram<R>::tableau,return_value_policy<copy_const_reference>())
    .def(self_ns::str(self))
  ;
}

template void export_linear_program<Rational>();
