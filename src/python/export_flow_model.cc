/***************************************************************************
 *            python/export_flow_model.cc
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

#include "function/flow_model.h"
#include "evaluation/integrator.code.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
FlowModel<R> 
integrate_(const TaylorDerivative< Interval<R> >& vf, const LinearAlgebra::Vector< Interval<R> > x,smoothness_type so) {
  return FlowModel<R>(Function::integrate(vf,Geometry::Point<Interval<R> >(x),so)); 
}

template<class R>
void export_flow_model() 
{
  class_< FlowModel<R> > flow_model_class("FlowModel", no_init);
  flow_model_class.def("evaluate",&FlowModel<R>::evaluate);
  flow_model_class.def(self_ns::str(self));

  def("integrate",&integrate_<R>);
}

template void export_flow_model<FloatPy>();
