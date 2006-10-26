/***************************************************************************
 *            python/export_lorenz_system.cc
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

#include "real_typedef.h"

#include "system/vector_field.h"
#include "system/lorenz_system.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;

#include <boost/python.hpp>
using namespace boost::python;

template<typename R>
void export_lorenz_system() 
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
  typedef Interval<R> I;
  typedef LorenzSystem<R> RLorenzSystem;
  
  class_< LorenzSystem<R>, bases< VectorField<R> > >("LorenzSystem",init<R,R,R>())
    .def("dimension", &LorenzSystem<R>::dimension)
    .def("smoothness", &LorenzSystem<R>::smoothness)
    //.def("__call__", (Vector<F>(LorenzSystem<R>::*)(const Point<R>&)const)(&LorenzSystem<R>::image))
    .def("__call__", (Vector<I>(LorenzSystem<R>::*)(const Rectangle<R>&)const)(&RLorenzSystem::image))
    //.def("jacobian", (Matrix<F>(LorenzSystem<R>::*)(const Point<R>&)const)(&RLorenzSystem::jacobian))
    .def("jacobian", (Matrix<I>(LorenzSystem<R>::*)(const Rectangle<R>&)const)(&RLorenzSystem::jacobian))
    .def(self_ns::str(self))
  ;
}

template void export_lorenz_system<Real>();
