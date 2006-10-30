/***************************************************************************
 *            python/export_hybrid_evolver.cc
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

#include "geometry/hybrid_set.h"
#include "evaluation/applicator.h"
#include "evaluation/integrator.h"
#include "evaluation/hybrid_evolver.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_hybrid_evolver() 
{
  class_< HybridEvolver<R> >("HybridEvolver",init<Applicator<R>&,Integrator<R>&>())
    .def("lower_reach",&HybridEvolver<R>::lower_reach)
    .def("upper_reach",&HybridEvolver<R>::upper_reach)
  ;
}


template void export_hybrid_evolver<Real>();
