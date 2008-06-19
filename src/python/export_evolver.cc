/***************************************************************************
 *            python/export_evolver.cc
 *
 *  Copyright  2005-8  Pieter Collins
 *
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
#include "geometry/zonotope.h"

#include "system/map.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/approximator_interface.h"
#include "evaluation/satisfier_interface.h"
#include "evaluation/subdivider_interface.h"
#include "evaluation/applicator_interface.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/evolver_interface.h"
#include "evaluation/map_evolver.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/set_based_hybrid_evolver.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_evolver()
{
  typedef Zonotope<R> ZES;

  class_< MapEvolver<ZES> > map_evolver_class("MapEvolver",no_init);
  map_evolver_class.def(init<const EvolutionParameters<R>&,const ApplicatorInterface<ZES>&,
                             const SubdividerInterface<ZES>&,const ReducerInterface<ZES>&>());
  map_evolver_class.def("evolution",&VectorFieldEvolver<ZES>::evolution);

  class_< VectorFieldEvolver<ZES> > vector_field_evolver_class("VectorFieldEvolver",no_init);
  vector_field_evolver_class.def(init<const EvolutionParameters<R>&,const IntegratorInterface<ZES>&,
                                      const SubdividerInterface<ZES>&,const ReducerInterface<ZES>&>());
  vector_field_evolver_class.def("evolution",&VectorFieldEvolver<ZES>::evolution);

  class_< SetBasedHybridEvolver<ZES> > set_based_hybrid_evolver_class("SetBasedHybridEvolver",no_init);
  set_based_hybrid_evolver_class.def(init<const EvolutionParameters<R>&,
                                          const ApplicatorInterface<ZES>&,const IntegratorInterface<ZES>&,const SatisfierInterface<ZES>&,
                                          const SubdividerInterface<ZES>&,const ReducerInterface<ZES>&>());
  set_based_hybrid_evolver_class.def("evolution",&SetBasedHybridEvolver<ZES>::evolution);


}

template void export_evolver<FloatPy>();
