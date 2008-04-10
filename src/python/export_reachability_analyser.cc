/***************************************************************************
 *            python/export_reachability_analyser.cc
 *
 *  Copyright  2008  Pieter Collins
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

#include "geometry/grid_approximation_scheme.h"
#include "geometry/hybrid_grid_approximation_scheme.h"

#include "system/map.h"
#include "system/vector_field.h"
#include "system/hybrid_automaton.h"

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
#include "evaluation/reachability_analyser.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

namespace Ariadne { namespace System {
template<class Sys> struct Name;
template<class R> struct Name< Map<R> > { static std::string string() { return "Map"; } };
template<class R> struct Name< VectorField<R> > { static std::string string() { return "VectorField"; } };
template<class R> struct Name< HybridAutomaton<R> > { static std::string string() { return "HybridAutomaton"; } };
} }

template<class Sys, class Aprx>
void export_system_reachability_analyser()
{
  typedef typename Aprx::real_type R;
  std::string name=System::Name<Sys>::string()+"ReachabilityAnalyser";

  class_< ReachabilityAnalyser<Sys,Aprx> > reachability_analyser_class(name.c_str(),no_init);
  //reachability_analyser_class.def(init<const EvolutionParameters<R>&());
  //reachability_analyser_class.def("parameters",(EvolutionParameters<R>&(ReachabilityAnalyser<Sys,Aprx>::*)())&ReachabilityAnalyser<Sys,Aprx>::parameters,reference_existing_object);
  reachability_analyser_class.def("lower_evolve",&ReachabilityAnalyser<Sys,Aprx>::lower_evolve,return_value_policy<manage_new_object>());
  reachability_analyser_class.def("lower_reach",&ReachabilityAnalyser<Sys,Aprx>::lower_reach,return_value_policy<manage_new_object>());
  reachability_analyser_class.def("upper_evolve",&ReachabilityAnalyser<Sys,Aprx>::upper_evolve,return_value_policy<manage_new_object>());
  reachability_analyser_class.def("upper_reach",&ReachabilityAnalyser<Sys,Aprx>::upper_reach,return_value_policy<manage_new_object>());
  reachability_analyser_class.def("chain_reach",&ReachabilityAnalyser<Sys,Aprx>::chain_reach,return_value_policy<manage_new_object>());
  reachability_analyser_class.def("viable",&ReachabilityAnalyser<Sys,Aprx>::viable,return_value_policy<manage_new_object>());
}

template<class R>
void export_reachability_analyser()
{
  export_system_reachability_analyser< Map<R>,GridApproximationScheme<R> >();
  export_system_reachability_analyser< VectorField<R>,GridApproximationScheme<R> >();
  export_system_reachability_analyser< HybridAutomaton<R>,HybridGridApproximationScheme<R> >();
}

template void export_reachability_analyser<FloatPy>();
