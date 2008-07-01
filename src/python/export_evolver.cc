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

#include "system/impact_system.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/approximator_interface.h"
#include "evaluation/satisfier_interface.h"
#include "evaluation/subdivider_interface.h"
#include "evaluation/applicator_interface.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/evolver_interface.h"
#include "evaluation/map_evolver.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/impact_system_evolver.h"
#include "evaluation/set_based_hybrid_evolver.h"

#include "evaluation/default_evolver.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class Sys, class ES>
class EvolverWrapper
  : public EvolverBase<Sys,ES>,
    public wrapper< EvolverInterface<Sys,ES> >
{
  typedef typename Sys::time_type T;
 public:
  EvolverWrapper<Sys,ES>* clone() const { return this->get_override("clone")(); }
 protected:
  void _evolution(ListSet<ES>&,ListSet<ES>&,const Sys&,const ES&,const T&,Semantics,bool) const  { this->get_override("_evolution")(); }
};


template<class R>
void export_evolver()
{
  class Rational;
  typedef Zonotope<R> ZES;
  typedef ListSet< Zonotope<R> > ZESL;
  typedef typename Map<R>::time_type MapT;

  class_< EvolverWrapper<Map<R>,ZES>, boost::noncopyable >("MapZonotopeEvolverInterface",init< >());
  class_< EvolverWrapper<VectorField<R>,ZES>, boost::noncopyable >("VectorFieldZonotopeEvolverInterface",init< >());
  class_< EvolverWrapper<ImpactSystem<R>,ZES>, boost::noncopyable >("ImpactSystemZonotopeEvolverInterface",init< >());
  class_< EvolverWrapper<HybridAutomaton<R>,ZES>, boost::noncopyable >("HybridAutomatonZonotopeEvolverInterface",init< >());

  class_< Evolver<Map<R>,ZES>, bases< EvolverInterface<Map<R>,ZES> > > 
    map_evolver_class("MapZonotopeEvolver",init<const EvolutionParameters<R>&,const ApplicatorInterface<ZES>&, const SubdividerInterface<ZES>&,const ReducerInterface<ZES>&>());
  map_evolver_class.def("evolve", (ZESL(Evolver<Map<R>,ZES>::*)(const Map<R>&,const ZES&,const MapT&)const) &Evolver<Map<R>,ZES>::evolve);
  map_evolver_class.def("reach", (ZESL(Evolver<Map<R>,ZES>::*)(const Map<R>&,const ZES&,const MapT&)const) &Evolver<Map<R>,ZES>::reach);
  //map_evolver_class.def("evolve",&Evolver<Map<R>,ZES>::evolve);
  //map_evolver_class.def("reach",&Evolver<Map<R>,ZES>::reach);
  map_evolver_class.def(self_ns::str(self));

  class_< Evolver<VectorField<R>,ZES>, bases< EvolverInterface<VectorField<R>,ZES> > > 
    vector_field_evolver_class("VectorFieldEvolver",no_init);
  vector_field_evolver_class.def(init<const EvolutionParameters<R>&,const IntegratorInterface<ZES>&,
                                      const SubdividerInterface<ZES>&,const ReducerInterface<ZES>&>());
  vector_field_evolver_class.def("evolve",&Evolver<VectorField<R>,ZES>::evolve);
  vector_field_evolver_class.def("reach",&Evolver<VectorField<R>,ZES>::reach);
  vector_field_evolver_class.def(self_ns::str(self));

  class_< Evolver<ImpactSystem<R>,ZES>, bases< EvolverInterface<ImpactSystem<R>,ZES> > > 
    impact_system_evolver_class("ImpactSystemEvolver",no_init);
  impact_system_evolver_class.def(init<const EvolutionParameters<R>&>());
  impact_system_evolver_class.def("evolve",&Evolver<ImpactSystem<R>,ZES>::evolve);
  impact_system_evolver_class.def("reach",&Evolver<ImpactSystem<R>,ZES>::reach);
  impact_system_evolver_class.def(self_ns::str(self));

  class_< Evolver<HybridAutomaton<R>,ZES> > set_based_hybrid_evolver_class("SetBasedHybridEvolver",no_init);
  set_based_hybrid_evolver_class.def(init<const EvolutionParameters<R>&,
                                          const ApplicatorInterface<ZES>&,const IntegratorInterface<ZES>&,const SatisfierInterface<ZES>&,
                                          const SubdividerInterface<ZES>&,const ReducerInterface<ZES>&>());
  set_based_hybrid_evolver_class.def("evolve",&Evolver<HybridAutomaton<R>,ZES>::evolve);
  set_based_hybrid_evolver_class.def("reach",&Evolver<HybridAutomaton<R>,ZES>::reach);
  set_based_hybrid_evolver_class.def(self_ns::str(self));

  def("default_evolver", (EvolverInterface< Map<R>, Zonotope<R> >*(*)(const Map<R>&, const Zonotope<R>&, const EvolutionParameters<R>&)) &make_default_evolver, return_value_policy<manage_new_object>());
  def("default_evolver",(EvolverInterface< VectorField<R>, Zonotope<R> >*(*)(const VectorField<R>&, const Zonotope<R>&, const EvolutionParameters<R>&)) &make_default_evolver, return_value_policy<manage_new_object>());
  def("default_evolver",(EvolverInterface< ImpactSystem<R>, Zonotope<R> >*(*)(const ImpactSystem<R>&, const Zonotope<R>&, const EvolutionParameters<R>&)) &make_default_evolver, return_value_policy<manage_new_object>());

  def("default_evolver", (EvolverInterface< Map<R>, Zonotope<R> >*(*)(const Map<R>&, const Zonotope<R>&)) &make_default_evolver, return_value_policy<manage_new_object>());
  def("default_evolver",(EvolverInterface< VectorField<R>, Zonotope<R> >*(*)(const VectorField<R>&, const Zonotope<R>&)) &make_default_evolver, return_value_policy<manage_new_object>());
  def("default_evolver",(EvolverInterface< ImpactSystem<R>, Zonotope<R> >*(*)(const ImpactSystem<R>&, const Zonotope<R>&)) &make_default_evolver, return_value_policy<manage_new_object>());

}



template void export_evolver<FloatPy>();
