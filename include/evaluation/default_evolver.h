/***************************************************************************
 *            default_evolver.h
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
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file default_evolver.h
 *  \brief Default evolver classes.
 */

#ifndef ARIADNE_DEFAULT_EVOLVER_H
#define ARIADNE_DEFAULT_EVOLVER_H

#include "evaluation/standard_applicator.h"
#include "evaluation/standard_integrator.h"
#include "evaluation/standard_satisfier.h"
#include "evaluation/standard_subdivider.h"
#include "evaluation/cascade_reducer.h"

#include "evaluation/map_evolver.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/impact_system_evolver.h"
#include "evaluation/hybrid_evolver.h"

#include "evaluation/standard_approximator.h"

#include "evaluation/evolution_parameters.h"

namespace Ariadne {
  



template<class R, class ES>
EvolverInterface<Map<R>,ES>*
make_default_evolver(const Map<R>& sys, const ES& encl, const EvolutionParameters<R>& params) {
  StandardApplicator< ES > applicator;
  StandardSubdivider< ES > subdivider;
  CascadeReducer< ES > reducer(3);
  return new Evolver< Map<R>, ES >(params,applicator,subdivider,reducer);
}


template<class R, class ES>
EvolverInterface<VectorField<R>,ES>*
make_default_evolver(const VectorField<R>& sys, const ES& encl, const EvolutionParameters<R>& params) {
  StandardIntegrator< ES > integrator;
  StandardSubdivider< ES > subdivider;
  CascadeReducer< ES > reducer(3);
  return new Evolver< VectorField<R>, ES >(params,integrator,subdivider,reducer);
}


template<class R, class ES>
EvolverInterface<ImpactSystem<R>,ES>*
make_default_evolver(const ImpactSystem<R>& sys, const ES& encl, const EvolutionParameters<R>& params) {
  StandardApplicator< ES > applicator;
  StandardIntegrator< ES > integrator;
  StandardSatisfier< ES > satisfier;
  StandardSubdivider< ES > subdivider;
  CascadeReducer< ES > reducer(3);
  return new Evolver< ImpactSystem<R>, ES >(params,applicator,integrator,satisfier,subdivider,reducer);
}


template<class Sys, class ES>
EvolverInterface<Sys,ES>*
make_default_evolver(const Sys& sys, const ES& encl) {
  EvolutionParameters<typename Sys::real_type> params;
  return make_default_evolver(sys,encl,params);
}




template<class Sys, class ES> 
EvolverInterface<Sys,ES>*
default_evolver() {
  Sys* sys=0; ES* encl=0;
  return make_default_evolver(*sys,*encl);
}


template<class Sys, class ES, class R> 
EvolverInterface<Sys,ES>*
default_evolver(const EvolutionParameters<R>& params) {
  Sys* sys=0; ES* encl=0;
  return make_default_evolver(*sys,*encl,params);
}




} // namespace Ariadne



#endif /* ARIADNE_DEFAULT_EVOLVER_H */
