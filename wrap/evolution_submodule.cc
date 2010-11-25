/***************************************************************************
 *            evolution_submodule.cc
 *
 *  Copyright 2008  Pieter Collins
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

#include <boost/python.hpp>

#include "orbit.h"
#include "vector_field_evolver.h"
#include "taylor_set.h"

#include "hybrid_evolver.h"
#include "hybrid_automaton.h"
#include "hybrid_set.h"
#include "hybrid_time.h"

#include "utilities.h"

using namespace boost::python;
using namespace Ariadne;


template<class Orb>
void export_orbit(const char* name)
{
    class_<Orb> orbit_class(name,no_init);
    orbit_class.def("reach", &Orb::reach,return_value_policy<copy_const_reference>());
    orbit_class.def("evolve", &Orb::final,return_value_policy<copy_const_reference>());
    orbit_class.def("final", &Orb::final,return_value_policy<copy_const_reference>());
    orbit_class.def(self_ns::str(self));
}


template<class Ev, class Init>
void export_evolver(const char* name)
{
    typedef typename Ev::SystemType Sys;
    typedef typename Ev::EnclosureType ES;
    typedef typename Ev::TimeType Tm;
    typedef typename Ev::OrbitType Orb;

    class_<Ev> evolver_class(name,Init());
    evolver_class.def("orbit",(Orb(Ev::*)(const Sys&,const ES&,const Tm&)) &Ev::orbit);
    evolver_class.def(self_ns::str(self));
}

void evolution_submodule()
{
    export_orbit< Orbit<TaylorConstrainedImageSet> >("ContinuousOrbit");
    export_orbit< Orbit<HybridTaylorConstrainedImageSet> >("HybridOrbit");
    //export_evolver<VectorFieldEvolver, init<ContinuousEvolutionParameters> >("VectorFieldEvolver");
    export_evolver<GeneralHybridEvolver, init<> >("GeneralHybridEvolver");
}
