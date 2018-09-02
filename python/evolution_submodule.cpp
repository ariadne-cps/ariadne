/***************************************************************************
 *            evolution_submodule.cpp
 *
 *  Copyright 2008--17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "boost_python.hpp"
#include "utilities.hpp"

#include "dynamics/orbit.hpp"
#include "dynamics/vector_field_evolver.hpp"

#include "hybrid/hybrid_evolver.hpp"
#include "hybrid/hybrid_automata.hpp"
#include "hybrid/hybrid_enclosure.hpp"
#include "hybrid/hybrid_time.hpp"

using namespace boost::python;
using namespace Ariadne;


template<class Orb>
Void export_orbit(const char* name)
{
    class_<Orb> orbit_class(name,no_init);
    orbit_class.def("reach", &Orb::reach,return_value_policy<copy_const_reference>());
    orbit_class.def("evolve", &Orb::final,return_value_policy<copy_const_reference>());
    orbit_class.def("final", &Orb::final,return_value_policy<copy_const_reference>());
    orbit_class.def(self_ns::str(self));
}


template<class Ev, class Init>
Void export_evolver(const char* name)
{
    typedef typename Ev::EnclosureType ES;
    typedef typename Ev::TerminationType Tm;
    typedef typename Ev::OrbitType Orb;

    class_<Ev> evolver_class(name,Init());
    evolver_class.def("orbit",(Orb(Ev::*)(const ES&,const Tm&,Semantics)const) &Ev::orbit);
    evolver_class.def(self_ns::str(self));
}

Void evolution_submodule()
{
    export_orbit< Orbit<HybridEnclosure> >("HybridOrbit");
    //export_evolver<VectorFieldEvolver, init<ContinuousEvolutionParameters> >("VectorFieldEvolver");
    export_evolver<GeneralHybridEvolver, init<GeneralHybridEvolver::SystemType const&> >("GeneralHybridEvolver");
}
