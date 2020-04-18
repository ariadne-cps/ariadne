/***************************************************************************
 *            evolution_submodule.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include "pybind11.hpp"
#include "utilities.hpp"

#include "dynamics/orbit.hpp"
#include "dynamics/map_evolver.hpp"
#include "dynamics/vector_field_evolver.hpp"

#include "hybrid/hybrid_evolver.hpp"
#include "hybrid/hybrid_automata.hpp"
#include "hybrid/hybrid_enclosure.hpp"
#include "hybrid/hybrid_time.hpp"

using namespace Ariadne;


template<class Orb>
Void export_orbit(pybind11::module& module, const char* name)
{
    pybind11::class_<Orb> orbit_class(module,name);
    orbit_class.def("reach", &Orb::reach);
    orbit_class.def("evolve", &Orb::final);
    orbit_class.def("final", &Orb::final);
    orbit_class.def("__str__", &__cstr__<Orb>);
}


template<class Ev, class... Params>
Void export_evolver(pybind11::module& module, const char* name)
{
    typedef typename Ev::EnclosureType ES;
    typedef typename Ev::TerminationType Tm;
    typedef typename Ev::OrbitType Orb;

    pybind11::class_<Ev> evolver_class(module,name);
    evolver_class.def(pybind11::init<Params...>());
    evolver_class.def("orbit",(Orb(Ev::*)(const ES&,const Tm&,Semantics)const) &Ev::orbit);
    evolver_class.def("__str__",&__cstr__<Ev>);
}

Void evolution_submodule(pybind11::module& module)
{
    export_orbit< Orbit<HybridEnclosure> >(module, "HybridOrbit");
    export_evolver<MapEvolver, IteratedMap>(module,"MapEvolver");
    export_evolver<VectorFieldEvolver, VectorField, IntegratorInterface const&>(module,"VectorFieldEvolver");
    export_evolver<GeneralHybridEvolver, GeneralHybridEvolver::SystemType const&>(module,"GeneralHybridEvolver");
}
