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

#include "python/python_float.h"

#include "geometry/hybrid_set.h"
#include "system/constraint_based_hybrid_automaton.h"
#include "evaluation/applicator.h"
#include "evaluation/integrator.h"
#include "evaluation/constraint_based_hybrid_evolver.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_constraint_based_hybrid_evolver() 
{
  class_< ConstraintBasedHybridEvolver<R> >("ConstraintBasedHybridEvolver",init<Applicator<R>&,Integrator<R>&>()) 
    .def("discrete_step",(HybridSet<R>(ConstraintBasedHybridEvolver<R>::*)(const ConstraintBasedHybridAutomaton<R>&,const HybridSet<R>&)const)&ConstraintBasedHybridEvolver<R>::discrete_step)
    .def("continuous_chainreach",(HybridSet<R>(ConstraintBasedHybridEvolver<R>::*)(const ConstraintBasedHybridAutomaton<R>&,const HybridSet<R>&,const HybridSet<R>&)const)&ConstraintBasedHybridEvolver<R>::continuous_chainreach)

    .def("lower_evolve",(HybridSet<R>(ConstraintBasedHybridEvolver<R>::*)(const ConstraintBasedHybridAutomaton<R>&,const HybridSet<R>&,time_type,size_type)const)&ConstraintBasedHybridEvolver<R>::lower_evolve)
    .def("upper_evolve",(HybridSet<R>(ConstraintBasedHybridEvolver<R>::*)(const ConstraintBasedHybridAutomaton<R>&,const HybridSet<R>&,time_type,size_type)const)&ConstraintBasedHybridEvolver<R>::upper_evolve)
    .def("lower_reach",(HybridSet<R>(ConstraintBasedHybridEvolver<R>::*)(const ConstraintBasedHybridAutomaton<R>&,const HybridSet<R>&,time_type,time_type,size_type)const)&ConstraintBasedHybridEvolver<R>::lower_reach)
    .def("upper_reach",(HybridSet<R>(ConstraintBasedHybridEvolver<R>::*)(const ConstraintBasedHybridAutomaton<R>&,const HybridSet<R>&,time_type,time_type,size_type)const)&ConstraintBasedHybridEvolver<R>::upper_reach)

    .def("chainreach",(HybridSet<R>(ConstraintBasedHybridEvolver<R>::*)(const ConstraintBasedHybridAutomaton<R>&,const HybridSet<R>&,const HybridSet<R>&)const)&ConstraintBasedHybridEvolver<R>::chainreach)
    .def("viable",(HybridSet<R>(ConstraintBasedHybridEvolver<R>::*)(const ConstraintBasedHybridAutomaton<R>&,const HybridSet<R>&)const)&ConstraintBasedHybridEvolver<R>::viable)
    .def("verify",(tribool(ConstraintBasedHybridEvolver<R>::*)(const ConstraintBasedHybridAutomaton<R>&,const HybridSet<R>&,const HybridSet<R>&)const)&ConstraintBasedHybridEvolver<R>::verify)
  ;
}


template void export_constraint_based_hybrid_evolver<Float>();
