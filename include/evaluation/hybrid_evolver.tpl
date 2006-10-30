/***************************************************************************
 *            hybrid_evolver.tpl
 *
 *  Copyright  2006  Alberto Casagrande,  Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
#include "../exceptions.h"
#include "../geometry/hybrid_set.h"
#include "../system/hybrid_automaton.h"
#include "hybrid_evolver.h"

template<class R>
Ariadne::Evaluation::HybridEvolver<R>::HybridEvolver(Applicator<R>& a, Integrator<R>& i)
  : _applicator(&a), _integrator(&i) 
{
}

template<class R>
Ariadne::Geometry::HybridGridMaskSet<R> 
Ariadne::Evaluation::HybridEvolver<R>::lower_reach(const System::HybridAutomaton<R>&, 
                                                   const Geometry::HybridGridMaskSet<R>&, 
                                                   time_type t1, 
                                                   time_type t2, 
                                                   size_type nmax)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
Ariadne::Geometry::HybridGridMaskSet<R> 
Ariadne::Evaluation::HybridEvolver<R>::upper_reach(const System::HybridAutomaton<R>&, 
                                                   const Geometry::HybridGridMaskSet<R>&, 
                                                   time_type t1, 
                                                   time_type t2, 
                                                   size_type nmax)
{
      throw NotImplemented(__PRETTY_FUNCTION__);
}
