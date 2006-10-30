/***************************************************************************
 *            hybrid_evolver.h
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
 
#ifndef _ARIADNE_HYBRID_EVOLVER_H
#define _ARIADNE_HYBRID_EVOLVER_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "../declarations.h"
#include "../system/hybrid_automaton.h"

namespace Ariadne {  
namespace Evaluation {
  

  /*! 
 *  \brief A class for computing the evolution of a hybrid system.
 */
template< class R >
class HybridEvolver
{
 public:
  HybridEvolver(Applicator<R>& applicator, Integrator<R>& integrator);
  
  Geometry::HybridGridMaskSet<R> lower_reach(const System::HybridAutomaton<R>&, 
                                             const Geometry::HybridGridMaskSet<R>&, 
                                             time_type, time_type, size_type);
  
  Geometry::HybridGridMaskSet<R> upper_reach(const System::HybridAutomaton<R>&, 
                                             const Geometry::HybridGridMaskSet<R>&, 
                                             time_type, time_type, size_type);
 private:
  Applicator<R>* _applicator;
  Integrator<R>* _integrator;
};


}
}

#endif /* _ARIADNE_HYBRID_EVOLVER_H */
