/***************************************************************************
 *            hybrid_automaton.cc
 *
 *  Copyright  2007  Pieter Collins
 *  pieter.collins@cwi.nl
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

#include "numeric/float.h"

#include "system/hybrid_automaton.h"
#include "system/hybrid_automaton.code.h"


namespace Ariadne {
  namespace System {
    using namespace Numeric;
    
#ifdef ENABLE_FLOAT64
    template class DiscreteMode<Float64>;
    template class DiscreteTransition<Float64>;
    template class HybridAutomaton<Float64>;
#endif
    
#ifdef ENABLE_FLOATMP
    template class DiscreteMode<FloatMP>;
    template class DiscreteTransition<FloatMP>;
    template class HybridAutomaton<FloatMP>;
#endif
      
  }
}
