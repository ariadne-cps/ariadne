/***************************************************************************
 *            discretiser.code.h
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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
 
#include <iosfwd>
#include "discretiser.h"

#include "output/logging.h"

namespace Ariadne {
  
template<class T, class Aprx, class ES>
Evaluation::Discretiser<System::NumericalSystemInterface<T,ES>,Aprx,ES>::
Discretiser(const EvolutionParameters<R>& parameters, 
            const ApproximatorInterface<Aprx,ES>& approximator)
  : _parameters(parameters.clone()), 
    _approximator(approximator.clone()) 
{ 
}
      

template<class T, class Aprx, class ES>
typename Aprx::CoverListSet 
Evaluation::Discretiser<System::NumericalSystemInterface<T,ES>,Aprx,ES>::
lower_evolve(const System& system, const BasicSet& initial_set, const Time& time) const 
{
  return this->_lower_approximation(system.evolve(this->_over_approximation(initial_set),time,lower_semantics)); 
}
    
template<class T, class Aprx, class ES>
typename Aprx::CoverListSet 
Evaluation::Discretiser<System::NumericalSystemInterface<T,ES>,Aprx,ES>::
lower_reach(const System& system, const BasicSet& initial_set, const Time& time) const 
{
  return this->_lower_approximation(system.reach(this->_over_approximation(initial_set),time,lower_semantics)); 
}
    
template<class T, class Aprx, class ES>
std::pair<typename Aprx::CoverListSet,typename Aprx::CoverListSet>
Evaluation::Discretiser<System::NumericalSystemInterface<T,ES>,Aprx,ES>::
lower_reach_evolve(const System& system, const BasicSet& initial_set, const Time& time) const 
{
  return this->_lower_approximation(system.reach_evolve(this->_over_approximation(initial_set),time,lower_semantics)); 
}


template<class T, class Aprx, class ES>
typename Aprx::PartitionListSet 
Evaluation::Discretiser<System::NumericalSystemInterface<T,ES>,Aprx,ES>::
upper_evolve(const System& system, const BasicSet& initial_set, const Time& time) const 
{
  return this->_outer_approximation(system.evolve(this->_over_approximation(initial_set),time,upper_semantics)); 
}
    
template<class T, class Aprx, class ES>
typename Aprx::PartitionListSet 
Evaluation::Discretiser<System::NumericalSystemInterface<T,ES>,Aprx,ES>::
upper_reach(const System& system, const BasicSet& initial_set, const Time& time) const 
{
  return this->_outer_approximation(system.reach(this->_over_approximation(initial_set),time,upper_semantics)); 
}
    
template<class T, class Aprx, class ES>
std::pair<typename Aprx::PartitionListSet,typename Aprx::PartitionListSet>
Evaluation::Discretiser<System::NumericalSystemInterface<T,ES>,Aprx,ES>::
upper_reach_evolve(const System& system, const BasicSet& initial_set, const Time& time) const 
{
  return this->_outer_approximation(system.reach_evolve(this->_over_approximation(initial_set),time,upper_semantics)); 
}





}
