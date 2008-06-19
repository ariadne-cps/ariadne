/***************************************************************************
 *            hybrid_timed_set.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
#ifndef ARIADNE_HYBRID_TIMED_SET_H
#define ARIADNE_HYBRID_TIMED_SET_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/tribool.h"
#include "linear_algebra/vector.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"
#include "evaluation/time_model.h"

#include "geometry/discrete_state.h"

namespace Ariadne {  
  
  



    template<class SetInterface> 
    class HybridTimedSet 
    {
     public:
      HybridTimedSet(const time_type& t, const DiscreteState& id, const SetInterface& s)
        : _time(t), _discrete_state(id), _continuous_state_set(s) { }
      const time_type& time() const { return _time; }
      const DiscreteState& discrete_state() const { return _discrete_state; }
      const SetInterface& continuous_state_set() const { return _continuous_state_set; } 
      
      bool operator==(const HybridTimedSet& other) const { 
        return this->_time == other._time 
          && this->_discrete_state==other._discrete_state 
          && this->_continuous_state_set==other._continuous_state_set; }
      bool operator!=(const HybridTimedSet& other) const { return !(*this==other); }
      bool operator<=(const HybridTimedSet& other) const { return this->_time <= other._time; }
     private:
      time_type _time;
      DiscreteState _discrete_state;
      SetInterface _continuous_state_set;
    };

    /*! \brief A class representing a hybrid time and a hybrid basic set. */
    template<class BS>
    class TimeModelHybridBasicSet 
      : public HybridBasicSet<BS>
    {                                           
      typedef typename BS::real_type R;
      typedef Interval<R> I;
     public:
      TimeModelHybridBasicSet(const DiscreteState& q, const BS& bs)
        : HybridBasicSet<BS>(q,bs), _time(I(0),Vector<I>(bs.number_of_generators())), _steps(0) { }
      TimeModelHybridBasicSet(const Rational& t, const Integer& n, const DiscreteState& q, const BS& bs)
        : HybridBasicSet<BS>(q,bs), _time(I(t),Vector<I>(bs.number_of_generators())), _steps(n) { }
      TimeModelHybridBasicSet(const TimeModel<R>& t, const Integer& n, const DiscreteState& q, const BS& bs)
        : HybridBasicSet<BS>(q,bs), _time(t), _steps(n) { 
        if(t.number_of_generators()!=bs.number_of_generators()) { std::cerr << "t=" << t << "\nbs=" << bs << std::endl; }
        assert(t.number_of_generators()==bs.number_of_generators()); }
      TimeModelHybridBasicSet(const TimeModelHybridBasicSet<BS>& thbs) 
        : HybridBasicSet<BS>(thbs), _time(thbs._time), _steps(thbs._steps) { }
      const TimeModel<R>& time() const { return this->_time; }
      const Integer& steps() const { return this->_steps; }
     private:
      TimeModel<R> _time;
      Integer _steps;
    };
  
    template<class BS> inline 
    std::ostream& 
    operator<<(std::ostream& os, const TimeModelHybridBasicSet<BS>& thbs) 
    {
      return os << "{ t=" << thbs.time() << ", n=" << thbs.steps() << ", q=" << thbs.discrete_state() << ", s=" << thbs.continuous_state_set() << "}";
    }

  }
}

#endif /* ARIADNE_HYBRID_TIME_H */
