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
  namespace Evaluation {
  



    template<class SetInterface> 
    class HybridTimedSet 
    {
     public:
      HybridTimedSet(const time_type& t, const Geometry::DiscreteState& id, const SetInterface& s)
        : _time(t), _discrete_state(id), _continuous_state_set(s) { }
      const time_type& time() const { return _time; }
      const Geometry::DiscreteState& discrete_state() const { return _discrete_state; }
      const SetInterface& continuous_state_set() const { return _continuous_state_set; } 
      
      bool operator==(const HybridTimedSet& other) const { 
        return this->_time == other._time 
          && this->_discrete_state==other._discrete_state 
          && this->_continuous_state_set==other._continuous_state_set; }
      bool operator!=(const HybridTimedSet& other) const { return !(*this==other); }
      bool operator<=(const HybridTimedSet& other) const { return this->_time <= other._time; }
     private:
      time_type _time;
      Geometry::DiscreteState _discrete_state;
      SetInterface _continuous_state_set;
    };

    /*! \brief A class representing a hybrid time and a hybrid basic set. */
    template<class BS>
    class TimeModelHybridBasicSet 
      : public Geometry::HybridBasicSet<BS>
    {                                           
      typedef typename BS::real_type R;
      typedef Numeric::Interval<R> I;
     public:
      TimeModelHybridBasicSet(const Geometry::DiscreteState& q, const BS& bs)
        : Geometry::HybridBasicSet<BS>(q,bs), _time(I(0),LinearAlgebra::Vector<I>(bs.number_of_generators())), _steps(0) { }
      TimeModelHybridBasicSet(const Numeric::Rational& t, const Numeric::Integer& n, const Geometry::DiscreteState& q, const BS& bs)
        : Geometry::HybridBasicSet<BS>(q,bs), _time(I(t),LinearAlgebra::Vector<I>(bs.number_of_generators())), _steps(n) { }
      TimeModelHybridBasicSet(const TimeModel<R>& t, const Numeric::Integer& n, const Geometry::DiscreteState& q, const BS& bs)
        : Geometry::HybridBasicSet<BS>(q,bs), _time(t), _steps(n) { 
        if(t.number_of_generators()!=bs.number_of_generators()) { std::cerr << "t=" << t << "\nbs=" << bs << std::endl; }
        assert(t.number_of_generators()==bs.number_of_generators()); }
      TimeModelHybridBasicSet(const TimeModelHybridBasicSet<BS>& thbs) 
        : Geometry::HybridBasicSet<BS>(thbs), _time(thbs._time), _steps(thbs._steps) { }
      const TimeModel<R>& time() const { return this->_time; }
      const Numeric::Integer& steps() const { return this->_steps; }
     private:
      TimeModel<R> _time;
      Numeric::Integer _steps;
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
