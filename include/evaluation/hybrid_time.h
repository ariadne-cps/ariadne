/***************************************************************************
 *            constraint_hybrid_evolver_plugin.h
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
 
#ifndef ARIADNE_HYBRID_TIME_H
#define ARIADNE_HYBRID_TIME_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "../base/types.h"
#include "../base/tribool.h"
#include "../linear_algebra/vector.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"
#include "../evaluation/declarations.h"
#include "../evaluation/integrator.h"
#include "../evaluation/time_model.h"

namespace Ariadne {  
  namespace Evaluation {
  
    /*! \brief A class containing a continuous time \c time() and a discrete time \c steps() . */
    class HybridTime {
     public:
      HybridTime(Numeric::Rational t) : _time(t), _steps(0) { }
      HybridTime(Numeric::Rational t, Numeric::Integer s) : _time(t), _steps(s) { }
      Numeric::Rational& time() { return this->_time; }
      Numeric::Integer& steps() { return this->_steps; }
      const Numeric::Rational& time() const { return this->_time; }
      const Numeric::Integer& steps() const { return this->_steps; }
      bool operator==(const HybridTime& other) const { return this->_time==other._time && this->_steps==other._steps; }
      bool operator!=(const HybridTime& other) const { return !(*this==other); }
      bool operator<(const HybridTime& other) const { return this->_time<other._time; }
      HybridTime operator+(const HybridTime& other) const { return HybridTime(this->_time+other._time,this->_steps+other._steps); }
      HybridTime operator+(const Numeric::Rational& other) const { return HybridTime(this->_time+other,this->_steps); }
     private:
      Numeric::Rational _time;
      Numeric::Integer _steps;
    };

    typedef HybridTime hybrid_time_type;  



    /*! \brief A class representing a hybrid time and a hybrid basic set. */
    template<class BS>
    class TimeModelHybridBasicSet 
      : public Geometry::HybridBasicSet<BS>
    {                                           
      typedef typename BS::real_type R;
      typedef Numeric::Interval<R> I;
     public:
      TimeModelHybridBasicSet(const id_type& q, const BS& bs)
        : Geometry::HybridBasicSet<BS>(q,bs), _time(I(0),LinearAlgebra::Vector<I>(bs.number_of_generators())), _steps(0) { }
      TimeModelHybridBasicSet(const Numeric::Rational& t, const Numeric::Integer& n, const id_type& q, const BS& bs)
        : Geometry::HybridBasicSet<BS>(q,bs), _time(I(t),LinearAlgebra::Vector<I>(bs.number_of_generators())), _steps(n) { }
      TimeModelHybridBasicSet(const TimeModel<R>& t, const Numeric::Integer& n, const id_type& q, const BS& bs)
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
