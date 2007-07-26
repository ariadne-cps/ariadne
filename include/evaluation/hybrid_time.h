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
#include "../geometry/declarations.h"
#include "../system/declarations.h"
#include "../evaluation/declarations.h"
#include "../evaluation/integrator.h"

namespace Ariadne {  
  namespace Evaluation {
  
    /*! \brief */
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


    template<class R>
    class TimeModel {
      typedef Numeric::Interval<R> I;
     public:
      TimeModel(const I& t, const LinearAlgebra::Vector<I>& dt)
        : _time(t), _time_gradient(dt) { }
      const I& time() const { return this->_time; }
      const LinearAlgebra::Vector<I>& time_gradient() const { return this->_time_gradient; }
      I time_interval() const { return this->_time+inner_product(this->_time_gradient,LinearAlgebra::Vector<I>(this->_time_gradient.size(),I(-1,1))); }
     private:
      I _time;
      LinearAlgebra::Vector<I> _time_gradient;
    };

    template<class R> inline
    TimeModel<R> operator+(const TimeModel<R>& t1, const TimeModel<R>& t2) {
      return TimeModel<R>(t1.time()+t2.time(),t1.time_gradient()+t2.time_gradient());
    }

    template<class R> inline
    TimeModel<R> operator+(const Numeric::Rational& t1, const TimeModel<R>& t2) {
      return TimeModel<R>(Numeric::Interval<R>(t1)+t2.time(),t2.time_gradient());
    }

    template<class R> inline
    TimeModel<R> operator+(const TimeModel<R>& t1, const Numeric::Rational& t2) {
      return TimeModel<R>(t1.time()+Numeric::Interval<R>(t2),t1.time_gradient());
    }

    template<class R> inline
    TimeModel<R> operator-(const Numeric::Rational& t1, const TimeModel<R>& t2) {
      return TimeModel<R>(Numeric::Interval<R>(t1)-t2.time(),-t2.time_gradient());
    }

    template<class R> inline
    bool operator<=(const TimeModel<R>& t1, const Numeric::Rational& t2) {
      return t1.time_interval().upper()<=t2;
    }

    template<class R> inline
    bool operator>=(const TimeModel<R>& t1, const Numeric::Rational& t2) {
      return t1.time_interval().lover() >=t2;
    }

    template<class R> inline
    bool operator==(const TimeModel<R>& t1, const Numeric::Rational& t2) {
      Numeric::Interval<R> t=t1.time_interval();
      return t1.time_interval().lower() <=t2 && t1.time_interval().upper()>=t2;
    }

    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const TimeModel<R>& t1) {
      return os << t1.time_interval();
    }


    template<class BS>
    class TimeModelHybridBasicSet 
      : public Geometry::HybridBasicSet<BS>
    {                                           
      typedef typename BS::real_type R;
      typedef Numeric::Interval<R> I;
     public:
      TimeModelHybridBasicSet(const Numeric::Rational& t, const Numeric::Integer& n, const id_type& q, const BS& bs)
        : Geometry::HybridBasicSet<BS>(q,bs), _time(I(t),LinearAlgebra::Vector<I>(bs.number_of_generators())), _steps(n) { }
      TimeModelHybridBasicSet(const TimeModel<R>& t, const Numeric::Integer& n, const id_type& q, const BS& bs)
        : Geometry::HybridBasicSet<BS>(q,bs), _time(t), _steps(n) { }
      const TimeModel<R>& time() const { return this->_time; }
      const Numeric::Integer& steps() const { return this->_steps; }
     private:
      TimeModel<R> _time;
      Numeric::Integer _steps;
    };
  


  }
}

#endif /* ARIADNE_HYBRID_TIME_H */
