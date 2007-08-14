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


    /*! \brief A class representing a time as a fuzzy affine function of independent variables.
     *
     * \internal This class is a "hack" and should at some point be incorporated into an "affine model" framework,
     * if used at all. */
    template<class R>
    class TimeModel {
      typedef Numeric::Interval<R> I;
     public:
      TimeModel() : _average(), _gradient() { }
      TimeModel(uint d) : _average(), _gradient(d) { }
      TimeModel(const I& t, const LinearAlgebra::Vector<I>& dt)
        : _average(t), _gradient(dt) { }
      TimeModel(const I& t, const LinearAlgebra::Vector<I>& dt, const I& ndt)
        : _average(t), _gradient(dt.size()+1) { for(uint i=0; i!=dt.size(); ++i) { this->_gradient(i)=dt(i); } this->_gradient(dt.size())=ndt; }
      TimeModel(const TimeModel<R>& t)
        : _average(t._average), _gradient(t._gradient) { }
      size_type number_of_generators() const { return this->_gradient.size(); }
      const I& average() const { return this->_average; }
      const LinearAlgebra::Vector<I>& gradient() const { return this->_gradient; }
      I bound() const { return this->_average+inner_product(this->_gradient,LinearAlgebra::Vector<I>(this->_gradient.size(),I(-1,1))); }
      I value(const LinearAlgebra::Vector<I>& e) const { return this->_average+inner_product(this->_gradient,e); }
     private:
      I _average;
      LinearAlgebra::Vector<I> _gradient;
    };

    template<class R> inline
    TimeModel<R> operator+(const TimeModel<R>& t1, const TimeModel<R>& t2) {
      return TimeModel<R>(t1.average()+t2.average(),t1.gradient()+t2.gradient());
    }

    template<class R> inline
    TimeModel<R> operator+(const Numeric::Rational& t1, const TimeModel<R>& t2) {
      return TimeModel<R>(Numeric::Interval<R>(t1)+t2.average(),t2.gradient());
    }

    template<class R> inline
    TimeModel<R> operator+(const TimeModel<R>& t1, const Numeric::Rational& t2) {
      return TimeModel<R>(t1.average()+Numeric::Interval<R>(t2),t1.gradient());
    }

    template<class R> inline
    TimeModel<R> operator-(const TimeModel<R>& t1, const TimeModel<R>& t2) {
      return TimeModel<R>(t1.average()-t2.average(),t1.gradient()-t2.gradient());
    }

    template<class R> inline
    TimeModel<R> operator-(const Numeric::Rational& t1, const TimeModel<R>& t2) {
      return TimeModel<R>(Numeric::Interval<R>(t1)-t2.average(),-t2.gradient());
    }

    template<class R> inline
    TimeModel<R> operator-(const TimeModel<R>& t1, const Numeric::Rational& t2) {
      return TimeModel<R>(t1.average()-Numeric::Interval<R>(t2),t1.gradient());
    }

    template<class R> inline
    TimeModel<R> operator*(const TimeModel<R>& t1, const int& s) {
      return TimeModel<R>(t1.average()*s,t1.gradient()*s);
    }

    template<class R> inline
    TimeModel<R> operator/(const TimeModel<R>& t1, const int& s) {
      return TimeModel<R>(t1.average()/s,t1.gradient()/s);
    }

    template<class R> inline
    tribool operator<=(const TimeModel<R>& t1, const Numeric::Rational& t2) {
      return t1.bound() <= t2;
    }

    template<class R> inline
    tribool operator>(const TimeModel<R>& t1, const Numeric::Rational& t2) {
      return t1.bound()>t2;
    }

    template<class R> inline
    tribool operator>=(const TimeModel<R>& t1, const Numeric::Rational& t2) {
      return t1.bound() >=t2;
    }

    template<class R> inline
    tribool operator>=(const TimeModel<R>& t1, const TimeModel<R>& t2) {
      return (t1-t2).bound() >= 0;
    }

    template<class R> inline
    tribool operator<=(const TimeModel<R>& t1, const TimeModel<R>& t2) {
      return (t1-t2).bound() <= 0;
    }

    template<class R> inline
    tribool operator==(const TimeModel<R>& t1, const Numeric::Rational& t2) {
      return t1.bound() == t2;
    }

    template<class R> inline
    TimeModel<R> upper_bound(const TimeModel<R>& t1, const TimeModel<R>& t2) {
      if(t1>=t2) { 
        return t1;
      } else if(t2>=t1) {
        return t2;
      } else {
        std::cerr << "Warning: upper_bound(TimeModel,TimeModel) not correctly implemented.\n";
        return t1;
      }
    }


    template<class R> inline
    TimeModel<R> lower_bound(const TimeModel<R>& t1, const TimeModel<R>& t2) {
      if(t1<=t2) { 
        return t1;
      } else if(t2<=t1) {
        return t2;
      } else {
        std::cerr << "Warning: lower_bound(TimeModel,TimeModel) not correctly implemented.\n";
        return t1;
      }
    }


    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const TimeModel<R>& t1) {
      //return os << t1.average() << "+" << LinearAlgebra::midpoint(t1.gradient()) << "xe";
      return os << t1.average() << "+" << t1.gradient() << "xe";
    }


    template<class R> inline
    TimeModel<R> lower_bound(TimeModel<R>& t) {
      typedef Numeric::Interval<R> I;
      LinearAlgebra::Vector<I> g=LinearAlgebra::midpoint(t.gradient());
      I a=t.average()+LinearAlgebra::inner_product(g-t.gradient(),LinearAlgebra::Vector<I>(g.size(),Numeric::Interval<I>(-1,1)));
      return TimeModel<R>(a.lower(),t);
    }
      
    template<class R> inline
    TimeModel<R> upper_bound(TimeModel<R>& t) {
      typedef Numeric::Interval<R> I;
      LinearAlgebra::Vector<I> g=LinearAlgebra::midpoint(t.gradient());
      I a=t.average()+LinearAlgebra::inner_product(g-t.gradient(),LinearAlgebra::Vector<I>(g.size(),Numeric::Interval<I>(-1,1)));
      return TimeModel<R>(a.upper(),t);
    }
      

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
