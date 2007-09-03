/***************************************************************************
 *            time_model.h
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
 
#ifndef ARIADNE_TIME_MODEL_H
#define ARIADNE_TIME_MODEL_H

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
      


  }
}

#endif /* ARIADNE_TIME_MODEL_H */
