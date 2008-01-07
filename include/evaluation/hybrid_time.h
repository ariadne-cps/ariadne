/***************************************************************************
 *            hybrid_time.h
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
      HybridTime& operator+=(const HybridTime& other) { this->_time+=other._time; this->_steps+=other._steps; return *this; }
      HybridTime& operator+=(const Numeric::Rational& time) { this->_time+=time; return *this; }
      HybridTime& operator++() { ++this->_steps; return *this; }
      HybridTime operator+(const HybridTime& other) const { return HybridTime(*this)+=other; }
      HybridTime operator+(const Numeric::Rational& other) const { return HybridTime(*this)+=other; }
     private:
      Numeric::Rational _time;
      Numeric::Integer _steps;
    };

    typedef HybridTime hybrid_time_type;  


  }
}

#endif /* ARIADNE_HYBRID_TIME_H */
