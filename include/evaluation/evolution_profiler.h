/***************************************************************************
 *            evolution_profiler.h
 *
 *  Copyright  2007  Pieter Collins
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
 
/*! \file evolution_profiler.h
 *  \brief Profile for controlling the accuracy of evaluation methods.
 */

#ifndef ARIADNE_EVOLUTION_PROFILE_H
#define ARIADNE_EVOLUTION_PROFILE_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "geometry/declarations.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief Profiler for evolution classes.
     */
    class EvolutionProfiler {
     public:
      size_type subdivisions;
      size_type transitions;
      size_type time_steps;
      Numeric::Rational total_stepping_time;
      Numeric::Rational minimum_time_step;
     public:
      EvolutionProfiler() { this->reset(); }
      void reset() { subdivisions=0; transitions=0; time_steps=0; total_stepping_time=0; minimum_time_step=1; }
    };

    inline std::ostream& operator<<(std::ostream& os, const EvolutionProfiler& p) { 
      return os << "subdivisions=" << p.subdivisions
                << ", transitions=" << p.transitions
                << ", time_steps=" << p.time_steps
                << ", average_time_steps=" << Numeric::Rational(p.total_stepping_time/Numeric::Integer(p.time_steps)).get_d()
                << ", minimum_time_step=" << p.minimum_time_step.get_d()
                << "\n";
    }

  }
}

#endif /* ARIADNE_EVOLUTION_PROFILE_H */
