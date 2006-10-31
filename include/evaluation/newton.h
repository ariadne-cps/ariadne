/***************************************************************************
 *            newton.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file newton.h
 *  \brief Newton and interval Newton methods.
 */

#ifndef _ARIADNE_NEWTON_H
#define _ARIADNE_NEWTON_H

#include <exception>
#include <stdexcept>
#include <string>

#include "../declarations.h"
#include "solver.h"

namespace Ariadne {
  namespace Evaluation {
      
    /*! \ingroup Solve
     *  \brief Interval Newton solver.
     */
    template<class R>
    class IntervalNewtonSolver : public Solver<R>
    {
     public:
      /*! \brief Constructor. */
      IntervalNewtonSolver(R max_error, uint max_steps) : Solver<R>(max_error,max_steps) { }
      
      /*! \brief Solve \f$f(x)=0\f$, using the interval Newton method. */
      Geometry::Rectangle<R>
      solve(const System::VectorField<R>& f, 
            const Geometry::Rectangle<R>& r); 
    };            
      
    /*! \ingroup Solve
     *  \brief Interval Newton solver.
     */
    template<class R>
    Geometry::Rectangle<R>
    interval_newton(const System::VectorField<R>& f, 
                    const Geometry::Rectangle<R>& r, 
                    const R& e,
                    uint max_steps=64);

    
  }
}

#endif /* _ARIADNE_NEWTON_H */
