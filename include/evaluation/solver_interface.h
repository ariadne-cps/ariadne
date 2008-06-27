/***************************************************************************
 *            solver.h
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
 
/*! \file solver.h
 *  \brief Base class for solving equations.
 */

#ifndef ARIADNE_SOLVER_H
#define ARIADNE_SOLVER_H

#include <exception>
#include <stdexcept>
#include <string>

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "function/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "function/function_interface.h"

namespace Ariadne {

    
    /*! \ingroup EvaluatorInterfaces \ingroup Solvers
     *  \brief %Interface for solving (nonlinear) equations. 
     */
    template<class R>
    class SolverInterface
    {
      typedef typename traits<R>::interval_type I;
     public:
      /*! \brief Constructor. */
      SolverInterface(R max_error, uint max_steps);
      
      /*! \brief Virtual destructor. */
      virtual ~SolverInterface();
        
      /*! \brief The maximum permissible error of the solution. */
      const R& maximum_error() const;
      /*! \brief The maximum number of steps allowed before the method must quit. */
      const uint& maximum_number_of_steps() const;
      
      /*! \brief Set the maximum error. */
      void set_maximum_error(R max_error);
      /*! \brief Set the maximum number of steps. */
      void set_maximum_number_of_steps(uint max_steps);
      
      /*! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt. */
      virtual Point<I> solve(const FunctionInterface<R>& f,const Point<I>& pt) = 0;
      /*! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt. */
      virtual Point<I> fixed_point(const Map<R>& f,const Point<I>& pt);
      
     private:
      R _max_error;
      uint _max_steps;
    };
    

} // namespace Ariadne

#include "solver.inline.h"

#endif /* ARIADNE_SOLVER_H */
