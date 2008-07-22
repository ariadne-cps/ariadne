/***************************************************************************
 *            solver_base.h
 *
 *  Copyright  2006-8  Pieter Collins
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
 
/*! \file solver_base.h
 *  \brief Common functionality for solver classes.
 */

#ifndef ARIADNE_SOLVER_BASE_H
#define ARIADNE_SOLVER_BASE_H

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

#include "evaluation/solver_interface.h"
#include "function/function_interface.h"
#include "function/difference_function.h"
#include "system/map.h"

namespace Ariadne {

    
/*! \ingroup \ingroup Solvers
 *  \brief %Common functionality for solving (nonlinear) equations. 
 */
template<class R>
class SolverBase
  : public SolverInterface<R>
{
  typedef typename traits<R>::interval_type I;
 public:
  /*! \brief Constructor. */
  SolverBase(R max_error, uint max_steps);
  
  /*! \brief The maximum permissible error of the solution. */
  R maximum_error() const { return this->_max_error; }
  /*! \brief Set the maximum error. */
  void set_maximum_error(R max_error) { this->_max_error=max_error; };
  
  /*! \brief The maximum number of steps allowed before the method must quit. */
  uint maximum_number_of_steps() const { return this->_max_steps; }
  /*! \brief Set the maximum number of steps. */
  void set_maximum_number_of_steps(uint max_steps) { this->_max_steps=max_steps; };

  /*! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt. */
  virtual Point<I> fixed_point(const Map<R>& f,const Point<I>& pt);
 private:
  R _max_error;
  uint _max_steps;
};

template<class R>
inline
SolverBase<R>::SolverBase(R max_error, uint max_steps)
  : _max_error(max_error), _max_steps(max_steps) 
{
}


template<class R> 
Point<typename SolverBase<R>::I> 
SolverBase<R>::fixed_point(const Map<R>& f,const Point<I>& pt) 
{
  return Point<I>(this->solve(DifferenceFunction<R>(f.function()),pt.position_vector())); 
}


} // namespace Ariadne

#endif /* ARIADNE_SOLVER_INTERFACE_H */
