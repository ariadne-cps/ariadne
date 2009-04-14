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

#include "solver_interface.h"
#include "function_interface.h"

#include "numeric.h"
#include "vector.h"
#include "pointer.h"
#include "function.h"

namespace Ariadne {

    
/*! \ingroup \ingroup Solvers
 *  \brief %Common functionality for solving (nonlinear) equations. 
 */
class SolverBase
  : public SolverInterface
{
 public:
  /*! \brief Constructor. */
  SolverBase(double max_error, uint max_steps);
  
  /*! \brief The maximum permissible error of the solution. */
  double maximum_error() const { return this->_max_error; }
  /*! \brief Set the maximum error. */
  void set_maximum_error(double max_error) { this->_max_error=max_error; };
  
  /*! \brief The maximum number of steps allowed before the method must quit. */
  uint maximum_number_of_steps() const { return this->_max_steps; }
  /*! \brief Set the maximum number of steps. */
  void set_maximum_number_of_steps(uint max_steps) { this->_max_steps=max_steps; };

  /*! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt. */
  virtual Vector<Interval> fixed_point(const FunctionInterface& f,const Vector<Interval>& pt);
 private:
  double _max_error;
  uint _max_steps;
};


class DifferenceFunction
    : public FunctionTemplate<DifferenceFunction>
{
  public:
    DifferenceFunction(const FunctionInterface& f) : fptr(f.clone()) { }
    virtual DifferenceFunction* clone() const { return new DifferenceFunction(*this); }
    virtual uint result_size() const { return fptr->result_size(); }
    virtual uint argument_size() const { return fptr->argument_size(); }
    virtual ushort smoothness() const { return fptr->smoothness(); }
    template<class Res, class Args> void _compute(Res& r, const Args& a) const { r=fptr->evaluate(a)-a; }
    template<class Res, class Args> void _compute_approx(Res& r, const Args& a) const { _compute(r,a); }
  private:
    boost::shared_ptr<FunctionInterface> fptr;
};


inline
SolverBase::SolverBase(double max_error, uint max_steps)
  : _max_error(max_error), _max_steps(max_steps) 
{
}



Vector<Interval> 
SolverBase::fixed_point(const FunctionInterface& f, const Vector<Interval>& pt) 
{
  return Vector<Interval>(this->solve(DifferenceFunction(f),pt)); 
}


} // namespace Ariadne

#endif /* ARIADNE_SOLVER_BASE_H */
