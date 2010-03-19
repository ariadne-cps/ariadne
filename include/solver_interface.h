/***************************************************************************
 *            solver_interface.h
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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

/*! \file solver_interface.h
 *  \brief Interface class for solving equations.
 */

#ifndef ARIADNE_SOLVER_INTERFACE_H
#define ARIADNE_SOLVER_INTERFACE_H

#include <exception>
#include <stdexcept>
#include <string>

#include "logging.h"

namespace Ariadne {

typedef unsigned int uint;

class Interval;
template<class T> class Set;
template<class T> class List;
template<class X> class Vector;
class VectorFunction;
class VectorTaylorFunction;
class ScalarFunction;
class ScalarTaylorFunction;

class SolverException : public std::runtime_error
{
  public:
    SolverException(const char* what) : std::runtime_error(what) { }
    SolverException(const std::string& what) : std::runtime_error(what) { }
};

class NoSolutionException : public SolverException
{
  public:
    NoSolutionException(const char* what) : SolverException(what) { }
    NoSolutionException(const std::string& what) : SolverException(what) { }
};



/*! \ingroup EvaluatorInterfaces \ingroup Solvers
 *  \brief %Interface for solving (nonlinear) equations.
 */
class SolverInterface
    : public Loggable
{
  public:
    /*! \brief Virtual destructor. */
    virtual ~SolverInterface() { };
    /*! \brief Make a dynamically-allocated copy. */
    virtual SolverInterface* clone() const = 0;

    /*! \brief The maximum permissible error of the solution. */
    virtual double maximum_error() const = 0;
    /*! \brief Set the maximum error. */
    virtual void set_maximum_error(double max_error) = 0;

    /*! \brief The maximum number of steps allowed before the method must quit. */
    virtual uint maximum_number_of_steps() const = 0;
    /*! \brief Set the maximum number of steps. */
    virtual void set_maximum_number_of_steps(uint max_steps) = 0;

    /*! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt. */
    virtual Vector<Interval> zero(const VectorFunction& f,const Vector<Interval>& pt) const = 0;
    /*! \brief Solve \f$f(x)=x\f$, starting in the interval point \a pt. */
    virtual Vector<Interval> fixed_point(const VectorFunction& f,const Vector<Interval>& pt) const = 0;

    /*! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt. Throws a SolverException if there is not a unique solution. */
    virtual Vector<Interval> solve(const VectorFunction& f,const Vector<Interval>& pt) const = 0;
    /*! \brief Solve \f$f(a,x)=0\f$ for a in \a par, looking for a solution with x in \a ix. */
    virtual VectorTaylorFunction implicit(const VectorFunction& f, const Vector<Interval>& par, const Vector<Interval>& ix) const = 0;
    /*! \brief Solve \f$f(a,x)=0\f$ for a in \a par, looking for a solution with x in \a ix. */
    virtual ScalarTaylorFunction implicit(const ScalarFunction& f, const Vector<Interval>& par, const Interval& ix) const = 0;

    /*! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt. Returns a set of boxes for which it can be <em>proved</em> that
     *  a solution exists in the box. This means that some solutions may be omitted if they are not sufficiently robust.
     */
    virtual Set< Vector<Interval> > solve_all(const VectorFunction& f,const Vector<Interval>& pt) const = 0;
};


} // namespace Ariadne

#endif /* ARIADNE_SOLVER_INTERFACE_H */
