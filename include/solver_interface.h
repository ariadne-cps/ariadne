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

typedef unsigned int uint;

namespace Ariadne {

class Float;
class Interval;
class Real;

template<class T> class Set;
template<class T> class List;
template<class X> class Vector;

template<class X> class ScalarFunction;
typedef ScalarFunction<Real> RealScalarFunction;
typedef ScalarFunction<Interval> IntervalScalarFunction;
template<class X> class VectorFunction;
typedef VectorFunction<Real> RealVectorFunction;
typedef VectorFunction<Interval> IntervalVectorFunction;

template<class X> class ScalarFunctionModel;
typedef ScalarFunctionModel<Interval> IntervalScalarFunctionModel;
template<class X> class VectorFunctionModel;
typedef VectorFunctionModel<Interval> IntervalVectorFunctionModel;

class SolverInterface;

//! \relates SolverInterface \brief An exception occurring in the solution of a (parameterised) algebraic equation.
class SolverException : public std::runtime_error
{
  public:
    SolverException(const char* what) : std::runtime_error(what) { }
    SolverException(const std::string& what) : std::runtime_error(what) { }
};

//! \relates SolverInterface \brief It cannot be shown that there is no solution to a (parameterised) algebraic equation,
//! \brief but no solution can be found by the solver.
class UnknownSolutionException : public SolverException
{
  public:
    UnknownSolutionException(const char* what) : SolverException(what) { }
    UnknownSolutionException(const std::string& what) : SolverException(what) { }
};

//! \relates SolverInterface \brief There are no solutions to a (parameterised) algebraic equation in the domain.
class NoSolutionException : public SolverException
{
  public:
    NoSolutionException(const char* what) : SolverException(what) { }
    NoSolutionException(const std::string& what) : SolverException(what) { }
};

//! \relates SolverInterface \brief The solutions to a (parameterised) algebraic equation are degenerate due to the Jacobian being singular.
class SingularJacobianException : public SolverException
{
  public:
    SingularJacobianException(const char* what) : SolverException(what) { }
    SingularJacobianException(const std::string& what) : SolverException(what) { }
};


//! \ingroup SolverModule EvaluationModule
//! \brief %Interface for solving (nonlinear) equations.
//! \sa SolverException
class SolverInterface
    : public Loggable
{
  public:
    //! \brief Virtual destructor.
    virtual ~SolverInterface() { };
    //! \brief Make a dynamically-allocated copy.
    virtual SolverInterface* clone() const = 0;
    //! \brief Write to an output stream.
    virtual void write(std::ostream& os) const = 0;


    //! \brief The maximum permissible error of the solution.
    virtual double maximum_error() const = 0;
    //! \brief Set the maximum error.
    virtual void set_maximum_error(double max_error) = 0;

    //! \brief The maximum number of steps allowed before the method must quit.
    virtual uint maximum_number_of_steps() const = 0;
    //! \brief Set the maximum number of steps.
    virtual void set_maximum_number_of_steps(uint max_steps) = 0;

    //! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt.
    virtual Vector<Interval> zero(const IntervalVectorFunction& f,const Vector<Interval>& pt) const = 0;
    //! \brief Solve \f$f(x)=x\f$, starting in the interval point \a pt.
    virtual Vector<Interval> fixed_point(const IntervalVectorFunction& f,const Vector<Interval>& pt) const = 0;

    //! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt. Throws a SolverException if there is not a unique solution.
    virtual Vector<Interval> solve(const IntervalVectorFunction& f,const Vector<Interval>& pt) const = 0;
    //! \brief Solve \f$f(a,x)=0\f$ for \f$a\f$ in \f$A\f$, looking for a solution with \f$x\f$ in \f$X\f$. The result is a function \f$h\f$ with domain \f$A\f$ such that \f$f(a,h(a))=0\f$ for all \f$a\in A\f$.
    //! \details A <em>strong solution</em> to the problem exists if for every \f$a\in A\f$, there is a unique solution \f$x\in X\f$ to the equation \f$f(a,x)=0\f$, and this solution varies continuously in \f$a\f$. A <em>weak solution</em> is a function \f$h\f$ with domain \f$A\f$ such that \f$f(a,h(a))=0\f$ for all \f$a\in A\f$, if \f$f(a,x)=0\f$ for \f$x\in X\f$ then \f$h(a)=x\f$, and there exists \f$(a,x)\in A\times X\f$ such that \f$f(a,x)=0\f$.
    //! If \f$D_2f(a,x)\f$ is nonsingular for all \f$(a,x)\in A\times X\f$, then any solution is guaranteed to be unique.
    //! May throw a NoSolutionException, but \em only if there are no solutions to \f$f(a,x)=0\f$ in \f$A\times X\f$.
    //! If there is a continuous branch of solutions \f$x=h(a)\f$ such that \f$h(a)\in X\f$ for some parameter values,
    //! but \f$h(a)\not\in X\f$ for others, then \f$h\f$ is a valid result for the function.
    virtual IntervalVectorFunctionModel implicit(const IntervalVectorFunction& f, const Vector<Interval>& A, const Vector<Interval>& X) const = 0;
    //! \brief Solve \f$f(a,x)=0\f$ for a in \a A, looking for a solution with x in \a X.
    virtual IntervalScalarFunctionModel implicit(const IntervalScalarFunction& f, const Vector<Interval>& A, const Interval& X) const = 0;

    //! \brief Solve \f$f(a,x)=0\f$ for x in \a X, and continue to a function \f$h\f$ solving \f$f(a,h(a))=0\f$ over \f$A\f$.
    virtual IntervalVectorFunctionModel continuation(const IntervalVectorFunction& f, const Vector<Float>& a, const Vector<Interval>& X,  const Vector<Interval>& A) const = 0;

    //! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt. Returns a set of boxes for which it can be <em>proved</em> that
    //! a solution exists in the box. This means that some solutions may be omitted if they are not sufficiently robust.
    //!
    virtual Set< Vector<Interval> > solve_all(const IntervalVectorFunction& f,const Vector<Interval>& pt) const = 0;
};

inline std::ostream& operator<<(std::ostream& os, const SolverInterface& solver) {
    solver.write(os); return os;
}

} // namespace Ariadne

#endif /* ARIADNE_SOLVER_INTERFACE_H */
