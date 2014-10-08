/***************************************************************************
 *            solver.h
 *
 *  Copyright  2006-9  Pieter Collins
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

/*! \file solver.h
 *  \brief Solver classes for algebraic equations.
 */

#ifndef ARIADNE_SOLVER_H
#define ARIADNE_SOLVER_H

#include <exception>
#include <stdexcept>
#include <string>

#include "solvers/solver_interface.h"
#include "function/function_interface.h"

#include "utility/logging.h"
#include "utility/attribute.h"
#include "utility/pointer.h"
#include "utility/container.h"
#include "numeric/numeric.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class FunctionModelFactoryInterface;

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

    /*! \brief The class which constructs functions for the implicit function solver. */
    const FunctionModelFactoryInterface<ValidatedTag>& function_factory() const;
    /*! \brief Set the class which constructs functions for the implicit function solver. */
    void set_function_factory(const FunctionModelFactoryInterface<ValidatedTag>& factory);

    /*! \brief Solve \f$f(x)=0\f$, starting in the box \a bx. */
    virtual Vector<ValidatedNumber> zero(const ValidatedVectorFunction& f,const ExactBox& bx) const;
    /*! \brief Solve \f$f(x)=0\f$, starting in the box \a bx. */
    virtual Vector<ValidatedNumber> fixed_point(const ValidatedVectorFunction& f,const ExactBox& bx) const;

    /*! \brief Solve \f$f(x)=0\f$, starting in the box \a bx. */
    virtual Vector<ValidatedNumber> solve(const ValidatedVectorFunction& f,const ExactBox& bx) const;
    virtual Vector<ValidatedNumber> solve(const ValidatedVectorFunction& f,const Vector<ValidatedNumber>& ipt) const;
    /*! \brief Solve \f$f(a,x)=0\f$ for a in \a par, looking for a solution with x in \a ix. */
    virtual ValidatedVectorFunctionModel implicit(const ValidatedVectorFunction& f, const ExactBox& par, const ExactBox& ix) const;
    /*! \brief Solve \f$f(a,x)=0\f$ for a in \a par, looking for a solution with x in \a ix. */
    virtual ValidatedScalarFunctionModel implicit(const ValidatedScalarFunction& f, const ExactBox& par, const ExactInterval& ix) const;
    //! \brief Solve \f$f(a,x)=0\f$ yielding a function \f$x=h(a)\f$ for a in \a A, looking for a solution with \f$h(A) \subset X\f$ and $h(a)\in x\f$.
    virtual ValidatedVectorFunctionModel continuation(const ValidatedVectorFunction& f, const Vector<ApproximateNumber>& a, const ExactBox& X, const ExactBox& A) const;


    /*! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt. */
    virtual Set< Vector<ValidatedNumber> > solve_all(const ValidatedVectorFunction& f,const ExactBox& bx) const;
  protected:
    /*! \brief Perform one iterative step of the contractor. */
    virtual Vector<ValidatedNumber> step(const ValidatedVectorFunction& f,const Vector<ValidatedNumber>& pt) const = 0;
    /*! \brief Perform one iterative step of the contractor. */
    virtual ValidatedVectorFunctionModel implicit_step(const ValidatedVectorFunction& f,const ValidatedVectorFunctionModel& p,const ValidatedVectorFunctionModel& x) const = 0;
  private:
    double _max_error;
    uint _max_steps;
    std::shared_ptr< FunctionModelFactoryInterface<ValidatedTag> > _function_factory_ptr;
};


/*! \ingroup Solvers
 *  \brief ExactInterval Newton solver. Uses the contractor \f$[x']=x_0-Df^{-1}([x])f(x_0)\f$.
 */
class IntervalNewtonSolver
    : public SolverBase
{
  public:
    /*! \brief Constructor. */
    IntervalNewtonSolver(double max_error, uint max_steps) : SolverBase(max_error,max_steps) { }
    IntervalNewtonSolver(MaximumError max_error, MaximumNumberOfSteps max_steps) : SolverBase(max_error,max_steps) { }
    /*! \brief Cloning operator. */
    virtual IntervalNewtonSolver* clone() const { return new IntervalNewtonSolver(*this); }
    /*! \brief Write to an output stream. */
    virtual void write(std::ostream& os) const;

    using SolverBase::implicit;
  public:
    virtual ValidatedVectorFunctionModel implicit_step(const ValidatedVectorFunction& f, const ValidatedVectorFunctionModel& p, const ValidatedVectorFunctionModel& x) const;

    virtual Vector<ValidatedNumber> step(const ValidatedVectorFunction& f, const Vector<ValidatedNumber>& pt) const;
};


/*! \ingroup Solvers
 *  \brief Krawczyk solver. Uses the contractor \f$[x']=x_0-Mf(x_0)+(I-MDf([x])([x]-x_0)\f$
 *  where \f$M\f$ is typically \f$Df^{-1}(x_0)\f$.
 */
class KrawczykSolver
    : public SolverBase
{
  public:
    /*! \brief Constructor. */
    KrawczykSolver(double max_error, uint max_steps) : SolverBase(max_error,max_steps) { }
    KrawczykSolver(MaximumError max_error, MaximumNumberOfSteps max_steps) : SolverBase(max_error,max_steps) { }
    /*! \brief Cloning operator. */
    virtual KrawczykSolver* clone() const { return new KrawczykSolver(*this); }
    /*! \brief Write to an output stream. */
    virtual void write(std::ostream& os) const;

    /*! \brief Solve \f$f(a,x)=0\f$ for a in \a par, looking for solutions with x in \a ix. */
    virtual ValidatedVectorFunctionModel implicit_step(const ValidatedVectorFunction& f, const ValidatedVectorFunctionModel& p, const ValidatedVectorFunctionModel& x) const;

  public:
    /*! \brief A single step of the Krawczyk contractor. */
    virtual Vector<ValidatedNumber>
    step(const ValidatedVectorFunction& f,
          const Vector<ValidatedNumber>& pt) const;
};


/*! \ingroup Solvers
 *  \brief Modified Krawczyk solver. Uses the contractor \f$[x']=x_0-J^{-1}(f(x_0)+(J-Df([x])([x]-x_0)\f$
 *  where \f$J\f$ is typically \f$Df(x_0)\f$.
 */
class FactoredKrawczykSolver
    : public KrawczykSolver
{
  public:
    /*! \brief Constructor. */
    FactoredKrawczykSolver(double max_error, uint max_steps) : KrawczykSolver(max_error,max_steps) { }
    FactoredKrawczykSolver(MaximumError max_error, MaximumNumberOfSteps max_steps) : KrawczykSolver(max_error,max_steps) { }
    /*! \brief Cloning operator. */
    virtual FactoredKrawczykSolver* clone() const { return new FactoredKrawczykSolver(*this); }
    /*! \brief Write to an output stream. */
    virtual void write(std::ostream& os) const;
  public:
    /*! \brief A single step of the modified Krawczyk contractor. */
    virtual Vector<ValidatedNumber>
    step(const ValidatedVectorFunction& f,
          const Vector<ValidatedNumber>& pt) const;
};



} // namespace Ariadne

#endif /* ARIADNE_SOLVER_H */
