/***************************************************************************
 *            nonlinear_programming.h
 *
 *  Copyright 2009  Pieter Collins
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

/*! \file nonlinear_programming.h
 *  \brief Nonlinear programming.
 */

#ifndef ARIADNE_NONLINEAR_PROGRAMMING_H
#define ARIADNE_NONLINEAR_PROGRAMMING_H

#include "logging.h"
#include "numeric.h"
#include "tuple.h"

using namespace boost::numeric;

namespace Ariadne {

template<class X> class Vector;
template<class X> class ScalarFunction;
typedef ScalarFunction<Real> RealScalarFunction;
template<class X> class VectorFunction;
typedef VectorFunction<Real> RealVectorFunction;

class NonlinearConstraint;

class InfeasibleProblemException : public std::exception { };
class DegenerateFeasibilityProblemException : public std::exception { };
class NearBoundaryOfFeasibleDomainException : public std::exception { };

//! \ingroup OptimisationModule EvaluationModule
//! Interface for nonlinear programming solvers.
class OptimiserInterface {
  public:
    typedef tribool Tribool;
    typedef Vector<Float> FloatVector;
    typedef Vector<Interval> IntervalVector;
    typedef Vector<Interval> BoxType;
  public:
    //! \brief Virtual destructor.
    virtual ~OptimiserInterface() { }
    //! \brief Create a dynamically-allocated copy.
    virtual OptimiserInterface* clone() const = 0;
    //! \brief Solve the linear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    virtual IntervalVector optimise(RealScalarFunction f, IntervalVector D, RealVectorFunction g, IntervalVector C) const = 0;

    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual tribool feasible(IntervalVector D, RealVectorFunction g, IntervalVector C) const = 0;

    //! \brief Tests if the point \a y is feasible.
    virtual bool is_feasible_point(IntervalVector D, RealVectorFunction g, IntervalVector C,
                                   FloatVector y) const = 0;
    //! \brief Tests if the Lagrange multipliers \a x are a certificate of infeasiblity.
    virtual bool is_infeasibility_certificate(IntervalVector D, RealVectorFunction g, IntervalVector C,
                                              FloatVector x) const = 0;
};

//! \ingroup OptimisationModule
//! Common routines for nonlinear optimisation
class OptimiserBase
    : public OptimiserInterface
    , public Loggable
{
  public:
    //! \brief Tests if the point \a y is feasible.
    virtual bool is_feasible_point(IntervalVector D, RealVectorFunction g, IntervalVector C,
                                   FloatVector y) const;
    //! \brief Tests if the Lagrange multipliers \a x are a certificate of infeasiblity.
    virtual bool is_infeasibility_certificate(IntervalVector D, RealVectorFunction g, IntervalVector C,
                                              FloatVector x) const;
};

    //! \ingroup OptimisationModule
//! Solver for linear programming problems using interior point methods.
class NonlinearInteriorPointOptimiser
    : public OptimiserBase
{
  public:
    virtual NonlinearInteriorPointOptimiser* clone() const { return new NonlinearInteriorPointOptimiser(*this); }

    //! \brief Solve the linear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    virtual IntervalVector optimise(RealScalarFunction f, IntervalVector D, RealVectorFunction g, IntervalVector C) const;
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual tribool feasible(IntervalVector D, RealVectorFunction g, IntervalVector C) const;

    //! \brief Test feasibility of the dual problem \f$f(y)\leq c;  y\in D\f$.
    //! Returns the pair (r,y,g) where r is the result, and y the (potential) feasible point,
    //! and g is a positive function on such that \f$g(y)\cdot f(y) \leq k\f$.
    tuple<tribool,FloatVector, RealVectorFunction >
    constrained_dual_feasible(const RealVectorFunction& f, const FloatVector& c, const IntervalVector& D) const;

    //! \brief Test if the constraints \f$g(y)\in C\f$ are solvable for \f$y\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the primal, dual and slack variables.
    Pair<Tribool,FloatVector> feasible(IntervalVector D, RealVectorFunction g, IntervalVector C,
                                       const Float& t0, const FloatVector& x0, const FloatVector& y0, const FloatVector& z0) const;

    void feasibility_step(const RealVectorFunction& g, FloatVector& x, FloatVector& y, FloatVector& z, Float& t) const;
    void optimization_step(const RealScalarFunction& f, const RealVectorFunction& g, FloatVector& x, FloatVector& y, FloatVector& z) const;
    void setup_feasibility(const IntervalVector& d, const RealVectorFunction& g, const IntervalVector& c,
                           FloatVector& x, FloatVector& y, FloatVector& z, Float& t) const;
    void feasibility_step(const IntervalVector& d, const RealVectorFunction& g, const IntervalVector& c,
                                         FloatVector& x, FloatVector& y, FloatVector& z, Float& t) const;
    void linearised_feasibility_step(const IntervalVector& d, const RealVectorFunction& g, const IntervalVector& c,
                                         FloatVector& x, FloatVector& y, FloatVector& z, Float& t) const;

  public:
    void compute_tz(const IntervalVector& d, const RealVectorFunction& g, const IntervalVector& b, const FloatVector& y, Float& t, FloatVector& z) const;
    void compute_z(const IntervalVector& d, const RealVectorFunction& g, const IntervalVector& b, const FloatVector& y, const Float& t, FloatVector& z) const;
  protected:
/*
    //! \brief Find approximate optimal solution of \f$\min c^T x \text{ s.t. } Ax=b; x\geq0\f$.
    //! Returns the triple (x,y,z) where x is the optimal point, and y the corresponding dual feasible point.
    tuple< FloatVector, FloatVector, FloatVector >
    _optimize(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
              Vector<Float>& x, Vector<Float>& y, Vector<Float>& z) const;
*/

};

//! \ingroup OptimisationModule
//! Solver for linear programming problems using interior point methods.
//! WARNING: This class currently does not work; maybe there is a problem with the algorithms.
class KrawczykOptimiser
    : public OptimiserBase
{

  public:
    virtual KrawczykOptimiser* clone() const { return new KrawczykOptimiser(*this); }

    //! \brief Solve the linear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    virtual IntervalVector optimise(RealScalarFunction f, IntervalVector D, RealVectorFunction g, IntervalVector C) const;
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual tribool feasible(IntervalVector D, RealVectorFunction g, IntervalVector C) const;

  public:
    //! \brief Try to solve the nonlinear constraint problem by applying the Krawczyk contractor to the Kuhn-Tucker conditions,
    //! hotstarting the iteration with the primal and dual variables.
    tribool optimise(RealScalarFunction f, IntervalVector D, RealVectorFunction g, IntervalVector C,
                     const Interval& t0, const IntervalVector& x0, const IntervalVector& y0, const IntervalVector& z0) const;

    //! \brief A primal-dual feasibility step for the problem \f$g(y)\in C;\ y\in D\f$.
    void optimisation_step(const IntervalVector& D, const RealVectorFunction& g, const IntervalVector& C,
                           IntervalVector& x, IntervalVector& y, IntervalVector& z, Interval& t) const;
    //! \brief A primal-dual feasibility step for the problem \f$g(y)\in C;\ y\in D\f$.
    void feasibility_step(const IntervalVector& D, const RealVectorFunction& g, const IntervalVector& C,
                          IntervalVector& x, IntervalVector& y, IntervalVector& z, Interval& t) const;

    //! \brief A primal feasibility step for the problem \f$g(y)\in C;\ y\in D\f$. \deprecated
    void feasibility_step(const IntervalVector& D, const RealVectorFunction& g, const IntervalVector& C,
                          IntervalVector& y, Interval& t) const;
    //! \brief A feasibility step for the problem \f$g(y)\leq 0\f$. \deprecated
    void feasibility_step(const RealVectorFunction& g,
                          IntervalVector& x, IntervalVector& y, IntervalVector& z, Interval& t) const;
    //! \brief An optimization step for the problem \f$\max f(y) \text{ s.t. } g(y)\leq 0\f$. \deprecated
    void optimisation_step(const RealScalarFunction& f, const RealVectorFunction& g,
                           IntervalVector& x, IntervalVector& y, IntervalVector& z) const;
  protected:
    void setup_feasibility(const IntervalVector& d, const RealVectorFunction& g, const IntervalVector& c,
                           IntervalVector& x, IntervalVector& y, IntervalVector& z, Interval& t) const;
    protected:
    void compute_tz(const IntervalVector& d, const RealVectorFunction& g, const IntervalVector& b, const IntervalVector& y, Interval& t, IntervalVector& z) const;
};


} // namespace Ariadne

#endif
