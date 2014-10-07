/***************************************************************************
 *            constraint_solver.h
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

/*! \file constraint_solver.h
 *  \brief Class for solving systems of constraints over a box.
 */

#ifndef ARIADNE_CONSTRAINT_SOLVER_H
#define ARIADNE_CONSTRAINT_SOLVER_H

#include "logging.h"
#include "container.h"

#include "tribool.h"
#include "numeric.h"
#include "constraint.h"

namespace Ariadne {

class Real;
class Interval;
typedef tribool Tribool;

template<class X> class Point;
typedef Point<ExactNumberType> ExactPoint;
class Box;
class UpperBox;
class GridTreeSet;

template<class X, class R> class Constraint;
typedef Constraint<EffectiveScalarFunction,EffectiveNumberType> EffectiveConstraint;
typedef Constraint<ValidatedScalarFunction,ValidatedNumberType> ValidatedConstraint;

template<class X> class Procedure;
typedef Procedure<ValidatedTag> ValidatedProcedure;

template<class X> class ScalarFunctionInterface;
typedef ScalarFunctionInterface<ValidatedTag> ValidatedScalarFunctionInterface;
template<class X> class ScalarFunction;
typedef ValidatedScalarFunction ValidatedScalarFunction;
template<class X> class VectorFunctionInterface;
typedef VectorFunctionInterface<ValidatedTag> ValidatedVectorFunctionInterface;
template<class X> class VectorFunction;
typedef ValidatedVectorFunction ValidatedVectorFunction;

class VectorTaylorFunction;

template<class X> struct FeasibilityState {
    X t;
    Vector<X> x;
    Vector<X> y;
    Vector<X> z;
};


//! \ingroup EvaluationModule OptimisationModule
//! \brief A class for finding solutions of systems of constraints of the form \f$g(y) \leq c\f$.
class ConstraintSolverInterface {
  public:
    //! \brief Test if the image of the box \a domain under the function \a function intersects \a codomain.
    virtual Pair<Tribool,ExactPoint> feasible(const Box& domain, const ValidatedVectorFunction& function, const Box& codomain) const = 0;
    //! \brief Test if \a point is in \a domain and the image of \a point under the function \a function lies in \a codomain.
    virtual Tribool check_feasibility(const Box& domain, const ValidatedVectorFunction& function, const Box& codomain, const ExactPoint& point) const = 0;
    //! \brief Try to reduce the size of the domain by propagating interval constraints. Returns \c true if the reduced domain is empty.
    virtual bool reduce(UpperBox& domain, const ValidatedVectorFunction& function, const Box& codomain) const = 0;

};



//! \ingroup OptimisationModule
//! \brief A class for finding solutions of systems of constraints of the form \f$g(y) \leq c\f$.
class ConstraintSolver
    : public ConstraintSolverInterface, public Loggable
{
  public:
    //! \brief Test if the image of the box \a domain under the function \a function intersects \a codomain.
    virtual Pair<Tribool,ExactPoint> feasible(const Box& domain, const ValidatedVectorFunction& function, const Box& codomain) const;
    //! \brief Test if \a point is in \a domain and the image of \a point under the function \a function lies in \a codomain.
    virtual Tribool check_feasibility(const Box& domain, const ValidatedVectorFunction& function, const Box& codomain, const ExactPoint& point) const;
    //! \brief Try to reduce the size of the domain by propagating interval constraints.
    virtual bool reduce(UpperBox& domain, const ValidatedVectorFunction& function, const Box& codomain) const;


    //! \brief Test if the constraints are solvable using a nonlinear feasibility test. Returns an approximate feasible point if the result is true. (Deprecated)
    virtual Pair<Tribool,ExactPoint> feasible(const Box& domain, const List<ValidatedConstraint>& constraints) const;
    //! \brief Try to reduce the size of the domain by propagating interval constraints. (Deprecated)
    virtual bool reduce(UpperBox& domain, const List<ValidatedConstraint>& constraints) const;

    //! \brief Try to enforce hull consistency by propagating several interval constraints at once.
    //! This method is sharp if each variable occurs at most once in the constraint.
    bool hull_reduce(UpperBox& bx, const ValidatedVectorFunctionInterface& function, const Box& codomain) const;
    bool hull_reduce(UpperBox& bx, const Vector<ValidatedProcedure>& procedure, const Box& codomain) const;
    //! \brief Try to enforce hull consistency by propagating an interval constraint.
    //! This method is sharp if each variable occurs at most once in the constraint.
    bool hull_reduce(UpperBox& bx, const ValidatedScalarFunctionInterface& function, const Interval& codomain) const;
    bool hull_reduce(UpperBox& bx, const ValidatedProcedure& procedure, const Interval& codomain) const;

    //! \brief Reduce the \a domain by testing intersection of \a multipliers inner product \a function(\a domain)
    //! with \a multipliers innner product \a codomain, centering at \a centre.
    //! Reduces \f$(\lambda\cdot f)(X) \cap (\lambda\cdot C)\f$, evaluating \f$g(x)=g(x^*)+Dg(X) (X-x^*)\f$.
    bool lyapunov_reduce(UpperBox& domain, const VectorTaylorFunction& function, const Box& codomain,
                         Vector<ExactFloatType> centre, Vector<ExactFloatType> multpliers) const;
    bool lyapunov_reduce(UpperBox& domain, const VectorTaylorFunction& function, const Box& codomain,
                         Vector<ApproximateNumberType> centre, Vector<ApproximateNumberType> multpliers) const;
    //! \brief Try to enforce hull consistency by reducing a constraint with respect to one variable.
    bool box_reduce(UpperBox& bx, const ValidatedScalarFunctionInterface& function, const Interval&, uint j) const;
    //! \brief Try to enforce hull consistency by reducing an a monotone dimension.
    //! This method is sharp if each variable occurs at most once in the constraint.
    bool monotone_reduce(UpperBox& bx, const ValidatedScalarFunctionInterface& function, const Interval&, uint j) const;

    //! Split the domain into two pieces to help try to solve the constraints.
    Pair<UpperBox,UpperBox> split(const UpperBox& domain, const ValidatedVectorFunction& function, const Box& codomain) const;

    // Deprecated functions.
    bool hull_reduce(UpperBox& bx, const ValidatedConstraint& constraint) const {
        return this->hull_reduce(bx,constraint.function(),constraint.bounds()); }
    bool box_reduce(UpperBox& bx, const ValidatedConstraint& constraint, uint j) const {
        return this->box_reduce(bx,constraint.function(),constraint.bounds(),j); }
    bool monotone_reduce(UpperBox& bx, const ValidatedConstraint& constraint, uint j) const {
        return this->monotone_reduce(bx,constraint.function(),constraint.bounds(),j); }

};


} //namespace Ariadne

#endif
