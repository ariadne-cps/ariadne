/***************************************************************************
 *            linear_programming.h
 *
 *  Copyright 2005-8  Alberto Casagrande, Pieter Collins
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

/*! \file linear_programming.h
 *  \brief Linear programming.
 */

#ifndef ARIADNE_LINEAR_PROGRAMMING_H
#define ARIADNE_LINEAR_PROGRAMMING_H

#include "vector.h"
#include "matrix.h"
#include "numeric.h"
#include "tuple.h"

using namespace boost::numeric;

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Affine;

struct StandardLinearProgram {
    StandardLinearProgram(uint m, uint n);
    Matrix<Float> A;
    Vector<Float> b;
    Vector<Float> c;
    Vector<Float> x;
    Vector<Float> y;
    Vector<Float> z;
};

StandardLinearProgram
feasibility_problem(const Matrix<Float>& B, const Vector<Float>& bl, const Vector<Float>& bu,
                    const Matrix<Float>& C, const Vector<Float>& c,
                    const Vector<Float>& dl, const Vector<Float>& du);

StandardLinearProgram
feasibility_problem(const Vector< Affine<Float> >& f, const Vector<Interval>& b,
                    const Vector< Affine<Float> >& g, const Vector<Float>& c,
                    const Vector<Interval>& d);



//! \ingroup OptimisationModule
//! Solver for linear programming problems using interior point methods.
class InteriorPointSolver {
  public:
    //! \brief Test feasibility of the problem \f$b\leq A^Ty\leq c; l\leq y\leq u\f$.
    //! Returns the pair (r,y) where r is the result, and y the (potential) feasible point.
    //!
    //! \internal We do not also allow constraints of the form Ay=b since these should be removed
    //! before starting the problem.
    tuple<tribool,Vector<Float> >
    constrained_dual_feasible(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
                              const Vector<Float>& l, const Vector<Float>& u) const;


    //! \brief Find approximate optimal solution of \f$\min c^T x \text{ s.t. } Ax=b; x\geq0\f$.
    //! Returns the pair (x,y) where x is the optimal point, and y the corresponding dual feasible point.
    tuple< Vector<Float>, Vector<Float>, Vector<Float> >
    optimize(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c) const;

    //! \brief Determine whether the problem \f$c_l \leq A^Ty \leq c_u;\ l\leq y\leq u\f$ is feasible.
    tribool
    feasible(const Matrix<Float>& A, const Vector<Float>& cl, const Vector<Float>& cu,
             const Vector<Float>& l, const Vector<Float>& u) const;

    //! \brief Determine whether the problem \f$b_l \leq B^Ty \leq b_u;\ C^Ty \leq c; \ d_l\leq y\leq d_u\f$ is feasible.
    tribool
    feasible(const Matrix<Float>& B, const Vector<Float>& bl, const Vector<Float>& bu,
             const Matrix<Float>& C, const Vector<Float>& c,
             const Vector<Float>& dl, const Vector<Float>& du) const;

    //! \brief Validate that \a x is (approximately) primal feasible and \a y is dual feasible.
    //! Returns the interval of possible optimal values.
    Interval validate(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c, const Vector<Float>& x, const Vector<Float>& y) const;
  public:
    tribool
    _feasible(const Matrix<Float>& A, const Vector<Float>& cl, const Vector<Float>& cu,
              const Vector<Float>& l, const Vector<Float>& u,
              Float& t, Vector<Float>& x, Vector<Float>& y, Vector<Float>& z) const;

    //! \brief Find approximate optimal solution of \f$\min c^T x \text{ s.t. } Ax=b; x\geq0\f$.
    //! Returns the triple (x,y,z) where x is the optimal point, and y the corresponding dual feasible point.
    tuple< Vector<Float>, Vector<Float>, Vector<Float> >
    _optimize(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
              Vector<Float>& x, Vector<Float>& y, Vector<Float>& z) const;

    //! \brief Perform a step of the optimization of \f$\min c^T x \text{ s.t. } Ax=b; x\geq0\f$.
    void
    _optimization_step(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
                       Vector<Float>& x, Vector<Float>& y, Vector<Float>& z) const;

};



//! \relates SimplexSolver \brief The type of variable; lower bounded, upper bounded, or basic.
enum VariableType { LOWER=-1, BASIS=0, UPPER=+1 };
std::ostream& operator<<(std::ostream& os, VariableType t);

//! \ingroup OptimisationModule
//! Solver for linear programming problems using the simplex algorithm.
template<class X>
class SimplexSolver {
  public:

    //! \ingroup LinearProgrammingModule
    //! Solve the linear programming problem \f$\min cx \text{ s.t. } Ax=b;\ x\geq0\f$.
    //! Returns the optimal vector. Throws an error if the problem is infeasible.
    Vector<X>
    optimize(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c);

    //! \ingroup LinearProgrammingModule
    //! Solve the linear programming problem \f$\min cx \text{ s.t. } Ax=b;\ l\leq x\leq u\f$.
    //! Returns the optimal vector. Throws an error if the problem is infeasible.
    Vector<X>
    optimize(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u);

    //! \ingroup LinearProgrammingModule
    //! Solve the linear programming problem \f$\min cx \text{ s.t. } Ax=b;\ l\leq x\leq u\f$.
    //! Returns the optimal vector. Throws an error if the problem is infeasible.
    //! Uses starting variable types \a vt.
    Vector<X>
    optimize(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt);

    //! \ingroup LinearProgrammingModule
    //! Solve the linear programming problem \f$\min cx \text{ s.t. } Ax=b;\ l\leq x\leq u\f$.
    //! Returns the optimal vector. Throws an error if the problem is infeasible.
    //! Uses starting variable types vt, permutation p and inverse basic matrix B.
    Vector<X>
    optimize(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, array<size_t>& p, Matrix<X>& B);

    //! \ingroup LinearProgrammingModule
    //! Test if there exists a point \f$x\f$ with \f$0 \leq x\f$ and \f$Ax=b\f$.
    tribool
    primal_feasible(const Matrix<X>& A, const Vector<X>& b);

    //! \ingroup LinearProgrammingModule
    //! Test if there exists a point \f$y\f$ with \f$yA\leq c\f$.
    tribool
    dual_feasible(const Matrix<X>& A, const Vector<X>& c);

    //! \ingroup LinearProgrammingModule
    //! Test if there exists a point \f$x\f$ with \f$l \leq x \leq u\f$ and \f$Ax=b\f$.
    tribool
    constrained_feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u);

    //! \ingroup LinearProgrammingModule
    //! Test if there exists a point \f$x\f$ with \f$l \leq x \leq u\f$ and \f$Ax=b\f$.
    Vector<X>
    feasible_point(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, array<size_t>& p, Matrix<X>& B);

    //! \ingroup LinearProgrammingModule
    //! Test if there exists a point \f$x\f$ with \f$l \leq x \leq u\f$ and \f$Ax=b\f$.
    //! If the initial array \a vt has size zero, then an initial basis is computed, otherwise the given variable types are assumed.
    //! The values of \a p and \a B corresponding to \a vt may be given.
    //! The values of \a x and \a y are output parameters, giving access to the final primal and dual variables.
    tribool
    constrained_feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x, Vector<X>& y);




    //! \ingroup LinearProgrammingModule
    //! Check whether the assignment of basis, lower and upper variables yields a certificate of feasibility or infeasibility.
    tribool
    verify_primal_feasibility(const Matrix<X>& A, const Vector<X>& b, const array<VariableType>& vt);

    //! \ingroup LinearProgrammingModule
    //! Check whether the assignment of basis, lower and upper variables yields a certificate of feasibility or infeasibility.
    tribool
    verify_dual_feasibility(const Matrix<X>& A, const Vector<X>& c, const array<VariableType>& vt);

    //! \ingroup LinearProgrammingModule
    //! Check whether the assignment of basis, lower and upper variables yields a certificate of feasibility or infeasibility.
    tribool
    verify_constrained_feasibility(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, const array<VariableType>& vt);

    tribool
    constrained_feasible_by_enumeration(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u);



    //! \ingroup LinearProgrammingModule
    //! Compute a permutation \f$p\f$ such that \f$p_{0},\ldots,p_{m-1}\f$ are the basis variables of \f$A\f$, and the inverse matrix \f$A_B^{-1}\f$.
    pair< array<size_t>, Matrix<X> >
    compute_basis(const Matrix<X>& A);

    //! \ingroup LinearProgrammingModule
    //! Perform a single step of the standard linear programming problem, updating the ordered variable array \a p, the inverse basis matrix \a B and the variables \a x.
    bool lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, array<size_t>& p, Matrix<X>& B, Vector<X>& x);

    //! \ingroup LinearProgrammingModule
    //! Perform a single step of the standard linear programming problem, updating the variable type array \a vt, the ordered variable array \a p, the inverse basis matrix \a B and the variables \a x.
    bool lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x);

    //! \ingroup LinearProgrammingModule
    //! Perform a step of the simplex algorithm, choosing basic variable \a s to exit the basis.
    //! The index s is with respect to the list of basic variables, and not the original variable list.
    //! The function to optimize is not needed for this procedure.
    //! The variable which left the basis is the new p[s]
    void lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x, size_t s);

    //! \ingroup LinearProgrammingModule
    //! Compute the inverse of the matrix \a A<sub>B</sub> for basic variables given by the first m items of a p.
    template<class XX> Matrix<XX> compute_B(const Matrix<X>& A, const array<size_t>& p);

    //! Compute the primal variables \a x for the standard linear programming problem.
    template<class XX> Vector<XX> compute_x(const Matrix<X>& A, const Vector<X>& b, const array<size_t>& p, const Matrix<XX>& B);

    //! \ingroup LinearProgrammingModule
    //! Compute primal variables \a x for the constrained linear programming problem.
    template<class XX> Vector<XX> compute_x(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, const array<VariableType>& vt, const array<size_t>& p, const Matrix<XX>& B);


  private:
    tribool _primal_feasible(const Matrix<X>& A, const Vector<X>& b, array<size_t>& p, Matrix<X>& B);
    tribool _dual_feasible(const Matrix<X>& A, const Vector<X>& c, array<size_t>& p, Matrix<X>& B);
    tribool _constrained_feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x);


};


} // namespace Ariadne

#endif
