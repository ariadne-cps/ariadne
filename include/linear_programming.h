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

enum VariableType { BASIS, LOWER, UPPER };
std::ostream& operator<<(std::ostream& os, VariableType t);

//! \ingroup LinearProgrammingModule 
//! Solve the linear programming problem \f$\min cx \text{ s.t. } Ax=b;\ x\geq0\f$.
//! Returns the optimal vector. Throws an error if the problem is infeasible.
template<class X> Vector<X> 
optimize(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c);

//! \ingroup LinearProgrammingModule 
//! Solve the linear programming problem \f$\min cx \text{ s.t. } Ax=b;\ l\leq x\leq u\f$.
//! Returns the optimal vector. Throws an error if the problem is infeasible.
template<class X> Vector<X> 
optimize(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u);

//! \ingroup LinearProgrammingModule 
//! Test if there exists a point \f$x\f$ with \f$0 \leq x\f$ and \f$Ax=b\f$.
template<class X> tribool 
primal_feasible(const Matrix<X>& A, const Vector<X>& b);

//! \ingroup LinearProgrammingModule 
//! Test if there exists a point \f$y\f$ with \f$yA\leq c\f$.
template<class X> tribool 
dual_feasible(const Matrix<X>& A, const Vector<X>& c);

//! \ingroup LinearProgrammingModule 
//! Test if there exists a point \f$x\f$ with \f$l \leq x \leq u\f$ and \f$Ax=b\f$.
template<class X> tribool 
constrained_feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u);




//! \ingroup LinearProgrammingModule 
//! Check whether the assignment of basis, lower and upper variables yields a certificate of feasibility or infeasibility.
template<class X> tribool 
verify_primal_feasibility(const Matrix<X>& A, const Vector<X>& b, const array<VariableType>& vt);

//! \ingroup LinearProgrammingModule 
//! Check whether the assignment of basis, lower and upper variables yields a certificate of feasibility or infeasibility.
template<class X> tribool 
verify_dual_feasibility(const Matrix<X>& A, const Vector<X>& c, const array<VariableType>& vt);

//! \ingroup LinearProgrammingModule 
//! Check whether the assignment of basis, lower and upper variables yields a certificate of feasibility or infeasibility.
template<class X> tribool 
verify_constrained_feasibility(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, const array<VariableType>& vt);

template<class X> tribool 
constrained_feasible_by_enumeration(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u);



//! \ingroup LinearProgrammingModule 
//! Compute a permutation \f$p\f$ such that \f$p_{0},\ldots,p_{m-1}\f$ are the basis variables of \f$A\f$, and the inverse matrix \f$A_B^{-1}\f$.
template<class X> pair< array<size_t>, Matrix<X> > 
compute_basis(const Matrix<X>& A);

//! \ingroup LinearProgrammingModule 
//! Perform a single step of the standard linear programming problem, updating the ordered variable array \a p, the inverse basis matrix \a B and the variables \a x.
template<class X> bool lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, array<size_t>& p, Matrix<X>& B, Vector<X>& x);

//! \ingroup LinearProgrammingModule 
//! Perform a single step of the standard linear programming problem, updating the variable type array \a vt, the ordered variable array \a p, the inverse basis matrix \a B and the variables \a x.
template<class X> bool lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x);

} // namespace Ariadne

#endif
