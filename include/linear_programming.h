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

enum VariableType { BASIS, NONBASIS, LOWER, UPPER, FEASIBLE, INFEASIBLE };
std::ostream& operator<<(std::ostream& os, VariableType t);


//! Test if there exists a point \f$x\f$ with \f$0 \leq x\f$ and \f$Ax=b\f$.
template<class X> tribool 
primal_feasible(const Matrix<X>& A, const Vector<X>& b);

//! Test if there exists a point \f$y\f$ with \f$yA\leq c\f$.
template<class X> tribool 
dual_feasible(const Matrix<X>& A, const Vector<X>& c);

//! Test if there exists a point \f$x\f$ with \f$l \leq x \leq u\f$ and \f$Ax=b\f$.
template<class X> tribool 
constrained_feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u);

//! Tests if \f$yb>\f$ and $\fyA\leq 0\f$.
template<class X> void 
verify_infeasibility(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& y);

template<class X> std::pair< array<size_t>, Matrix<X> > 
compute_basis(const Matrix<X>& A);
   
    
template<class X>
bool lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, array<size_t>& p, Matrix<X>& B, Vector<X>& x);

} // namespace Ariadne

#endif
