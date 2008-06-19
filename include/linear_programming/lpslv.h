/***************************************************************************
 *            lpslv.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
 ****************************************************************************/
/*
 * Based on the linear programming algorithms in PPL-0.8
 *   Copyright (C) 2001-2006 Roberto Bagnara <bagnara@cs.unipr.it>
 */

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
 
/*! \file lpslv.h
 *  \brief Linear programming solver.
 */

#ifndef ARIADNE_LPSLV_H
#define ARIADNE_LPSLV_H

#include "linear_algebra/declarations.h"

namespace Ariadne { 

  

    /*! \ingroup LinearProgramming
     *  \brief Solver the linear programming problem \f$\max c^Tx\text{ s.t. } Ax=b; \ x\geq0\f$.
     *.
     * This problem is equivalent to the dual problem \f$\min b^Ty\text{ s.t. } A^Ty\leq c\f$ with slack variables \f$z=c-A^Ty\f$.
     *
     * \return d is the optimum value (if the problem has a finite feasible solution).
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param c is an \f$n\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p. 
     * \param x is an output parameter giving the optimum point (if the primal problem is feasible and finite).
     * \param y is an output parameter giving the optimum solution of the dual problem (if the dual problem is feasible and finite).
     *
     * For a given basis \f$B\f$ with nonbasic variables \f$N\f$, the matrix \f$A_B\f$ is the matrix formed by 
     * the columns of \f$A\f$ with indices in \f$B\f$. 
     * The basic solution is given by:
     * \f[ d=c^TA_B^{-1}b;\ x_B=A_N^{-1}b; \ y=A_B^{-1}c_B; \ z_N = c_N - A_NA_{B}^{-1}c_B  \f]
     */
    template<class R, class AP>
    AP 
    lpslv(const Matrix<R>& A, 
          const Vector<R>& b, 
          const Vector<R>& c, 
          Permutation& p,
          Vector<AP>& x,
          Vector<AP>& y);
  

    /*! \ingroup LinearProgramming
     *  \brief Solver the constrained linear programming problem \f$\max c^Tx\text{ s.t. } Ax=b; \ l\leq x\leq0\f$.
     *
     * \return d is the optimum value (if the problem has a finite feasible solution).
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param c is an \f$n\f$-vector.
     * \param l is an \f$n\f$-vector.
     * \param u is an \f$n\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p. 
     * \param x is an output parameter giving the optimum point (if the primal problem is feasible and finite).
     * \param y is an output parameter giving the optimum solution of the dual problem (if the dual problem is feasible and finite).
     *
     * For a given basis \f$B\f$, the non-basic variables \f$N\f$ are partitioned into a set \f$L\f$ taking value \f$l\f$ and a set \f$U\f$ taking value \f$u\f$.
     * The basic solution is given by:
     * \f[ d=c^Tx; \quad x_L=l_L; \ x_U=u_U;\ x_B=A_N^{-1}(b-A_Nx_N) . \f]
     */
    template<class R, class AP>
    AP 
    lpslvc(const Matrix<R>& A, 
           const Vector<R>& b, 
           const Vector<R>& c, 
           const Vector<R>& l, 
           const Vector<R>& u,
           Permutation& p,
           Vector<AP>& x, 
           Vector<AP>& y);
 
    
     /*! \ingroup LinearProgramming
     *  \brief Solver for linear programming problems. (Deprecated) \deprecated
     *
     *  \param m Number of free constraints
     *  \param n Number of free variables
     *  \param (A,rincA,cincA) An m-by-n matrix
     *  \param (B,incB) An m-element vector
     *  \param (C,incC) An n-element vector
     *  \param d A scalar
     *  \param piv A consecutive (m+n)-element vector of integers.
     *
     * Solve the linear programming problem 
     * \f$ text{minimize}\ c^Tx \text{ subject to } Ax+y=b,\ x,y\geq0\f$ where the
     * current point is given by \f$x=0,\ y=b\f$ and the current value is \f$-d\f$. 
     */
    template<class R>
    void 
    lpslv(int m, int n, 
          R* A, int rincA, int cincA, 
          R* B, int incB, 
          R* C, int incC, 
          R& d, 
          int* piv);


} // namespace Ariadne

#include "lpslv.template.h"

#endif /* ARIADNE_LPSLV_H */
