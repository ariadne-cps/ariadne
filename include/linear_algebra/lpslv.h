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

#include <iosfwd>
#include <cassert>
#include <map>


#include "../linear_algebra/matrix.h"
#include "../output/logging.h"

namespace Ariadne { 
  namespace LinearAlgebra {
  
    // Forward declarations
    class Permutation;
    template<class X> class Vector;
    template<class X> class Matrix;


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
    template<class R, class A>
    A lpslv(const Matrix<R>& A, const Vector<R>& b, const Vector<R>& c, 
            Permutation& p,
            Vector<A>& x, Vector<A>& y);



    /*! \ingroup LinearProgramming
     *  \brief Solver the constrained linear programming problem \f$\max c^Tx\text{ s.t. } Ax=b; \ l\leq x\leq0\f$.
     *
     * For a given basis \f$B\f$, the non-basic variables \f$N\f$ are partitioned into a set \f$L\f$ taking value \f$l\f$ and a set \f$U\f$ taking value \f$u\f$.
     * The basic solution is given by:
     * \f[ d=c^Tx; \quad x_L=l_L; \ x_U=u_U;\ x_B=A_N^{-1}(b-A_Nx_N) . \f]
     */
    template<class R, class A>
    A lpslvc(const Matrix<R>& A, const Vector<R>& b, const Vector<R>& c, 
             const Vector<R>& l, const Vector<R>& u,
             Permutation& p,
             Vector<A>& x, Vector<A>& y);




    /*! \ingroup LinearProgramming
     *  \brief Solver the primal feasibility problem \f$Ax=b; \ x\geq0\f$ using approximate arithmetic. (Only implement if useful)
     *
     * \return True if a solution is found; false otherwise.
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p. 
     * \param x is an output parameter giving a fesible point (if the problem is solvable).
     * \param y is an output parameter giving a dual feasible point (if the problem is unsolvable).
     *
     */
    template<class AP>
    bool lpfsp(const Matrix<AP>& A, const Vector<AP>& b, 
               Permutation& p,
               Vector<AP>& x, Vector<AP>& y);



    /*! \ingroup LinearProgramming
     *  \brief Solver the constrained primal feasibility problem \f$Ax=b; \ l\leq x\leq u\f$ using approximate arithmetic. (Only implement if useful)
     *
     * \return True if a solution is found; false otherwise.
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param l is an \f$n\f$-vector.
     * \param u is an \f$n\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p. 
     * \param x is an output parameter giving a fesible point (if the problem is solvable).
     * \param y is an output parameter giving a dual feasible point (if the problem is unsolvable).
     *
     * \internal Only implement this if necessary!
     */
    template<class AP>
    bool lpfsc(const Matrix<AP>& A, const Vector<AP>& b, 
               Vector<AP>& l, Vector<AP>& u,
               Permutation& p,
               Vector<AP>& x, Vector<AP>& y);


    /*! \ingroup LinearProgramming
     *  \brief Solver the dual feasibility problem \f$A^Ty\leq c\f$ using approximate arithmetic. (Only implement if useful)
     *
     * \return True if a solution is found; false otherwise.
     * \param A is an \f$m\times n\f$ matrix.
     * \param c is an \f$nm\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p. 
     * \param x is an output parameter giving a dual fesible point (if the problem is solvable).
     * \param y is an output parameter giving a primal feasible point (if the problem is unsolvable).
     *
     * \internal Only implement this if necessary!
     */
    template<class AP>
    bool lpfsd(const Matrix<AP>& A, const Vector<AP>& c, 
               Permutation& p,
               Vector<AP>& x, Vector<AP>& y);




    /*! \ingroup LinearProgramming
     *  \brief Solver the robust primal feasibility problem \f$Ax=b; \ x>0\f$.
     *
     * \return If a solution is found, return true; if a certificate of unsolvability \f$b^Ty>0;\ A^Ty<0\f$ is found, return false; otherwise return indeterminate.
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p. 
     * \param x is an output parameter giving a fesible point (if the problem is solvable).
     * \param y is an output parameter giving a robust dual feasible point  optimum solution of the dual problem (if the dual problem is unsolvable).
     *
     */
    template<class R, class A>
    tribool lprfsp(const Matrix<R>& A, const Vector<R>& b, 
                   Permutation& p,
                   Vector<A>& x, Vector<A>& y);



    /*! \ingroup LinearProgramming
     *  \brief Solve the robust constrained feasibility problem \f$Ax=b; \ l< x< u\f$.
     *
     * \return If a solution is found, return true; if a certificate of unsolvability is found, return false; otherwise return indeterminate.
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param l is an \f$n\f$-vector.
     * \param u is an \f$n\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p. 
     * \param x is an output parameter giving a fesible point (if the problem is robustly solvable).
     * \param y is an output parameter giving a robust dual feasible point (if the problem is robustly unsolvable).
     *
     */
    template<class R, class A>
    tribool lprfsc(const Matrix<R>& A, const Vector<R>& b, 
                   Vector<R>& l, Vector<R>& u,
                   Permutation& p,
                   Vector<A>& x, Vector<A>& y);



  /*! \ingroup LinearProgramming
     *  \brief Solver the robust dual feasibility problem \f$A^T y<c\f$.
     *
     * \return If a solution is found, return true; if a certificate of unsolvability \f$c^Tx<0;\ Ax=0\f$ is found, return false; otherwise return indeterminate.
     * \param A is an \f$m\times n\f$ matrix.
     * \param c is an \f$n\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p. 
     * \param x is an output parameter giving a fesible point (if the problem is solvable).
     * \param y is an output parameter giving a robust dual feasible point  optimum solution of the dual problem (if the dual problem is unsolvable).
     *
     */
    template<class R, class A>
    tribool lprfsd(const Matrix<R>& A, const Vector<R>& c, 
                   Permutation& p,
                   Vector<A>& x, Vector<A>& y);
   

    /*!\ingroup LinearProgramming
     * \brief A step of the simplex algorithm for the linear programming problem \f$\min c^Tx \text{ s.t. } Ax=b; \ x\geq0\f$, or
     * equivalently, of the dual problem \f$\max b^Ty \text{ s.t. } A^Tx\leq c\f$. (Only implement if useful)
     * 
     * \return True if the current basic solution is optimal.
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param c is an \f$n\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p. 
     * \param B is an input-output parameter giving the inverse of the basis matrix; \f$B=A_B^{-1}\f$.
     * \param x is an \f$n\f$-vector giving the current variables.(Optional)
     * \param y is an \f$m\f$-vector giving the current dual variables.(Optional)
     * \param z is an \f$n\f$-vector giving the current slack variables. (Optional)
     *
     * The basic primal, dual and slack variables are
     * \f[ x_B=A_B^{-1}b; \ x_N=0; y=A_B^{-1}c_B;\ z_B=0; z_N=c_N-A_NA_B^{-1}c_B . \f]
     * The algorithm proceeds as follows.
     *  -# Find an entering variable \f$x_j\f$ such that satisfying \f$z_j<0\f$.
     *  -# Find the leaving variable \f$x_k\f$ such that \f$i=\arg\min\{ t_k = x_k/d_k \mid d_k>0  \} \f$where \f$d = A_B^{-1} a_j\f$ and \f$a_j\f$ is the \f$j^\mathrm{th}\f$ column of \f$A\f$.
     *  -# Update the inverse basis matrix
     *     \f$ A_{B'}^{-1} = A_B^{-1} - (d-e_i) r^T / d_j . \f$
     *    Here, \f$e_n\f$ is the n-th unit vector, \f$d:=A_B^{-1}a_j\f$ is the direction of change of the state, and \f$r^T:= e_i^TA_B^{-1}\f$ is the i-th row of \f$A_B^{-1}\f$.
     *
     * \internal Use x and y (and maybe z) only if it is more efficient. This routine might not be so useful; maybe use matrix pointers instead.
     */
    template<class AP>
    bool lpstp(const Matrix<AP>& A, const Vector<AP>& b, const Vector<AP>& c, 
               Permutation& p, Matrix<AP>& B
               // Use these variables too if it's more efficient
               , Vector<AP>& x, Vector<AP>& y, Vector<AP>& z
               );





    /*!\ingroup LinearProgramming
     * \brief A step of the simplex algorithm for the constrained linear programming problem \f$\min c^Tx \text{ s.t. } Ax=b; \ l\leq x\leq u\f$, or the equivalent dual problem, using approximate (i.e. double or mpfr_class) arithmetic. 
     * 
     * \return True if the current basic solution is optimal.
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param c is an \f$n\f$-vector.
     * \param l is an \f$n\f$-vector.
     * \param u is an \f$n\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p. 
     * \param B is an input-output parameter giving the inverse of the basis matrix; \f$B=A_B^{-1}\f$.
     * \param x is the current point.
     *
     * \internal Use appropriate dual and slack variables if this makes other routines more efficient. This routine might not be so useful; maybe use matrix pointers instead.
     */
    template<class AP>
    bool lpstpc(const Matrix<AP>& A, const Vector<AP>& b, const Vector<AP>& c, 
                const Vector<AP>& l, const Vector<AP>& u,
                Permutation& p, Matrix<AP>& B, Vector<AP>& x);

    /*!\ingroup LinearProgramming
     * \brief A step of the simplex algorithm for the linear programming problem \f$\min c^Tx \text{ s.t. } Ax=b; \ x\geq0\f$, or
     * equivalently, of the dual problem \f$\max b^Ty \text{ s.t. } A^Tx\leq c\f$, using pointers to matrix data and approximate (i.e. double or mpfr_class) arithmetic.
     * 
     * \return True if the current basic solution is optimal.
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param c is an \f$n\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p. 
     * \param B is an input-output parameter giving the inverse of the basis matrix; \f$B=A_B^{-1}\f$.
     * \param x is an \f$n\f$-vector giving the current variables.(Optional)
     * \param y is an \f$m\f$-vector giving the current dual variables.(Optional)
     * \param z is an \f$n\f$-vector giving the current slack variables. (Optional)
     *
     * The basic primal, dual and slack variables are
     * \f[ x_B=A_B^{-1}b; \ x_N=0; y=A_B^{-1}c_B;\ z_B=0; z_N=c_N-A_NA_B^{-1}c_B . \f]
     * The algorithm proceeds as follows.
     *  -# Find an entering variable \f$x_j\f$ such that satisfying \f$z_j<0\f$.
     *  -# Find the leaving variable \f$x_k\f$ such that \f$i=\arg\min\{ t_k = x_k/d_k \mid d_k>0  \} \f$where \f$d = A_B^{-1} a_j\f$ and \f$a_j\f$ is the \f$j^\mathrm{th}\f$ column of \f$A\f$.
     *  -# Update the inverse basis matrix
     *     \f$ A_{B'}^{-1} = A_B^{-1} - (d-e_i) r^T / d_j . \f$
     *    Here, \f$e_n\f$ is the n-th unit vector, \f$d:=A_B^{-1}a_j\f$ is the direction of change of the state, and \f$r^T:= e_i^TA_B^{-1}\f$ is the i-th row of \f$A_B^{-1}\f$.
     *
     * \internal Use x, y and z only if it is more efficient for other routines.
     */
    template<class AP>
    bool lpstp(uint m, uint n, 
               const AP* Aptr, uint Arinc, uint Acinc,
               const AP* bptr, uint binc,
               const AP* cptr, uint cinc,
               uint* pptr,
               AP* Bptr, uint Brinc, uint Bcinc
               // use only if it's more efficient
               , AP* xptr, uint xinc
               , AP* yptr, uint yinc
               , AP* zptr, uint zinc
              );



    /*!\ingroup LinearProgramming
     * \brief A step of the simplex algorithm for the constrained linear programming problem \f$\min c^Tx \text{ s.t. } Ax=b; \ l\leq x\leq u\f$, or the equivalent dual problem, using pointers to matrix data and approximate (i.e. double or mpfr_class) arithmetic.
     * 
     * \return True if the current basic solution is optimal.
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param c is an \f$n\f$-vector.
     * \param l is an \f$n\f$-vector.
     * \param u is an \f$n\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p. 
     * \param B is an input-output parameter giving the inverse of the basis matrix; \f$B=A_B^{-1}\f$.
     * \param x is the current point.
     */
    template<class AP>
    bool lpstpc(uint m, uint n, 
                const AP* Aptr, uint Arinc, uint Acinc,
                const AP* bptr, uint binc,
                const AP* cptr, uint cinc,
                const AP* lptr, uint linc,
                const AP* uptr, uint uinc,
                uint* pptr,
                AP* Bptr, uint Brinc, uint Bcinc
                , AP* xptr, uint xinc
                // input/output dual variables if it's more efficient
                , AP* yptr, uint yinc
                // input/output slack variables if it's more efficient
                , AP* zptr, uint zinc
               );

     /*! \ingroup LinearAlgebra
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

}}

#include "lpslv.template.h"

#endif /* ARIADNE_LPSLV_H */
