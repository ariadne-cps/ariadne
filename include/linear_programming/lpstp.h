/***************************************************************************
 *            lpstp.h
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

/*! \file lpstp.h
 *  \brief Linear programming solver step.
 */

#ifndef ARIADNE_LPSLV_H
#define ARIADNE_LPSLV_H

namespace Ariadne {
  namespace LinearAlgebra {
    
    
    /*!\ingroup LinearProgramming
     * \brief A step of the simplex algorithm for the linear programming problem \f$\min c^Tx \text{ s.t. } Ax=b; \ x\geq0\f$, or
     * equivalently, of the dual problem \f$\max b^Ty \text{ s.t. } A^Tx\leq c\f$. (Only implement if useful)
     *
     * \return True if the current basic solution is optimal; Indeterminate if problem is unbounded.
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
    tribool lpstp(const Matrix<AP>& A, const Vector<AP>& b, const Vector<AP>& c,
    Permutation& p, Matrix<AP>& B
    // Use these variables too if it's more efficient
    , Vector<AP>& x, Vector<AP>& y, Vector<AP>& z
    );
    
    
    
    /*!\ingroup LinearProgramming
     * \brief A step of the simplex algorithm for the constrained linear programming problem \f$\min c^Tx \text{ s.t. } Ax=b; \ l\leq x\leq u\f$, or the equivalent dual problem, using approximate (i.e. double or mpfr_class) arithmetic.
     *
     * \return True if the current basic solution is optimal; Indeterminate if problem is unbounded.
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
    tribool lpstpc(const Matrix<AP>& A, const Vector<AP>& b, const Vector<AP>& c,
    const Vector<AP>& l, const Vector<AP>& u,
    Permutation& p, Matrix<AP>& B, Vector<AP>& x);
    
    
    
    /*!\ingroup LinearProgramming
     * \brief A step of the simplex algorithm for the linear programming problem \f$\min c^Tx \text{ s.t. } Ax=b; \ x\geq0\f$, or
     * equivalently, of the dual problem \f$\max b^Ty \text{ s.t. } A^Tx\leq c\f$, using pointers to matrix data and approximate (i.e. double or mpfr_class) arithmetic.
     *
     * \return True if the current basic solution is optimal; Indeterminate if problem is unbounded.
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
    tribool lpstp(uint m, uint n,
    const AP* Aptr, uint Arinc, uint Acinc,
    const AP* bptr, uint binc,
    const AP* cptr, uint cinc,
    const AP* dptr,
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
     * \return True if the current basic solution is optimal; Indeterminate if problem is unbounded.
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param c is an \f$n\f$-vector.
     * \param d is an \f$n\f$-vector.
     * \param l is an \f$n\f$-vector.
     * \param u is an \f$n\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p.
     * \param B is an input-output parameter giving the inverse of the basis matrix; \f$B=A_B^{-1}\f$.
     * \param x is the current point.
     */
    template<class AP>
    tribool lpstpc(uint m, uint n,
    const AP* Aptr, uint Arinc, uint Acinc,
    const AP* bptr, uint binc,
    const AP* cptr, uint cinc,
    const AP* lptr, uint linc,
    const AP* uptr, uint uinc,
    const AP* dptr,
    uint* pptr,
    AP* Bptr, uint Brinc, uint Bcinc
    , AP* xptr, uint xinc
    // input/output dual variables if it's more efficient
    , AP* yptr, uint yinc
    // input/output slack variables if it's more efficient
    , AP* zptr, uint zinc
    );
    
    
  } // namespace LinearAlgebra
} // namespace Ariadne


#endif /* ARIADNE_LPSLV_H */
