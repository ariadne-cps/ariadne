/***************************************************************************
 *            lpfsp.h
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

/*! \file lpfsp.h
 *  \brief Linear programming feasibility solver.
 */

#ifndef ARIADNE_LPSLV_H
#define ARIADNE_LPSLV_H

namespace Ariadne {
  namespace LinearAlgebra {
    
    
    /*! \ingroup LinearProgramming
     *  \brief Solver the primal feasibility problem \f$Ax=b; \ x\geq0\f$ using approximate arithmetic. (Only implement if useful)
     *
     * \return True if a solution is found; false otherwise.
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p.
     * \param x is an output parameter giving a feasible point (if the problem is solvable).
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
     * \param x is an output parameter giving a feasible point (if the problem is solvable).
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
     * \param x is an output parameter giving a dual feasible point (if the problem is solvable).
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
     * \param x is an output parameter giving a feasible point (if the problem is solvable).
     * \param y is an output parameter giving a robust dual feasible point  optimum solution of the dual problem (if the dual problem is unsolvable).
     *
     */
    template<class R, class AP>
    tribool lprfsp(const Matrix<R>& A, const Vector<R>& b,
    Permutation& p,
    Vector<AP>& x, Vector<AP>& y);
    
    
    /*! \ingroup LinearProgramming
     *  \brief Solve the robust constrained feasibility problem \f$Ax=b; \ l< x< u\f$.
     *
     * \return If a solution is found, return true; if a certificate of unsolvability is found, return false; otherwise return indeterminate.
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param l is an \f$n\f$-vector.
     * \param u is an \f$n\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p.
     * \param x is an output parameter giving a feasible point (if the problem is robustly solvable).
     * \param y is an output parameter giving a robust dual feasible point (if the problem is robustly unsolvable).
     *
     */
    template<class R, class AP>
    tribool lprfsc(const Matrix<R>& A, const Vector<R>& b,
    Vector<R>& l, Vector<R>& u,
    Permutation& p,
    Vector<AP>& x, Vector<AP>& y);
    
    
    /*! \ingroup LinearProgramming
     *  \brief Solver the robust dual feasibility problem \f$A^T y<c\f$.
     *
     * \return If a solution is found, return true; if a certificate of unsolvability \f$c^Tx<0;\ Ax=0\f$ is found, return false; otherwise return indeterminate.
     * \param A is an \f$m\times n\f$ matrix.
     * \param c is an \f$n\f$-vector.
     * \param p is an input-output parameter giving a permutation of the rows of \a A. The initial basis elements are in the first \a m elements of \a p.
     * \param x is an output parameter giving a feasible point (if the problem is solvable).
     * \param y is an output parameter giving a robust dual feasible point  optimum solution of the dual problem (if the dual problem is unsolvable).
     *
     */
    template<class R, class AP>
    tribool lprfsd(const Matrix<R>& A, const Vector<R>& c,
    Permutation& p,
    Vector<AP>& x, Vector<AP>& y);
    

  }//namespace LinearAlgebra
}//namespace Ariadne


#endif /* ARIADNE_LPSLV_H */
