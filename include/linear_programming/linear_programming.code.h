/***************************************************************************
 *            linear_programming.code.h
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
 
#include "lp.h"
#include "linear_program.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "base/stlio.h"

namespace Ariadne {
  

template<class R> void instantiate();

template<class AP> tribool lpstp(uint, uint, const AP*, uint, uint, const AP*, uint, const AP*, uint, AP*, uint*, AP*, uint, uint, int);



template<class R>
void
instantiate()
{
  Matrix<R>* M=0;
  Vector<R>* v=0;
  
  solve(*M,*v,*v);
  solve(*M,*v,*v,*v,*v);
  feasible(*M,*v);
  feasible(*M,*v,*v,*v);
  dual_feasible(*M,*v);
}


template<class R>
R
solve(const Matrix<R>& A, const Vector<R>& b, const Vector<R>& c) 
{
  LinearProgram<R> lp(A,b,c);
  return lp.optimal_value();
}


template<class R>
R
solve(const Matrix<R>& A, const Vector<R>& b, const Vector<R>& c,
                         const Vector<R>& l, const Vector<R>& u) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R>
tribool
feasible(const Matrix<R>& A, const Vector<R>& b) 
{
  Vector<R> c(A.number_of_columns(),0);
  LinearProgram<R> lp(A,b,c);
  return lp.is_feasible();
}


template<class R>
tribool
feasible(const Matrix<R>& A, const Vector<R>& b,
                            const Vector<R>& l, const Vector<R>& u) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  /*
  Vector<R> c(A.number_of_columns(),0);
  LinearProgram<R> lp(A,b,c);
  return lp.is_feasible();
  */
}


template<class R>
tribool
dual_feasible(const Matrix<R>& A, const Vector<R>& c) 
{
  // The problem A'y<=c is feasible if and only if the problem min c'x st Ax=0; x>=0 has a negative solution
  Vector<R> b(A.number_of_rows(),0);
  LinearProgram<R> lp(A,b,c);
  R d=lp.optimal_value();
  return (d<0);
}




template<class AP>
tribool 
lpstp(uint m, uint n,
                         const AP* Aptr, uint Arinc, uint Acinc,
                         const AP* bptr, uint binc,
                         const AP* cptr, uint cinc,
                         AP* dptr,
                         uint* pptr,
                         AP* Bptr, uint Brinc, uint Bcinc,
                         int aux // index of additional variable to solve the auxiliary problem
                        )
{
  // A step of the linear programming problem min c'x st Ax=b; x>=0 
  //   A is an mxn-matrix with m<=n
  //   b is an n-vector
  //   c is an m-vector
  //   d is a scalar
  //   pptr is a list of the m pivot variables
  //   B in an mxm-matrix giving the inverse of A restricted to the pivot variables
}


} // namespace Ariadne 
