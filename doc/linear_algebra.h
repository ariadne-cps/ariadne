/***************************************************************************
 *            documentation/linear_algebra.h
 *
 *  Copyright  2004-7  Pieter Collins
 *  Pieter.Collins@cwi.nl
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


/*! 

\file documentation/linear_programming.h
\brief Documentation on linear algebra



\page linear_algebra Linear Algebra

\section matrix_inverse Inverting Matrices

\subsection schulz_matrix_inverse Schulz's method 
Schulz's method is an iterative method for computing matrix inverses:
\f[ B_{k+1} = B_k \sum_{i=0}^{n} (I-AB_k)^i .  \f]
The interval version is given by:
\f[ [B_{k+1}] = m[B_k] \,  \Bigl( \sum_{i=0}^{n-1} \bigl(I-A\,m[B_k]\bigr)^i \Bigr) \, + \, [B_k]\,\bigl(1-A\,m[B_k]\bigr)^n .  \f]

\subsection gauss_seidell_matrix_inverse Gauss-Seidel Iteration

Gauss-Seidel iteration is a method for solving linear equations \f$Ax=b\f$, which can also be used to compute inverses using \f$AB=I\f$. The scheme is
\f[   x^{(k+1)} = \left( {D + L} \right)^{ - 1} \left( {b-U x^{(k)} } \right) \f]
where \f$A=L+D+U\f$; the matrices \f$D\f$, \f$L\f$ and \f$U\f$ represent the diagonal, negative strictly lower triangular, and negative strictly upper triangular parts of \f$A\f$. When implementing Gauss Seidel, an explicit entry-by-entry approach is commonly used:
\f[  x^{(k+1)}_i = \frac{1}{a_{ii}} \left(b_i - \sum_{j<i}a_{ij}x^{(k+1)}_j-\sum_{j>i}a_{ij}x^{(k)}_j\right),\, i=1,2,\ldots,n. \f] 

*/
