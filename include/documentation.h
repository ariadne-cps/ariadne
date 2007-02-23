/***************************************************************************
 *            documentation.h
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

#ifndef _ARIADNE_DOCUMENTATION_H
#define _ARIADNE_DOCUMENTATION_H

/*! \file documentation.h
 * \brief Miscellaneous documentation pages
 */

/*! 

\mainpage

\section Introduction

%Ariadne is a C++ package for set-based analysis of dynamical and control systems, including reachability analysis and verification.

\section Download

The homepage of %Ariadne is <a href="http://fsv.dimi.uniud.it/ariadne/">http://fsv.dimi.uniud.it/ariadne/</a>.

You can check out the latest version on the Subversion repository by typing:

  \code svn checkout https://fsv.dimi.uniud.it/svn/ariadne ariadne \endcode

To make the code documentation, change to the ariadne/trunk/ directory and type:

  \code doxygen \endcode

\section Requirements

To compile %Ariadne, you will need a C++ compiler (we recommend g++ version 4.0.2 or higher, which can be downloaded from <a href="http://gcc.gnu.org/">http://gcc.gnu.org</a>).

You will also need the following libraries:
  - The GNU Multiple-Precision Library (version 4.1.2 or higher)  <a href="http://www.swox.com/gmp/">http://www.swox.com/gmp/</a>.
  - MPFR Library (version 2.2.1 or higher) <a href="http://www.mpfr.org/">http://www.mpfr.org/</a>.
  - The Boost C++ Libraries (version 1.33.1) <a href="http://www.boost.org/">http://www.boost.org/</a>.
  - The Parma Polyhedra Library (version 0.9 or higher) <a href="http://www.cs.unipr.it/ppl/">http://www.cs.unipr.it/ppl/</a>.
  - TBLAS (version 0.4.1 or higher) <a href="http://homepages.cwi.nl/~collins/software/">http://homepages.cwi.nl/~collins/software/</a>.
  - TLAPACK (version 0.4.0 or higher) <a href="http://homepages.cwi.nl/~collins/software/">http://homepages.cwi.nl/~collins/software/</a>.

For the Python interface, you will also need:
  - Python (version 2.4) <a href="http://www.python.org/">http://www.python.org/</a>.

To make the source code documentation, you will need:
  - Doxygen (version 1.4.6 or higher is recommended) 
       <a href="http://www.stack.nl/~dimitri/doxygen/">http://www.stack.nl/~dimitri/doxygen/</a>

\section Installation

From main source directory, type
\code 
  ./configure
  make
  make install
\endcode
The default installation directory for the library is $/usr/local/lib/, and for the Python interface is /usr/local/python/.
These defaults can be changed by using the --prefix flag to ./configure, as in
\code 
  ./configure --prefix=$HOME
\endcode

If installing from the Subversion repository, change to the ariadne/trunk/ directory and type 
\code 
  ./bootstrap
\endcode
Then follow the directions above.




 
\page real Real Number Types

\section Introduction

Ariadne currently supports three different real number types, \a \ref Float64, \a MPFloat and \a Rational.
The \a Float64 type is a \em finite-precision type, the \a MPFloat type is a \em multiple-precision type
and \a Rational is an \em exact number type. Future realeases will also support an \em arbitrary-precision
real number type \a Real.

\section floatingpoint Floating-Point Numbers

Real numbers are traditionally described by floating-point types \a float and \a double.
The set of elements which can be represented by these typesis <em>finite</em>.
Arithmetic operations can only be performed approximately to maximal precision determined by the data type.
However, these types have the advantage of requiring a known amount of memory, which means they can be statically
allocated, and having hardware-supported arithmetical approximations. This makes them especially suitable in
situations where execution speed and memory usage are more important than a knowledge of the computational accuracy.

\section finiteprecision Finite-Precision Real Numbers

For a package such as %Ariadne, in which we are concerned with keeping track of numerical errors,
we cannot use these types directly. Instead, we provide four different arithmetic operations for
the same type. We can specify that the result can be rounded up or down, or that the result
should be rounded to the nearest representable value. However, the default option is to
<em>return an interval</em> containing the exact value of the operation. To summarize, we have the following 
operations.
\code
// Exact arithmetic where possible
FPReal neg_exact(FPReal);
FPReal neg(FPReal);
FPReal operator-(FPReal);

// Approximate and rounded arithmetic
FPReal add_approx(FPReal);
FPReal add_down(FPReal);
FPReal add_up(FPReal);

// Interval arithmetic
Interval<FPReal> operator+(Interval<FPReal>,Interval<FPReal>);

// Automatic conversion to interval arithmetic
Interval<FPReal> operator+(FPReal,FPReal);
\endcode
For types based on floating point numbers, the arithmetical operators should follow IEEE standards.

\code
// Approximate and rounded algebraic transcendental functions
FPReal exp_approx(FPReal);
FPReal exp_down(FPReal);
FPReal exp_up(FPReal);

// Interval algebraic and transcendental functions
Interval<FPReal> exp(Interval<FPReal>);
\endcode
The precision of the algebraic and transcendental functions is not guarenteed,
except that the functions must give values in the correct range.

The upper and lower rounded versions satisfy the mathematical postcondition 
\f$ \underline{f}(x) \leq f(x) \leq \overline{f}(x) \f$
and the implementation postcondition 
\code f_down(x) <= f_approx(x) <= f_up(x) \endcode

\note
No other conditions are required of the implementation.
Hence a valid (but useless) implementation of \f$ \sin(x) \f$ is
\code
FPReal sin_approx(FPReal x) { return 0.0; }
FPReal sin_down(FPReal x) { return -1.0; }
FPReal sin_up(FPReal x) { return 1.0; }
\endcode


\section multipleprecision Multiple-precision types.

If we have a problem for which a fixed-precision type is not sufficient to obtain an accurate answer,
we can switch to a \em multiple-precision type. The semantics of arithmetic on multiple-precision types 
is mostly the same as that of a fixed-precision type; arithmetic is approximate, and the default is to
return an interval. However, a multiple-precision type has a 
\code MPReal::set_precision(unsigned int) \endcode
method, which sets the precision to which the type can store its result.
Using a higher precision yields a more accurate answer. Further, the result of an
arithmetic operation is guarenteed to converge as the precision is increased.

For efficiency, all elements of an array of a multiple precision type 

\section arbitraryprecision Arbitrary-precision types.

An arbitrary-precision type stores a number in a form so that it can be
recovered to any desired precision. This is typically acheived by expressing
the number as a formula in terms of other arbitrary-precision or exact number
types. However, the high computational overhead of such numbers makes them
impractical for describing higher-order types, such as matrices or sets.

Arbitrary-precision numbers may be used to store constants occurring in a 
system definition.

\code
// Exact arithmetical operators.
APReal operator-(APReal);
APReal operator+(APReal,APReal);

// Exact algebraic and transcendental functions
APReal sqrt(APReal);
APReal exp(APReal);
APReal sin(APReal);
\endcode
\section exactarithmetic Exact arithmetic types.

Elements of a countable set of numbers, such as integer, dyadic, rational or algebraic numbers,
may be stored using a finite amount of data, though the size of the data depends on the 
element used. Hence in an array of an exact arithmetic type, each element needs to be
dynamically allocated, which is expensive in terms of spacial overhead. 
Since these sets are typically closed under arithmetical operations, 
arithmetic can be performed exactly for these types. Hence these types are
appropriate where time and space overhead are not at a premium.

Since these types require arbitrarily
large amounts of memory which is typically dynamically allocated, and hardware support for arithmetic does not
exist, they are typically less efficient than the fixed-precision and multiple-precision types.

In %Ariadne, we support arithmetic on the Rational number type, but no algebraic or transcendental functions.
This type is primarily useful for testing.

\code
// Exact arithmetical operators.
ExReal operator-(ExReal);
ExReal operator+(ExReal,ExReal);

// No algebraic or transcendental functions
\endcode

Unlike a finite-precision type, arbitrary-precision types are intended for use when precise error specifications are required.
Hence, whenever a function returns an approximation, the suffix _approx is \em always added to the function name.
This ensures that the user is always aware of the use of approximations.

The Dyadic type does not in general support exact division.
The exception is that division by a power of 2 yields an exact answer.
Generic code intended for use with all exact arithmetic types should \em not use division, except by a power of two.



\section interval Interval Arithmetic

Interval arithmetic is commonly used in finite-precision computations to obtain guarenteed bounds for the result.
However, it can also be used in arbtitrary-precision computations to compute the image of basic sets.
The semantics is as follows:

\subsection finite_precision_interval Finite-precision interval arithmetic

Finite precision functions are of the form
\code
Interval<FPReal> f(Interval<FPReal>);
\endcode
and satisfy the mathematical postcondition \f$ \forall x\in I,\ f(x)\in\f$\c f(I),
and the implementation postcondition <tt>I.contains(x)</tt> implies <tt>f(I).contains(f_approx(x))</tt>.

\subsection multiple_precision_interval Multiple-precision interval arithmetic

Multiple-precision interval functions follow the same conditions as fixed-precision
interval functions, together with the convergence criterion
that the length of \c f(I) approaches 0 as the length of \c I approaches 0.




\page pivot Pivoting and Permutations

Re-ordering of matrix rows/columns and vector elements can be described by <em>permutations</em> and <em>pivoting</em>.
Essentially, a permutation gives the re-ordering directly, and is most useful if the re-ordered elements are to be used in-place or copied, whereas pivoting describes the re-ordering as a product of transpositions, and is most useful if the elements are to be re-ordered in place.

We write a permuation
\f[ \pi=\left(\begin{array}{cccccc}0&1&\cdots&n\!-\!1\\\pi_0&\pi_1&\cdots&\pi_{n-1}\end{array}\right); \quad \pi_i\neq\pi_j \text{ if } i\neq j \f]
and a pivot
\f[ p = (0\;p_0)\ (1\;p_1)\ \cdots\ (n\!-\!1,p_{n-1}); \quad p_i\geq i . \f]
where the transpositions are applied from left to right.

For example, the permutation
\f[ \pi=\left(\begin{array}{ccccc}0&1&2&3&4\\1&3&4&2&0\end{array}\right)  \f]
can be writing in pivot form
\f[ p=(0\;1)\ (1\;3)\ (2\;4)\ (3\;4)\ (4\;4) . \f]
The inverses are
\f[ \pi^{-1}=\left(\begin{array}{ccccc}0&1&2&3&4\\4&0&3&1&2\end{array}\right); \qquad 
     p^{-1}=(4\;4)\ (3\;4)\ (2\;4)\ (1\;3)\ (0\;1) \f]

Applying the permuation to the vector \f$v=(v_0,v_1,v_2,v_3,v_4)^T\f$ yields the vector \f$\pi v=(v_1,v_3,v_4,v_2,v_0)^T\f$.
Applying the pivots to the vector yield successively \f$(v_1,v_0,v_2,v_3,v_4)^T\f$, \f$(v_1,v_3,v_2,v_0,v_4)^T\f$, \f$(v_1,v_3,v_4,v_0,v_2)^T\f$ and \f$Pv=(v_1,v_3,v_4,v_2,v_0)^T\f$.

The pivot form is used by linear algebra packages such as LAPACK.
For the triangular factorization of a matrix \f$A\f$ with partial pivoting yields matrices \f$P,L,U\f$ such that \f$PA=LU\f$, where \f$P\f$ is a pivot matrix.
The matrix
\f[ A=\left(\begin{matrix}0&0&0&0&1\\1&x&x&x&x\\0&0&0&1&x\\0&1&x&x&x\\0&0&1&x&x\end{matrix}\right)\f]
yields a pivot matrix \f$P\f$ described by the pivot \f$p=(0\;1)\ (1\;3)\ (2\;4)\ (3\;4)\ (4\;4)\f$. In other words, during the factorisation, we first interchange row 0 and row 1, then row 1 (which corresponds to the old row 0) with row 3, and so on.

The permuation form is mostly used by linear programming packages to select basis variables.
<br>
Consider the \f$m\times n\f$ matrix \f$\mathrm{A}\f$ with \f$n>m\f$, and suppose \f$\mathrm{A}\f$ has full row rank.
Let \f$B=\{ j_0,j_1,\ldots,j_{m-1}\}\f$ be a set of \f$m\f$ values from \f$\{0,\ldots,n\!-\!1\}\f$, and define the matrix \f$\mathrm{B}=\mathrm{A}_B=\bigl( a_{i_0}\ a_{i_1}\ \cdots\ a_{i_{n-1}}\bigr)\f$.
Here, a (partial) permutation is used to select the pivot elements.
<br>
In the reduced simplex algorithm, we store an intermediate form of the problem as \f$x_B+\mathrm{A}x_N=b\f$
We now store both \f$B=(j_0,j_1,\ldots,j_{m-1})\f$ giving the columns of the basic variables, and \f$N=(j_{m},j_{m+1},\ldots,j_{n-1})\f$ giving the columns of the non-basic variables.






\page linearprogramming Linear Programming

\section standardprimaldual Standard primal and dual problems

The standard linear programming problem is
\f[ \text{(P)} \qquad \min c^Tx \text{ s.t. } Ax=b;\ x\geq0.  \f]
Without loss of generality we can take \f$b\geq0\f$.
We let \f$x^*\f$ be an optimal point.

The dual to the standard problem is
\f[ \text{(D)} \qquad \max b^T y \text{ s.t. } A^Ty\leq c \f]
or alternatively
\f[ \text{(D)} \qquad \max b^T y \text{ s.t. } A^Ty+z=c;\ z\geq 0. \f]
The variables \f$y\f$ are called <em>dual variables</em> and the variables \f$z\f$ are <em>slack variables</em>.

<b>Theorem</b> If \f$x\f$ is feasible for (P) and \f$y\f$ is feasible for (D), then \f[c^Tx \geq b^Ty.\f]
Further, if \f$x^*\f$ is optimal for (P) and \f$(y^*,z^*)\f$ is optimal for (D), then \f[c^Tx^* = b^Ty^* \qquad \text{and} \qquad x^*\!\cdot\! z^*=0 .\f]

The second condition is called <em>complementary slackness</em>.


\section optimalbasis Optimal basic solutions 

A basic solution is given by a set of column indices \f$B\f$ such that the square matrix \f$A_B\f$ formed by the \f$B\f$ columns of \f$A\f$ is nonsingular.
Then
\f[x_B=A_B^{-1}b, \ x_N=0; \quad y=(A_B^T)^{-1}c_B; \quad z_B=0, \ z_N=c_N^T-c_B^TA_B^{-1}A_N; \qquad c^Tx = c_B^TA_B^{-1}b = b^T(A_B^T)^{-1}c_B = b^Ty. \f]

Suppose \f$x^*\f$ is optimal. Then \f$x_B\geq0\f$ and \f[c^Tx = c_B^T x_B + c_N^T x_N = c_B^TA_B^{-1}(b-A_Nx_N)+c_N^Tx_N = c_B^TA_B^{-1}b + (c_N^T-c_B^TA_B^{-1}A_N)x_N = c^Tx^* + (c_N^T-c_B^TA_B^{-1}A_N)x_N,\f] and hence \f$z_N=c_N^T-c_B^TA_B^{-1}A_N\geq0\f$.


\section robustprimaldual Robust primal and dual problems

The robust standard linear programming problem is
\f[ \text{(RP)} \qquad \min c^Tx \text{ s.t. } Ax=b;\ x>0.  \f]
and additionally find column indices \f$B\f$ such that \f$A_B\f$ is nonsingular.
The dual to the standard problem is
\f[ \text{(RD)} \qquad \max b^T y \text{ s.t. } A^Ty<c \f]
or alternatively
\f[ \text{(RD)} \qquad \max b^T y \text{ s.t. } A^Ty+z=c;\ z>0. \f]


\section feasibleprimaldual Primal and dual feasibility problems

The primal feasibility problem is
\f[ \text{(PF)} \qquad Ax=b;\ x\geq0.  \f]
with dual problems
\f[ \text{(DPF)} \qquad A^Ty \leq0;\ -b^Ty\leq-\!1 \quad \text{ or } \quad \max b^Ty \text{ s.t. } A^Ty\leq 0.  \f]
Problem (PF) is solvable iff the first form of (DPF) is unsolvable or the optimal value in the second form is strictly positive.

The dual feasibility problem is
\f[ \text{(DF)} \qquad A^Ty\leq c \f]
which has dual
\f[ \text{(DDF)} \qquad Ax=0,\ c^Tx=-\!1;\ x\geq0 \quad \text{ or } \quad \min c^Tx \text{ s.t. } Ax=0,\ x\geq 0. \f]
Problem (DF) is solvable iff the first form of (DDF) is unsolvable or the optimal value in the second form is negative.

\section robustprimaldual Robust feasibility certificates

A <em>robust feasibility certificate</em> for the primal problem is a base \f$B\f$ such that \f$A_B\f$ is nonsingular, and a vector \f$x_N>0\f$ such that \f$-A_B^{-1}A_Nx_N>0\f$, and an infeasibility certificate is a vector \f$y>0\f$ such that \f$A^Ty<c\f$ and \f$b^Ty>0\f$.

A <em>robust feasibility certificate</em> for the dual problem is a vector \f$y\f$ such that \f$A^Ty<c\f$, and an infeasibility certificate is a vector \f$x_N\f$ such that \f$-A_B^{-1}A_Nx_N>0\f$ and \f$c^Tx<0\f$.


\section dualcertificate Certificates of infeasibility / Farka's Lemma

A certificate of infeasibility of \f$Ax=b;\ x\geq0\f$ is a vector \f$y\f$ such that \f$A^Ty\leq0\f$ and \f$b^Ty>0\f$.
For then \f$0\geq y^TAx = y^Tb > 0\f$, a contradiction.

A certificate of infeasibility of \f$A^Ty\leq c\f$ is a vector \f$x\geq0\f$ such that \f$Ax=0\f$ and \f$c^Tx<0\f$.
For then \f$0=y^TAx \leq c^Tx < 0\f$, a contradiction.

\section robustcertificate Certificates of robust (in)feasibility

A robust certificate for the primal feasibility problem \f$Ax=b;\ x\geq0\f$ is a base \f$B\f$ such that \f$A_B\f$ is nonsingular, and an \f$x_N>0\f$ such that \f$-A_B^{-1}A_Nx_N>0\f$; this is equivalent to \f$Ax=b;\ x>0\f$.
A robust certificate of infeasibility is a point \f$y\f$ such that \f$A^Ty<0\f$ and \f$b^Tx>0\f$.

<b>Theorem</b>
Suppose \f$Ax=b;\ x\geq0\f$ is robust. Then either \f$A\f$ has full row rank and there exists \f$x>0\f$ such that \f$Ax=b\f$, or there exists \f$y\f$ such that \f$A^T<0\f$ and \f$b^Ty>0\f$.

A robust certificate for the dual feasibility problem \f$A^Ty\leq c\f$ is a point \f$y\f$ such that \f$A^Ty<c\f$.
A robust certificate of infeasibility is a base \f$B\f$ such that \f$A_B\f$ is nonsingular, and a vector \f$x>0\f$ such that \f$Ax=0\f$ and \f$c^Tx<0\f$.
We prove \f$Ax=0\f$ by setting \f$x_B=-A_B^{-1}A_Nx_N\f$, so that \f$Ax=A_Bx_B+A_Nx_N=-A_BA_B^{-1}A_Bx_N+A_Nx_n=0\f$.

<b>Theorem</b>
Suppose \f$A^Ty\leq c\f$ is robust. Then either there exists \f$y\f$ such that \f$A^Ty<c\f$, or \f$A\f$ has full row rank and there exists \f$x>0\f$ such that \f$Ax=0\f$ and \f$c^Tx<0\f$.

<i>Proof</i> 
If \f$A^Ty<c\f$, then this holds also for perturbations of \f$y\f$.
If \f$Ax=0\f$ and \f$A_B\f$ is nonsingular, then \f$x_B=-x_NA_NA_B^{-1}\f$.
Perturbing \f$A,b\f$ and keeping \f$x_N\f$ constant, we obtain a perturbation of \f$x_B\f$, and hence a certificate for the perturbed problem.
<br>
Conversely, suppose the problem is robustly solvable.
Then the problem \f$A^Ty\leq c-\epsilon p\f$ is solvable for \f$p>0\f$ and \f$\epsilon\f$ sufficiently small. Hence there exists \f$y\f$ such that \f$A^Ty<c\f$.
<br>
Suppose the problem is robustly unsolvable
Let \f$P\f$ be a matrix with all positive entries. Since \f$A^Ty\leq c\f$ is robustly unsolvable, \f$(I+\epsilon P)^T A^Ty\leq (I+\epsilon P)^T c\f$ is unsolvable for some \f$\epsilon>0\f$. Then there exists \f$x_\epsilon\f$ such that \f$A(I+\epsilon P)x_\epsilon=0\f$, \f$x_\epsilon\geq0\f$ and \f$(I+\epsilon P)x_\epsilon c<0\f$. Then if \f$x=(I+\epsilon P)x_\epsilon\f$, then \f$Ax=0\f$, \f$c^Tx<0\f$ and \f$x>0\f$ since \f$x_\epsilon\geq0;\ x_\epsilon\neq0\f$ and \f$(I+\epsilon P)>0\f$.


\section robustsolve Converting robust feasibility problems to feasibility problems

To solve the robust primal feasibility problem,
\f[ \text{(RPF)} \qquad Ax=b;\ x>0 \f]
we choose \f$p>0\f$ and consider the problem
\f[ \min -s \text{ s.t. } Ax + Ap\,s = b;\ x,s\geq0 .\f]
Let \f$\hat{x}^T=(x\;s)\f$, \f$\hat{A}=(A\;Ap)\f$ and \f$\hat{c}^T=(0\;\mbox{}-\!1)\f$.
Then we obtain the standard primal optimisation problem 
\f$ \min \hat{c}^T\hat{x} \text{ s.t. } \hat{A}\hat{x} = b;\ \hat{x}\geq0 . \f$
If the optimal value is negative, then we have found \f$x^*,s^*\f$ such that \f$A(x^*+ps^*)=b; \ x^*\geq0,\ s^*>0\f$, so taking \f$\tilde{x}=x^*+ps^*\f$, we have \f$A\tilde{x}=b;\ \tilde{x} = x^*+s^*p \geq s^*p > 0\f$.
If the optimal value is non-negative, then we can attempt to solve the dual robust optimisation problem
\f$ \max b^Ty \text{ s.t. } Ay<0 \f$.
Since if this problem is feasible, it has unbounded solutions, we can instead choose \f$q>0\f$ look for a positive optimal value of
\f[ \max b^Ty \text{ s.t. } Ay\leq-q . \f]

To solve the robust dual feasibility problem, 
\f[ \text{(RDF)} \qquad A^T y + z = c; \ z>0 \quad \text{or} \quad  A^Ty < c, \f]
we choose \f$q>0\f$ and consider the problem
\f[ \max t \text{ s.t. } A^T y + qt \leq c . \f]
Let \f$ \hat{A}^T = (A^T \; q)\f$, \f$\hat{b}^T = (0^T \; 1)\f$ and \f$\hat{y}^T = ( y^T\;t )\f$.
Then we obtain the standard dual optimisation problem 
\f$ \max \hat{b}^T \hat{y} \text{ s.t. } \hat{A}^T \hat{y} \leq c; \f$
If the optimal value is positive, then we have found \f$y^*,t^*\f$ such that \f$A^Ty^*\leq c-t^*q < c\f$.
If the optimal value is zero or negative, then we can attempt to solve the primal robust feasibility problem 
\f$ \min b^T x \text{ s.t. } A x = 0, \ q^T x = 1, \ x\geq 0 . \f$
However, even if the dual feasibility problem is unsolvable (i.e. \f$t^* < 0\f$), the primal may become solvable by a perturbation of \f$A\f$.
We therefore consider a robust version
\f$ \min b^T x \text{ s.t. } A x = 0,\ x>0 . \f$
and introduce \f$p>0\f$ to make a problem
\f$ \min b^T x \text{ s.t. } A x = 0,\ x-p\geq0 . \f$
Note that if this problem has negative value, then we have found \f$x^*\f$ such that \f$b^Tx^*<0\f$, \f$Ax^*=0\f$ and \f$x^*>0\f$, which implies that the original dual problem has no solution.
Taking \f$\tilde{x} = x-p\f$, we obtain
\f[ \min b^T \tilde{x} + b^T p \text{ s.t. } A\tilde{x} = -Ap,\ \tilde{x}\geq 0 . \f]



\section reducedform Reduced linear programming problem

The reduced form of the standard linear programming problem is
\f[ \min c^Tx \text{ s.t. } x_B+\tilde{A}x_N=b;\ x\geq0 . \f]
Here, \f$x_B\f$ are the basic variables, and \f$x_N\f$ the non-basic variables.

If \f$A\f$ is of full row rank, by choosing a basis such that the basis matrix \f$A_B\f$ is invertible, we can put any linear programming problem in standard form.


\section boundedoptimisation  Lower and upper bounds on variables

The constrained primal linear programming problem is
\f[ \text{(CP)} \qquad \min c^Tx \text{ s.t. } Ax=b;\ l\leq x\leq u . \f]
Given lower bounded variables \f$x_L\f$ and upper bounded variables \f$x_L\f$, the problem becomes 
\f[ \min c^Tx \text{ s.t. } Ax=b;\ x_L\geq l_L;\ x_U\leq u_U\f]
and the dual problem is
\f[ \text{(DCP)} \qquad \max\ (b^T - l_L^TA_L^T - u_U^T A_U^T)y+l_L^Tc_L + u_U^Tc_U \ \text{ s.t. }\  A_L^Ty\leq c_L;\ A_U^Ty\geq c_U . \f]
Note that the objective function can be written \f$b^T y +l_L^T (c_L - A_L^Ty) + u_U^T(c_U - A_U^Ty)\f$

<b>Theorem</b>
Suppose \f$Ax=b;\ l\leq x\leq u\f$ is robust. Then either \f$A\f$ is nonsingular, and there exists \f$x\f$ such that \f$Ax=b\f$ and \f$l<x<u\f$, or there exists column index sets \f$L,U\f$ and \f$y\f$ such that \f$ (b^T-l_L^TA_L^T-u_U^TA_U^T)y>0,\ A_L^Ty<0,\ A_U^Ty>0 . \f$

To solve the robust problem, we first find a solution to the standard problem, and then try to force saturated constraints to be positive as before.


\section simplexalgorithm The simplex algorithm

  Suppose we wish to update a basis of the standard linear programming problem.
   - The current point \f$v=\mathrm{A}_B^{-1} (b - \mathrm{A}_N x_N)\f$ (typically, \f$x_N=0\f$).
   - The reduced costs are \f$c_N-(c_B\mathrm{A}_B^{-1})\mathrm{A}_N\f$.
   - If the reduced costs are all positive, the algorithm terminates. Otherwise, select \f$j\f$ such that \f$c_j<0\f$.
   - The direction to move is \f$d=\mathrm{A}_B^{-1}a_j\f$.
   - Choose \f$t\f$ maximal so that \f$l_B \leq v-td\leq u_B\f$; if the update is being used for feasibility, constraints violated by \f$v\f$ may be violated by \f$v-td\f$. Choose \f$k\f$ such that \f$i=\pi_k\f$ corresponds to a saturated constraint.
   - Replace \f$x_j\f$ by \f$x_k\f$ in the basis and update \f$\mathrm{B}:=\mathrm{A}_B{-1}\f$.
     - We have \f$\mathrm{B} \mathrm{A}_{\pi_i}=e_i\f$ for \f$i\neq k\f$ and we want \f$\mathrm{B} \mathrm{A}_{\pi_k}=e_k\f$.
     - Let \f$\mathrm{B}\mathrm{A}_{\pi_k} = a\f$.
     - For \f$i\neq k\f$, subtract \f$\mathrm{B}_{kj}\,a_i/a_k\f$ from \f$\mathrm{B}_{ij}\f$ for all \f$j\f$.
     - Then divide \f$\mathrm{B}_{kj}\f$ by \f$a_k\f$ for all \f$j\f$.


\section simplexefficiency Efficiency of the simplex algorithm

For a linear programming problem of standard form, with \f$A\f$ an \f$m\times n\f$ matrix, the number of iterations of the simplex algorithm for practical problems grows approximately as \f$m\log n\f$.

\section feasibilityalgorithms Algorithms for feasibility 

 - Constrained feasibility problem with equalities  \f$ Ax=b;\ l\leq x\leq u\ (m\leq n)\f$.<br>
   Find a set of basic variables \f$B\f$ so that \f$ A_B\f$ is nonsingular, where \f$ A_B\f$ is the matrix formed from the columns of \f$A\f$ in \f$B\f$.
   Initialise \f$x_N\f$ to \f$l\f$ for non-basic variables, and set \f$x_B=A_B^{-1}(b-A_Nx_N)\f$. Then \f$Ax=b\f$, but possibly not \f$l_B\leq x_B\leq u_B\f$.
   <br>
   Let \f$c_i=-1\f$ if \f$x_i<l_i\f$ and \f$c_i=-1\f$ if \f$x_i>u_i\f$.
   Now minimise \f$c^Tx\f$, but relax the currently violated constraints.

   \b Remark: Since we do not assume the existence of \f$\pm\infty\f$ in our number types, we use \f$l=0,\ u=-1\f$ for the constraint \f$x\geq0\f$; this is the only unbounded constraint we allow.

   See Chvatal [Chapter 8, pp 129] for more details.
 
 - Unconstrained feasibility problem with inequalities  \f$ Ax\leq b;\ l\leq x\leq u\ (m\geq n)\f$.<br>
   Let \f$I\f$ be a set of basic indices such that \f$A_I\f$ is invertible, where \f$A_I\f$ is the matrix formed from the \em rows corresponding to the \f$I\f$.
   Set \f$v=A_I^{-1}b_I\f$, and let \f$S\f$ be the constraints satisfied by \f$v\f$.
   Find \f$i\f$ such that \f$a_iv>v\f$ where \f$a_i\f$ is the \f$i^\mathrm{th}\f$ row of \f$A\f$ and let \f$c=a_i\f$.
   Change basis until either \f$a_i v\leq b\f$ or it is impossible to reduce \f$a_i x\f$ without violating \f$A_Sx\leq b_S\f$.





\section geometricfeasibility Feasibility problems for geometric operations

In the Geometry module, we need to solve the following linear programming problems to test intersection.
\f[ \begin{array}{|l||c|c|c|c|}\hline
      &\text{Polyhedron}&\text{Polytope}&\text{Zonotope}\\\hline\hline
      \text{Point} & Ap\leq b & p=Vs;\ 1\!\cdot\!s=1;\ s\geq0 & p=c+Ge;\ -1\leq e\leq1 \\\hline
      \text{Rectangle} & Ax\leq b;\ l\leq x\leq u & x=Vs;\ 1\!\cdot\!s=1;\ l\leq x\leq u;\ s\geq0 & x=c+Ge;\ l\leq x\leq u; \ -1\leq e\leq1 \\\hline
      \text{Zonotope} & A(c+Ge)\leq b;\ -1\leq e\leq 1 & Vs=c+Ge;\ 1\!\cdot s=1;\ -1\leq e\leq1;\ s\geq0 & c_1+G_1e_1=c_2+G_2e_2;\ -1\leq e_1,e_2\leq1 \\\cline{0-3}
      \text{Polytope} & AVs\leq b;\ 1\!\cdot\!s=1;\ s\geq0 & V_1s_1=V_2s_2;\ 1\!\cdot s_1=1;\ 1\cdot s_2=1;\ s_1,s_2\geq0 \\\cline{0-2}
      \text{Polyhedron} & A_1x\leq b_1;\ A_2x\leq b_2 \\\cline{0-1}
    \end{array}
\f]
We notice that by introducing slack variables, we can convert all problems into a standard linear programming problem with constraints.
 - Standard primal feasibility problem \f$ Ax=b;\ x\geq 0\f$ 
 - Constrained primal feasibility problem \f$ Ax=b;\ l\leq x\leq u\f$ 

 - Standard dual feasibility problem \f$ Ax\leq b\f$
 - Constrained dual feasibility problem \f$ Ax\leq b;\ l\leq x\leq u\f$ 

We can convert the standard dual feasibility problem into a primal linear programming problem \f$\min b^Ty\text{ s.t. } A^Ty=0\f$, but it is not so straightforward to convert a constrained dual feasibility problem into its dual. Instead we add slack variables and solve
\f$ Ax+z=b;\ l\leq x\leq u\f$. We can use the reduced simplex algorithm to take advantage of sparseness.

\section ariadnelpsolvers Linear programming solvers provided by Ariadne.

 - lpstp() Perform one step of the standard linear programming problem.
      Input: \f$\mathrm{A},b,c\f$, InOut: \f$\pi,\mathrm{A}_B^{-1}\f$.

 - lpcstp() Perform one step of the standard linear programming problem with constraints \f$l\leq x\leq u\f$.
     Any constraints which are violated are assumed to remain violated; this allows for constraints with infinities.
     Input: \f$\mathrm{A},b,c\f$, InOut: \f$\pi,\mathrm{A}_B^{-1}\f$.

 - lpupd() Update the matrix \f$\mathrm{A}_B^{-1}\f$ so that \f$\mathrm{A}_B^{-1}a=e_i\f$ by pivoting on the \f$i^\textrm{th}\f$ row.
      Input: \f$a,i\f$; InOut\f$\mathrm{A}_B^{-1}\f$.

 - lpslv() Solve the standard linear programming problem \f$\min c^Tx \text{ s.t. }Ax=b;\ x\geq0\f$.
      Input: \f$A,b,c\f$; Output: \f$\pi,\ \mathrm{A}_B^{-1},\ x^*,\ y^*,\ z^*\f$.

 - lpcslv() Solve the standard linear programming problem with constraints \f$\min c^Tx \text{ s.t. } Ax=b;\ l\leq x\leq u\f$.
      Input: \f$A,b,c\f$; Output: \f$\pi,\ \mathrm{A}_B^{-1},\ x^*,\ z^*\f$.

 - lprslv() Solve the standard linear programming problem with constraints given in reduced form \f$\min \tilde{c}^Tx_N \text{ s.t. } \tilde{A}x_N+x_B=b;\ l\leq x\leq u\f$ \f$\tilde{A}x_N+x_B=b\f$. (This is useful if a starting basis can easily be found, and there are almost as many constraints as variables.)
      Input: \f$\tilde{A},b,c,\pi\f$; Output: \f$\pi, x^*, z^*, y^*\f$.

 - lpfeas() Solve the feasibility problem \f$Ax=b;\ x\geq0\f$.

 - lpcfeas() Solve the feasibility problem \f$Ax=b\f$ with \f$l\leq x\leq u\f$.

 - lpdfeas() Solve the dual feasibility problem \f$Ax\leq b\f$ directly.




\page geometric Geometric Representation

\deprecated The material in this section is out-of-date

\section Introduction

One of the main concerns of Ariadne is the computation of points and sets in Euclidean space.
Unfortunately, the set of points in Euclidean space is uncountable (it has continuum cardinality) and so there
is no way of modelling all points in Euclidean space using binary words. The set of subsets of Euclidean space
is also uncountable, and has an even greater cardinality!

To represent an uncountable set, we can use approximations.
We take a set of elements which we model exactly, such as \a double reals or \a rational reals, which we call <em>denotable</em> elements.
This denotable set may be finite or countably infinite, the former being denoted by words of a fixed length, and the latter by arbitrary length words.
An <em>approximation</em> to an element is then a pair consisting of a denotable element and an <em>error</em>.
Ideally, we wish to be able to approximate to arbitrary precision;
this means that the set of denotable elements must be <em>dense</em> in the set of all elements.
We can then represent arbitrary elements be a <em>convergent sequence</em> of approximations.
In terms of a concrete representation on a computer, we can think of a such a sequence as a neverending <em>stream</em> of data.

The material in this section is heavily influenced by the book "Computational Analysis" by Klaus Weihrauch.

\section state Representation of points in Euclidean space

Points in Euclidean space can be represented by the templated class \a DenotablePoint<T>, where \a T is a numeric type such as double or rational. From this, we can define a PointApproximation class, which gives an approximation of a state, plus an error bound in terms of the sup norm or Euclidean norm.

\code
concept State
{
  typename real_type; 

  State(const State &);
  State& operator=(const State &);
 
  real_type operator[] (size_type) const;
};
\endcode 

\section basicset Basic sets.

A \c BasicSet provides a building block for representing more complicated sets.
Mathematically, a \c BasicSet type represents elements of a countable base of a topological space.
More precisely, a \c BasicSet represents the \em closure of a basic set for the topology.

In many cases, it is not possible to evaluate geometric predicates concerning two sets using a given real number type. For this reason, the results of a geometric predicate returns an object of type \a tribool, which may be \a true, \a false or \a indeterminate (unknown).
 


\code
// The basic set concept.
concept BasicSet
{
  type real_type; // The type of denotable real number used for the representation.
  type state_type; // The type of denotable point the set contains.

  BasicSet(const std::string &); // Construct from a string literal (optional).
  BasicSet(const BasicSet &); // Copy constructor.
  BasicSet & operator=(const BasicSet &); // Assignment operator.

  dimension_type dimension() const; // The dimension of the set.
  state_type centre() const; // A point in the set (typically, the "centre" point, if this makes sense).
  real_type radius() const; // The maximum distance from the centre to another point in the set in an appropriate metric. 
  real_type volume() const; // An approximation to the volume of the set. (Optional).

  tribool empty() const; // Tests if the set is empty.
  tribool bounded() const; // Tests if the set is bounded.

  tribool contains(const State &) const; // Tests if the set contains a point.

  Rectangle<real_type> bounding_box() const; // A rectangle containing the set.
  tribool disjoint(const Rectangle &); // Tests if the set is disjoint from a Rectangle
  tribool superset(const Rectangle &); // Tests if the set contains a Rectangle.

};

 tribool equal(const BasicSet1 &, const BasicSet2 &); // Returns indeterminate if the two sets are equal; if false, then they are not equal.
 tribool disjoint(const BasicSet1 &, const BasicSet2 &); // If true, then the sets are disjoint; if false, then they robustly intersect.
 tribool subset(const BasicSet1 &, const BasicSet2 &); // If true, then the first set is a subset of the second interior of the second; 
                                                     // if false, then the first set is not a subset of the second.

 // Optional, depending on whether the operation yields a basic set of the same type.
 BasicSet open_intersection(const BasicSet &, const BasicSet &); // The closure of the intersection of the interiors of the two sets.
 BasicSet closed_intersection(const BasicSet &, const BasicSet &); // The intersection of the two (closed) sets.
 BasicSet convex_hull(const BasicSet &, const BasicSet &); // The convex hull of the two sets.
 BasicSet minkowski_sum(const BasicSet &, const BasicSet &); // The Minkowski (pointwise) sum of the two sets.
 BasicSet minkowski_difference(const BasicSet &, const BasicSet &); // The Minkowski (pointwise) difference of the two sets.

\endcode

Classes fulfilling the BasicSet concept are \ref Rectangle (or Cuboid), Simplex, Parallelotope, Zonotope, Polytope, Polyhedron, Sphere and Ellipsoid.
Actually, these are templates, parameterised by the real number type real_type.


\section denotable_set Denotable Sets

A DenotableSet implements a set as a union of basic sets type \c DenotableSet::basic_set_type.
\code
concept DenotableSet
{
  type real_type;
  type state_type;
  type basic_set_type;

  type const_iterator; // Must satisfy the requirements of a ForwardIterator.

  // No default constructor required.

  DenotableSet(const DenotableSet &);
  DenotableSet & operator=(const DenotableSet &);

  // No equality operator required.

  // Set-theoretic operations
  dimension_type dimension() const;
  tribool empty() const;
  tribool contains(const state_type &) const;

  Rectangle<real_type> bounding_box() const; // Optional.

  void adjoin(const basic_set_type &);
  void adjoin(const DenotableSet &);

  // List operations
  const_iterator begin() const;
  const_iterator end() const;

  size_type size() const; // Only required if the iterator is a RandomAccessIterator.
  basic_set_type operator[] (size_type) const; // Only required if the iterator is a RandomAccessIterator.

  void push_back(const basic_set_type &); // Only used if the DenotableSet is an ordered list. (Optional)
  basic_set_type pop_back(); // Only used if the DenotableSet is an ordered list. (Optional)

  void insert(const basic_set_type &); // Only used if the DenotableSet is an unordered or sorted list. (Optional)
  void remove(const basic_set_type &); // Only used if the DenotableSet is an unordered or sorted list. (Optional)
};

tribool subset(const BasicSet &, const DenotableSet &); // Optional, but highly recommended.

tribool disjoint(const DenotableSet &, const DenotableSet &);
tribool subset(const DenotableSet &, const DenotableSet &);

DenotableSet join(const DenotableSet &, const DenotableSet &);
DenotableSet open_intersection(const BasicSet &, const DenotableSet &); // Optional.
DenotableSet difference(const DenotableSet&, const DenotableSet&); // Optional

\endcode

\section set_approximation Approximating Sets

%Ariadne provides operators for approximating sets. All the operators have one
of the following forms.

\code
Result outer_approximation(Argument,Error); // postcondition: inner_subset(Argument,Result)
Result over_approximation(Argument,Error);  // postcondition: subset(Argument,Result)
Result lower_approximation(Argument,Error); // postcondition: 
Result under_approximation(Argument,Error); // postcondition: subset(Result,Argument)
Result inner_approximation(Argument,Error); // postcondition: inner_subset(Result,Argument)
\endcode

The error specification depends on the type of approximation used. 
 - When approximating by a GridSet, the error is determined by the grid.
 - When approximating by a PartitionTreeSet, the error is determined by the 
      partition scheme and the depth.
 - When approximating by a ListSet, the error is a metric error bound.








\page GeometricOps Geometric Operations on Basic Sets.

The core geometric types used by %Ariadne to represent are Rectangle, Zonotope, Polytope (described by generators)
and Polyhedron (described by constraints). 
The core geometric operations are contains(A,p), subset(A,B) and disjoint(A,B) (equivalent to intersects(A,B) ) .

  - Rectangle representation: \f$l\leq x\leq u\f$.
  - Zonotope representation: \f$x=c+Ge,\ -1\leq e\leq1\f$.
  - Polytope representation: \f$x=Vs,\ \sum s_i=1,\ s\geq0\f$.
  - Polyhedron representation: \f$Ax\leq b\f$.


Alternative representations for polytopes and polyhedron are given by taking \f$\hat{x}=(x,1)\f$.
  - Rectangle representation: \f$x_i\in[l_i,u_i]\f$.
  - Zonotope representation: \f$x = \hat{G}\hat{e},\ e_0=1,\ -1\leq e\leq1\f$.
  - Polytope representation: \f$\hat{x} = Rs\f$, \f$s\geq0\f$.
  - Polyhedron representation: \f$C\hat{x}\geq0\f$.

We sometimes use the based polytopic representation
  - Polytope representation: \f$x=v+Ws,\ \sum s_i\leq1,\ s\geq0\f$.

Henceforth, we shall always assume \f$-1\leq e\leq1\f$ and \f$1\!\cdot\!s=1;\ s\geq0\f$.


\section geometricpreprocessing Preprocessing sets

To simplify certain computations, it may be useful to pre-process the matrices describing a Zonotope, Polytope or Polyhedron.

Given an \f$m\times n\f$ matrix \f$A\f$ of rank \f$k\f$, find an index set \f$I\subset\{0,1,\ldots,n\!-\!1\}\f$ of cardinality \f$k\f$, and an \f$k\times m\f$ matrix \f$B^{-1}\f$ such that
\f[ B^{-1} A_I = \mathbf{I} \f]
where \f$ A_I \f$ is the \f$m\times k\f$ matrix consisting of the columns of \f$A\f$ with index \f$j\in I\f$.


\section conversion Converting between basic sets

A %Rectangle can be directly converted to a %Polytope or %Polyhedron without using arithmetic. 
A %Rectangle can be easily converted to a %Zonotope, and a %Zonotope to a %Polytope, but these conversions require arithmetic.
Conversion between a %Polytope to a %Polyhedron can be performed using the double description algorithm.

The conversion from a %Rectangle or %Zonotope to a %Polytope or conversion between %Polytope and %Polyhedron may be of exponential complexity.

 \section contains Testing inclusion

  - contains(Rectangle,Point) : Check \f$l\leq x\leq u\f$.
  - contains(Zonotope,Point) : Solve \f$x=c+Ge;\ -1\leq e\leq1\f$.
  - contains(Polytope,Point) : Solve \f$x=Vs;\ 1\!\cdot\!s=1;\ s\geq 0\f$.
  - contains(Polyhedron,Point) : Check \f$Ax\leq b\f$.

 \section Intersection Testing intersection/disjointness

  - intersects(Rectangle,Rectangle) : Check \f$l_1\leq u_2;\ l_2\leq u_1\f$.
  - intersects(Rectangle,Zonotope) : Solve \f$x=c+Ge;\ l\leq x\leq u;\ -1\leq e\leq1\f$ or \f$l\leq c+Ge\leq u;\ -1\leq e\leq 1\f$.
  - intersects(Zonotope,Zonotope) : Solve \f$c_1+G_1e_1 = c_2+G_2e_2;\ -1\leq e_1\leq1;\ -1\leq e_2\leq1\f$.
  - intersects(Rectangle,Polytope) : Solve \f$x=Vs;\ 1\!\cdot\!s=1;\ l\leq x\leq u;\ s\geq 0\f$ or \f$l\leq Vs\leq u;\ 1\!\cdot\!s=1;\ s\leq 1\f$.
  - intersects(Zonotope,Polytope) : Solve \f$c+Ge = Vx;\ \ 1\!\cdot\!s=1;\ -1\leq e\leq 1;\ s\geq0\f$.
  - intersects(Polytope,Polytope) : Solve \f$V_1s_1=V_2s_2;\ 1\!\cdot\!s_1=1;\ 1\!\cdot\!s_2=1;\ \cdot s_1\geq0;\ s_2\geq0\f$.

  - intersects(Rectangle,Polyhedron) : Solve \f$Ax\leq b;\ l\leq x\leq u\f$.
  - intersects(Zonotope,Polyhedron) : Solve \f$Ac+AGe \leq b;\ -1\leq e\leq 1\f$.
  - intersects(Polytope,Polyhedron) : Solve \f$AVs\leq b;\ 1\!\cdot\!s=1;\ s\geq0\f$.
  - intersects(Polyhedron,Polyhedron) : Solve \f$A_1x\leq b_1;\ A_2x\leq b_2\f$.

All these problems can be solved using (a variant of) the simplex algorithm. 
When computing over-approximations of a set on a grid using disjoint(Rectangle R, ConvexSet S), it may be useful to pre-process S to simplify computation of an initial vertex for the feasibility problem.


\section Subset Testing subset

Testing whether a convex set is a subset of a polyhedron, or whether a polytope is a subset of a convex set are easy:
  - subset(ConvexSet S, Polyhedron P) : Check disjoint(S,H) for all complementary halfspaces H of P.
  - subset(Polytope P, ConvexSet S) : Check contains(S,v) for all vertices v of P.

Checking whether a Rectangle, Zonotope or Polytope is a subset of a Rectangle can be performed directly:
  - subset(Rectangle, Rectangle) : Check \f$l_1\geq l_2\ \textrm{and}\ u_1\leq u_2\f$.
  - subset(Zonotope, Rectangle) : Check \f$\sum_{j} |G_{ij}|\leq \min\{c_i-l_i,u_i-c_i\} \forall\,i\f$.
  - subset(Polytope, Rectangle) : Check \f$l\leq v\leq u\f$ for all vertices \f$v\f$ of \f$P\f$.

The following operations cannot be performed efficiently, but can be easily implemented by testing the vertices of the Rectangle, and are provided for conformance to the Set interface.
  - subset(Rectangle, Zonotope)
  - subset(Rectangle, Polytope)

The following operation is provided to help with checking zonotopic over-approximations.
  - subset(Zonotope, Zonotope)

While a Rectangle can be directly converted to both a Polytope and a Polyhedron, the conversion to a Polytope yields \f$2^d\f$ vertices.
Hence testing subset(Rectangle, ConvexSet) is in general of exponential complexity, as is subset(Zonotope, ConvexSet).
When testing subset(Rectangle R, ConvexSet S) for under-approximationg S on a grid, we need to check contains(S,v) for each grid vertex v once only, considerably reducing the amount of work.

Note that for a single test of subset(Rectangle R, Polyhedron P), we prefer to use disjoint(R,H), but for under-approximating a polyhedron P on a grid, we prefer to check contains(P,v) for the grid vertices.


\section PolyhedralConversion Converting between a polyhedron and a polytope.
 
Augment the state by \f$\hat{x}=(x,1)\f$.

  - Define a Polytope by generators \f$x=\lambda g\f$, given as columns of the augmented generator matrix \f$G'\f$.
  - Define a Polyhedron by constraints \f$a\cdot x\geq 0\f$, given as rows of the augmented constraints matric \f$A'\f$.

Construct the <em>saturation matrix</em> by 
  - Generator \f$g\f$ \e violates constraint \f$a\f$ if \f$a\cdot g<0\f$, 
  - Generator \f$g\f$ \e saturates constraint \f$a\f$ if \f$a\cdot g=0\f$.
  - Generator \f$g\f$ \e satisfies constraint \f$a\f$ if \f$a\cdot g=0\f$.
Saturation matrix \f$S=\mathrm{sgn}(AG)\f$.

Generators are \em adjacent if the corresponding columns of the saturation matrix
differ only in one row.

\section zonotope Zonotopic reduction methods

Throughout this sections, we use the supremum norm on \f$R^n\f$, and the corresponding operator norm on \f$\mathbb{R}^{m\times n}\f$.

Given a zonotope \f$ Z=\{ c+Ae \mid ||e||\leq 1 \}\subset \mathbb{R}^n\f$, where \f$A\in \mathbb{R}^{n\times p}\f$, 
we wish to compute a zonotope \f$Z' = \{ c + A' e' \mid ||e'||\leq 1\}\f$ with fewer generators 
i.e. \f$A'\in \mathbb{R}^{n\times p'}\f$ with \f$p'<p\f$.
The general reduction method is to choose \f$ A'\f$ such that \f$ A = A' B\f$ with \f$||B||\leq 1\f$.
The key to zonotopic reduction is to choose a method with good properties.

A simple criterion to note is that if \f$\sum_{j=1}^{p} |b_{ij}|<1\f$ for some \f$i\f$, then we can improve the approximation by taking
\f$D=\mathrm{diag}(d_{i})\f$ with \f$d_{i}=\sum_{j=1}^{p} |b_{ij}|\f$, 
and \f$B'= D^{-1}B\f$ which has \f$b'_{ij}=b_{ij}/\sum_{k=1}^{p} |b_{ik}|\f$.

It is clear that if the rows of \f$B\f$ are close to a set of mutually orthogonal coordinate vectors, then the approximation is good, 
since the image of \f$B\f$ is close to the unit ball. 


\subsection interval_zonotope Interval zonotopic reduction

An <em>interval zonotope</em> is a set of the form \f$ \{ y = c + A e \mid c\in R,\ A\in\mathcal{A} \text{ and } ||e||\leq1 \} 
= R + \mathcal{A} B\f$.
To reduce an interval zonotope, we first write \f$ R = \{ c + Be\mid ||e||\leq 1 \}\f$ and combine this in \f$\mathcal{A}\f$.
To reduce \f$ \mathcal{A} \f$, write \f$ \mathcal{A} = A \mathcal{B} \mathcal{A} \f$ where \f$ A\in\mathcal{A}\f$ and \f$A\mathcal{B}\ni I\f$.
Then take \f$ \mathcal{C} = \mathcal{B} \mathcal{A} \f$ and \f$ || \mathcal{C} || 
  = \sup_{i]1}^{n} \sum_{j=1}^{p} \max |\mathcal{C}_{ij}| \f$, 
where \f$ \max |\mathcal{C}_{ij}| = \max\{ |x| \mid x\in \mathcal{C}_{ij}\f$.


 
\page function Function Evaluation

%Ariadne is primarily a module for set-based computation.
For this reason, functions are best defined by their actions on sets.
However, we sometimes also want to compute function values on points, and to evaluate real-valued functions.
For this reason, we also allow definition of functions on points.

We distinguish between computations on \em fixed-precision and \em multiple-precision types, 
between computations on \em points and \em sets, and between \em exact and \em approximate computations.

The basic computation on sets is to compute \em over-approximations to the image of <em>basic sets</em>
From these computations, arbitrarily-accurate computations can be performed.
The basic computation on points can be \em exact and \em approximate, as appropriate.

\section set_functions  Computations on sets.

A valid arbitrary-precision computation on sets is defined as follows
If \f$f\f$ is a mathematical function, and \c BS is a basic set type,
then a valid representation \c f of \f$f\f$ is a function which, 
for any basic set \f$A\f$, returns a set \f$B\f$ such that \f$f(A)\subset B\f$,
and such that whenever \f$A_n\f$ is a decreasing sequence of sets with \f$\bigcap_{n=1}^{\infty} A_n=\{x\}\f$,
then \f$\bigcap_{n=1}^{\infty} B_n=\{y\}\f$, where \f$y=f(x)\f$.

A valid fixed-precision computation on sets is an over-approximation. 
In other words, a valid representation \c f of \f$f\f$ is a function which,
for any basic set \f$A\f$, returns a set \f$B\f$ such that \f$f(A)\subset B\f$.
No guarentees on the accuracy are required.
Note that it does not make sense to consider a sequence \f$A_n\f$ converging to a point for fixed-precision types. 

\section point_functions Computations on points.

If a denotable state type is invariant under a class of functions (e.g. polynomial functions on a ring),
then the image of a point is given exactly. Otherwise, a \a fuzzy \a point is given,
which is a point defined using interval coefficients. A Point< Interval<R> > can be automatically 
converted to a Rectangle<R>.

Note that even if \c f is exact, it is impossible to compute 
the action of \f$f\f$ on a set just from the action of \c f on points,
unless a modulus of continuity for \f$f\f$ is known.


\section real_functions Computations on real numbers.

Most continuous functions used in science and engineering are built up from elementary real-valued functions,
either arithmetical or defined using Taylor series expansions. 
For this reason, Ariadne provides extended operations for computation on real-valued functions.

Arbitrary-precision computations may be exact or approximate. 
The function \c f(x) computes \f$f(x)\f$ exactly, if possible. 
If exact computation is not possible, then \c f(x) either returns an interval,
or gives an error at compile time.

Although it is not, in general, possible to perform evaluation of functions on sets from their definitions on points,
in many cases such a computation can be extracted. 
In particular, we can construct interval computations from pointwise computations in the following cases:
<ul>
   <li>The function is Lipschitz with a known Lipschitz constant e.g. \f$\sin,\ \cos\f$.</li>
   <li>The function is monotone e.g. \f$\exp\f$.</li>
   <li>The function is piecewise-monotone with known branches e.g. arithmetic.</li>
</ul>


Fixed-precision computations are inherantly approximate, 
though in some cases they may happen to give exact result. 
The function \c f_approx(x) computes \f$f(x)\f$ approximately, with no control on the error bound.
The function \c f_down(x) computes a lower-approximation to \f$f(x)\f$, 
and \c f_up(x) computes an upper-approximation.

Without an error specification, \c f(x) either gives an eFor consistency with existing practise, we also allow \c fnc(x) as a valid alternative to \c fnc_approx(x).

Note that in many cases, including arithmetic and simple functions, it is possible to compute an interval \f$J\f$ containing \f$f(I)\f$ 
using \c f_down and \c f_upp. This allows an implementation of the standard set-based function \c f(I).

\section function_syntax Syntax for continuous functions and function objects.

Elementary and arithmetical functions are represented as ordinary functions
with the same syntax.

\code
// Action of functions on points
Point<ExReal> f(Point<ExReal> x); // Compute f(x) exactly.

Point< Interval<FPReal> > f(Point<ExReal> x); // Compute a fuzzy point approximation to f(x)
BasicSet<FPReal> f(Point<ExReal> x); // Compute a set containing f(x)

// Action of functions on sets
BasicSet<FPReal> f(BasicSet<FPReal> A); // Compute an over-approximation to f(A).
BasicSet<MPReal> f(BasicSet<APReal> A); // Compute a convergent over-approximation to f(A).
\endcode

 
\page integration Integration methods

\section taylor Taylor methods 
All integration methods are based on Taylor expansion of solutions curves.

\f[ \begin{array}{rl} \displaystyle
     \frac{dx}{dt} &= f(x(t)) \\
 \frac{d^2x}{dt^2} &= Df(x(t))f(x(t)) \\
 \frac{d^3x}{dt^3} &= D^2f(x(t))f(x(t))f(x(t))+Df(x(t))Df(x(t))f(x(t)) \\
 \frac{d^4x}{dt^4} &= D^3f(x(t))f(x(t))f(x(t))f(x(t)) + 4D^2f(x(t))Df(x(t))f(x(t))f(x(t)) + Df(x(t))Df(x(t))Df(x(t))f(x(t)) \\
 \frac{d^5x}{dt^5} &= D^4f(x(t))f(x(t))f(x(t))f(x(t))f(x(t)) + 7D^3f(x(t))Df(x(t))f(x(t))f(x(t))f(x(t)) \\
                   &\ \qquad + 4D^2f(x(t))D^2f(x(t))f(x(t))f(x(t))f(x(t)) + 11D^2f(x(t))Df(x(t))Df(x(t))f(x(t))f(x(t)) \\
                   &\ \qquad + Df(x(t))Df(x(t))Df(x(t))Df(x(t))f(x(t)) \end{array} \f]


\subsection euler Euler method
\f[ x_1 = x_0 + hf(x_0) \approx x(h) \f]
One-step error
\f[ ||x_1-x(h)|| = h ||f(x_0)-f(\xi)|| = \frac{h^2}{2} || Df(\xi)f(\xi) || \f]

\subsection second_order_taylor 2nd Order Taylor Method
\f[ x_1 = x_0 + hf(x_0) + \frac{h^2}{2} Df(x_0)f(x_0) \approx x(h) \f]
One-step error
\f[ ||x_1-x(h)|| = \frac{h^2}{2} || Df(x_0)f(x_0)-Df(\xi)f(\xi)|| = \frac{h^3}{6} || D^2f(\xi)f(\xi)f(\xi) + Df(\xi)Df(\xi)f(\xi) || \f]

\subsection second_order_rk 2nd Order Runge-Kutta Method

\f[ x_1 = x_0 + \frac{h}{2}\left( f(x_0)+f(x_0+hf(x_0))\right) \approx x_0 + hf(x_0) + \frac{h^2}{2} Df(x_0)f(x_0) 
     + \frac{h^3}{4} D^2f(\xi)f(\xi)f(\xi) \f]



\page references References

\section computable_analysis_references Computable Analysis

Ker-I Ko <em>Complexity Theory of Real Functions</em>, Birkh\"aser, 1991, ISBN 3-7643-3586-6.

Klaus Weihrauch, <em>Computable Analysis</em>, Springer, 2000.

\section interval_references Interval Arithmetic
Ramon E. Moore, <em>Methods and applications of interval analysis</em>,
  SIAM Studies in Applied Mathematics, 2.
Society for Industrial and Applied Mathematics (SIAM), Philadelphia, Pa., 1979. xi+190 pp. ISBN 0-89871-161-4 

Baker R. Kearfott, "Interval computations: introduction, uses, and resources",
  <em>Euromath Bull.</em> <b>2</b> (1996), no. 1, 95--112. 

Marcel Gavriliu, "Towards more efficient interval analysis: corner forms and a remainder Newton method", 
  Ph.D. Thesis, California Institute of Technology, 2005. <br>

R. Krawczyk, "A class of interval-Newton-operators",
  <em>Computing</em> <b>37</b> (1986), no. 2, 179--183.

Arnold Neumaier, <em>Interval methods for systems of equations</em>,
  Encyclopedia of Mathematics and its Applications, 37.
  Cambridge University Press, Cambridge, 1990. xvi+255 pp. ISBN 0-521-33196-X 

A. Neumaier, "The wrapping effect, ellipsoid arithmetic, stability and confidence regions",
  <em>Computing Supplementum</em> <b>9</b> (1993), 175-190.

\section linear_programming_references Linear Programming

V. Chv\'atal, <em>Linear Programming</em>, Freeman, 1983, ISBN 0-7167-1195-8.

\section zonotope_references Zonotopes

K. Fukuda, "From the zonotope construction to the Minkowski addition of convex polytopes", Preprint, 2003. 

Leonidas J. Guibas, An Nguyen and Li Zhang, "Zonotopes as bounding volumes", Preprint.

Ari Ingimundarson, Jose Manuel Bravo, Vicenc Puig and Teodoro Alama, 
  "Robust Fault Diagnosis using Parallelotope-based Set-membership Consistency Tests",
  In <em>Proceedings of CDC-ECC 2005</em>.

\section integration_references Integration

Rudolf J. Lohner, "Enclosing the solutions of ordinary initial and boundary value problems",
  <em>Computer Arithmetic</em>, 255--286, Teubner, Stuttgart, 1987. 

Rudolf J. Lohner, "Computation of guaranteed enclosures for the solutions of ordinary initial and boundary value problems",
  <em>Computational ordinary differential equations (London, 1989)</em>,  425--435.

N. S. Nedialkov, K. R. Jackson and G. F. Corliss,
  "Validated solutions of initial value problems for ordinary differential equations",
  <em>Appl. Math. Comput.</em> <b>105</b> (1999), no. 1, 21--68.

Piotr Zgliczynski, "C^1 Lohner algorithm", <i>Found. Comput. Math.</i> <b>2</b> (2002), no. 4, 429--465.

\section reachability_references Reachability Analysis and Control

Antoine Girard, Colas Le Guernic and Oded Maler, 
"Efficient Computation of Reachable Sets of Linear Time-Invariant Systems with Inputs",

Antoine Girard, "Reachability of Uncertain Linear Systems using Zonotopes," 
  in <em>Proceedings of HSCC 2005</em>, LNCS 3414, pp 291--305, 2005.

Alex Kurzhanskiy and Pravin Varaiya, "Ellipsoidal Techniques for Reachability Analysis of Discrete-Time Linear Systems",

Alexander Kurzhanski and Pravin Varaiya, "On ellipsoidal techniques for reachability analysis", 
  <em>Optim. Methods Softw.</em> <b>17</b> (2002), no. 2, 207--237

S. V. Rakovic and D. Q. Mayne, "Set Robust Control Invariance for Linear Discrete Time Systems", 
  in <em>Proceedings of CDC-ECC 2005</em>.

F. Lydoire and P. Poignet, "Nonlinear Model Predictive Control via Interval Analysis", 
  in <em>Proceedings of CDC-ECC 2005</em>.

Dietmar Szolnoki, "Set oriented methods for computing reachable sets and control sets",
  <em>Discrete Contin. Dyn. Syst. Ser. B</em> <b>3</b> (2003), no. 3, 361--382.

\section spacial_data_structure_references Spacial Data Structures

Hanan Samet, <em>The Design and Analysis of Spacial Data Structures</em>, Addison-Wesley, 1990, ISBN 0-201-50255-0.

\section algebraic_topology_references Algebraic Topology

Tomasz Kaczynski, Konstantin Mischaikow, Marian Mrozek, <em>Computational Homology</em>, Springer-Verlag, 2004, ISBN 0-387-40853-3.

Afra J. Zomorodian, <em>Topology for Computing</em>, Cambridge University Press, 2005, ISBN 0-521-83666-2.

\section automatic_differentiation_references Automatic Differentiation

Andreas Giewank, <em>Evaluating Derivatives</em>, SIAM, 2000, ISBN 0-89871-451-6.
*/
 
#endif /* _ARIADNE_DOCUMENTATION_H/
