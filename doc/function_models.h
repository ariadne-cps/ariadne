/***************************************************************************
 *            function_models.h
 *
 *  Copyright  2007  Pieter Collins
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

/*! 

\file function_models.h
\brief Documentation on polynomial models of functions



\page function_models Function Models

A \em model for a function \f$f:\mathbb{R}^m\rightarrow \mathbb{R}^n\f$ is a polynomial approximation of \f$f\f$ on some domain \f$D\f$.
A model is described by two parameters, the \em degree and the \em smoothness.

Hence the general form of a model of \f$f\f$ with degree \f$d\f$ and smoothness \f$s\f$ is
\f[
  f^i(x) \in \sum_{0\leq|\alpha|\leq s} c_{\alpha}\,\interval{a}^{\;\!i}_{\,\alpha}\,x^\alpha
     + \sum_{s<|\alpha|\leq d} c_\alpha\,a^{i}_{\,\alpha} x^\alpha .
\f]
Here, \f$\alpha\f$ is a multi-index \f$(\alpha_1,\ldots,\alpha_m)\f$ with degree \f$|\alpha|=\alpha_1+\cdots+\alpha_m\f$, and \f$x^\alpha := x_1^{\alpha_1}\,x_2^{\alpha_2}\,\cdots\,x_m^{\alpha_m}\f$, and \f$\interval{a}\f$ denotes an interval coefficient, and \f$a\f$ a numerical coefficient. The coefficients \f$c_\alpha\f$ are constants depending on the representation of the model. Taking \f$c_\alpha=1\f$ gives a polynomial model, and taking \f$c_\alpha=\prod_{i=1}^{k} \alpha_i!\f$ gives a model where \f$a^i_{\,\alpha} = D_\alpha f^i(0)\f$ for \f$|\alpha|\leq s\f$.

The interpretation of \f$s\f$ is that the model contains sufficient information to compute \f$f\f$ and its derivatives up to order \f$s\f$.

If \f$s=-1\f$ then the model is assumed to be \em exact, i.e. a polynomial.

\section function_operation Function Operations

The following operations on functions are supported:
 - Evaluation: If \f$p\f$ is a \f$(d,s)\f$-model \f$\mathbb{R}^m\rightarrow\mathbb{R}^n\f$, then \f$p(c_1,\ldots,c_m)\f$ is an interval.
 - Partial evaluation: If \f$p\f$ is a \f$(d,s)\f$-model \f$\mathbb{R}^m\rightarrow\mathbb{R}^n\f$, then \f$p(\cdot,\ldots,c_j,\ldots,\cdot)\f$ is a \f$(d,s)\f$-model \f$\mathbb{R}^{m-1}\rightarrow\mathbb{R}^n\f$.

 - Addition/Subtraction: If \f$p\f$ and \f$q\f$ are \f$(d,s)\f$-models \f$\mathbb{R}^m\rightarrow\mathbb{R}^n\f$, then so are \f$p+q\f$ and \f$p-q\f$.
 - Multiplication: If \f$s\f$ is a scalar \f$(d,s)\f$-model and \f$p\f$ is a \f$(d,s)\f$-model, then \f$s\times p\f$ is a \f$(d,s)\f$-model.
 - Division: If \f$s\f$ is a scalar model and \f$s(x)\neq0\f$, then \f$1/s\f$ is a \f$(d,s)\f$-model.

 - Differentiation: If \f$p\f$ is an \f$(d,s)\f$-model, then the single-variable derivative \f$\frac{\partial p}{\partial x_i}\f$ is a \f$(d-1,s-1)\f$-model.

 - Direct sum (combine): If \f$p_1:\mathbb{R}^{m_1}\rightarrow\mathbb{R}^{n_1}\f$ and \f$p_2:\mathbb{R}^{m_2}\rightarrow\mathbb{R}^{n_2}\f$ are \f$(d,s)\f$ models, then so is \f$p_1\oplus p_2 : \mathbb{R}^{m_1+m_2}\rightarrow\mathbb{R}^{n_1+n_2}\f$ given by \f$(p_1\oplus p_2)(x_1,x_2):=(p_1(x_1),p_2(x_2))\f$.
 - Cartesian product (join): If \f$p_1:\mathbb{R}^{m}\rightarrow\mathbb{R}^{n_1}\f$ and \f$p_2:\mathbb{R}^{m}\rightarrow\mathbb{R}^{n_2}\f$ are \f$(d,s)\f$ models, then so is \f$p_1 \times p_2 : \mathbb{R}^{m}\rightarrow\mathbb{R}^{n_1+n_2}\f$ given by \f$(p_1\times p_2)(x):= (p_1(x),p_2(x))\f$.

 - Composition: If \f$p_1:\mathbb{R}^{m}\rightarrow\mathbb{R}^{l}\f$ and \f$p_2:\mathbb{R}^{l}\rightarrow\mathbb{R}^{n}\f$ are \f$(d,s)\f$-models, then \f$p_2\circ p_1\f$ is a \f$(d,s)\f$-model \f$\mathbb{R}^{m}\rightarrow\mathbb{R}^{n}\f$.
 - Inverse: If \f$p:\mathbb{R}^{n}\rightarrow\mathbb{R}^{n}\f$ is a \f$(d,s)\f$-model with \f$s\geq1\f$, and if \f$Dp\f$ is nonsingular, then \f$p^{-1}\f$ is a \f$(d,s)\f$ model. The inverse should be centred around \f$y=p(x)\f$, where \f$x\f$ is given, since locally \f$p^{-1}(y)\f$ may have several branches.
 - Implicit: If \f$p:\mathbb{R}^{m}\times\mathbb{R}^{n}\rightarrow\mathbb{R}^{n}\f$ is a \f$(d,s)\f$-model with \f$s\geq1\f$, and if \f$D_2p\f$ is nonsingular near \f$(x,y)\f$, then there is a \f$(d,s)\f$-model \f$q\f$ such that \f$y\in q(x)\f$ and \f$p(x,q(x))=z\f$.

 - Reduce: If \f$p:\mathbb{R}^{m}\times\mathbb{R}^{n}\rightarrow\mathbb{R}^{n}\f$ is a \f$(d_1,s_1)\f$-model, \f$d_2\leq d_1\f$ and \f$s_2\leq\min\{d_2,s_1\}\f$, then \f$p\f$ can be reduced to a \f$(d_2,s_2)\f$-model.


\section affine_models Affine Models

An affine model represents a function \f$f\f$ on \f$X\f$ by \f$f(x) = b + A (x-c)\f$. 
Given \f$f\f$, an \f$s=1\f$ model can be computed by taking 
\f[ [b] \ni f(c); \quad [A] \ni Df(X) \f]


\section matrix_inverse_derivative Derivative of Matrix Inverse

If \f$B=A^{-1}\f$, then 
\f[ \fbox{$ \dot{B} = B\dot{A}B $} \f]


*/
