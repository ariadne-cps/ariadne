/***************************************************************************
 *            integration_methods.h
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

\file integration_methods.h
\brief Documentation on integration methods



\page integration Integration methods

Integration in Ariadne is performed using <em>integration steps</em> and <em>reachability steps</em> on <em>basic sets</em>, typically cuboids or zonotopes.
From time-to-time, the basic set may be <em>regularised</em> or <em>subdivided</em> to prevent loss of accuracy.

Most integration methods used in Ariadne are based on Taylor expansion of solutions curves.

\f[ \begin{aligned}
     \frac{dy}{dt} &= f(y) \\[\jot]
 \frac{d^2y}{dt^2} &= Df(y)\,f(y) \\[\jot]
 \frac{d^3y}{dt^3} &= D^2f(y)\,(f(y))^2 + (Df(y))^2\,f(y) \\[\jot]
 \frac{d^4y}{dt^4} &= D^3f(y)\,(f(y))^3 + 4\,D^2f(y)\,Df(y)\,(f(y))^2 + (Df(y))^3\,f(y) \\[\jot]
 \frac{d^5y}{dt^5} &= D^4f(y)\,(f(y))^4 + 7\,D^3f(y)\,Df(y)\,(f(y))^3 + 4\,(D^2f(y))^2\,(f(y))^3 \\
                   &\ \qquad + 11\,D^2f(y)\,(Df(y))^2\,(f(y))^2 + (Df(y))^4\,f(y) \end{aligned} \f]

The first spacial derivative \f$dy/dx\f$ of the solution satisfies
\f[ \begin{aligned}
 \frac{d}{dt}     \frac{dy}{dx} &= \bigl( Df(y) \bigr) \, \frac{dy}{dx} \\[\jot]
 \frac{d^2}{dt^2} \frac{dy}{dx} &= \bigl( D^2f(y)\,f(y)  + (Df(y))^2 \bigr) \, \frac{dy}{dx} \\[\jot]
 \frac{d^3}{dt^3} \frac{dy}{dx} &= \bigl( D^3f(y)\,(f(y))^2 + 4\,D^2f(x)\,Df(x)\,f(x) + (Df(x))^3 \bigr) \, \frac{dy}{dx} \end{aligned} \f]

If we define functions \f$g_i:\mathbb{R}^n\rightarrow\mathbb{R}^{n\times n}\f$ by \f$g_0(y)=f(y)\f$ and \f$g_{i+1}(y)={\displaystyle Dg_i(y)\,f(y)}\f$,
then 
\f[ \frac{d^ny}{dt^n}   = g_n(y); \qquad \frac{d^n}{dt^n} \frac{dy}{dx} = Dg_{n}(y) \frac{dy}{dx} \f]
 


\section boundingenclosure Bounding enclosures

For most integration methods, it is necessary to find a rough <em>bounding enclosure</em> \a B for the time-h integration of an initial set X.
A sufficient condition for B to be an enclosure is
\f[ X+[0,h]f(B) \subset B . \f]
This condition can always be satisfied for a continuous vector field by choosing the step size h sufficiently small.
If \f$B\f$ is an enclosure for \f$\Phi([0,h],X)\f$, then \f$X+[0,h]f(B_k)\f$ is a tighter enclosure.

To compute approximations to the spacial derivative \f$\frac{dy}{dx} = D\Phi(t,x)\f$, it is necessary to find a rough enclosure for \f$ \frac{dy}{dx}\f$.

Suppose \f$V(t)\f$ is a matrix-valued function satisfying \f$ \dot{V}(t) = A(t) V(t) \f$.
Then we have
\f[ ||V(t)-V(0)|| \leq ||V(0)|| \bigl( e^{L(t)} -1  \bigr) \text{ where } L(t) = \int_0^t ||A(\tau)|| \,d\tau \f]
Additionally, 
\f[ ||V(t)|| \leq e^{l(t)} \, ||V(0)|| \text{ where } l(t) = \int_0^t \lambda(A(\tau))\,d\tau \text{ and } \lambda(A) = \lim_{h\to0} \frac{||I+Ah||-1}{h} \f]
The quantity \f$\lambda(A)\f$ is the <em>logarithmic derivative</em>. In %Ariadne, we use the supremum norm, yielding
\f[ ||v||_\infty = \max_{i}|v_i|; \quad ||A||_\infty = \max_i\sum_j |a_{ij}|; \quad \lambda_\infty(A) = \max_{i} \bigl( a_{ii} + \sum_{j\neq i} |a_{ij} | \bigr)   . \f]

Hence if \f$B\f$ is a bounding enclosure for \f$\Phi([0,h],x)\f$, then 
\f[ ||D\Phi(t,x)-I||_\infty \leq e^{Lt} - 1 \quad  \text{ where } \quad L = ||Df(B)||_\infty . \f]


\section affineintegrator Affine Integrator

An affine vector field \f$ \dot{x}=Ax+b \f$ can be integrated directly using the formula
\f[ \begin{aligned} x(t) &=e^{At}x_0 + A^{-1}(e^{At}-I)b \\ 
                         &= x_0+ t \sum_{n=0}^{\infty} \frac{A^nt^n}{(n+1)!} \bigl(Ax_0+b\bigr) \end{aligned} \f]
The spacial derivative \f$Dx(t)\f$ is given by
\f[  Dx(t) = e^{At} \f]
When applied to a zonotope \f$c_0+G_0e\f$, the the centre and generators of the integrated zonotope become
\f[ c_1 = e^{Ah}c_0 + A^{-1}(e^{Ah}-I)b; \qquad  G_1 = e^{Ah} G_0, \f] 
which are easily computed using interval arithmetic.                         



\section eulerintegrator Euler method

The classical Euler method is 
\f$ x_1 = x_0 + hf(x_0) \approx x(h) \f$
with one-step error bound
\f$ ||x_1-x(h)|| = h ||f(x_0)-f(\xi)|| = \frac{h^2}{2} || Df(\xi)f(\xi) || . \f$
The method can be implemented directly using interval arithmetic by
\f[ R_1 = R_0 + hf(B) . \f]
where the initial set \f$R_0\f$ and the final set \f$R_1\f$ are cuboids.
The accuracy is usually too poor for this method to be useful.



\internal I'm not sure that this is the "standard" Lohner integrator.

\section lohnerintegrator Lohner Integrator

The classical Lohner integrator is a first-order in space integration scheme operating on cuboids.
Let \f[ P_n(x,h) = \sum_{j=0}^{n} \frac{h^j}{j!} f^{(j)}(x), \quad  R_n(B,h)= \frac{h^{n+1}}{(n+1)!} f^{(n+1)}(B), \qquad
\text{where} \quad f^{(n)} = Df^{(n-1)}\cdot f, \quad f^{(0)}(x) = x . \f]
If \f$X_0\f$ is the initial set with centre \f$x_0\f$, and \f$B\f$ is an enclosure, then the Lohner scheme is
\f[ X_1 = P_n(x_0,h) + DP_n(x_0,h)\cdot(X_0-x_0) + R_n(B,h) \f]
Taking \f$n=1\f$, we obtain
\f[ X_1 = x_0 + h f(x_0) + (I+hDf(x_0)) \cdot (X_0-x_0) + \frac{h^2}{2} Df(B) f(B) \f]
Taking \f$n=2\f$, we obtain
\f[ X_1 = x_0 + h f(x_0) + \frac{h^2}{2} Df(x_0)f(x_0) \;+\; \Bigl(I+hDf(x_0)+\frac{h^2}{2} \bigl( D^2f(x_0) f(x_0) + (Df(x_0))^2 \bigr) \Bigr) \cdot (X_0-x_0) \;+\; \frac{h^3}{3} \bigl( D^2f(B) f(B) + (Df(B)^2) \bigr) f(B) \f]

This scheme cannot be used directly to compute the spacial variation derivative \f$V(t,x)=D\Phi(t,x)\f$.
To compute the spacial derivative \f$D\Phi(t,x)\f$, we need a bound \f$B\f$ for \f$\Phi([0,t],X)\f$ and a bound \f$W\f$ for \f$V([0,t],X)\f$.
Once this has been obtained, we can use the first-order update
\f[   V_1 = \bigl( I + t Df(B) W \bigr) \, V_0 \f]
For a second-order (in time) approximation, we can use
\f[   V_1 \in \bigl( I + t Df(X_0) + (D^2f(B)f(B)+ Df(B)^2) W \bigr) \, V_0 \f]

In Ariadne, we use a \f$C^1\f$ in space algorithm similar to that of Zgliczynski.

If the initial set is a zonotope with centre \f$c_0\f$ and generators \f$G_0\f$, 
then the new centre and generators can be computed by
\f[ c_1 = c_0 + h f(c_0) + \frac{h^2}{2} Df(B) f(B) ;  \qquad  G_1 = \bigl(I  + h Df(c_0) \bigr)\, G_0. \f]
or as a variant, taking \f$B_c\f$ as a rough bounding enclosure for \f$c\f$, by
\f[ c_1 = c_0 + h f(c_0) + \frac{h^2}{2} Df(B_c) f(B_c) ;  \qquad  G_1 = \bigl(I+hDf(B)\bigr) \, W \, G_0. \f]
This scheme may be easily implemented using interval arithmetic.

For a reachability step over \f$t\in[h_1,h_2]\f$, the centre is displaced to \f$c_1\f$ and a new generator \f$g^*\f$ is introduced.
\f[ c_1 = c_0 + \frac{h_1+h_2}{2} f(B) ; \qquad  g^* = \frac{h_2-h_1}{2} f(B) . \f]
For a more accurate reachability step over \f$[0,h]\f$, we first perform an integration step of size \f$h/2\f$, and then a reachability step over \f$[-h/2,h/2]\f$.

To reduce the wrapping effect in the generator matrix, the interval matrix \f$G\f$ may be replaced by its average \f$\widetilde{G}\f$, and the error absorbed into the centre \f$\tilde{c} = c+(G-\tilde{G})e\f$, where \f$e\f$ is the interval vector with elements \f$[-1,1]\f$.

The Lohner integration step is usually used in conjunction with an orthogonal regularisation scheme.
The generator matrix is factorised \f$ G = QR \f$ where \f$Q\f$ is orthonormal and \f$R\f$ is upper-triangular, 
and the new generators \f$G'\f$ are given by \f$ G' = QD \f$, where \f$D\f$ is diagonal with \f$D_{ii}=\sum_{j} |R_{ij}|\f$.


\section taylorintegrator Taylor methods

Methods based on higher-order Taylor expansions can be devised. A second-order in time and first-order in space integrator for zonotopes is given by
\f[ \begin{aligned} 
      c_1 &= c_0 \;+\; h\,f(c_0) \;+\; \frac{h^2}{2}\,Df(c_0)\,f(c_0) \;+\; \frac{h^3}{6}\,\bigl( D^2f(B_c)\,f(B_c) + (Df(B_c))^2 \bigr)\,f(B_c); \\[\jot]
      G_1 &= G_0 \;+\; h\,Df(c_0) \;+\; \frac{h^2}{2}\,\bigl( D^2f(B)\,f(B) + (Df(B))^2 \bigr) 
    \end{aligned} \f]
For higher order in space, Taylor sets must be used.


*/
