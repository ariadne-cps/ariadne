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



*/
