/***************************************************************************
 *            hybrid_evolution_methods.h
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

\file hybrid_evolution_methods.h
\brief Documentation on methods for evolution of hybrid systems



\page hybrid_evolution Hybrid Evolution Methods


\section force_transitions Forced transitions

Let \f$\tau(x)\f$ be the time needed to flow from \f$x\f$ to the guard set \f$g(x)=0\f$.
Then the transition is given by 
\f[ \Psi(x,t) = \Phi_2(R(\Phi_1(x,\tau(x))),t-\tau(x)) . \f]
and the Jacobian derivative is
\f[ D\Psi(x,t) \in D\Phi_2(B_2)\,DR(B_1)\,D\Phi_1(B_1) \, + \, \bigl( D\Phi_2(B_2) \, DR(B_1) \, F_1(B_1) - F_2(B_2) \bigr) \nabla\tau(x) . \f]
where
\f[ - \nabla\tau(x) = \frac{\nabla g(B)\cdot D\Phi_1(B)}{\nabla g(B)\cdot F_1(B)} . \f]
Suppose \f$\tau(c)\f$ is known. Then
\f[ \Psi(x,t) = \Psi(c,t) + D\Phi_2(B_2) DR(B_1) D\Phi_1(B_1) + (F_1(B_1) - F_2(B_2)) \nabla\tau(x) R(\Phi_1(x,\tau(x))),t-\tau(x)) . \f]


*/
