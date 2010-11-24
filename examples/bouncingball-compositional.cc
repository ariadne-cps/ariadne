/***************************************************************************
 *      bouncingball.cc
 *
 *  Copyright  2008-9  Davide Bresolin, Pieter Collins
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

#include <cstdarg>
#include "numeric.h"
#include "expression.h"
#include "hybrid_automaton.h"

using namespace Ariadne;
using std::cout;

int main()
{
    /// Set the system parameters
    Constant<Real> a("a",0.5);
    Constant<Real> g("g",9.8);

    AtomicHybridAutomaton bouncingball("Bouncing Ball");
    AtomicDiscreteLocation falling("-");
    DiscreteEvent bounce("bounce");

    RealVariable x("x");
    RealVariable v("v");

    bouncingball.new_mode(falling, (dot(x)=v, dot(v)=-g) );
    bouncingball.new_invariant(falling, bounce, x>=0);
    bouncingball.new_transition(falling, bounce, x<=0 && v<0, falling, (next(x)=x, next(v)=-a*v) );

    cout << "BouncingBall = \n" << std::boolalpha << bouncingball;
}
