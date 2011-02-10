/***************************************************************************
 *            watertank-nonlinear-monolithic-hysteresis.h
 *
 *  Copyright  2011  Luca Geretti
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

#ifndef WATERTANK_NONLINEAR_MONOLITHIC_HYSTERESIS_H_
#define WATERTANK_NONLINEAR_MONOLITHIC_HYSTERESIS_H_

#include "ariadne.h"

namespace Ariadne {

HybridAutomaton getWatertankNonlinearMonolithicHysteresis()
{
    /// Create a HybridAutomaton object
    HybridAutomaton system("watertank-nl-mono-hy");

    /// Set the system parameters
	RealConstant a("a",0.065);
	RealConstant b("b",Interval(0.3,0.32863));
	RealConstant T("T",4.0);
	RealConstant hmin("hmin",5.5);
	RealConstant hmax("hmax",8.0);
	RealConstant Delta("Delta",0.1);
  
    // Accessible constants
    system.register_accessible_constant(Delta);
    system.register_accessible_constant(T);
    system.register_accessible_constant(b);
    system.register_accessible_constant(hmin);
    system.register_accessible_constant(hmax);

    /// Create four discrete states
    DiscreteState opened("opened");
    DiscreteState closed("closed");
    DiscreteState opening("opening");
    DiscreteState closing("closing");
  
    /// Create the discrete events
    DiscreteEvent b_opening("b_opening");
    DiscreteEvent e_opening("e_opening");
    DiscreteEvent b_closing("b_closing");
    DiscreteEvent e_closing("e_closing");
    
    // System variables
    RealVariable x("x");    // water level
    RealVariable y("y");    // valve aperture
    List<RealVariable> varlist;
    varlist.append(x);
    varlist.append(y);

    // Water level dynamics
    RealExpression x_opening_closing = -a*sqrt(x) + b*y;
    RealExpression x_opened = -a*sqrt(x) + b;
    RealExpression x_closed = -a*sqrt(x);
    
    // Valve Aperture dynamics
    RealExpression y_opening = 1.0/T;
    RealExpression y_closing = -1.0/T;
    RealExpression y_opened_closed = 0.0;
    
    // Dynamics at the different modes
    List<RealExpression> exprlist;
    exprlist.append(x_opened);
    exprlist.append(y_opened_closed);
    VectorFunction dyn_opened(exprlist, varlist);
    exprlist[0] = x_closed;
    VectorFunction dyn_closed(exprlist, varlist);
    exprlist[0] = x_opening_closing;
    exprlist[1] = y_opening;
    VectorFunction dyn_opening(exprlist, varlist);
    exprlist[1] = y_closing;
    VectorFunction dyn_closing(exprlist, varlist);    
      
    // Reset functions
    RealExpression idx = x;
    RealExpression zero = 0.0;
    RealExpression one = 1.0;
    exprlist[0] = idx;
    exprlist[1] = zero;
    VectorFunction reset_y_zero(exprlist, varlist);
    exprlist[1] = one;
    VectorFunction reset_y_one(exprlist, varlist);

    // Create the guards.
    // Guards are true when f(x) >= 0
    RealExpression x_leq_min = -x + hmin + Delta;       // x <= hmin + Delta
    ScalarFunction guard_b_opening(x_leq_min, varlist);
    RealExpression y_geq_one = y - 1.0;                 // y >= 1
    ScalarFunction guard_e_opening(y_geq_one, varlist);
    RealExpression x_geq_max = x - hmax + Delta;        // x >= hmax - Delta
    ScalarFunction guard_b_closing(x_geq_max, varlist);
    RealExpression y_leq_zero = -y;                     // y <= 0
    ScalarFunction guard_e_closing(y_leq_zero, varlist);

    // Create the invariants.
    // Invariants are true when f(x) = Ax + b < 0
    // forced transitions do not need an explicit invariant, 
    // we need only the invariants for location open and closed
    RealExpression x_leq_max = x - hmax - Delta;    // x <= hmax + Delta
    ScalarFunction inv_opened(x_leq_max, varlist);
    RealExpression x_geq_min = -x + hmin - Delta;   // x >= hmin - Delta
    ScalarFunction inv_closed(x_geq_min, varlist);
  
    /// Build the automaton
    system.new_mode(opened,dyn_opened);
    system.new_mode(closing,dyn_closing);
    system.new_mode(closed,dyn_closed);
    system.new_mode(opening,dyn_opening);

    system.new_invariant(opened,inv_opened);
    system.new_invariant(closed,inv_closed);

    system.new_unforced_transition(b_closing,opened,closing,reset_y_one,guard_b_closing);
    system.new_forced_transition(e_closing,closing,closed,reset_y_zero,guard_e_closing);
    system.new_unforced_transition(b_opening,closed,opening,reset_y_zero,guard_b_opening);
    system.new_forced_transition(e_opening,opening,opened,reset_y_one,guard_e_opening);

    return system;
}

}

#endif /* WATERTANK_NONLINEAR_MONOLITHIC_HYSTERESIS_H_ */
