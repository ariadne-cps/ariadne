/***************************************************************************
 *            watertank.cc
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
#include "real.h"
#include "expression.h"
#include "hybrid_system.h"
#include "hybrid_automaton.h"

using namespace Ariadne;
using std::cout; using std::endl;

int main()
{
    // Create the system object
    HybridSystem watertank;

    // Declare a discrete variable
    StringVariable valve("valve"); // Values "closed", "opening", "open", "closing"

    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0);
    RealConstant hmin("hmin",5.5);
    RealConstant hmax("hmax",8.0);
    RealConstant delta("delta",0.05);
    RealConstant lambda("lambda",0.02);
    RealConstant b("b",-0.3);

    // Declare the system variables
    RealVariable h("h");
    RealVariable alpha("alpha");

    // Declare the events we use
    Event start_opening("start_opening");
    Event start_closing("start_closing");
    Event finished_opening("finished_opening");
    Event finished_closing("finished_closing");

    // The water level is always given by the same dynamic
    watertank.new_dynamic(dot(h)=-lambda*h+b*alpha);

    // Specify the equation for how the valve opens/closes
    watertank.new_dynamic(valve=="opening", dot(alpha)=+1.0/T);
    watertank.new_dynamic(valve=="closing", dot(alpha)=-1.0/T);

    // When the valve is open or closed, alpha is constant.
    // Note that since we know alpha=0.0 or alpha=1.0, we should not need to consider alpha as a state variable.
    // This requires some cleverness on the part of the symbolic system analyser.
    // It would be possible to model the system with alpha explicitly a constant, but this would require some
    // cleverness on the part of the modeler.
    watertank.new_dynamic(valve=="closed" || valve=="open", dot(alpha)=0.0);

    // Specify the condition that the valve starts opening when hmax <= h <= hmax+delta
    // using an invariant and guard.
    watertank.new_invariant(h<=hmax+delta);
    watertank.new_guard(start_opening,valve=="closed" || valve=="closing", h>=hmax);

    // Specify the condition that the valve starts closing when hmin <= h <= hmin+delta
    // using a combined 'invariant and activation'. The event may occur when h<=hmin, and
    // must occur while h>=hmin-delta.
    watertank.new_guard(start_closing,valve=="open" || valve=="opening", h<=hmin,h>=hmin-delta);

    // Specify the guards for when the valve reaches the desired position
    watertank.new_guard(finished_opening, valve=="opening", alpha>=1.0, alpha<=1.0);
    watertank.new_guard(finished_closing, valve=="closing", alpha<=1.0, alpha>=1.0);

    // Explicitly disallow events when they don't make sense
    watertank.new_guard(finished_opening, !(valve=="opening"), false);
    watertank.new_guard(finished_closing, !(valve=="closing"), false);
    watertank.new_guard(start_closing, !(valve=="open" || valve=="opening"), false);
    watertank.new_guard(start_opening, !(valve=="closed" || valve=="closing"), false);

    // Specify the discrete resets
    watertank.new_transition(finished_opening, next(valve)="open");
    watertank.new_transition(finished_closing, next(valve)="closed");
    watertank.new_transition(start_opening, next(valve)="opening");
    watertank.new_transition(start_closing, next(valve)="closing");

    // For any event occurring in any location, the value of h and alpha are not updated.
    watertank.new_reset(next(h)=h);
    watertank.new_reset(next(alpha)=alpha);


    /// Finished building the automaton

    cout << "Watertank = " << std::boolalpha << watertank << endl << endl;


{
    MonolithicHybridAutomaton watertank;

    const bool urgent=true;
    const bool permissive=true;

    RealConstant o("1.0",1.0);
    RealConstant hmax("hmax",8.0);
    RealConstant hmin("hmin","8.0");
    RealVariable h("height");
    RealVariable a("alpha");

    DiscreteEvent start_closing("start_closing");
    DiscreteEvent start_opening("start_opening");
    DiscreteEvent finish_closing("finish_closing");
    DiscreteEvent finish_opening("finish_opening");
    DiscreteState open("open");
    DiscreteState closed("closed");
    DiscreteState opening("opening");
    DiscreteState closing("closing");

    RealExpression dh=6.0*alpha-h;
    watertank.new_mode(open,alpha=1.0,dot(h)=dh);
    watertank.new_mode(closed,alpha=0.0,dot(h)=dh);
    watertank.new_mode(opening,(dot(h)=dh,dot(alpha)=1.0));
    watertank.new_mode(closing,(dot(h)=dh,dot(alpha)=-1.0));

    watertank.new_transition(start_closing,open,closing,(next(alpha)=alpha,next(h)=h),h>=hmax,urgent);
    watertank.new_transition(start_closing,opening,closing,(next(alpha)=alpha,next(h)=h),h>=hmax,urgent);
    watertank.new_transition(start_opening,closed,opening,(next(alpha)=alpha,next(h)=h),h<=hmin,urgent);
    watertank.new_transition(start_opening,closing,opening,(next(alpha)=alpha,next(h)=h),h<=hmin,urgent);
    watertank.new_transition(finish_closing,closing,closed,(next(h)=h),alpha>=1.0,urgent);
    watertank.new_transition(finish_opening,opening,closed,(next(h)=h),alpha<=1.0,urgent);

    std::cout << "WatertankAutomaton = " << watertank << endl;

}

}