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
#include "hybrid_automaton.h"
#include "hybrid_evolver.h"
#include "hybrid_set.h"


using namespace Ariadne;
using std::cout; using std::endl;

int main()
{
    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0);
    RealConstant hmin("hmin",5.5);
    RealConstant hmax("hmax",8.0);
    RealConstant delta("delta",0.05);
    RealConstant lambda("lambda",0.02);
    RealConstant rate("rate",0.3);




    // Declare the system variables
    RealVariable height("height");
    RealVariable alpha("alpha");

    // Create the tank object
    AtomicHybridAutomaton tank("tank");

    // Declare a trivial discrete mode.
    AtomicDiscreteLocation trivial("");

    // The water level is always given by the same dynamic
    // The inflow is controlled by the valve alpha, the outflow depends on the
    // pressure, which is proportional to the water height.
    tank.new_mode(trivial,(dot(height)=-lambda*height+rate*alpha));



    // Describe the valve model

    // Declare the events we use
    DiscreteEvent start_opening("start_opening");
    DiscreteEvent start_closing("start_closing");
    DiscreteEvent finished_opening("finished_opening");
    DiscreteEvent finished_closing("finished_closing");

    // Declare the locations we use
    AtomicDiscreteLocation open("open");
    AtomicDiscreteLocation opening("opening");
    AtomicDiscreteLocation closed("closed");
    AtomicDiscreteLocation closing("closing");

    AtomicHybridAutomaton valve("valve");

    // Since alpha is a known constant when the valve is open or closed,
    // specify alpha by an algebraic equation.
    valve.new_mode(open,(alpha=+1.0));
    valve.new_mode(closed,(alpha=-1.0));
    // Specify the differential equation for how the valve opens/closes.
    valve.new_mode(opening,(dot(alpha)=+1.0/T));
    valve.new_mode(closing,(dot(alpha)=-1.0/T));

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    valve.new_invariant(open,start_closing,height<=hmax);
    valve.new_invariant(opening,start_closing,height<=hmax);
    valve.new_invariant(opening,finished_opening,alpha<=1.0);
    valve.new_invariant(closed,start_opening,height>=hmin);
    valve.new_invariant(closing,start_opening,height>=hmin);
    valve.new_invariant(closing,finished_closing,alpha>=0.0);

    valve.new_transition(closed,start_opening,opening,(next(alpha)=alpha),height<=hmin);
    valve.new_transition(closing,start_opening,opening,(next(alpha)=alpha),height<=hmin);
    valve.new_transition(open,start_closing,closing,(next(alpha)=alpha),height>=hmax);
    valve.new_transition(opening,start_closing,closing,(next(alpha)=alpha),height>=hmax);

    // Set the transitions for when the valve finished opening.
    // Since alpha is defined by an algebraic equation in the new mode,
    // it may not be specified in the reset.
    valve.new_transition(opening,finished_opening,open,alpha>=1.0);
    valve.new_transition(closing,finished_closing,closed,alpha<=0.0);


    CompositeHybridAutomaton watertank_system((tank,valve));
    std::cout << "watertank_system:\n" << watertank_system << "\n";




/*
    // Compute the system evolution

    // Create a HybridEvolver object
    HybridEvolver evolver;
    evolver.verbosity = 1;

    // Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = 0.25;
    evolver.parameters().maximum_step_size = 0.125;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution starting from location l2, x = 0.0, y = 1.0" << std::endl;

    DiscreteLocation initial_location=(trivial,opening);
    Box initial_box(2, 0.0,0.00, 0.0,0.00);
    HybridEnclosureType initial_enclosure(DiscreteLocation((trivial,opening)),initial_box);
    Box bounding_box(2, -0.1,9.1, -0.1,1.1);

    HybridTime evolution_time(80.0,5);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(watertank_system,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.final size="<<orbit.final().size()<<std::endl;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    //std::cout << "Plotting orbit... "<<std::flush;
    //plot("watertank_compositional-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit);
    std::cout << "done." << std::endl;
*/
}
