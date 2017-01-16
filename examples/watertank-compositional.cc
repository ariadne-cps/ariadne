/***************************************************************************
 *            watertank-compositional.cc
 *
 *  Copyright  2008-16  Davide Bresolin, Pieter Collins
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
#include "ariadne.h"


using namespace Ariadne;
using std::cout; using std::endl;

Int main()
{
    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0_decimal);
    RealConstant hmin("hmin",5.5_decimal);
    RealConstant hmax("hmax",8.0_decimal);
    RealConstant delta("delta",0.05_decimal);
    RealConstant lambda("lambda",0.02_decimal);
    RealConstant rate("rate",0.3_decimal);

    // Declare the shared system variables
    TimeVariable time;
    RealVariable aperture("aperture");

    // Create the tank object
    HybridAutomaton tank_automaton("tank_automaton");
    RealVariable height("height");

    // Declare a trivial discrete mode.
    DiscreteLocation draining;
    // The water level is always given by the same dynamic
    // The inflow is controlled by the valve aperture, the outflow depends on the
    // pressure, which is proportional to the water height.
    tank_automaton.new_mode(draining,{dot(height)=-lambda*height+rate*aperture});

    // Describe the valve model

    // Declare the events we use
    DiscreteEvent start_opening("start_opening");
    DiscreteEvent start_closing("start_closing");
    DiscreteEvent finished_opening("finished_opening");
    DiscreteEvent finished_closing("finished_closing");


    HybridAutomaton valve_automaton("valve_automaton");
    StringVariable valve("valve");
    // Declare the values the valve can variable can have
    StringConstant open("open");
    StringConstant opening("opening");
    StringConstant closed("closed");
    StringConstant closing("closing");

    // Since aperture is a known constant when the valve is open or closed,
    // specify aperture by an algebraic equation.
    valve_automaton.new_mode(valve|open,{aperture=+1.0_decimal});
    valve_automaton.new_mode(valve|closed,{aperture=-1.0_decimal});
    // Specify the differential equation for how the valve opens/closes.
    valve_automaton.new_mode(valve|opening,{dot(aperture)=+1/T});
    valve_automaton.new_mode(valve|closing,{dot(aperture)=-1/T});

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    //valve.new_invariant(open,height<=hmax,start_closing);
    //valve.new_invariant(opening,height<=hmax,start_closing);
    //valve.new_invariant(opening,aperture<=1.0,finished_opening);
    //valve.new_invariant(closed,height>=hmin,start_opening);
    //valve.new_invariant(closing,height>=hmin,start_opening);
    //valve.new_invariant(closing,aperture>=0.0,finished_closing);

    valve_automaton.new_transition(valve|closed,start_opening,valve|opening,{next(aperture)=aperture},height<=hmin);
    valve_automaton.new_transition(valve|closing,start_opening,valve|opening,{next(aperture)=aperture},height<=hmin);
    valve_automaton.new_transition(valve|open,start_closing,valve|closing,{next(aperture)=aperture},height>=hmax);
    valve_automaton.new_transition(valve|opening,start_closing,valve|closing,{next(aperture)=aperture},height>=hmax);

    // Set the transitions for when the valve finished opening.
    // Since aperture is defined by an algebraic equation in the new mode,
    // it may not be specified in the reset.
    valve_automaton.new_transition(valve|opening,finished_opening,valve|open,aperture>=1);
    valve_automaton.new_transition(valve|closing,finished_closing,valve|closed,aperture<=0);

    CompositeHybridAutomaton watertank_system({tank_automaton,valve_automaton});
    std::cout << "watertank_system:\n" << watertank_system << "\n";

    // Compute the system evolution

    // Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver(watertank_system);
    evolver.verbosity = 0;

    // Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(0.25);
    evolver.configuration().set_maximum_step_size(0.125);
    std::cout << evolver.configuration() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolver::EnclosureType HybridEnclosureType;
    typedef GeneralHybridEvolver::OrbitType OrbitType;
    typedef GeneralHybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution starting from location l2, x = 0.0, y = 0.0" << std::endl;
    DiscreteLocation initial_location={valve|opening};
    HybridSet initial_set({valve|opening},{height==0,aperture==0});

    HybridTime evolution_time(80.0,5);

    std::cout << "Computing orbit... " << std::flush;

    OrbitType orbit = evolver.orbit(initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << orbit << std::endl;

    std::cout << "Orbit.final size="<<orbit.final().size()<<std::endl;

    // FIXME: Aperture variable not a state variable in all modes, so cannot be plotted
    // Axes2d axes(-0.1<=height<=9.1, -0.1<=aperture<=1.1);
    Axes2d axes(0<=time<=80,-0.1<=height<=9.1);
    std::cout << "Plotting orbit... "<<std::flush;
    plot("watertank_compositional-orbit",axes, Colour(0.0,0.5,1.0), orbit);
    std::cout << "done." << std::endl;

    Axes2d height_aperture_axes(-0.1,height,9.1, -0.1,aperture,1.3);

    std::cout << "Discretising orbit" << std::flush;
    HybridGrid grid(watertank_system.state_space());
    HybridGridTreeSet hgts(grid);

    for (ListSet<HybridEnclosure>::ConstIterator it = orbit.reach().begin(); it != orbit.reach().end(); it++)
    {
        std::cout<<"."<<std::flush;
        it->adjoin_outer_approximation_to(hgts,4);
    }
    std::cout << "done." << std::endl;

    // The following currently fails since auxiliary variables are not tracked
    plot("watertank_compositional-reach", height_aperture_axes, Colour(0.0,0.5,1.0), hgts);

}
