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
#include "ariadne.h"


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
    RealVariable aperture("aperture");
    TimeVariable time;
    // Create the tank object
    AtomicHybridAutomaton tank("tank");

    // Declare a trivial discrete mode.
    AtomicDiscreteLocation draining("draining");
    // The water level is always given by the same dynamic
    // The inflow is controlled by the valve aperture, the outflow depends on the
    // pressure, which is proportional to the water height.
    tank.new_mode(draining,(dot(height)=-lambda*height+rate*aperture));

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

    // Since aperture is a known constant when the valve is open or closed,
    // specify aperture by an algebraic equation.
    valve.new_mode(open,(aperture=+1.0));
    valve.new_mode(closed,(aperture=-1.0));
    // Specify the differential equation for how the valve opens/closes.
    valve.new_mode(opening,(dot(aperture)=+1/T));
    valve.new_mode(closing,(dot(aperture)=-1/T));

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    //valve.new_invariant(open,height<=hmax,start_closing);
    //valve.new_invariant(opening,height<=hmax,start_closing);
    //valve.new_invariant(opening,aperture<=1.0,finished_opening);
    //valve.new_invariant(closed,height>=hmin,start_opening);
    //valve.new_invariant(closing,height>=hmin,start_opening);
    //valve.new_invariant(closing,aperture>=0.0,finished_closing);

    valve.new_transition(closed,start_opening,opening,(next(aperture)=aperture),height<=hmin);
    valve.new_transition(closing,start_opening,opening,(next(aperture)=aperture),height<=hmin);
    valve.new_transition(open,start_closing,closing,(next(aperture)=aperture),height>=hmax);
    valve.new_transition(opening,start_closing,closing,(next(aperture)=aperture),height>=hmax);

    // Set the transitions for when the valve finished opening.
    // Since aperture is defined by an algebraic equation in the new mode,
    // it may not be specified in the reset.
    valve.new_transition(opening,finished_opening,open,aperture>=1);
    valve.new_transition(closing,finished_closing,closed,aperture<=0);

    CompositeHybridAutomaton watertank_system((tank,valve));
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
    DiscreteLocation initial_location=(tank|draining,valve|opening);
    HybridSet initial_set((tank|draining,valve|opening),(height==0,aperture==0));

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
    for (ListSet<HybridEnclosure>::const_iterator it = orbit.reach().begin(); it != orbit.reach().end(); it++)
    {
        std::cout<<"."<<std::flush;
        it->adjoin_outer_approximation_to(hgts,4);
    }
    std::cout << "done." << std::endl;

    plot("watertank_compositional-reach", height_aperture_axes, Colour(0.0,0.5,1.0), hgts);

}
