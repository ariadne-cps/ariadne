/***************************************************************************
 *            watertank.cc
 *
 *  Copyright  2008  Davide Bresolin
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


int main() 
{
  
    /// Set the system parameters
    double a = 0.02;
    double b = 0.3;
    double T = 4.0;
    double hmin = 5.5;
    double Delta = 0.05;
    double hmax = 8.0;

    /// Create a HybridAutomton object
    HybridAutomaton watertank_system("Watertank");
  
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
    RealVariable x("x");            // water level
    RealVariable y("y");    // valve aperture
    List<RealVariable> varlist;
    varlist.append(x);
    varlist.append(y);    
    std::cout << "Variables: " << varlist << std::endl << std::endl;
    
    // Water level dynamics
    RealExpression x_opening_closing = -a*x + b*y;
    RealExpression x_opened = -a*x + b;
    RealExpression x_closed = -a*x;
    
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
    std::cout << "Dynamics" << std::endl;
    std::cout << "   Opened: " << dyn_opened << std::endl;
    std::cout << "   Closed: " << dyn_closed << std::endl;
    std::cout << "  Opening: " << dyn_opening << std::endl;
    std::cout << "  Closing: " << dyn_closing << std::endl << std::endl;
      
      
    // Reset functions
    RealExpression idx = x;
    RealExpression zero = 0.0;
    RealExpression one = 1.0;
    exprlist[0] = idx;
    exprlist[1] = zero;
    VectorFunction reset_y_zero(exprlist, varlist);
    exprlist[1] = one;
    VectorFunction reset_y_one(exprlist, varlist);
    cout << "Resets" << endl;
    cout << "  reset_y_zero:" << reset_y_zero << endl;
    cout << "   reset_y_one:" << reset_y_one << endl << endl;

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
    cout << "Guards" << endl;
    cout << "  x <= hmin + Delta: " << guard_b_opening << endl;
    cout << "             y >= 1: " << guard_e_opening << endl;
    cout << "  x >= hmax - Delta: " << guard_b_closing << endl;
    cout << "             y <= 0: " << guard_e_closing << endl << endl;

    // Create the invariants.
    // Invariants are true when f(x) = Ax + b < 0
    // forced transitions do not need an explicit invariant, 
    // we need only the invariants for location open and closed
    RealExpression x_leq_max = x - hmax - Delta;    // x <= hmax + Delta
    ScalarFunction inv_opened(x_leq_max, varlist);
    RealExpression x_geq_min = -x + hmin - Delta;   // x >= hmin - Delta
    ScalarFunction inv_closed(x_geq_min, varlist);
    cout << "Invariants" << endl;
    cout << "  Opened: " << inv_opened << endl;
    cout << "  Closed: " << inv_closed << endl << endl;
  
    /// Build the automaton
    watertank_system.new_mode(opened,dyn_opened);
    watertank_system.new_mode(closing,dyn_closing);
    watertank_system.new_mode(closed,dyn_closed);
    watertank_system.new_mode(opening,dyn_opening);

    watertank_system.new_invariant(opened,inv_opened);
    watertank_system.new_invariant(closed,inv_closed);

    watertank_system.new_unforced_transition(b_closing,opened,closing,reset_y_one,guard_b_closing);
    watertank_system.new_forced_transition(e_closing,closing,closed,reset_y_zero,guard_e_closing);
    watertank_system.new_unforced_transition(b_opening,closed,opening,reset_y_zero,guard_b_opening);
    watertank_system.new_forced_transition(e_opening,opening,opened,reset_y_one,guard_e_opening);


    /// Finished building the automaton

    cout << "Automaton = " << watertank_system << endl << endl;
       
    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver;
    evolver.verbosity = 1;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_cell = Vector<Float>(2,0.25);
    evolver.parameters().maximum_step_size = 0.125;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution starting from location closed, x = 7.0, y = 1.0" << std::endl;

    Box initial_box(2, 7.0,7.00, 0.0,0.00);
    HybridEnclosureType initial_enclosure(closed,initial_box);
    Box bounding_box(2, -0.1,9.1, -0.1,1.1);
  
    HybridTime evolution_time(80.0,5);
  
    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(watertank_system,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.final size="<<orbit.final().size()<<std::endl;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    std::cout << "Plotting orbit... "<<std::flush;
    plot("watertank-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit);
    std::cout << "done." << std::endl;

    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    evolver.verbosity=0;
    analyser.parameters().lock_to_grid_time = 30.0;
	analyser.parameters().maximum_grid_depth = 0;
    watertank_system.set_grid(Grid(Vector<Float>(2),Vector<Float>(2,0.2)));
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[closed]=initial_box;

    HybridTime reach_time(64.0,2);

    plot("watertank-initial_set1",bounding_box, Colour(0.0,0.5,1.0), initial_set);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)
    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet lower_reach_set_ptr = analyser.lower_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-lower_reach1",bounding_box, Colour(0.0,0.5,1.0), lower_reach_set_ptr);

    // Compute evolved sets and reach sets using upper semantics.
    // These functions compute over-approximations to the evolved and reachabe sets. Subdivision is used
    // as necessary to keep the local errors reasonable. The accumulated global error may be very large.
    std::cout << "Computing upper reach set... " << std::flush;
    HybridGridTreeSet upper_reach_set_ptr = analyser.upper_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-upper_reach1",bounding_box, Colour(0.0,0.5,1.0), upper_reach_set_ptr);

/*
   
	std::cout << "Computing evolution starting from location l1, x = 0.0, y = 0.0" << std::endl;

    Box initial_box2(2, 0.0,0.001, 0.0,0.001);
    HybridImageSet initial_set2;
    initial_set2[l1]=initial_box2;

    plot("watertank-initial_set2",bounding_box, Colour(0.0,0.5,1.0), initial_set2);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)
    std::cout << "Computing lower reach set... " << std::flush;
    lower_reach_set_ptr = analyser.lower_reach(watertank_system,initial_set2,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-lower_reach2",bounding_box, Colour(0.0,0.5,1.0), *lower_reach_set_ptr);

    // Compute evolved sets and reach sets using upper semantics.
    // These functions compute over-approximations to the evolved and reachabe sets. Subdivision is used
    // as necessary to keep the local errors reasonable. The accumulated global error may be very large.
    std::cout << "Computing upper reach set... " << std::flush;
    upper_reach_set_ptr = analyser.upper_reach(watertank_system,initial_set2,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-upper_reach2",bounding_box, Colour(0.0,0.5,1.0), *upper_reach_set_ptr);
*/

}
