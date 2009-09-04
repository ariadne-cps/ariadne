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
    double a = -0.02;
    double b = 0.3;
    double T = 1.5;
    double hmin = 5.55;
    double Delta = 0.05;
    double hmax = 5.70;
    double tmax = 20.0;
    double dmax = 1.0;
    
    double A1[9]={a,b,0,
                  0,0,0,
                  0,0,0};
    double b1[3]={0,1.0/T,1.0};

    double A2[9]={a,  0.0,0.0,
                  0.0,0.0,0.0,
                  0.0,0.0,0.0};
    double b2[3]={b,0.0,1.0};

    double A3[9]={a,b,0,
                  0,0,0,
                  0,0,0};
    double b3[3]={0,-1.0/T,1.0};

    double A4[9]={a,  0.0,0.0,
                  0.0,0.0,0.0,
                  0.0,0.0,0.0};
    double b4[3]={0.0,0.0,1.0};

    /// Build the Hybrid System
  
    /// Create a HybridAutomton object
    HybridAutomaton watertank_system;
  
    /// Create four discrete states
    DiscreteState l1(1);
    DiscreteState l2(2);
    DiscreteState l3(3);
    DiscreteState l4(4);
  
    /// Create the discrete events
    DiscreteEvent e12(12);
    DiscreteEvent e13(13);
    DiscreteEvent e23(23);
    DiscreteEvent e34(34);
    DiscreteEvent e31(31);
    DiscreteEvent e41(41);
  
    /// Create the dynamics
    AffineFunction dynamic1(Matrix<Float>(3,3,A1),Vector<Float>(3,b1));
    AffineFunction dynamic2(Matrix<Float>(3,3,A2),Vector<Float>(3,b2));
    AffineFunction dynamic3(Matrix<Float>(3,3,A3),Vector<Float>(3,b3));
    AffineFunction dynamic4(Matrix<Float>(3,3,A4),Vector<Float>(3,b4));
    
    cout << "dynamic1 = " << dynamic1 << endl << endl;
    cout << "dynamic2 = " << dynamic2 << endl << endl;
    cout << "dynamic3 = " << dynamic3 << endl << endl;
    cout << "dynamic4 = " << dynamic4 << endl << endl;

    /// Create the resets
    IdentityFunction reset_id(3);
    AffineFunction reset_y_zero(Matrix<Float>(3,3,
                        1.0,0.0,0.0,
                        0.0,0.0,0.0,
                        0.0,0.0,1.0),
                        Vector<Float>(3,0.0,0.0,0.0));
    cout << "reset_y_zero=" << reset_y_zero << endl << endl;
    AffineFunction reset_y_one(Matrix<Float>(3,3,
                        1.0,0.0,0.0,
                        0.0,0.0,0.0,
                        0.0,0.0,1.0),
                        Vector<Float>(3,0.0,1.0,0.0));
    cout << "reset_y_one=" << reset_y_one << endl << endl;

    /// Create the guards.
    /// Guards are true when f(x) = Ax + b > 0
    AffineFunction guard12(Matrix<Float>(1,3,0.0,1.0,0.0),Vector<Float>(1,-1.0));
    cout << "guard12=" << guard12 << endl << endl;
    AffineFunction guard23(Matrix<Float>(1,3,1.0,0.0,0.0),Vector<Float>(1, - hmax + Delta));
    cout << "guard23=" << guard23 << endl << endl;
    AffineFunction guard34(Matrix<Float>(1,3,0.0,-1.0,0.0),Vector<Float>(1,0.0));
    cout << "guard34=" << guard34 << endl << endl;
    AffineFunction guard41(Matrix<Float>(1,3,-1.0,0.0,0.0),Vector<Float>(1,hmin + Delta));
    cout << "guard41=" << guard41 << endl << endl;

    /// Create the invariants.
    /// Invariants are true when f(x) = Ax + b < 0
    /// forced transitions do not need an explicit invariant
    /// x < hmax + Delta
    AffineFunction inv2(Matrix<Float>(1,3,1.0,0.0,0.0),Vector<Float>(1, - hmax - Delta));//
    cout << "inv2=" << inv2 << endl << endl;
    /// x > hmin - Delta
    AffineFunction inv4(Matrix<Float>(1,3,-1.0,0.0,0.0),Vector<Float>(1, hmin - Delta));
    cout << "inv4=" << inv4 << endl << endl;
  
    /// Build the automaton
    watertank_system.new_mode(l1,dynamic1);
    watertank_system.new_mode(l2,dynamic2);
    watertank_system.new_mode(l3,dynamic3);
    watertank_system.new_mode(l4,dynamic4);

    watertank_system.new_invariant(l1,inv2);
    watertank_system.new_invariant(l2,inv2);
    watertank_system.new_invariant(l3,inv4);
    watertank_system.new_invariant(l4,inv4);

    watertank_system.new_forced_transition(e12,l1,l2,reset_y_one,guard12);
    watertank_system.new_unforced_transition(e13,l1,l3,reset_id,guard23);    
    watertank_system.new_unforced_transition(e23,l2,l3,reset_y_one,guard23);
    watertank_system.new_forced_transition(e34,l3,l4,reset_y_zero,guard34);
    watertank_system.new_unforced_transition(e31,l3,l1,reset_id,guard41);
    watertank_system.new_unforced_transition(e41,l4,l1,reset_y_zero,guard41);


    /// Finished building the automaton

    cout << "Automaton = " << watertank_system << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = 0.1;
    evolver.parameters().maximum_step_size = 0.1;
    evolver.verbosity=1;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution starting from location l2, x = 5.0, y = 1.0, t = 0.0" << std::endl;

    Box initial_box(3, 5.0,5.001, 1.0,1.001, 0.0,0.001);
    HybridEnclosureType initial_enclosure(l2,initial_box);
    Box bounding_box(3, -0.1,9.1, -0.1,1.1, -0.1,tmax+0.1);
  
    HybridTime evolution_time(tmax,dmax);
  
    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(watertank_system,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit="<<orbit<<std::endl;
    Figure g;
    Box graphic_box(2, -0.1,tmax+0.1, 4.1,6.1);
    g.set_bounding_box(graphic_box);
    array<uint> p(2,2,0);
    g.set_projection_map(ProjectionFunction(p,3));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << orbit;
    g.write("watertank-time-orbit");
/*
    std::cout << "Computing reach set using HybridEvolver... " << std::flush;
    EnclosureListType reach = evolver.reach(watertank_system,initial_enclosure,evolution_time);
    std::cout << "done." << std::endl;

    std::cout << "Orbit="<<reach<<std::endl;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    plot("watertank-reach-evolver",bounding_box, Colour(0.0,0.5,1.0), reach);


    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time = 10.0;
    analyser.parameters().maximum_grid_depth = 10;
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[l2]=initial_box;

    HybridTime reach_time(tmax,dmax);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)
    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet lower_reach_set = analyser.lower_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;

    Figure g;
    Box graphic_box(2, -0.1,tmax+0.1, 4.1,6.1);
    g.set_bounding_box(graphic_box);
    array<uint> p(2,2,0);
    g.set_projection_map(ProjectionFunction(2,3,p));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << lower_reach_set;
    g.write("watertank-time-lower");
*/

}
