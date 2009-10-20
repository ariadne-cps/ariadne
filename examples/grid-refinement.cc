/***************************************************************************
 *            grid-refinement.cc
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
    double T = 4.0;
    double hmin = 5.5;
    double Delta = 0.05;
    double hmax = 8.0;
    
    // Initial grid depth and maximum grid depth
    int min_grid_depth = 2;
    int max_grid_depth = 7;
    
    // System is verified safe at the 4th iteration
    double minsafe=5.25;
    double maxsafe=8.25;

    // System is verified unsafe at the 2nd iteration
    //double maxsafe=8.1;

    // System's safety cannot be verified
    //double maxsafe=8.2;

    
    double A1[4]={a,b,0,0};
    double b1[2]={0,1.0/T};

    double A2[4]={a,0.0,0.0,0.0};
    double b2[2]={b,0.0};

    double A3[4]={a,b,0,0};
    double b3[2]={0,-1.0/T};

    double A4[4]={a,0.0,0.0,0.0};
    double b4[2]={0.0,0.0};

    /// Build the Hybrid System
  
    /// Create a HybridAutomton object
    MonolithicHybridAutomaton watertank_system;
  
    /// Create four discrete states
    AtomicDiscreteLocation l1(1);
    AtomicDiscreteLocation l2(2);
    AtomicDiscreteLocation l3(3);
    AtomicDiscreteLocation l4(4);
  
    /// Create the discrete events
    DiscreteEvent e12(12);
    DiscreteEvent e23(23);
    DiscreteEvent e34(34);
    DiscreteEvent e41(41);
  
    /// Create the dynamics
    VectorAffineFunction dynamic1(Matrix<Float>(2,2,A1),Vector<Float>(2,b1));
    VectorAffineFunction dynamic2(Matrix<Float>(2,2,A2),Vector<Float>(2,b2));
    VectorAffineFunction dynamic3(Matrix<Float>(2,2,A3),Vector<Float>(2,b3));
    VectorAffineFunction dynamic4(Matrix<Float>(2,2,A4),Vector<Float>(2,b4));
    
    cout << "dynamic1 = " << dynamic1 << endl << endl;
    cout << "dynamic2 = " << dynamic2 << endl << endl;
    cout << "dynamic3 = " << dynamic3 << endl << endl;
    cout << "dynamic4 = " << dynamic4 << endl << endl;

    /// Create the resets
    VectorAffineFunction reset_y_zero(Matrix<Float>(2,2,1.0,0.0,0.0,0.0),Vector<Float>(2,0.0,0.0));
    cout << "reset_y_zero=" << reset_y_zero << endl << endl;
    VectorAffineFunction reset_y_one(Matrix<Float>(2,2,1.0,0.0,0.0,0.0),Vector<Float>(2,0.0,1.0));
    cout << "reset_y_one=" << reset_y_one << endl << endl;

    /// Create the guards.
    /// Guards are true when f(x) = Ax + b > 0
    VectorAffineFunction guard12(Matrix<Float>(1,2,0.0,1.0),Vector<Float>(1,-1.0));
    cout << "guard12=" << guard12 << endl << endl;
    VectorAffineFunction guard23(Matrix<Float>(1,2,1.0,0.0),Vector<Float>(1, - hmax + Delta));
    cout << "guard23=" << guard23 << endl << endl;
    VectorAffineFunction guard34(Matrix<Float>(1,2,0.0,-1.0),Vector<Float>(1,0.0));
    cout << "guard34=" << guard34 << endl << endl;
    VectorAffineFunction guard41(Matrix<Float>(1,2,-1.0,0.0),Vector<Float>(1,hmin + Delta));
    cout << "guard41=" << guard41 << endl << endl;

    /// Create the invariants.
    /// Invariants are true when f(x) = Ax + b < 0
    /// forced transitions do not need an explicit invariant, 
    /// we need only the invariants for location 2 and 4
    VectorAffineFunction inv2(Matrix<Float>(1,2,1.0,0.0),Vector<Float>(1, - hmax - Delta));//
    cout << "inv2=" << inv2 << endl << endl;
    VectorAffineFunction inv4(Matrix<Float>(1,2,-1.0,0.0),Vector<Float>(1, hmin - Delta));
    cout << "inv4=" << inv4 << endl << endl;
  
    /// Build the automaton
    watertank_system.new_mode(l1,dynamic1);
    watertank_system.new_mode(l2,dynamic2);
    watertank_system.new_mode(l3,dynamic3);
    watertank_system.new_mode(l4,dynamic4);

    watertank_system.new_invariant(l2,inv2);
    watertank_system.new_invariant(l4,inv4);

    watertank_system.new_forced_transition(e12,l1,l2,reset_y_one,guard12);
    watertank_system.new_unforced_transition(e23,l2,l3,reset_y_one,guard23);
    watertank_system.new_forced_transition(e34,l3,l4,reset_y_zero,guard34);
    watertank_system.new_unforced_transition(e41,l4,l1,reset_y_zero,guard41);


    /// Finished building the automaton

    cout << "Automaton = " << watertank_system << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver;
    evolver.verbosity = 0;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = 0.25;
    evolver.parameters().maximum_step_size = 0.125;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    Box initial_continuous_box(2, 6.0,6.00, 1.0,1.00);
    HybridImageSet initial_set;
    initial_set[l2]=initial_continuous_box;
    Box bounding_box(2, -0.1,9.1, -0.1,1.1);
  
    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time = 32.0;
    analyser.parameters().maximum_grid_depth=2;
    std::cout <<  analyser.parameters() << std::endl;

    // Set the grid
    watertank_system.set_grid(Grid(Vector<Float>(2, 1.0, 0.1)));

    HybridTime reach_time(80.0,5);

    Float eps = 1.0/4.0;
    // GRID REFINEMENT LOOP
    //
    // This verification loop starts with a grid_depth of 2 (that is, cells of size 0.25), and proceeds as follows:
    //      1. Compute an over-approximation of the reachability set of the automaton
    //      2. Check if the result is safe (i.e., water level between minsafe and maxsafe)
    //      3. If the result is safe, exit with success.
    //      4. If the result is not safe, compute a lower-approximation
    //      5. Check if the lower-approx is not safe (i.e. water lever below minsafe-eps or above maxsafe+eps)
    //      6. If the lower-approx is not safe, exit with false, 
    //      7. Otherwise, increade grid_depth by 1 and repeat.
    //
    std::cout << "Starting verification loop. Water level should be kept between "<< minsafe<< " and "<< maxsafe<<std::endl<<std::endl;
    int grid_depth;
    for(grid_depth = min_grid_depth; grid_depth <= max_grid_depth ; grid_depth+=1) {
        std::cout << "Computing upper reach set with grid depth " << grid_depth <<" and cell size "<<eps<<"..." << std::flush;
        analyser.parameters().maximum_grid_depth=grid_depth;
        HybridGridTreeSet upper_reach_set = analyser.upper_reach(watertank_system,initial_set,reach_time);
        std::cout << "done." << std::endl;
        char filename[30];
        sprintf(filename, "watertank-upper_reach-%d", grid_depth);
        plot(filename,bounding_box, Colour(0.0,0.5,1.0), upper_reach_set);
        textplot(filename, upper_reach_set);
        
        Box check_box(2, minsafe,maxsafe, -1.0,2.0);
        HybridBoxes hcheck_box;
        hcheck_box[l1]=check_box;
        hcheck_box[l2]=check_box;
        hcheck_box[l3]=check_box;
        hcheck_box[l4]=check_box;
        
        Box lcheck_box(2, minsafe-eps,maxsafe+eps, -1.0,2.0);
        HybridBoxes hlcheck_box;
        hlcheck_box[l1]=lcheck_box;
        hlcheck_box[l2]=lcheck_box;
        hlcheck_box[l3]=lcheck_box;
        hlcheck_box[l4]=lcheck_box;
        
            
        if(upper_reach_set.subset(hcheck_box)) {
            std::cout << "Result is safe, exiting refinement loop." << std::endl;
            break;
        } else {
            std::cout << "Upper reach set is not safe, computing lower reach set... " << std::flush;
            HybridGridTreeSet lower_reach_set = analyser.lower_reach(watertank_system,initial_set,reach_time);
            std::cout << "done." << std::endl;
            sprintf(filename, "watertank-lower_reach-%d", grid_depth);
            plot(filename,bounding_box, Colour(0.0,0.5,1.0), lower_reach_set);
            textplot(filename,lower_reach_set);
            
            if(!lower_reach_set.subset(hlcheck_box)) {
                std::cout << "Lower reach set is not safe, exiting." << std::endl;
                break;
            }
            std::cout << "Lower reach set is safe, refining the grid." << std::endl;
        } 
        eps = eps/2.0;
    }

    if(grid_depth > max_grid_depth) {
        std::cout << "Maximum allowed grid depth reached, cannot verify the system."<<std::endl;
    }
}
