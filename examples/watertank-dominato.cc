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

//
// Definition of the nonlinear dynamics for locations Opening and Closing
// the parameters are a and T
//
struct OpenValve : FunctionData<4,4,2> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
        r[0] = - p[0] * x[0] + x[2]; // xdot = -ax + b
        r[1] = 0.0;                    // ydot = 0
        r[2] = 0.0;                         // b and Delta are constant
        r[3] = 0.0;
    }
};

struct Opening : FunctionData<4,4,2> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
        r[0] = - p[0] * x[0] + x[2] * x[1]; // xdot = -ax + by
        r[1] = 1.0/p[2];                    // ydot = 1/T
        r[2] = 0.0;                         // b and Delta are constant
        r[3] = 0.0;
    }
};

struct Closing : FunctionData<4,4,2> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
        r[0] = - p[0] * x[0] + x[2] * x[1]; // xdot = -ax + by
        r[1] = -1.0/p[2];                   // ydot = -1/T
        r[2] = 0.0;                         // b and Delta are constant
        r[3] = 0.0;
    }
};

struct ClosedValve : FunctionData<4,4,2> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
        r[0] = - p[0] * x[0]; // xdot = -ax
        r[1] = 0.0;                    // ydot = 0
        r[2] = 0.0;                         // b and Delta are constant
        r[3] = 0.0;
    }
};

int main() 
{
     
    /// Set the system parameters
    double a = 0.02;
    double bmin = 0.3;
    double bmax = 0.34;
    double T = 1.5;
    double hmin = 5.55;
    double Delta = 0.05;
    double hmax = 5.70;

    Vector<Interval> system_parameters(2);
    system_parameters[0] = a;
    system_parameters[1] = T;

    // System variables:
    // x: water level
    // y: valve level
    // b: input pressure
    // Delta: sensor error

    // xdot = -a x + b, ydot = 0
    double A1[16]={-a,0,1,0,
                  0,0,0,0,
                  0,0,0,0,
                  0,0,0,0};
    double b1[4]={0,0,0,0};

  
    // xdot = -a x + b, ydot = -1/T
    double A2[16]={-a,0,1,0,
                  0,0,0,0,
                  0,0,0,0,
                  0,0,0,0};
    double b2[4]={0,-1.0/T,0,0};

    // xdot = -a x, ydot = 0
    double A3[16]={-a,0,0,0,
                  0,0,0,0,
                  0,0,0,0,
                  0,0,0,0};
    double b3[4]={0,0,0,0};

    // xdot = -a x, ydot = 1/T
    double A4[16]={-a,0,0,0,
                  0,0,0,0,
                  0,0,0,0,
                  0,0,0,0};
    double b4[4]={0,1.0/T,0,0};


    /// Build the Hybrid System
  
    /// Create a HybridAutomton object
    HybridAutomaton watertank_system;
  
    /// Create four discrete states
    DiscreteState l1(1);      // Valve open
    DiscreteState l2(2);      // Closing
    DiscreteState l3(3);      // Valve closed
    DiscreteState l4(4);      // Opening
  
    /// Create the discrete events
    DiscreteEvent e12(12);
    DiscreteEvent e23(23);
    DiscreteEvent e24(24);
    DiscreteEvent e34(34);
    DiscreteEvent e41(41);
    DiscreteEvent e42(42);
  
    /// Create the dynamics
    AffineFunction dynamic1(Matrix<Float>(4,4,A1),Vector<Float>(4,b1));
    AffineFunction dynamic2(Matrix<Float>(4,4,A2),Vector<Float>(4,b2));
    AffineFunction dynamic3(Matrix<Float>(4,4,A3),Vector<Float>(4,b3));
    AffineFunction dynamic4(Matrix<Float>(4,4,A4),Vector<Float>(4,b4));
//    Function<OpenValve> dynamic1(system_parameters);
//    Function<Closing> dynamic2(system_parameters);
//    Function<ClosedValve> dynamic3(system_parameters);
//    Function<Opening> dynamic4(system_parameters);
    
    cout << "dynamic1 = " << dynamic1 << endl << endl;
    cout << "dynamic2 = " << dynamic2 << endl << endl;
    cout << "dynamic3 = " << dynamic3 << endl << endl;
    cout << "dynamic4 = " << dynamic4 << endl << endl;

    /// Create the resets
    AffineFunction reset_y_zero(Matrix<Float>(4,4,
                                1.0,0.0,0.0,0.0,
                                0.0,0.0,0.0,0.0,
                                0.0,0.0,1.0,0.0,
                                0.0,0.0,0.0,1.0),
                                Vector<Float>(4,0.0,0.0,0.0,0.0));
    cout << "reset_y_zero=" << reset_y_zero << endl << endl;
    AffineFunction reset_y_one(Matrix<Float>(4,4,
                                1.0,0.0,0.0,0.0,
                                0.0,0.0,0.0,0.0,
                                0.0,0.0,1.0,0.0,
                                0.0,0.0,0.0,1.0),
                                Vector<Float>(4,0.0,1.0,0.0,0.0));
    cout << "reset_y_one=" << reset_y_one << endl << endl;
    IdentityFunction reset_id(4);
    cout << "reset_id="<< reset_id << endl << endl;

    /// Create the guards.
    /// Guards are true when f(x) = Ax + b > 0
    /// x >= hmax - Delta
    AffineFunction guard1(Matrix<Float>(1,4,1.0,0.0,0.0,1.0),
                          Vector<Float>(1,-hmax));
    cout << "guard1=" << guard1 << endl << endl;
    /// y < 0
    AffineFunction guard2(Matrix<Float>(1,4,0.0,-1.0,0.0,0.0),
                          Vector<Float>(1,0.0));
    cout << "guard2=" << guard2 << endl << endl;
    /// x <= hmin - delta
    AffineFunction guard3(Matrix<Float>(1,4,-1.0,0.0,0.0,-1.0),
                          Vector<Float>(1,hmin));
    cout << "guard3=" << guard3 << endl << endl;
    /// y >= 1
    AffineFunction guard4(Matrix<Float>(1,4,0.0,1.0,0.0,0.0),
                          Vector<Float>(1,-1.0));
    cout << "guard4=" << guard4 << endl << endl;

    /// Create the invariants.
    /// Invariants are true when f(x) = Ax + b < 0
    /// forced transitions do not need an explicit invariant, 
    /// hence we do not need invariants
  
    /// Build the automaton
    watertank_system.new_mode(l1,dynamic1);
    watertank_system.new_mode(l2,dynamic2);
    watertank_system.new_mode(l3,dynamic3);
    watertank_system.new_mode(l4,dynamic4);

    watertank_system.new_forced_transition(e12,l1,l2,reset_y_one,guard1);
    watertank_system.new_forced_transition(e23,l2,l3,reset_y_zero,guard2);
    watertank_system.new_forced_transition(e24,l2,l4,reset_id,guard3);    
    watertank_system.new_forced_transition(e34,l3,l4,reset_y_zero,guard3);
    watertank_system.new_forced_transition(e41,l4,l1,reset_y_one,guard4);
    watertank_system.new_forced_transition(e42,l4,l2,reset_id,guard1);

    /// Finished building the automaton

    cout << "Automaton = " << watertank_system << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver;

    /// Set the evolution parameters
    double maximum_step_size= 0.05;
    evolver.parameters().maximum_enclosure_radius = 0.05;
    evolver.parameters().maximum_step_size = maximum_step_size;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;
    typedef ListSet<HybridEvolver::ContinuousEnclosureType> ListSetType;
    typedef HybridEvolver::TimedEnclosureListType TimedEnclosureListType;

    double time_step = 0.25;
    double total_time = 19.5;
    double skip_time = 25.0;
    HybridTime evolution_time(skip_time,1);

    std::cout << "Computing timed evolution starting from location l1, x = 0.0, y = 1.0 for " << skip_time << " seconds" << std::endl;

    Box initial_box(4, 0.0,0.0, 1.0,1.0, bmin,bmax, -Delta,Delta);
    HybridEnclosureType initial_enclosure(l1,initial_box);

    TimedEnclosureListType result = evolver.timed_evolution(watertank_system,initial_enclosure,evolution_time,UPPER_SEMANTICS,true);

    Box graphic_box(2, 18.0,32.0 , 5.0,6.0);
    Figure g1, g2;
    g1.set_bounding_box(graphic_box);
    g1 << fill_colour(Colour(1.0,1.0,0.0));
    g1 << Box(2, 18,32, hmax - Delta, hmax + Delta);
    g1 << fill_colour(Colour(0.0,1.0,1.0));
    g1 << Box(2, 18,32, hmin - Delta, hmin + Delta);
    g1 << fill_colour(Colour(0.0,0.5,1.0));

    g2.set_bounding_box(graphic_box);
    g2 << fill_colour(Colour(1.0,1.0,0.0));
    g2 << Box(2, 18,32, hmax - Delta, hmax + Delta);
    g2 << fill_colour(Colour(0.0,1.0,1.0));
    g2 << Box(2, 18,32, hmin - Delta, hmin + Delta);
    g2 << fill_colour(Colour(0.0,0.5,1.0));

    cout << "Result size: " << result.size() << endl;
    for( int i = 0; i < result.size(); i++ ) {
        Interval t = result[i].first;
        Box b = result[i].second.second.bounding_box();
        Interval x = b[0];
        std::cout << "  t = " << t << ", loc = " << result[i].second.first << ", x = " << x << std::endl;
        if(result[i].second.first == l1) {
          g1 << Box(2, t.lower(), t.lower()+maximum_step_size, x.lower(), x.upper());        
        } else {
          g2 << Box(2, t.lower(), t.lower()+maximum_step_size, x.lower(), x.upper());        
        }
    }

    g1.write("watertank-dominato-time-l1");
    g2.write("watertank-dominato-time-l2");

/*  
    OrbitType orbit = evolver.timed_evolution(watertank_system,initial_enclosure,evolution_time,UPPER_SEMANTICS,true);
    EnclosureListType final = orbit.final();        
    typedef EnclosureListType::const_iterator const_iterator;
    double xmin=100.0;
    double xmax=-100.0;
    for(const_iterator iter=final.begin(); iter != final.end(); ++iter) {
        Interval x = iter->second.bounding_box()[0];
        if(x.lower() < xmin) xmin = x.lower();
        if(x.upper() > xmax) xmax = x.upper();
    }
    cout << " xmin = " << xmin << ", xmax = " << xmax << endl << endl;
      
    evolution_time.continuous_time=time_step;

    Box graphic_box(2, skip_time-0.1,skip_time+total_time+0.1, -0.1,6.1);
    Figure g;
    g.set_bounding_box(graphic_box);
    g << fill_colour(Colour(0.0,1.0,1.0));
    g << Box(2, skip_time,skip_time+time_step, xmin,xmax);
    g << fill_colour(Colour(0.0,0.5,1.0));

    std::cout << "Computing upper reach timed set... " << std::flush;
 
    std::cout << "t = " << skip_time << std::endl;
    std::cout << "final set = " << final << endl << endl;
    
    final = orbit.initial();
    skip_time = 0.0;
 
    for(double t = time_step; t <= total_time; t = t + time_step) {
        std::cout << "t = " << skip_time + t << flush;
        EnclosureListType reach;
        for(const_iterator iter=final.begin(); iter != final.end(); ++iter) {
        orbit = evolver.orbit(watertank_system,*iter,evolution_time,UPPER_SEMANTICS);
        reach.adjoin(orbit.final());
        std::cout << "." << flush;
        }
//        std::cout << endl << "final set before clearing = " << reach  << endl; 
        final.clear();
        // Remove redundant elements from final set
        int i = 0;
        for(const_iterator iter1=reach.begin(); iter1 != reach.end(); ++iter1) {
            bool flag = 0;
            for(const_iterator iter2=iter1; iter2 != reach.end() && flag == 0; ++iter2) {
                if( iter1!= iter2 && iter1->first == iter2->first ) {
                    Box b1 = iter1->second.bounding_box();
                    Box b2 = iter2->second.bounding_box();
                    if(subset(b1,b2)) {
                        // cout << "iter1 = " << iter1->second << endl ;
                        // cout << "is a subset of iter2 = " << iter1->second << endl;
                        flag = 1;
                        i++;
                    }                 
                }
             }  
            if (flag == 0) final.adjoin(*iter1);
        }        
//        cout << i << " elements removed. " << endl;
        cout << " final set after clearing = " << final << endl << endl;
        reach.clear();
        //     std::cout << "final set[l1] = " << final[l1] << endl;
        // std::cout << "final set[l1][0] = " << final[l1][0].bounding_box() << endl  << endl;
        xmin=100.0;
        xmax=-100.0;
        for(const_iterator iter=final.begin(); iter != final.end(); ++iter) {
            Interval x = iter->second.bounding_box()[0];
            if(x.lower() < xmin) xmin = x.lower();
            if(x.upper() > xmax) xmax = x.upper();
        }
        cout << " xmin = " << xmin << ", xmax = " << xmax << endl << endl;
        g << Box(2, skip_time+t,skip_time+t+time_step, xmin,xmax);
    }
      
    std::cout << "done." << std::endl;

    g.write("watertank-dominato-orbit");

    std::cout << "Orbit="<<orbit<<std::endl;
    Box bounding_box(2, -0.1,9.1, -0.1,1.1);
    Figure g;
    g.set_bounding_box(bounding_box);
    array<uint> p(2,0,1);
    g.set_projection_map(ProjectionFunction(p,4));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << orbit;
    g.write("watertank-dominato-orbit");


    std::cout << "Computing reach set using HybridEvolver... " << std::flush;
    EnclosureListType reach = evolver.reach(watertank_system,initial_enclosure,evolution_time);
    std::cout << "done." << std::endl;

    std::cout << "Orbit="<<reach<<std::endl;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    plot("watertank-reach-evolver",bounding_box, Colour(0.0,0.5,1.0), reach);


    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time = 32.0;
    analyser.parameters().grid_lengths = 0.05;
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[l2]=initial_box;

    HybridTime reach_time(64.0,2);

    plot("watertank-initial_set1",bounding_box, Colour(0.0,0.5,1.0), initial_set);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)
    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet* lower_reach_set_ptr = analyser.lower_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-lower_reach1",bounding_box, Colour(0.0,0.5,1.0), *lower_reach_set_ptr);

    // Compute evolved sets and reach sets using upper semantics.
    // These functions compute over-approximations to the evolved and reachabe sets. Subdivision is used
    // as necessary to keep the local errors reasonable. The accumulated global error may be very large.
    std::cout << "Computing upper reach set... " << std::flush;
    HybridGridTreeSet* upper_reach_set_ptr = analyser.upper_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-upper_reach1",bounding_box, Colour(0.0,0.5,1.0), *upper_reach_set_ptr);

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
