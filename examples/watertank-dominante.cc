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

typedef HybridEvolver::EnclosureListType EnclosureListType;
typedef HybridEvolver::ContinuousEnclosureType ContinuousEnclosureType;


HybridGridTreeSet 
outer_approximation(const EnclosureListType& hls,
                    const HybridGrid& hgr,
                    const int accuracy)
{
    HybridGridTreeSet result;
    for(EnclosureListType::const_iterator 
            iter=hls.begin(); iter!=hls.end(); ++iter)
        {
            AtomicDiscreteLocation loc=iter->first;
            const ContinuousEnclosureType& es=iter->second;
            if(result.find(loc)==result.locations_end()) {
                result.insert(make_pair(loc,GridTreeSet(hgr[loc])));
            }
            GridTreeSet& gts=result[loc];
            gts.adjoin_outer_approximation(ImageSet(es.range()),accuracy);
            //gts.adjoin_outer_approximation(ModelSet<ES>(es),accuracy);
        }
    return result;
}


//
// Definition of the nonlinear dynamics for locations Opening and Closing
// the parameters are a and T
//
struct ZeroSaturated : VectorFunctionData<5,5,4> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
        r[0] = - p[0] * x[0] + x[2] * x[1];   // xdot = -ax + b y
        r[1] = - x[1] / p[1];                 // ydot = -1/r y
        r[2] = 0.0;                           // b and Delta are constant
        r[3] = 0.0;
        r[4] = 1.0;
    }
};

struct NotSaturated : VectorFunctionData<5,5,4> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
        r[0] = - p[0] * x[0] + x[2] * x[1];   // xdot = -ax + b y
        r[1] = ( p[3] * ( p[2] - x[0] - x[3] ) - x[1] ) / p[1];  // ydot = 1/r (Kp(R-x-Delta)-y)
        r[2] = 0.0;                         // b and Delta are constant
        r[3] = 0.0;
        r[4] = 1.0;
    }
};

struct OneSaturated : VectorFunctionData<5,5,4> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
        r[0] = - p[0] * x[0] + x[2] * x[1];   // xdot = -ax + b y
        r[1] = (1.0 - x[1]) / p[1];                 // ydot = (1 - y)/r
        r[2] = 0.0;                         // b and Delta are constant
        r[3] = 0.0;
        r[4] = 1.0;
    }
};


int main(int argc,char *argv[]) 
{
    if(argc != 3) {
      std::cerr << "Usage: watertank-dominante bmin bmax" <<std::endl;
      return 1;
    }
    
    /// Set the system parameters
    double a = 0.02;
    double bmin = atoi(argv[1])*0.0001;
    double bmax = atoi(argv[2])*0.0001;
//    double bmax = 0.29999;
    double r = 1.25;
    double Rif = 5.67;
//    double Delta = 0.001;
    double Delta = 0.05;
    double Kp = 15;
    double bstep = 0.0025;
    double dstep = 0.01;
    
    std::cout << "bmin = " << bmin <<", bmax = "<< bmax << std::endl << std::flush;

    Vector<Interval> system_parameters(4);
    system_parameters[0] = a;
    system_parameters[1] = r;
    system_parameters[2] = Rif;
    system_parameters[3] = Kp;
    

    // System variables:
    // x: water level
    // y: valve level
    // b: input pressure
    // Delta: sensor error

    /// Build the Hybrid System
  
    /// Create a HybridAutomton object
    MonolithicHybridAutomaton watertank_system;
  
    /// Create four discrete states
    AtomicDiscreteLocation l1(1);      // Zero saturated
    AtomicDiscreteLocation l2(2);      // Not saturated
    AtomicDiscreteLocation l3(3);      // One saturated
  
    /// Create the discrete events
    DiscreteEvent e12(12);
    DiscreteEvent e21(21);
    DiscreteEvent e23(23);
    DiscreteEvent e32(32);
  
    /// Create the dynamics
    VectorUserFunction<ZeroSaturated> dynamic1(system_parameters);
    VectorUserFunction<NotSaturated> dynamic2(system_parameters);
    VectorUserFunction<OneSaturated> dynamic3(system_parameters);
    
    cout << "dynamic1 = " << dynamic1 << endl << endl;
    cout << "dynamic2 = " << dynamic2 << endl << endl;
    cout << "dynamic3 = " << dynamic3 << endl << endl;

    /// Create the resets
    IdentityFunction reset_id(5);
    cout << "reset_id="<< reset_id << endl << endl;

    /// Create the guards.
    /// Guards are true when f(x) = Ax + b > 0
    /// x <= Rif - Delta
    VectorAffineFunction guard12(Matrix<Float>(1,5, -1.0,0.0,0.0,-1.0,0.0),
                          Vector<Float>(1,Rif));
    cout << "guard12=" << guard12 << endl << endl;
    /// x >= Rif - Delta
    VectorAffineFunction guard21(Matrix<Float>(1,5, 1.0,0.0,0.0,1.0,0.0),
                          Vector<Float>(1,-Rif));
    cout << "guard21=" << guard21 << endl << endl;
    /// x <= Rif - 1/Kp - Delta
    VectorAffineFunction guard23(Matrix<Float>(1,5, -1.0,0.0,0.0,-1.0,0.0),
                          Vector<Float>(1,(Rif-1.0/Kp)));
    cout << "guard23=" << guard23 << endl << endl;
    /// x >= Rif - 1/Kp - Delta
    VectorAffineFunction guard32(Matrix<Float>(1,5, 1.0,0.0,0.0,1.0,0.0),
                          Vector<Float>(1,(1.0/Kp - Rif)));
    cout << "guard32=" << guard32 << endl << endl;

    /// Create the invariants.
    /// Invariants are true when f(x) = Ax + b < 0
    /// forced transitions do not need an explicit invariant, 
    /// hence we do not need invariants
  
    /// Build the automaton
    watertank_system.new_mode(l1,dynamic1);
    watertank_system.new_mode(l2,dynamic2);
    watertank_system.new_mode(l3,dynamic3);

    watertank_system.new_forced_transition(e12,l1,l2,reset_id,guard12);
    watertank_system.new_forced_transition(e21,l2,l1,reset_id,guard21);
    watertank_system.new_forced_transition(e23,l2,l3,reset_id,guard23);    
    watertank_system.new_forced_transition(e32,l3,l2,reset_id,guard32);

    /// Finished building the automaton

    cout << "Automaton = " << watertank_system << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver;
    evolver.verbosity = 1;

    /// Set the evolution parameters
    double maximum_step_size= 0.05;
    evolver.parameters().maximum_enclosure_radius = 0.5;
    evolver.parameters().maximum_step_size = maximum_step_size;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;
    typedef ListSet<HybridEvolver::ContinuousEnclosureType> ListSetType;
    typedef HybridEvolver::TimedEnclosureListType TimedEnclosureListType;

    double time_step = 0.25;
    double total_time = 35.0;
    double skip_time = 35.0;
    HybridTime evolution_time(skip_time,6);

    Box graphic_box(2, 18.0,skip_time, 5.0,6.0);
    array<uint> tx(2,4,0);

    Vector<Float> lengths(5, 0.25, 1.0, 1.0, 1.0, 1.0);
    HybridGridTreeSet hgts(watertank_system.state_space(), lengths);
    uint grid_depth = 9;
    uint grid_height = 8;
    
    std::cout << "Computing timed evolution starting from location l3, x = 0.0, y = 1.0 for " << skip_time << " seconds" << std::endl;
    for(double b=bmin ; b < bmax+bstep ; b += bstep) {
        for(double d=-Delta ; d < Delta+dstep ; d += dstep) {
            cout << "b = "<< b <<", Delta = "<<d<<std::endl;
            Box initial_box(5, 0.0,0.0, 1.0,1.0, b,b, d,d, 0.0,0.0);
            HybridEnclosureType initial_enclosure(l3,initial_box);
            OrbitType result = evolver.orbit(watertank_system,initial_enclosure,evolution_time,UPPER_SEMANTICS);
            cout<<"Orbit.final=" << result.final() << endl;
            cout<<"Adjoining result to the grid..."<<std::flush;
            hgts.adjoin(outer_approximation(result.reach(),hgts.grid(),grid_depth));
            cout<<"done:"<<hgts.size()<<" total cells."<<std::endl;
            char filename[30];
            sprintf(filename,"wt-best-%d-%d",int(b*10000),int(d*10000));
            cout<<"Saving result to "<<filename<<"..."<<std::flush;
            Figure g2;            
            g2.set_bounding_box(graphic_box);
            g2.set_projection_map(ProjectionFunction(tx,5));        
            g2 << fill_colour(Colour(0.9,0.9,0.0));
            g2 << hgts;
            g2.write(filename);
            g2.clear();
            cout<<"done."<<endl<<std::flush;
        }
    }

/*
    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.verbosity = 4;
    analyser.parameters().lock_to_grid_time = total_time;
    analyser.parameters().maximum_grid_depth= 7;
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[l1]=result.final()[l1][0].range();

    HybridTime reach_time((total_time-skip_time),4);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)
    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet* lower_reach_set_ptr = analyser.lower_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;


//    Box graphic_box(2, 18.0,32.0 , 5.0,6.0);
    Box graphic_box(2, skip_time,total_time , 5.0,7.0);
    Figure g1;
    array<uint> tx(2,4,0);
    g1.set_bounding_box(graphic_box);
    g1.set_projection_map(ProjectionFunction(tx,5));    
    g1 << Box(2, 18,32, hmax - Delta, hmax + Delta);
    g1 << fill_colour(Colour(0.0,1.0,1.0));
    g1 << Box(2, 18,32, hmin - Delta, hmin + Delta);

    g1 << fill_colour(Colour(0.0,0.5,1.0));
    g1 << result;
    
    g1 << fill_colour(Colour(1.0,1.0,0.0));
    g1 << result.final();
    
    g1 << fill_colour(Colour(0.0,1.0,1.0));
    g1 << *lower_reach_set_ptr;

    g1.write("watertank-dominato-time");
//    g2.write("watertank-dominato-time-l2");

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
    return 0;
}
