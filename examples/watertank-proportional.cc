/***************************************************************************
 *            watertank-proportional.cc
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
#include <include/hybrid_set.h>

using namespace Ariadne;

typedef GeneralHybridEvolver::EnclosureListType EnclosureListType;
typedef GeneralHybridEvolver::EnclosureType EnclosureType;
typedef EnclosureType::ContinuousStateSetType ContinuousEnclosureType;


HybridGridTreeSet
outer_approximation(const EnclosureListType& hls,
                    const HybridGrid& hgr,
                    const int accuracy)
{
    HybridGridTreeSet result(hgr);
    for(EnclosureListType::const_iterator
            iter=hls.begin(); iter!=hls.end(); ++iter)
        {
            DiscreteLocation loc=iter->location();
            const ContinuousEnclosureType& es=iter->continuous_state_set();
            GridTreeSet& gts=result[loc];
            gts.adjoin_outer_approximation(ImageSet(es.bounding_box()),accuracy);
            //gts.adjoin_outer_approximation(ModelSet<ES>(es),accuracy);
        }
    return result;
}


int main(int argc,char *argv[])
{
    if(argc != 3) {
      std::cerr << "Usage: watertank-proportional bmin bmax" <<std::endl;
      return 1;
    }

    /// Set the system parameters
    Real a("0.02");
    Real r("1.25");
    double Rif = 5.67;
    double Kp = 15;

    double bmin = atoi(argv[1])*0.0001;
    double bmax = atoi(argv[2])*0.0001;

    double bstep = 0.0025;
    double dstep = 0.01;
    double Delta = 0.05;

    std::cout << "bmin = " << bmin <<", bmax = "<< bmax << std::endl << std::flush;


    // System variables
    RealScalarFunction x=RealScalarFunction::coordinate(5,0); // water level
    RealScalarFunction y=RealScalarFunction::coordinate(5,1); // valve level
    RealScalarFunction b=RealScalarFunction::coordinate(5,2); // input pressure
    RealScalarFunction delta=RealScalarFunction::coordinate(5,3); // sensor error
    RealScalarFunction t=RealScalarFunction::coordinate(5,4); // time

    RealVariable water("water"); // Water level
    RealVariable aperture("aperture"); // Valve aperture
    RealVariable pressure("pressure"); // Input pressure
    RealVariable error("error"); // Sensor error
    RealVariable time("time"); // Time


    /// Build the Hybrid System

    /// Create a HybridAutomton object
    HybridAutomaton watertank_system;

    /// Create discrete states
    DiscreteLocation zero_saturated("zero_saturated");      // Zero saturated
    DiscreteLocation not_saturated("not_saturated");      // Not saturated
    DiscreteLocation one_saturated("one_saturated");      // One saturated

    /// Create the discrete events
    DiscreteEvent event_zero_to_not("e0?");
    DiscreteEvent event_not_to_zero("e?0");
    DiscreteEvent event_not_to_one("e?1");
    DiscreteEvent event_one_to_not("e1?");

    /// Create the dynamics

    DottedRealAssignment dwater(dot(water)=-a*water+pressure*aperture);
    DottedRealAssignment dpressure(dot(pressure)=0);
    DottedRealAssignment derror(dot(error)=0);
    DottedRealAssignment dtime(dot(time)=1);

    DottedRealAssignment daperture_zero( dot(aperture) = -aperture/r );
    DottedRealAssignment daperture_not( dot(aperture) = -aperture/r*(Kp*(Rif-water-error)-aperture) );
    DottedRealAssignment daperture_one( dot(aperture) = (1-aperture)/r );

    DottedRealAssignments zero_saturated_dynamic((dwater,daperture_zero,dpressure,derror,dtime));
    DottedRealAssignments not_saturated_dynamic((dwater,daperture_not,dpressure,derror,dtime));
    DottedRealAssignments one_saturated_dynamic((dwater,daperture_one,dpressure,derror,dtime));
    //RealVectorFunction notsaturated_d((-a*x+b*y,-y/r*(Kp*(Rif-x-delta)-y),zero,zero,one));
    //RealVectorFunction onesaturated_d((-a*x+b*y,(1-y)/r,zero,zero,one));

    cout << "zero-saturated dynamic = " << zero_saturated_dynamic << "\n\n";
    cout << "not-saturated dynamic = " << not_saturated_dynamic << "\n\n";
    cout << "one-saturated dynamic = " << one_saturated_dynamic << "\n\n";

    /// Create the resets
    PrimedRealAssignments reset_id((next(water)=water,next(aperture)=aperture,next(pressure)=pressure,next(error)=error,next(time)=time));
    cout << "reset_id="<< reset_id << endl << endl;

    /// Create the guards.
    /// x <= Rif - Delta
    ContinuousPredicate guard12(water<=Rif-error);
    cout << "guard12=" << guard12 << endl << endl;
    /// x >= Rif - Delta
    ContinuousPredicate guard21(water>=Rif-error);
    cout << "guard21=" << guard21 << endl << endl;
    /// x <= Rif - 1/Kp - Delta
    ContinuousPredicate guard23(water<=Rif-1/Kp-error);
    cout << "guard23=" << guard23 << endl << endl;
    /// x >= Rif - 1/Kp - Delta
    ContinuousPredicate guard32(water>=Rif-1/Kp-error);
    cout << "guard32=" << guard32 << endl << endl;

    /// Create the invariants.
    /// Invariants are true when f(x) = Ax + b < 0
    /// forced transitions do not need an explicit invariant,
    /// hence we do not need invariants

    /// Build the automaton
    watertank_system.new_mode(zero_saturated,zero_saturated_dynamic);
    watertank_system.new_mode(not_saturated,not_saturated_dynamic);
    watertank_system.new_mode(one_saturated,one_saturated_dynamic);

    watertank_system.new_transition(zero_saturated,event_zero_to_not,not_saturated,reset_id,guard12,urgent);
    watertank_system.new_transition(not_saturated,event_not_to_zero,zero_saturated,reset_id,guard21,urgent);
    watertank_system.new_transition(not_saturated,event_not_to_one,one_saturated,reset_id,guard23,urgent);
    watertank_system.new_transition(one_saturated,event_one_to_not,not_saturated,reset_id,guard32,urgent);

    /// Finished building the automaton

    cout << "Automaton = " << watertank_system << endl << endl;

    /// Compute the system evolution

    /// Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver;
    evolver.verbosity = 1;

    /// Set the evolution parameters
    double maximum_step_size= 0.05;
    evolver.parameters().maximum_enclosure_radius = 0.5;
    evolver.parameters().maximum_step_size = maximum_step_size;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolver::EnclosureType HybridEnclosureType;
    typedef GeneralHybridEvolver::OrbitType OrbitType;
    typedef GeneralHybridEvolver::EnclosureListType EnclosureListType;

    double time_step = 0.25;
    double total_time = 35.0;
    double skip_time = 35.0;
    HybridTime evolution_time(skip_time,6);

    Box graphic_box(2, 18.0,skip_time, 5.0,6.0);
    Array<uint> tx(2,4,0);

    HybridScaling scaling( (water|0.25, aperture|1.0, pressure|1.0, error|1.0, time|1.0) );
    HybridGrid hg(watertank_system.state_space(),scaling);
    HybridGridTreeSet hgts(hg);
    uint grid_depth = 9;
    uint grid_height = 8;

    std::cout << "Computing timed evolution starting from location one_saturated, x = 0.0, y = 1.0 for " << skip_time << " seconds" << std::endl;
    for(double b=bmin ; b < bmax+bstep ; b += bstep) {
        for(double d=-Delta ; d < Delta+dstep ; d += dstep) {
            cout << "b = "<< b <<", Delta = "<<d<<std::endl;
            RealVariableBox initial_box((water==0.0, aperture==1.0, pressure==b, error==d, time==0.0));
            HybridSet initial_set(one_saturated,initial_box);
            OrbitType result = evolver.orbit(watertank_system,initial_set,evolution_time,UPPER_SEMANTICS);
            cout<<"Orbit.final=" << result.final() << endl;
            /*cout<<"Adjoining result to the grid..."<<std::flush;
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
            g2.clear();*/
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
    initial_set[zero_saturated]=result.final()[zero_saturated][0].range();

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
    Array<uint> tx(2,4,0);
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
//    g2.write("watertank-dominato-time-not_saturated");

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
        //     std::cout << "final set[zero_saturated] = " << final[zero_saturated] << endl;
        // std::cout << "final set[zero_saturated][0] = " << final[zero_saturated][0].bounding_box() << endl  << endl;
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
    Array<uint> p(2,0,1);
    g.set_projection_map(ProjectionFunction(p,4));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << orbit;
    g.write("watertank-dominato-orbit");


    std::cout << "Computing reach set using GeneralHybridEvolver... " << std::flush;
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

    std::cout << "Computing evolution starting from location zero_saturated, x = 0.0, y = 0.0" << std::endl;

    Box initial_box2(2, 0.0,0.001, 0.0,0.001);
    HybridImageSet initial_set2;
    initial_set2[zero_saturated]=initial_box2;

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
