/***************************************************************************
 *            tutorial.cc
 *
 *  Copyright  2008  Pieter Collins
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


//! \file tutorial.cc

#include "ariadne.h"

#include "hybrid_automaton-composite.h"
#include "hybrid_set.h"
#include "hybrid_evolver.h"
#include "hybrid_simulator.h"
#include "hybrid_graphics.h"
#include <include/hybrid_enclosure.h>

using namespace Ariadne;

typedef GeneralHybridEvolver HybridEvolverType;


void press_enter_to_continue() {
    std::cout << "Press ENTER to continue... " << flush;
    std::cin.ignore( std::numeric_limits <std::streamsize> ::max(), '\n' );
}

template<class SET1,class SET2,class SET3,class SET4>
void nolines_plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2,
          const Colour& fc3, const SET3& set3, const Colour& fc4, const SET4& set4) {
    HybridFigure g;  g.set_axes(axes); g.set_line_style(false); g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2);
    g.set_fill_colour(fc3); draw(g,set3); g.set_fill_colour(fc4); draw(g,set4); g.write(filename); }


template<class SET1,class SET2,class SET3,class SET4,class SET5>
void nolines_plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2,
          const Colour& fc3, const SET3& set3, const Colour& fc4, const SET4& set4, const Colour& fc5, const SET5& set5) {
    HybridFigure g;  g.set_axes(axes); g.set_line_style(false); g.set_line_width(0.0); g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2);
    g.set_fill_colour(fc3); draw(g,set3); g.set_fill_colour(fc4); draw(g,set4);
    g.set_fill_colour(fc5); draw(g,set5); g.write(filename); }

int main(int argc, const char* argv[])
{
    uint evolver_verbosity = 0;
    if(argc>1) { evolver_verbosity=atoi(argv[1]); }

    // Create the system
    // Set the system dynamic parameters
    RealConstant P("P",4.0);
    RealConstant K("K",1.0);
    RealConstant Tav("Tav",16.0);
    RealConstant Tamp("Tamp",8.0);

    // Set the system control parameters
    RealConstant Tmax("Tmax",23.0);
    RealConstant Tmin("Tmin",14.0);
    RealConstant Toff("Toff",21.0);
    RealConstant Ton("Ton",15.0);
    RealConstant Ton_upper("Ton_upper",15.25);
    RealConstant Ton_lower("Ton_lower",14.75);

    // Create the discrete states
    StringVariable heating("heating");
    StringConstant on("on");
    StringConstant off("off");

    // Create the discrete events
    DiscreteEvent must_switch_on("must_switch_on");
    DiscreteEvent switch_on("switch_on");
    DiscreteEvent switch_off("switch_off");
    DiscreteEvent midnight("midnight");

    // Declare the system variables.
    RealVariable T("T");
    RealVariable C("C");
    TimeVariable t;

    cerr<<"WARNING: Using different event labels for guard and invariant.\n";
    // Create the heater subsystem
    HybridAutomaton heater;
    heater.new_mode( heating|on, (dot(T)=P+K*(Tav-Tamp*cos(2.0*pi*C)-T)) );
    heater.new_mode( heating|off, (dot(T)=K*(Tav-Tamp*cos(2.0*pi*C)-T)) );
    heater.new_invariant( heating|off, T>=Ton_lower, must_switch_on );
    heater.new_transition( heating|off, switch_on, heating|on, (next(T)=T), T<=Ton_upper, permissive );
    // Comment out above two lines and uncomment the line below to make the switch_on transition urgent
    //heater.new_transition( heating|off, switch_on, heating|on, (next(T)=T), T<=Ton_upper, urgent );
    heater.new_transition( heating|on, switch_off, heating|off, (next(T)=T), T>=Toff, urgent );

    // Create the clock subsystem
    HybridAutomaton clock;
    clock.new_mode( (dot(C)=1.0) );
    clock.new_transition( midnight, next(C)=0.0, C>=1.0, urgent );

    CompositeHybridAutomaton heating_system((clock,heater));
    cout << "heating_system=" << heating_system << "\n" << "\n";

    // Create the analyser classes

    TaylorSeriesIntegrator series_integrator(1e-3);
    series_integrator.set_maximum_spacial_order(6);
    series_integrator.set_maximum_temporal_order(12);
    series_integrator.verbosity=0;
    TaylorPicardIntegrator picard_integrator(1e-5);
    IntervalNewtonSolver solver(1e-12,8);

    // Create a GeneralHybridEvolver object
    HybridEvolverType evolver(heating_system);
    evolver.set_solver(solver);

    // Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(0.25);
    evolver.configuration().set_maximum_step_size(7.0/16);
    evolver.verbosity=evolver_verbosity;
    cout <<  evolver.configuration() << endl << endl;

    evolver.configuration().set_enable_reconditioning(true);
    evolver.configuration().set_enable_subdivisions(true);



    // Compute the system evolution

    // Set the initial set.
    double r=1.0/1024; double Ti=17.0; 
    Real Tinitmin=Ti+r; Real Tinitmax=Ti+3*r; Real Cinitmin=0+r; Real Cinitmax=0+3*r; // Tinit=16.0;
    HybridExpressionSet initial_set(heating|off, (Tinitmin<=T<=Tinitmax,Cinitmin<=C<=Cinitmax) );
    cout << "initial_set=" << initial_set << endl;
    // Compute the initial set as a validated enclosure.
    HybridEnclosure initial_enclosure = evolver.enclosure(initial_set);
    cout << "initial_enclosure="<<initial_enclosure << endl << endl;

    HybridTime evolution_time(2.75,127);
    cout << "evolution_time=" << evolution_time << endl;


    cout << "\nComputing orbit using series integrator... \n" << flush;
    evolver.set_integrator(series_integrator);
    Orbit<HybridEnclosure> series_orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    cout << "    done." << endl;

    cout << "\nComputed " << series_orbit.reach().size() << " reach enclosures and " << series_orbit.final().size() << " final enclosures.\n";
    
    DRAWING_METHOD = AFFINE_DRAW;
    DRAWING_ACCURACY += 1;

    Colour guard_colour(0.5,0.5,0.5);
    Colour midnight_guard_colour(0.75,0.75,0.75);
    Colour picard_orbit_colour(0.0,1.0,1.0);
    Colour series_orbit_colour(0.0,0.0,1.0);

    Real tmax=evolution_time.continuous_time();
    double dTmin=Tmin.value().get_d(); double dTmax=Tmax.value().get_d();
    HybridBox guard(heating|off,(Ton_lower.value()<=T<=Ton_upper.value(),0<=C<=1,0<=t<=tmax));
    HybridBox midnight_guard(heating|off,(dTmin<=T<=dTmax,0.0<=C<=1.0,1.0<=t<=2.0));
    cout << "\nPlotting time trace of orbit... " << flush;
    plot("heating-orbit-time.png",Axes2d(0.0<=t<=tmax,dTmin<=T<=dTmax), midnight_guard_colour, midnight_guard, guard_colour, guard, series_orbit_colour, series_orbit);
    cout << "done." << endl << endl;

    HybridSet init_set(initial_set,RealSpace({C,T}));
    HybridReachabilityAnalyser analyser(heating_system,GeneralHybridEvolverFactory());
    analyser.configuration().set_lock_to_grid_time(1+1.0/1024);
    analyser.configuration().set_lock_to_grid_steps(1);
    std::cerr<<"max grid depth="<<analyser.configuration().maximum_grid_depth();
    std::cerr<<"transient_time="<<analyser.configuration().transient_time();
    std::cerr<<"transient_steps="<<analyser.configuration().transient_steps();
    analyser.configuration().set_maximum_grid_depth(4);
    analyser.verbosity=5;
    cout << "\nComputing chain-reachable set... \n" << flush;
    HybridGridTreeSet chain_reach_set = analyser.outer_chain_reach(init_set);
    cout << "done." << endl << endl;
    cout << "\nPlotting chain-reachable set... " << flush;
    plot("heating-chainreach-time.png",Axes2d(0.0<=C<=1.0,dTmin<=T<=dTmax), guard_colour, guard, series_orbit_colour, chain_reach_set);
    cout << "done." << endl << endl;

}
