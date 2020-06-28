/***************************************************************************
 *            LOVO20.hpp
 *
 *  Copyright  2020  Luca Geretti
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

#include "arch.hpp"

using namespace Ariadne;

void LOVO20()
{
    ArchBenchmark benchmark("LOVO20");

    ARIADNE_LOG_PRINTLN("Lotka-Volterra benchmark " << benchmark.name() << ":");

    RealVariable x("x");
    RealVariable y("y");
    RealVariable cnt("cnt");

    RealConstant alpha("alpha",3);
    RealConstant beta("beta",3);
    RealConstant gamma("gamma",1);
    RealConstant delta("delta",1);
    RealConstant radius("radius",0.15_dec);

    StringVariable lotkavolterra("lotkavolterra");
    StringConstant outside("outside");
    StringConstant inside("inside");
    HybridAutomaton automaton(lotkavolterra.name());

    DiscreteEvent enter("enter");
    DiscreteEvent exit("exit");

    Real cx = gamma/delta;
    Real cy = alpha/beta;

    automaton.new_mode(lotkavolterra|outside,{dot(x)=alpha*x-beta*x*y, dot(y)= delta*x*y-gamma*y, dot(cnt)=0});
    automaton.new_mode(lotkavolterra|inside,{dot(x)=alpha*x-beta*x*y, dot(y)= delta*x*y-gamma*y, dot(cnt)=1});
    automaton.new_transition(lotkavolterra|outside,enter,lotkavolterra|inside,{next(x)=x,next(y)=y,next(cnt)=cnt},sqr(x-cx)+sqr(y-cy)<=sqr(radius),EventKind::IMPACT);
    automaton.new_transition(lotkavolterra|inside,exit,lotkavolterra|outside,{next(x)=x,next(y)=y,next(cnt)=cnt},sqr(x-cx)+sqr(y-cy)>=sqr(radius),EventKind::IMPACT);

    MaximumError max_err=1e-5;
    TaylorPicardIntegrator integrator(max_err);

    GeneralHybridEvolver evolver(automaton);
    evolver.set_integrator(integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.09);
    evolver.configuration().set_maximum_spacial_error(1e-5);

    RealPoint ic({1.3_dec,1.0_dec});
    Real e(0.008_dec);
    HybridSet initial_set(lotkavolterra|outside,{ic[0]-e<=x<=ic[0]+e,y==ic[1],cnt==0});
    HybridTime evolution_time(3.64,3);

    StopWatch sw;

    ARIADNE_LOG_PRINTLN("Computing evolution... ");
    auto orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);

    ARIADNE_LOG_PRINTLN("Checking properties... ");

    ARIADNE_LOG_RUN_AT(2,auto final_bounds = orbit.final().bounding_box());

    Bool has_any_zero_transitions;
    Bool has_any_two_transitions;
    for (HybridEnclosure encl : orbit.final()) {
        if ((not has_any_zero_transitions) and encl.previous_events().size() == 0)
            has_any_zero_transitions = true;
        if ((not has_any_two_transitions) and encl.previous_events().size() == 2)
            has_any_two_transitions = true;
    }

    sw.click();

    if (has_any_zero_transitions) ARIADNE_LOG_PRINTLN("A final set with zero transitions has been found correctly.")
    else ARIADNE_LOG_PRINTLN("No final set with zero transitions has been found!");
    if (has_any_two_transitions) ARIADNE_LOG_PRINTLN("A final set with two transitions has been found correctly.")
    else ARIADNE_LOG_PRINTLN("No final set with two transitions has been found!");

    ARIADNE_LOG_PRINTLN("Trajectory stays within " << (radius*100).get_d() << "% of the equilibrium for at most " <<
                        final_bounds[cnt].upper() << " time units: " << ((definitely(final_bounds[cnt] <= 0.2_dec)) ?
                        "constraint satisfied." : "constraint not satisfied!"));

    ARIADNE_LOG_PRINTLN("Done in " << sw.elapsed() << " seconds.");
    ARIADNE_LOG_PRINTLN("# of final sets: " << orbit.final().size());
    ARIADNE_LOG_PRINTLN("Final set area: " << final_bounds[x].width()*final_bounds[y].width());

    auto instance = benchmark.create_instance();
    if (has_any_zero_transitions and has_any_two_transitions and (definitely(final_bounds[cnt] <= 0.2_dec))) {
        instance.set_verified(1)
                .set_execution_time(sw.elapsed())
                .add_loss((final_bounds[x].width()*final_bounds[y].width()).get_d())
                .add_loss(final_bounds[cnt].upper().get_d());
    }
    instance.write();

    HybridAutomaton circle;
    DiscreteLocation rotate;
    RealVariable t("t");
    circle.new_mode(rotate,{let(x)=cx+radius*cos(t),let(y)=cy-radius*sin(t)},{dot(t)=1});
    HybridSimulator simulator;
    simulator.set_step_size(0.1);
    HybridRealPoint circle_initial(rotate,{t=0});
    HybridTime circle_time(2*pi,1);
    ARIADNE_LOG_RUN_MUTED(auto circle_orbit = simulator.orbit(circle,circle_initial,circle_time));

    ARIADNE_LOG_RUN_AT(2,plot(benchmark.name().c_str(),Axes2d(0.6<=x<=1.4,0.6<=y<=1.4), ariadneorange, orbit, black, circle_orbit));
    ARIADNE_LOG_PRINTLN("File " << benchmark.name() << ".png written.");
}
