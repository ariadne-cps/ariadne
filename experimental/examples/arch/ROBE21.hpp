/***************************************************************************
 *            ROBE21.hpp
 *
 *  Copyright  2021  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "arch.hpp"

using namespace Ariadne;

void ROBE21() {

    ArchBenchmark benchmark("ROBE21");

    RealVariable x("x"), y("y"), z("z"), s("s");
    RealConstant alpha("alpha",0.4_dec);

    ARIADNE_LOG_PRINTLN("Robertson chemical reaction system:");

    StepMaximumError max_err = 1e-9;
    GradedTaylorSeriesIntegrator integrator(max_err);

    ListSet<LabelledEnclosure> reach1, reach2, reach3;

    Real evolution_time = 40;

    {
        RealConstant beta("beta",100);
        RealConstant gamma("gamma",1000);

        ARIADNE_LOG_PRINTLN("Instance 1:");

        VectorField dynamics({dot(x) = -alpha*x + beta*y*z, dot(y) = alpha*x - beta*y*z - gamma*sqr(y), dot(z) = gamma*sqr(y)},{let(s)=x+y+z});

        RealVariablesBox initial_set({x==1,y==0,z==0});

        Stopwatch<Milliseconds> sw;

        ARIADNE_LOG_PRINTLN_AT(1,"Computing orbit... ");
        VectorFieldEvolver evolver(dynamics, integrator);
        evolver.configuration().set_maximum_enclosure_radius(1.0);
        evolver.configuration().set_maximum_step_size(0.004);
        evolver.configuration().set_maximum_spacial_error(1e-6);
        ARIADNE_LOG_RUN_AT(1,auto orbit = evolver.orbit(initial_set, evolution_time, Semantics::UPPER));
        reach1 = orbit.reach();

        auto bb = orbit.final().bounding_box();
        auto width = bb[s].width();
        ARIADNE_LOG_PRINTLN_AT(1,"Reach size = " << orbit.reach().size());
        ARIADNE_LOG_PRINTLN_AT(1,"Final x+y+z width = " << width);
	
	    sw.click();
        ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.");

        auto instance = benchmark.create_instance("1");
        instance.set_verified(1).set_execution_time(sw.elapsed_seconds()).add_loss(orbit.reach().size()).add_loss(width.get_d());
        instance.write();
    }

    {
        RealConstant beta("beta",1000);
        RealConstant gamma("gamma",100000);

        ARIADNE_LOG_PRINTLN("Instance 2:");

        VectorField dynamics({dot(x) = -alpha*x + beta*y*z, dot(y) = alpha*x - beta*y*z - gamma*sqr(y), dot(z) = gamma*sqr(y)},{let(s)=x+y+z});

        RealVariablesBox initial_set({x==1,y==0,z==0});

        Stopwatch<Milliseconds> sw;

        ARIADNE_LOG_PRINTLN_AT(1,"Computing orbit... ");
        VectorFieldEvolver evolver(dynamics, integrator);
        evolver.configuration().set_maximum_enclosure_radius(1.0);
        evolver.configuration().set_maximum_step_size(0.004);
        evolver.configuration().set_maximum_spacial_error(1e-6);
        ARIADNE_LOG_RUN_AT(1,auto orbit = evolver.orbit(initial_set, evolution_time, Semantics::UPPER));
        reach2 = orbit.reach();

        auto bb = orbit.final().bounding_box();
        auto width = bb[s].width();
        ARIADNE_LOG_PRINTLN_AT(1,"Reach size = " << orbit.reach().size());
        ARIADNE_LOG_PRINTLN_AT(1,"Final x+y+z width = " << width);

        sw.click();
        ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.");

        auto instance = benchmark.create_instance("2");
        instance.set_verified(1).set_execution_time(sw.elapsed_seconds()).add_loss(orbit.reach().size()).add_loss(width.get_d());
        instance.write();
    }

    {
        RealConstant beta("beta",1000);
        RealConstant gamma("gamma",10000000);

        ARIADNE_LOG_PRINTLN("Instance 3:");
        
        VectorField dynamics({dot(x) = -alpha*x + beta*y*z, dot(y) = alpha*x - beta*y*z - gamma*sqr(y), dot(z) = gamma*sqr(y)},{let(s)=x+y+z});

        RealVariablesBox initial_set({x==1,y==0,z==0});

        Stopwatch<Milliseconds> sw;

        ARIADNE_LOG_PRINTLN_AT(1,"Computing orbit... ");
        VectorFieldEvolver evolver(dynamics, integrator);
        evolver.configuration().set_maximum_enclosure_radius(1.0);
        evolver.configuration().set_maximum_step_size(0.004);
        evolver.configuration().set_maximum_spacial_error(1e-6);
        ARIADNE_LOG_RUN_AT(1,auto orbit = evolver.orbit(initial_set, evolution_time, Semantics::UPPER));
        reach3 = orbit.reach();

        auto bb = orbit.final().bounding_box();
        auto width = bb[s].width();
        ARIADNE_LOG_PRINTLN_AT(1,"Reach size = " << orbit.reach().size());
        ARIADNE_LOG_PRINTLN_AT(1,"Final x+y+z width = " << width);

        sw.click();
        ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.");

        auto instance = benchmark.create_instance("3");
        instance.set_verified(1).set_execution_time(sw.elapsed_seconds()).add_loss(orbit.reach().size()).add_loss(width.get_d());
        instance.write();
    }

    ARIADNE_LOG_PRINTLN("Plotting...");
    LabelledFigure fig(Axes2d({0<=TimeVariable()<=evolution_time,0.999<=s<=1.001}));
    fig << line_style(false);
    fig << fill_colour(black);
    fig.draw(reach2);
    fig << fill_colour(white);
    fig.draw(reach1);
    fig << fill_colour(orange);
    fig.draw(reach3);
    ARIADNE_LOG_RUN_AT(2,fig.write(benchmark.name().c_str()))

    ARIADNE_LOG_PRINTLN("File " << benchmark.name() << ".png written.");
}
