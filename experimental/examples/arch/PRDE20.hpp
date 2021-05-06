/***************************************************************************
 *            PRDE20.hpp
 *
 *  Copyright  2020  Luca Geretti
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

void PRDE20() {

    ArchBenchmark benchmark("PRDE20");

    RealVariable x("x"), y("y"), z("z");
    RealVariable a("a");

    VectorField dynamics({dot(x) = -x * y / (x + 1),
                                 dot(y) = x * y / (x + 1) - a * y,
                                 dot(z) = a * y,
                                 dot(a) = 0
                         });

    ARIADNE_LOG_PRINTLN("Production-destruction system:");

    MaximumError max_err = 1e-6;
    TaylorPicardIntegrator integrator(max_err);

    VectorFieldEvolver evolver(dynamics, integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.04);
    evolver.configuration().set_maximum_spacial_error(1e-6);

    ListSet<LabelledEnclosure> reach1, reach2, reach3;

    Real evolution_time(100);

    {
        RealVariablesBox initial_set({9.5_dec<=x<=10,y==0.01_dec,z==0.01_dec,a==0.3_dec});

        Stopwatch<Milliseconds> sw;

        ARIADNE_LOG_PRINTLN_AT(1,"Computing orbit for 'I' setup... ");
        ARIADNE_LOG_RUN_AT(1,auto orbit = evolver.orbit(initial_set, evolution_time, Semantics::UPPER));
        reach1 = orbit.reach();
        ARIADNE_LOG_PRINTLN_AT(1,"Verifying properties...");

        auto bb = orbit.final().bounding_box();
        unsigned int num_failures = 0;
        if (not possibly(bb[x] >= 0)) { ++num_failures; ARIADNE_LOG_PRINTLN_AT(1,"x is not >= 0"); }
        if (not possibly(bb[y] >= 0)) { ++num_failures; ARIADNE_LOG_PRINTLN_AT(1,"y is not >= 0"); }
        if (not possibly(bb[z] >= 0)) { ++num_failures; ARIADNE_LOG_PRINTLN_AT(1,"z is not >= 0"); }
        if (definitely(not contains(bb[x]+bb[y]+bb[z],100.0))) {
            ++num_failures;
            ARIADNE_LOG_PRINTLN_AT(1,"x+y+z does not contain 100");
        }

        auto volume = bb[x].width()*bb[y].width()*bb[z].width();
        ARIADNE_LOG_PRINTLN_AT(1,"Final volume = " << volume);
        ARIADNE_LOG_PRINTLN_AT(1,"Range for final z = " << bb[z]);
	
	    sw.click();
        ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.");

        auto instance = benchmark.create_instance("I");
        if (num_failures==0)
            instance.set_verified(1).set_execution_time(sw.elapsed_seconds()).add_loss(volume.get_d());
        instance.write();
    }

    {
        RealVariablesBox initial_set({x==10,y==0.01_dec,z==0.01_dec,0.296_dec<=a<=0.304_dec});

        Stopwatch<Milliseconds> sw;

        ARIADNE_LOG_PRINTLN_AT(1,"Computing orbit for 'P' setup... ");
        ARIADNE_LOG_RUN_AT(1,auto orbit = evolver.orbit(initial_set, evolution_time, Semantics::UPPER));
        reach2 = orbit.reach();
        ARIADNE_LOG_PRINTLN_AT(1,"Verifying properties...");

        auto bb = orbit.final().bounding_box();
        unsigned int num_failures = 0;
        if (not possibly(bb[x] >= 0)) { ++num_failures; ARIADNE_LOG_PRINTLN_AT(1,"x is not >= 0"); }
        if (not possibly(bb[y] >= 0)) { ++num_failures; ARIADNE_LOG_PRINTLN_AT(1,"y is not >= 0"); }
        if (not possibly(bb[z] >= 0)) { ++num_failures; ARIADNE_LOG_PRINTLN_AT(1,"z is not >= 0"); }
        if (definitely(not contains(bb[x]+bb[y]+bb[z],100.0))) {
            ++num_failures;
            ARIADNE_LOG_PRINTLN_AT(1,"x+y+z does not contain 100");
        }

        auto volume = bb[x].width()*bb[y].width()*bb[z].width();
        ARIADNE_LOG_PRINTLN_AT(1,"Final volume = " << volume);
        ARIADNE_LOG_PRINTLN_AT(1,"Range for final z = " << bb[z]);

        sw.click();
        ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.");

        auto instance = benchmark.create_instance("P");
        if (num_failures==0)
            instance.set_verified(1).set_execution_time(sw.elapsed_seconds()).add_loss(volume.get_d());
        instance.write();
    }

    {
        RealVariablesBox initial_set({9.7_dec<=x<=10,y==0.01_dec,z==0.01_dec,0.298_dec<=a<=0.302_dec});

        Stopwatch<Milliseconds> sw;

        ARIADNE_LOG_PRINTLN_AT(1,"Computing orbit for 'I+P' setup... ");
        ARIADNE_LOG_RUN_AT(1,auto orbit = evolver.orbit(initial_set, evolution_time, Semantics::UPPER));
        reach3 = orbit.reach();
        ARIADNE_LOG_PRINTLN_AT(1,"Verifying properties...");

        auto bb = orbit.final().bounding_box();
        unsigned int num_failures = 0;
        if (not possibly(bb[x] >= 0)) { ++num_failures; ARIADNE_LOG_PRINTLN_AT(1,"x is not >= 0"); }
        if (not possibly(bb[y] >= 0)) { ++num_failures; ARIADNE_LOG_PRINTLN_AT(1,"y is not >= 0"); }
        if (not possibly(bb[z] >= 0)) { ++num_failures; ARIADNE_LOG_PRINTLN_AT(1,"z is not >= 0"); }
        if (definitely(not contains(bb[x]+bb[y]+bb[z],100.0))) {
            ++num_failures;
            ARIADNE_LOG_PRINTLN_AT(1,"x+y+z does not contain 100");
        }

        auto volume = bb[x].width()*bb[y].width()*bb[z].width();
        ARIADNE_LOG_PRINTLN_AT(1,"Final volume = " << volume);
        ARIADNE_LOG_PRINTLN_AT(1,"Range for final z = " << bb[z]);

        sw.click();
        ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.");

        auto instance = benchmark.create_instance("IP");
        if (num_failures==0)
            instance.set_verified(1).set_execution_time(sw.elapsed_seconds()).add_loss(volume.get_d());
        instance.write();
    }

    ARIADNE_LOG_PRINTLN("Plotting...");
    LabelledFigure fig(Axes2d({0<=TimeVariable()<=100,0<=z<=11}));
    fig << line_colour(0.0,0.0,0.0);
    fig << line_style(false);
    fig << fill_colour(0.6,0.6,0.6);
    fig.draw(reach2);
    fig << fill_colour(1.0,0.0,0.0);
    fig.draw(reach3);
    fig << fill_colour(1.0,0.75,0.5);
    fig.draw(reach1);
    fig.write(benchmark.name().c_str());
    ARIADNE_LOG_PRINTLN("File " << benchmark.name() << ".png written.");
}
