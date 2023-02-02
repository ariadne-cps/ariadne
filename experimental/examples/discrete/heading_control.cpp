/***************************************************************************
 *            heading_control.cpp
 *
 *  Copyright  2023  Luca Geretti
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

#include "ariadne_main.hpp"
#include "utility/stopwatch.hpp"

void ariadne_main()
{
    Real deltat=0.1_dec;
    Real v=0.5_dec;
    RealVariable x("x"), y("y"), theta("theta"), Kx("Kx"), Ky("Ky"), Kt("Kt"), b("b");
    IteratedMap heading({next(x)=x+deltat*v*cos(theta),next(y)=y+deltat*v*sin(theta),next(theta)=theta+deltat*(Kx*x+Ky*y+Kt*theta+b),
                         next(Kx)=Kx,next(Ky)=Ky,next(Kt)=Kt,next(b)=b});
    CONCLOG_PRINTLN_VAR(heading);

    // Set up the evaluators
    IteratedMapEvolver evolver(heading);

    Real d = 0.5_dec;
    // Set-up initial set and time for evolution
    RealExpressionBoundedConstraintSet initial_set={0.0_dec<=x<=0.5_dec,0.0_dec<=y<=0.5_dec,0.0_dec<=theta<=0.3926_dec,-d<=Kx<=d,-d<=Ky<=d,-d<=Kt<=d,-d<=b<=d};

    // Set up the evolution parameters and grid
    Integer evolve_time(1);

    Stopwatch<Milliseconds> sw;
    CONCLOG_PRINTLN("Computing evolution...");
    SizeType num_iterations = 552*160;

    auto function = heading.function();
    List<UpperBoxType> boxes;

    typedef BinaryWord SWord;
    typedef BinaryWord CWord;
    typedef GridTreePaving SPaving;
    typedef GridTreePaving CPaving;

    FloatDP eps(0.00390625_x,DoublePrecision());
    Grid sgrid({0.5,0.5,2*pi/8});
    ExactBoxType sdomain({{0+eps,5-eps},{0+eps,5-eps},{0+eps,2*pi-eps}});
    SPaving sdomain_paving(sgrid);
    sdomain_paving.adjoin_outer_approximation(sdomain,0);

    ExactBoxType graphics_box({{-5,10},{-5,10},{0,7}});
    Figure fig(graphics_box,0,1);
    fig << sdomain_paving;
    fig.write("sdomain_paving");

    SPaving targets_paving;
    SPaving obstacles_paving;
    Map<SWord,Map<CWord,SPaving>> forward_graph;
    Map<SWord,Map<CWord,SPaving>> backward_graph;

    auto initial_box =initial_set.euclidean_set(heading.state_space()).bounding_box();
    for (SizeType i=0; i<num_iterations; ++i) {
        boxes.push_back(apply(function,initial_box));
    }
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Using Function representation of the map: " << sw.elapsed_seconds()/num_iterations << " seconds on average.");

/*
    LabelledFigure fig(Axes2d(0<=x<=20,0<=y<=20));
    fig << orbit;
    fig.write("heading_control");*/
}