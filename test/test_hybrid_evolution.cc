/***************************************************************************
 *            test_hybrid_evolution.cc
 *
 *  Copyright  2006-9  Pieter Collins
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

#include <fstream>
#include <iostream>

#include "config.h"
#include "tuple.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "box.h"
#include "zonotope.h"
#include "list_set.h"
#include "integrator.h"
#include "orbit.h"
#include "hybrid_automaton.h"
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "hybrid_evolver.h"
#include "graphics_interface.h"
#include "graphics.h"
#include "hybrid_graphics.h"
#include "logging.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

int evolver_verbosity=0;


EffectiveScalarFunction c=EffectiveScalarFunction::constant(2,1);
EffectiveScalarFunction x0=EffectiveScalarFunction::coordinate(2,0);
EffectiveScalarFunction x1=EffectiveScalarFunction::coordinate(2,1);
DiscreteLocation q("q");
DiscreteEvent e("e");

Colour reach_set_colour(0.25,0.25,0.50);
Colour intermediate_set_colour(0.50,0.50,0.75);
Colour final_set_colour(0.75,0.75,1.00);
Colour initial_set_colour(0.75,0.75,1.00);
Colour guard_set_colour(0.75,0.75,0.75);

// Test evolution of realistic hybrid systems
class TestHybridEvolution
{
    int verbosity;
    mutable shared_ptr<HybridEvolverBase> evolver;
  public:
    TestHybridEvolution(const int verb);
  private:
    void _set_evolver(const HybridAutomatonInterface& system) const;
  public:
    void test() const;
    void test_bouncing_ball() const;
    void test_water_tank() const;
};

TestHybridEvolution::TestHybridEvolution(const int verb) : verbosity(verb) { }

void TestHybridEvolution::_set_evolver(const HybridAutomatonInterface& system) const
{
    evolver.reset(new GeneralHybridEvolver(system));
    //evolver->set_integrator(TaylorSeriesIntegrator(1e-5));
    evolver->set_integrator(TaylorPicardIntegrator(1e-5));
    evolver->verbosity=verbosity;
    evolver->configuration().set_maximum_step_size(1./4);
    evolver->configuration().set_maximum_enclosure_radius(1./8);
    evolver->configuration().set_maximum_enclosure_radius(1./2);
}

void TestHybridEvolution::test() const {
    ARIADNE_TEST_CALL(test_bouncing_ball());
    ARIADNE_TEST_CALL(test_water_tank());
};

void TestHybridEvolution::test_bouncing_ball() const {
    HybridAutomaton bouncing_ball;
    Real one(1);
    RealVariable x("x");
    RealVariable v("v");

    Real lambda(0.5);
    bouncing_ball.new_mode(q,(dot(x)=v,dot(v)=-one));
    bouncing_ball.new_transition(q,e,q,(next(x)=x,next(v)=-lambda*v),x<=0,impact);
    ARIADNE_TEST_PRINT(bouncing_ball);

    double height=2.0;
    double radius=1.0/64;
    HybridBox initial(q,bouncing_ball.continuous_state_space(q),Box{{height-radius,height+radius},{-radius,+radius}});
    HybridTime time(4.5,3);

    this->_set_evolver(bouncing_ball);

    Orbit<HybridEnclosure> orbit=evolver->orbit(initial,time,UPPER_SEMANTICS);
    ListSet<HybridEnclosure> const& orbit_final=orbit.final();

    //ARIADNE_TEST_PRINT(orbit);
    if(orbit_final.size()!=1u) {
        ARIADNE_TEST_WARN("orbit.final().size()="<<orbit_final.size()<<"; expected 1. "
                          "This may indicate over-zealous splitting, and/or errors in detecting the end conditions.");
    }


    Decimal exl(0.12), exu(+0.13), evl(-0.04), evu(+0.04); // Expected bounds
    HybridBox expected_orbit_final_bounding_box=HybridBox(q,{x.in(exl,exu),v.in(evl,evu)});
    for(ListSet<HybridEnclosure>::const_iterator iter=orbit_final.begin(); iter!=orbit_final.end(); ++iter) {
        const HybridEnclosure& orbit_final_set=*iter;
        ARIADNE_TEST_PRINT(orbit_final_set.bounding_box());
        ARIADNE_TEST_BINARY_PREDICATE(inside,orbit_final_set,expected_orbit_final_bounding_box);
    }
    ARIADNE_TEST_PRINT(orbit.final().size());
    ARIADNE_TEST_PRINT(expected_orbit_final_bounding_box);

    Dyadic xl(-0.5), xu(+2.5), vl(-4.0), vu(+4.0);
    Axes2d bounding_box(xl<=x<=xu,vl<=v<=vu);
    plot("test_hybrid_evolution-bouncing_ball",bounding_box,
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


void TestHybridEvolution::test_water_tank() const {
    // Declare some constants
    Real T(Rational(4));
    Real hmin(Rational(11,2));
    Real hmax(Rational(8));
    Real delta(Rational(1,20));
    Real lambda(Rational(1,50));
    Real rate(Rational(3,10));

    RealVariable height("height");
    RealVariable aperture("aperture");
    Real zero(0);
    Real one(1);

    DiscreteLocation open("open");
    DiscreteLocation opening("opening");
    DiscreteLocation closed("closed");
    DiscreteLocation closing("closing");

    DiscreteEvent start_opening("start_opening");
    DiscreteEvent start_closing("start_closing");
    DiscreteEvent finished_opening("finished_opening");
    DiscreteEvent finished_closing("finished_closing");

    // Create the tank object
    HybridAutomaton watertank;

    watertank.new_mode(open,dot((height,aperture))=(-lambda*height+rate*aperture,zero));
    watertank.new_mode(closed,dot((height,aperture))=(-lambda*height+rate*aperture,zero));
    watertank.new_mode(opening,dot((height,aperture))=(-lambda*height+rate*aperture,1/T));
    watertank.new_mode(closing,dot((height,aperture))=(-lambda*height+rate*aperture,-1/T));

    watertank.new_transition(closed,start_opening,opening,next((height,aperture))=(height,aperture),height<=hmin,urgent);
    watertank.new_transition(open,start_closing,closing,next((height,aperture))=(height,aperture),height>=hmax,urgent);
    watertank.new_transition(opening,finished_opening,open,next((height,aperture))=(height,aperture),aperture>=1,urgent);
    watertank.new_transition(closing,finished_closing,closed,next((height,aperture))=(height,aperture),aperture<=0,urgent);

    ARIADNE_TEST_PRINT(watertank);

    DiscreteLocation initial_location=opening;
    HybridBox initial_box(initial_location,(0<=height<=one/16,0<=aperture<=one/64));
    HybridSet initial(initial_location,(0<=height<=one/16,0<=aperture<=one/64));

    //HybridTime evolution_time(80.0,5);
    HybridTime evolution_time(80.0,8);

    _set_evolver(watertank);

    evolver->configuration().set_maximum_step_size(1.0);

    Orbit<HybridEnclosure> orbit = evolver->orbit(initial,evolution_time,UPPER_SEMANTICS);
    if(orbit.final().size()!=1u) {
        ARIADNE_TEST_WARN("orbit.final().size()="<<orbit.final().size()<<"; expected 1. "
                          "This may indicate over-zealous splitting, and/or errors in detecting the end conditions.");
    }
    HybridEnclosure final_enclosure=HybridEnclosure(*orbit.final().begin());
    ARIADNE_TEST_PRINT(final_enclosure.bounding_box());
    Dyadic ehl(7.6875), ehu(8.0); Decimal eal(0.999), eau(1.001); // Expected bounds
    ARIADNE_TEST_BINARY_PREDICATE(inside,final_enclosure,HybridBox(open,(height.in(ehl,ehu),aperture.in(eal,eau))));

    Decimal hl(-0.1), hu(+9.1), al(-0.3), au(+1.3);
    Axes2d bounding_box(hl<=height<=hu, al<=aperture<=au);
    plot("test_hybrid_evolution-water_tank",bounding_box,
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());

}


int main(int argc, const char* argv[])
{
    int evolver_verbosity=get_verbosity(argc,argv);

    DRAWING_METHOD = AFFINE_DRAW; DRAWING_ACCURACY = 1u;

    TestHybridEvolution(evolver_verbosity).test();
    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

