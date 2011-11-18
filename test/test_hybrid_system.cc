/***************************************************************************
 *            test_hybrid_system.cc
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

#include <iostream>
#include <fstream>

#include "config.h"
#include "test.h"

#include "vector.h"
#include "matrix.h"

#include "function.h"
#include "hybrid_system.h"

using namespace std;
using namespace Ariadne;


class TestHybridSystem {
  public:
    void test();
  private:
    void test_build_hybrid_system();
    void test_static_analysis();
  private:
    HybridSystem _system;
};

void
TestHybridSystem::test()
{
    ARIADNE_TEST_CALL(test_build_hybrid_system());
    ARIADNE_TEST_CALL(test_static_analysis());
}


void
TestHybridSystem::test_build_hybrid_system()
{
    // Declare the hyrbid system object
    HybridSystem& system=this->_system;

    // Declare the events, types and variables
    Event turn_on("turn_on");
    Event turn_off("turn_off");
    Event midnight("midnight");

    EnumeratedType swtch("Switch",(build_array,"on","off"));

    //EnumeratedVariable heater("heater",swtch);
    StringVariable heater("heater");
    RealVariable T("T");
    RealVariable t("t");
    RealVariable u("u");
    RealVariable y("y");

    // Define the continuous dynamics
    system.new_dynamic(heater=="on", dot(T)=20+cos(t)+u);
    system.new_dynamic(heater=="off", dot(T)=10+cos(t));
    system.new_dynamic(dot(t)=1);

    //system.new_equation(y=RealFormula(T));

    // Define the nontrivial update rules
    system.new_transition(turn_off,heater=="on",next(heater)="off");
    system.new_transition(turn_on,next(heater)=EnumeratedValue("on"));
    system.new_reset(midnight,next(t)=0.0);

    // Define the guard sets and invariants
    system.new_guard(turn_on,heater=="on",false);
    system.new_guard(turn_on,heater=="off",T<=16.0);
    system.new_guard(turn_off,heater=="off",false);
    system.new_guard(turn_off,heater=="on",T>=22.0);
    system.new_guard(midnight,t>=1.0);

    system.new_guard(!(turn_on,turn_off,midnight),false);

    system.new_invariant(heater=="off",T>=16.0);
    system.new_invariant(heater=="on",T<=22.0);
    system.new_invariant(t<=1.0);

    // Define the trivial resets for nonjumping variables.
    system.new_reset(!midnight,next(t)=t);
    system.new_reset(next(T)=T);
    system.new_transition(!(turn_on,turn_off),next(heater)=heater);

    ARIADNE_TEST_PRINT(system);


    // Temporary test code
    ARIADNE_TEST_PRINT(_system.events());
    ARIADNE_TEST_PRINT(_system.discrete_variables());

    DiscreteValuation location;

    //location.set(heater,EnumeratedValue("on"));
    location.set(heater,"on");
    ARIADNE_TEST_PRINT(location);
    ARIADNE_TEST_PRINT(_system.state_variables(location));
    ARIADNE_TEST_PRINT(_system.auxiliary_variables(location));
    ARIADNE_TEST_PRINT(_system.input_variables(location));

    location.set(heater,"off");
    ARIADNE_TEST_PRINT(location);
    ARIADNE_TEST_PRINT(_system.state_variables(location));

    RealVariable xa("xa"), xb("xb"), xc("xc"), xd("xd"), xe("xe"), xf("xf");
    Event e1("e1"), e2("e2"), e3("e3");
    StringVariable q1("q1");
    StringVariable q2("q2");
    HybridSystem sys;
    sys.new_dynamic(dot(xa)=xc-xb+xd);
    sys.new_equation(xc=xb+xd);
    sys.new_equation(xd=xb+xe);
    sys.new_dynamic(dot(xb)=xe+xb+xc);
    sys.new_equation(xe=xa+xb);
    DiscreteValuation loc;
    loc.set(q1,"on");
    loc.set(q2,"off");

    sys.new_transition(EventSet(e1),next(q1)=q2);
    sys.new_reset(e1,next(xa)=xc);
    sys.new_reset(e1,next(xb)=xa+xe);

    sys.new_guard(e1,xa<=0);
    sys.new_guard(e1,xc<=0);
    sys.new_guard(e2,xb<=0);
    sys.new_guard(!((e1,e2)),false);
    //sys.dynamic(loc);

    DiscreteValuation src=loc;
    DiscreteValuation trg=sys.target(e1,src);

    ARIADNE_TEST_PRINT(sys.check_dynamic(loc));
    ARIADNE_TEST_PRINT(sys.check_guards(src));
    ARIADNE_TEST_PRINT(sys.check_reset(e1,src,trg));

    ARIADNE_TEST_PRINT(sys.target(e1,loc));
    ARIADNE_TEST_PRINT(sys.dynamic(loc));
    ARIADNE_TEST_PRINT(sys.equations(loc));
    ARIADNE_TEST_PRINT(sys.guards(loc));
}



void
TestHybridSystem::test_static_analysis()
{
    ARIADNE_TEST_PRINT(_system.events());
    ARIADNE_TEST_PRINT(_system.discrete_variables());

    //ARIADNE_TEST_PRINT(_system.state_variables(location));
}


int main() {
    TestHybridSystem().test();
    return ARIADNE_TEST_FAILURES;
}

