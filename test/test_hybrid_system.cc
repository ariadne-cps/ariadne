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

#include "test.h"

#include "vector.h"
#include "matrix.h"

#include "hybrid_system.h"

using namespace std;
using namespace Ariadne;


class TestHybridSystem {
  public:
    void test();
  private:
    void test_build_hybrid_system();
};

void
TestHybridSystem::test()
{
    ARIADNE_TEST_CALL(test_build_hybrid_system());
}


void
TestHybridSystem::test_build_hybrid_system()
{
    // Declare the hyrbid system object
    HybridSystem system;

    // Declare the events, types and variables
    Event turn_on("turn_on");
    Event turn_off("turn_off");
    Event midnight("midnight");

    EnumeratedType swtch("Switch",(build_array,"on","off"));

    EnumeratedVariable heater("heater",swtch);
    RealVariable T("T");
    RealVariable t("t");

    // Define the continuous dynamics
    system.new_dynamic(heater=="on", dot(T)=20+cos(t));
    system.new_dynamic(heater=="off", dot(T)=10+cos(t));
    system.new_dynamic(dot(t)=1);

    // Define the nontrivial update rules
    system.new_reset(turn_off,heater=="on",next(heater)="off");
    system.new_reset(turn_on,next(heater)=DiscreteValue("on"));
    system.new_reset(midnight,next(t)=0.0);

    // Define the guard sets and invariants
    system.new_guard(turn_on,T<=16.0);
    system.new_guard(turn_off,T>=22.0);
    system.new_guard(midnight,t>=1.0);

    system.new_invariant(T>=16.0);
    system.new_invariant(T<=22.0);
    system.new_invariant(t<=1.0);

    // Define the trivial resets for nonjumping variables.
    system.new_reset(!midnight,next(t)=t);
    system.new_reset(next(T)=T);
    system.new_reset(!(turn_on,turn_off),next(heater)=heater);

    ARIADNE_TEST_PRINT(system);
}



int main() {
    TestHybridSystem().test();
    return ARIADNE_TEST_FAILURES;
}

