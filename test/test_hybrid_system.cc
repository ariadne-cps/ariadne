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

#include "discrete_automaton.h"
#include "formula.h"
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
    HybridSystem system;
    DiscreteEvent turn_on("turn_on");
    DiscreteEvent turn_off("turn_off");
    DiscreteEvent midnight("midnight");
    DiscreteType swtch("Switch",(build_array,"on","off"));
    DiscreteVariable heater("heater",swtch);
    RealVariable T("Te");
    RealVariable t("t");
    RealVariable c("c");

    system.new_dynamic(heater=="on", dot(T)=20+cos(t));
    system.new_dynamic(heater=="off", dot(T)=10+cos(t));
    system.new_equation(true, c=t+0.5);
    system.new_reset(turn_off,heater=="on",next(heater),DiscreteValue("off"));
    system.new_reset(turn_on,heater=="off",next(heater),DiscreteValue("on"));
    system.new_reset(midnight,DiscretePredicate(true),next(t)=RealFormula(0.0));
    system.new_guard(midnight,DiscretePredicate(true),t>RealFormula(1.0));
    system.new_invariant(DiscretePredicate(true),t<=RealFormula(1.0));

    ARIADNE_TEST_PRINT(system);

    //system.new_invariant();
}



int main() {
    TestHybridSystem().test();
    return ARIADNE_TEST_FAILURES;
}

