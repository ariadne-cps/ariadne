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

    DiscreteValue off("off");
    DiscreteValue on("on");
    DiscreteType swtch("switch",(ArrayBuilder(),on,off));
    //system.new_invariant();
}



int main() {
    TestHybridSystem().test();
    return ARIADNE_TEST_FAILURES;
}

