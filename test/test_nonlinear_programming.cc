/***************************************************************************
 *            test_nonlinear_programming.cc
 *
 *  Copyright  2010  Pieter Collins
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

#include "numeric.h"
#include "vector.h"
#include "function.h"
#include "nonlinear_programming.h"
#include "box.h"


using namespace std;
using namespace Ariadne;

class TestOptimiser
{
  private:
    scoped_ptr<OptimiserInterface> optimiser;
  public:
    TestOptimiser(const OptimiserInterface& opt)
        : optimiser(opt.clone()) { }

    void test() {
        ARIADNE_TEST_CALL(test_linear_feasibility());
    }

    void test_linear_feasibility() {
        // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
        List<ScalarFunction> x=ScalarFunction::coordinates(2);
        VectorFunction g=VectorFunction(1u, 2*x[0]+x[1]);
        ARIADNE_TEST_PRINT(g);
        Box D(2, 0.0,2.0, 0.0,2.0);
        Box C(1, -2.0,1.0);

        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        C=Box(1, 1.0,1.5);
        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        D=Box(2, 1.0,1.5,0.5,1.0);
        ARIADNE_TEST_ASSERT(!optimiser->feasible(D,g,C));
    }

};

int main(int argc, const char* argv[]) {
    NonlinearInteriorPointOptimiser nlo;
    nlo.verbosity=get_verbosity(argc,argv);
    TestOptimiser(nlo).test();

    //KrawczykOptimiser kro;
    //kro.verbosity=get_verbosity(argc,argv);
    //TestOptimiser(kro).test();
    std::cerr<<"WARNING: Interval-based Krawczyk optimiser does not work.\n";

    return ARIADNE_TEST_FAILURES;
}

