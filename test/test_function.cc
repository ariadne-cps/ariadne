/***************************************************************************
 *            test_function.cc
 *
 *  Copyright  2009  Pieter Collins
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

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <fenv.h>

#include "function.h"
#include "predicate.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

class TestFunction 
{
  public:
    void test();
  private:
    void test_concept();
};

void TestFunction::test()
{
}

void TestFunction::test_concept()
{
    Vector<ScalarFunctionInterface> vsf(0);
    FunctionElement fe(vsf,0u);

    ScalarFunctionInterface* e1=0;
    ScalarFunctionInterface* e2=0;
    ScalarFunctionInterface* e3=0;
    //Vector<ScalarFunctionInterface> ve=join(3,*e1,*e2,*e3);
    Vector<ScalarFunctionInterface> ve=join(3,e1,e2,e3);

    shared_ptr<const VectorFunctionInterface> fp,gp;
    ComposedFunction h(gp,fp);

    shared_ptr<const VectorFunctionInterface> f1p,f2p;
    JoinedFunction jf(f1p,f2p);
    CombinedFunction cf(f1p,f2p);

    Vector< Polynomial<Float> > p;
    VectorPolynomialFunction pf(p);

    Vector<Float> b; Matrix<Float> A;
    VectorAffineFunction aff(A,b);

}

int main() {
    TestFunction().test();

    return ARIADNE_TEST_FAILURES;
}

