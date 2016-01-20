/***************************************************************************
 *            test_formula.cc
 *
 *  Copyright  2010-4  Pieter Collins
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

#include "config.h"

#include "expression/formula.h"
#include "numeric/numeric.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

class TestFormula
{
    //static ApproximateFormula o;;
    //static ApproximateFormula x;
    //static ApproximateFormula y;
  public:
    TestFormula();
    void test();
  private:
    void test_constructors();
    void test_pointer_constructor();
};


TestFormula::TestFormula()
{
}

void TestFormula::test()
{
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_pointer_constructor());
}

void TestFormula::test_constructors()
{
    ApproximateFormula o(ApproximateFormula::constant(1.0));
    ApproximateFormula x(ApproximateFormula::coordinate(0));
    ApproximateFormula y(ApproximateFormula::coordinate(1));

    ARIADNE_TEST_CONSTRUCT( ApproximateFormula, r, (sqrt(pow(x,2)+pow(y,2))) );
    ARIADNE_TEST_EQUALS(evaluate(r,Vector<FloatApproximation>{6,8}),10);
}

void TestFormula::test_pointer_constructor()
{
    unsigned int zero=0u;
    ARIADNE_TEST_CONSTRUCT(ApproximateFormula, ncz, (zero));
    ARIADNE_TEST_ASSERT(ncz.node_ptr()!=nullptr);

    ARIADNE_TEST_WARN("Constructor Formula<X>(0) yields nullptr");
}

int main() {
    TestFormula().test();
    return ARIADNE_TEST_FAILURES;
}

