/***************************************************************************
 *            test_ltl_formula.cpp
 *
 *  Copyright  2020  Pieter Collins
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

#include <fstream>
#include <iostream>

#include "dynamics/ltl_formula.hpp"

#include "../test.hpp"

using namespace Ariadne;

class TestLTLFormula {
  public:
    Void test() const;
};

Void TestLTLFormula::test() const {
    True t;
    False f;

    AtomicProposition p("p");
    AtomicProposition q("q");

    PathFormula phi=next(p) and (q or f);
    PathFormula psi=release(eventually(phi),p);
    StateFormula spec=all(psi and t);
}

int main() {
    ARIADNE_TEST_CALL(TestLTLFormula().test());
    return ARIADNE_TEST_FAILURES;
}
