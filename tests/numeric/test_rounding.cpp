/***************************************************************************
 *            test_rounding.cpp
 *
 *  Copyright 2008--18 Luca Geretti
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

#include <iostream>

#include "numeric/rounding.hpp"

#include "config.hpp"
#include "../test.hpp"

namespace Ariadne { }

using namespace Ariadne;

class TestRounding
{
  public:
    TestRounding() { }
    void test() const;
};

void TestRounding::test() const {
    #if defined ARIADNE_C99_ROUNDING
        std::cout << "Using standard fenv.hpp C header file for setting the rounding mode." << std::endl;
    #elif defined ARIADNE_BOOST_ROUNDING
        #if defined BOOST_NUMERIC_INTERVAL_DETAIL_C99_ROUNDING_CONTROL_HPP
            std::cout << "Using Boost interval library standard fenv.hpp C header for setting the rounding mode." << std::endl;
        #else
            std::cout << "Using Boost interval library hardware rounding for setting the rounding mode." << std::endl;
        #endif
    #elif defined ARIADNE_GCC_ROUNDING
        std::cout << "Using ordinary GCC inline assembler for setting the rounding mode." << std::endl;
    #elif defined ARIADNE_EGCC_ROUNDING
        std::cout << "Using extended GCC inline assembler for setting the rounding mode." << std::endl;
    #elif defined ARIADNE_SSE_ROUNDING
        std::cout << "Using SSE <xmmintrin.h> header file for setting the rounding mode." << std::endl;
    #elif defined ARIADNE_MSVC_ROUNDING
        std::cout << "Using Microsoft Visual Studio inline assembler for setting the rounding mode." << std::endl;
    #else
        std::cout << "Error: no rounding mode available." << std::endl;
        ++ARIADNE_TEST_FAILURES;
    #endif
}


int main() {
   TestRounding().test();
   return ARIADNE_TEST_FAILURES;
}

