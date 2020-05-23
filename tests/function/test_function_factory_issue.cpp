/***************************************************************************
 *            test_function_factory_issue.cpp
 *
 *  Copyright  2018-20  Luca Geretti
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
#include <fstream>
#include <sstream>
#include <string>

#include "config.hpp"

#include "function/taylor_function.hpp"
#include "algebra/algebra.hpp"

#include "../test.hpp"

using namespace Ariadne;

Int main(Int argc, const char* argv[]) {

    TaylorFunctionFactory factory = TaylorFunctionFactory(Sweeper<FloatDP>());

    ExactBoxType mydom({{0,1},{2,3},{4,5}});
    UpperBoxType mybx({{0,1},{2,3}});
    std::cout << factory.create_constants(mydom,cast_singleton(mybx)) << std::endl;
    std::cout << factory.create_constants(mydom,cast_singleton(mybx)) << std::endl;
    std::cout << "done." << std::endl;

    return 0;
}
