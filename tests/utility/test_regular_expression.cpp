/***************************************************************************
 *            test_regular_expression.cpp
 *
 *  Copyright  2018  Pieter Collins
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

#include "utility/regular_expression.hpp"

using namespace Ariadne;

void test() {
    
    RegularExpression<Char> ab={'a','b'};
    RegularExpression<Char> c={'c'};
    RegularExpression<Char> ab_or_c = ab|c;
    RegularExpression<Char> r = {ab_or_c++,'a'};

    std::cout << r << "\n";

}
    
int main() { 
    test();
}
