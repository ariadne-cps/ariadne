/***************************************************************************
 *            test_discrete_automaton.cpp
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
#include "dynamics/discrete_automaton.hpp"

using namespace Ariadne;

using Char = char;

void test() {

    Map<Pair<Char,Char>,Set<Char>> trns({ {{'u','0'},{'0','1'}}, {{'u','1'},{'0'}} });
    Map<Char,Char> outp( {{'0','0'},{'1','1'}} );
    Set<Char> init({'0'});
    FiniteNondeterministicAutomaton<Char,Char,Char> a(trns,outp,init);
    std::cout << a.trajectories(Word<Char>({'u','u','u','u','u','u','u'})) << "\n";

    {
        Char a='a', b='b', c='c', d='d', y0='0', y1='1';   
        Set<Char> s0 = {a};
        Map<Pair<Char,Char>,RegularExpression<Char>> s = { {{a,a},y0},{{a,b},y0},{{a,c},y0} , {{b,c},y1},{{c,d},y1},{{d,b},y1},{{d,a},y1} };
        std::cerr<<"s="<<s<<"\n";
        s=reduce(s,s0);
        std::cerr<<"s="<<s<<"\n";
        s=reduce(s,s0);
        std::cerr<<"s="<<s<<"\n";
        s=reduce(s,s0);
        std::cerr<<"s="<<s<<"\n";
    }
}
    
int main() { 
    test();
}
