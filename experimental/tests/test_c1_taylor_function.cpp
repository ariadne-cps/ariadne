/***************************************************************************
 *            test_c1_taylor_function.cpp
 *
 *  Copyright  2008--17  Pieter Collins
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

#include "include/c1_taylor_function.hpp"

using namespace Ariadne;

Int test_taylor_series() {
    C1TaylorSeries f1=C1TaylorSeries::constant(2.0);
    f1._coefficients={2,3,5,7};
    f1._uniform_error=0.1;
    std::cout << "f1="<<f1<<"\n";

    C1TaylorSeries f2=C1TaylorSeries::coordinate();
    f2._coefficients={1,2,3,4,5};
    f2._uniform_error=0.01;
    std::cout << "f2="<<f2<<"\n";

    const FloatDPValue c1=3.0;
    const ExactIntervalType c2=ExactIntervalType(1.0)/ExactIntervalType(3.0);
    C1TaylorSeries f1pc1=f1; f1pc1+=c1;
    std::cout << "f1pc1="<<f1pc1<<"\n";
    C1TaylorSeries f1tc2=f1; f1tc2*=c2;
    std::cout << "f1tc2="<<f1tc2<<"\n";


    C1TaylorSeries f1pf2=f1+f2;
    std::cout << "f1pf2="<<f1pf2<<"\n";
    C1TaylorSeries f1tf2=f1*f2;
    std::cout << "f1tf2="<<f1tf2<<"\n";
    std::cout << "\n\n";
    //std::cout << "f2*f2="<<f2*f2<<"\n";

    ExactIntervalType x=2.0;
    std::cout << "f1(x)="<<evaluate(f1,x)<<"\n";
   
}

Int test_taylor_function() {
    C1TaylorFunction f1=C1TaylorFunction(2);
    f1._expansion={{{1,1},7},{{1,0},5},{{0,1},3},{{0,0},2}};
    f1._uniform_error=0.0;
    std::cout << "f1="<<f1<<"\n";
    std::cout << "(f1+=13)="<<(f1+=13)<<"\n";
    std::cout << "(f1+=1/5)="<<(f1+=FloatDPValue(0.2))<<"\n";
    std::cout << "(f1*=1/5)="<<(f1*=FloatDPValue(0.2))<<"\n";

    C1TaylorFunction f2=C1TaylorFunction::coordinate(2,1);
    f2._uniform_error=0.01;
    std::cout << "f2="<<f2<<"\n";

    C1TaylorFunction f3=C1TaylorFunction::constant(2,3.0);
    std::cout << "f3="<<f3<<"\n\n";

    f1._expansion={{{1,1},7},{{1,0},5},{{0,1},3},{{0,0},2}};
    f2._expansion={{{2,0},1.3},{{1,1},1.1},{{1,0},0.7},{{0,2},0.5},{{0,1},0.3},{{0,0},0.2}};
    f1._expansion.sort(ReverseLexicographicKeyLess());
    f2._expansion.sort(ReverseLexicographicKeyLess());
    std::cout << "f1="<<f1<<"\n";
    std::cout << "f2="<<f2<<"\n";
    C1TaylorFunction f1pf2=f1+f2;
    std::cout << "f1pf2="<<f1pf2<<"\n";
    C1TaylorFunction f1tf2=f1*f2;
    std::cout << "f1tf2="<<f1tf2<<"\n";
    std::cout << "\n\n";
    //std::cout << "f2*f2="<<f2*f2<<"\n";

    Vector<ExactIntervalType> x={2,3};
    std::cout << "x="<<x<<"\n";
    std::cout << "f1(x)="<<evaluate(f1,x)<<"\n";

    C1TaylorSeries s;
    s._coefficients={1.,1./2,1./3,1./4};
    std::cout << "s="<<s<<"\n";
    std::cout << "f="<<f1<<"\n";
    std::cout << "compose(s,f)="<<compose(s,f1)<<"\n";
    compose(f1,{f1,f2});
}


Int main() {
    test_taylor_series();
    test_taylor_function();
    std::cout << "Done\n";
}
