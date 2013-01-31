#include <iostream>

#include "include/c1_taylor_function.h"

using namespace Ariadne;

int main() {
    C1TaylorSeries f1=C1TaylorSeries::constant(2.0);
    f1._coefficients={2,3,5,7};
    f1._uniform_error=0.1;
    std::cout << "f1="<<f1<<"\n";

    C1TaylorSeries f2=C1TaylorSeries::coordinate();
    f2._coefficients={1,2,3,4,5};
    f2._uniform_error=0.01;
    std::cout << "f2="<<f2<<"\n";

    C1TaylorSeries f1pf2=f1+f2;
    std::cout << "f1pf2="<<f1pf2<<"\n";
    C1TaylorSeries f1tf2=f1*f2;
    std::cout << "f1tf2="<<f1tf2<<"\n";
    std::cout << "Done\n";
    //std::cout << "f2*f2="<<f2*f2<<"\n";
}