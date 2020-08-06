/***************************************************************************
 *            exact_real_computation.cpp
 *
 *  Copyright  2020  Pieter Collins
 *      (Based on joint work with S. Park, M. Ziegler, M. Konecny, and N. Mueller)
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
 *  You should have received a copy of the GNU G3c767e04cec413f9afb4c30b521ca71ceb5b0409eneral Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <vector>
#include <iostream>
#include <functional>
#include <initializer_list>

#include "utility/array.hpp"
#include "numeric/logical.hpp"
#include "numeric/accuracy.hpp"
#include "numeric/builtin.hpp"
#include "numeric/integer.hpp"
#include "numeric/twoexp.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/real.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/float_bounds.hpp"


using namespace Ariadne;

void trisection_sqrt() {
    Real x=3;
    auto f=[&x](Real y){return sqr(y)-x;};

    Real a=1; Real b=2;
    Dyadic tol=exp2(-24);

    while (nondeterministic_choose_index(Array<Kleenean>{tol>b-a,b-a>hlf(tol)})==1) {
        Real c=(2*a+b)/3;
        Real d=(a+2*b)/3;

        if (nondeterministic_choose_index(Array<Kleenean>{f(c)<0,f(d)>0})==0) {
            a=c;
        } else {
            b=d;
        }
        std::cerr<<"b-a="<<(b-a).compute(Accuracy(exp2(-48)))<<"\n";
    }

    Accuracy acc(exp2(-48));
//    std::cout<<"a="<<a<<", b="<<b<<"\n";
    DyadicBounds y_bnds(a.compute(acc).get().lower(),b.compute(acc).get().upper());
    std::cout<<"y_bnds="<<y_bnds<<"\n";
    std::cout<<"y_bnds.radius()="<<hlf(y_bnds.upper()-y_bnds.lower())<<"\n";
    assert(y_bnds.upper()-y_bnds.lower()<exp2(-20));
    std::cout<<"sqr(y_bnds)-x="<<y_bnds*y_bnds-3<<"\n";
}


Real compute_heron_sqrt(Accuracy acc, Real x) {
    // Algorithm from CID paper
    Natural steps=0u;
    Real y=x; Real z=x/y;
    Dyadic tol=acc.error();
    while (nondeterministic_choose_index(Array<Kleenean>{tol>y-z,y-z>hlf(tol)})==1) {
        y=(y+z)/2;
        z=x/y;
        ++steps;
    }
    std::cerr<<steps<<" steps\n";
    return y;
}


Real heron_sqrt(Real x) {
    return limit(FastCauchySequence<Real>([&x](Natural p){return compute_heron_sqrt(Accuracy(exp2(-p)),x);}));
}

Bounds<FloatMP> heron_sqrt(FloatMP x) {
    // Assume for simplicity that x>1, and tha.
    ARIADNE_ASSERT(x>1);
    FloatMP y=x;
    FloatMP z(1,x.precision());
    FloatMP e=shft(FloatMP::eps(x.precision()),3u); // Multiply by 2^3

    while(sub(up,y,z)>e) {
        y=hlf(add(up,y,z));
        z=div(down,x,y);
    }
    return Bounds<FloatMP>(y,z);
}


void test_heron_sqrt() {
    Real x=3;
    Real y=heron_sqrt(x);

    Dyadic err=exp2(-256);
    DyadicBounds y_bnds=y.compute(Accuracy(err)).get();
    std::cout << "y_bnds=" << y_bnds <<"\n";
    std::cout<<"sqr(y_bnds)-x="<<y_bnds*y_bnds-3<<"\n";
}


void try_heron_sqrt() {
    // Compute sqrt(x) by Heron's method to an accuracy of within eps
    Real x=Real(3);
    Dyadic tol=exp2(-128);
    Natural steps=0u;

    // Algorithm from CID paper
    Real y=x; Real z=x/y;
    while (nondeterministic_choose_index(Array<Kleenean>{tol>y-z,y-z>hlf(tol)})==1) {
        y=(y+z)/2;
        z=x/y;
        ++steps;
    }
    std::cout<<steps<<" steps\n";

    // Compute bounds on y to within an accuracy of 2^(-160)
    Dyadic err=exp2(-256);
    std::cout<<"y="<<y<<"\n";
    DyadicBounds y_bnds=y.compute(Accuracy(err)).get();
    std::cout<<"y_bnds="<<y_bnds<<"\n";
    std::cout<<"y_bnds.radius()="<<hlf(y_bnds.upper()-y_bnds.lower())<<"\n";
    DyadicBounds sqrt_x_bnds(y_bnds.lower()-tol,y_bnds.upper()+tol);
    std::cout<<"sqrt_x_bnds="<<sqrt_x_bnds<<"\n";

    std::cout<<"sqr(y_bnds)-x="<<y_bnds*y_bnds-3<<"\n";
    std::cout<<"sqr(sqrt_x_bnds)-x="<<sqrt_x_bnds*sqrt_x_bnds-3<<"\n";


    assert(sqrt_x_bnds.upper()-sqrt_x_bnds.lower()<2*tol+2*err);

    DyadicBounds mpfr_sqrt_x_bnds=sqrt(x).compute(Accuracy(exp2(-192)));
    assert(refines(mpfr_sqrt_x_bnds,sqrt_x_bnds));
}

int main() {

    Dyadic::set_default_writer(DecimalWriter());

    trisection_sqrt();
//    try_heron_sqrt();
//    test_heron_sqrt();

    return 0;

}
// It's entirely feasible to have higher-order constructs in the language.
// It seems impossible to do so using usual C++ syntax.
// We *need* functions, and since the arguments to a function must be *expressions*, not *statements*,
// we *need* to generate a list of "statements" in the code of a loop or conditional.
// We can't generate a list of C++ statements (without global variables) so need to generate a list of (assignment) expressions.
// This list must be separated by commas, not semicolons, and cannot declare any new variables.
// It can either go through an initializer_list, or use the sequencing operator.
// Alternatively, we can use a functional approach to compound statements.
