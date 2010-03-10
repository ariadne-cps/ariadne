/***************************************************************************
 *            test_misc.cc
 *
 *  Copyright 2008-9  Pieter Collins
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

#include <iostream>

#include "test.h"

#include "function.h"
#include "user_function.h"

namespace Ariadne {

struct Radius : ScalarFunctionData<2,0> {
    template<class R, class A, class P>
    static void compute(R& r, const A& x, const P& p) {
        r = sqrt( x[0]*x[0] + x[1]*x[1] );
    }
};

struct Henon : public VectorFunctionData<2,2,2> {
    template<class R, class A, class P>
    static void compute(R& r, const A& x, const P& p) {
        r[0] = p[0] - x[0]*x[0] - p[1]*x[1];
        r[1] = x[0];
    }
};

} // namespace Ariadne

using namespace Ariadne;

int main() {

    ScalarUserFunction<Radius> g=ScalarUserFunction<Radius>();
    std::cout << g(Vector<Float>(2,3.0,4.0)) << "\n\n\n";

    VectorUserFunction<Henon> h=VectorUserFunction<Henon>(Vector<Real>(2,1.5,0.375));
    Vector<Real> x(2,3.0,4.0);
    Vector<TaylorModel> t=TaylorModel::variables(2);

    std::cout << h.evaluate(t) << "\n";

    std::cout << h(Vector<Float>(x)) << "\n";
    std::cout << h(Vector<Interval>(x)) << "\n";
    std::cout << h.jacobian(Vector<Float>(x)) << "\n";
    std::cout << h.jacobian(Vector<Interval>(x)) << "\n";

}