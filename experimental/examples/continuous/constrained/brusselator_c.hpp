/***************************************************************************
 *            brusselator_c.hpp
 *
 *  Copyright  2023  Luca Geretti
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

#include "ariadne.hpp"
#include "constrained.hpp"

using namespace std;
using namespace Ariadne;

SystemSpecification BRU_c()
{
    RealVariable x("x"), y("y");
    VectorField dynamics({dot(x)=-y-1.5_dec*pow(x,2)-0.5_dec*pow(x,3)-0.5_dec,dot(y)=3*x-y});

    Real e1=5/100_q; Real e2=7/100_q;
    RealExpressionBoundedConstraintSet initial_set({1-e1<=x<=1+e1,1-e2<=y<=1+e2});

    Real evolution_time = 1.0_dec;

    return {"brusselator",dynamics,initial_set,evolution_time};
}
