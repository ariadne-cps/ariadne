/***************************************************************************
 *            lotka-volterra_c.hpp
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

SystemSpecification LOT_c()
{
    RealVariable x("x"), y("y");
    DottedRealAssignments dynamics={dot(x)=3*x*(1-y),dot(y)=y*(x-1)};

    RealExpressionBoundedConstraintSet initial_set={{x==1.2_dec},{y==1.1_dec}};

    Real evolution_time=10;

    return {"lotka-volterra",dynamics,initial_set,evolution_time};
}
