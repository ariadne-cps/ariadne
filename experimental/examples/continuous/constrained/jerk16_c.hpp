/***************************************************************************
 *            jerk16_c.hpp
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

SystemSpecification J16_c()
{
    RealVariable x("x"), y("y"), z("z");

    DottedRealAssignments dynamics={dot(x)=y,dot(y)=z,dot(z)=-y+pow(x,2)-0.03_dec};

    Real e=1/1024_q;
    RealExpressionBoundedConstraintSet initial_set={{-e<=x<=e},{-e<=y<=e},{-e<=z<=e}};

    Real evolution_time=10;

    return {"jerk16",dynamics,initial_set,evolution_time};
}
