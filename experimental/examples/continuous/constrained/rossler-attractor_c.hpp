/***************************************************************************
 *            rossler-attractor_c.hpp
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

SystemSpecification ROS_c()
{
    RealVariable x("x"), y("y"), z("z");
    DottedRealAssignments dynamics={dot(x)=-y-z,dot(y)=x+y*0.1_dec,dot(z)=z*(x-6)+0.1_dec};

    Real e=1/100_q;
    RealExpressionBoundedConstraintSet initial_set={{-9-e<=x<=-9+e},{-e<=y<=e},{0.01_dec-e<=z<=0.01_dec+e}};

    Real evolution_time=12;

    return {"rossler-attractor",dynamics,initial_set,evolution_time};
}
