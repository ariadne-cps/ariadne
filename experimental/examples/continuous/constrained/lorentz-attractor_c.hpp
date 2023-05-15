/***************************************************************************
 *            lorentz-attractor_c.hpp
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

SystemSpecification LOR_c()
{
    RealVariable x("x"), y("y"), z("z");
    RealConstant sigma("sigma",10);
    RealConstant beta("beta",8/3_q);
    DottedRealAssignments dynamics={dot(x)=sigma*(y-x), dot(y)=x*(28-z)-y,dot(z)=x*y-beta*z};

    Real e=1/100_q;
    RealExpressionBoundedConstraintSet initial_set={1-e<=x<=1+e,1-e<=y<=1+e,1-e<=z<=1+e};

    Real evolution_time=1;

    return {"lorentz-attractor",dynamics,initial_set,evolution_time};
}
