/***************************************************************************
 *            cvdp_c.hpp
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

SystemSpecification COU_c()
{
    RealConstant mu("mu",1);
    RealVariable x1("x1"), y1("y1"), x2("x2"), y2("y2");

    VectorField dynamics({dot(x1)=y1, dot(y1)=mu*(1-sqr(x1))*y1+x2-2*x1, dot(x2)=y2, dot(y2)=mu*(1-sqr(x2))*y2+x1-2*x2});

    RealExpressionBoundedConstraintSet initial_set({1.25_dec<=x1<=1.55_dec,2.35_dec<=y1<=2.45_dec,1.25_dec<=x2<=1.55_dec,2.35_dec<=y2<=2.45_dec});

    Real evolution_time = 7;

    return {"coupled-vdp",dynamics,initial_set,evolution_time};
}
