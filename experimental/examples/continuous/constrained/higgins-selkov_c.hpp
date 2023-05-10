/***************************************************************************
 *            higgins-selkov_c.hpp
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

SystemSpecification HIG_c()
{
    RealVariable S("S"), P("P");
    VectorField dynamics({dot(S)=-S*pow(P,2)+1,dot(P)=S*pow(P,2)-P});

    Real e=1/100_q;
    RealExpressionBoundedConstraintSet initial_set={{2-e<=S<=2+e},{1-e<=P<=1+e}};

    Real evolution_time=10;

    return {"higgins-selkov",dynamics,initial_set,evolution_time};
}
