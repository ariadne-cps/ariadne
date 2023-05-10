/***************************************************************************
 *            pi-controller_c.hpp
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

SystemSpecification PIC_c()
{
    RealVariable v("v"), x("x");
    RealExpression dynv = -0.101_dec*(v-20)+1.3203_dec*(x-0.1616_dec)-0.01_dec*pow(v,2);
    VectorField dynamics({dot(v)=dynv,dot(x)=-dynv + 3*(20-v)});

    Real e=1/1024_q;
    RealExpressionBoundedConstraintSet initial_set={{5<=v<=10},{-e<=x<=+e}};

    Real evolution_time=5;

    return {"pi-controller",dynamics,initial_set,evolution_time};
}
