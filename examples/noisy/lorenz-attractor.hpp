/***************************************************************************
 *            lorenz-attractor.hpp
 *
 *  Copyright  2008-18 Luca Geretti
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

using namespace Ariadne;


inline Tuple<String,DottedRealAssignments,RealVariablesBox,RealVariablesBox,Real,double> LA()
{
    RealVariable x("x"), y("y"), z("z"), u("u");
    DottedRealAssignments dynamics={dot(x)=10*(y-x),
                                         dot(y)=x*(28 - z) - y + x*u,
                                         dot(z)=x*y - z*8/3_q};
    RealVariablesBox inputs={-1/100_q<=u<=1/100_q};

    Real e=1/1024_q;
    RealVariablesBox initial={{1-e<=x<=1+e},{1-e<=y<=1+e},{1-e<=z<=1+e}};

    Real evolution_time=1;
    double step=1.0/256;

    return make_tuple("LA",dynamics,inputs,initial,evolution_time,step);
}
