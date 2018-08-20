/***************************************************************************
 *            lorenz-attractor.hpp
 *
 *  Copyright  2008-18 Luca Geretti
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

#include "ariadne.hpp"

using namespace Ariadne;


Tuple<String,DottedRealAssignments,RealVariablesBox,RealVariablesBox,Real,double> LA()
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
