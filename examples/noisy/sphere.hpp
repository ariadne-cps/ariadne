/***************************************************************************
 *            sphere.hpp
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


Tuple<String,DottedRealAssignments,RealVariablesBox,RealVariablesBox,Real,double> SP()
{
    RealVariable x("x"), y("y"), z("z"), u1("u1"), u2("u2"), u3("u3");
    RealExpression cuberadius(pow(x,3)+pow(y,3)+pow(z,3));
    DottedRealAssignments dynamics={dot(x)=pow(x,2) - u1*x*cuberadius,
                                    dot(y)=pow(y,2) - u2*y*cuberadius,
                                    dot(z)=pow(z,2) - u3*z*cuberadius};
    RealVariablesBox inputs={-1/100_q+1<=u1<=1/100_q+1,-1/100_q+1<=u2<=1/100_q+1,-1/100_q+1<=u3<=1/100_q+1};

    Real e=1/1000000000_q;
    RealVariablesBox initial={{1/4_q-e<=x<=1/4_q+e},{1/8_q-e<=y<=1/8_q+e},{1/10_q-e<=z<=1/10_q+e}};

    Real evolution_time=9;
    double step=1.0/16;

    return make_tuple("SP",dynamics,inputs,initial,evolution_time,step);
}
