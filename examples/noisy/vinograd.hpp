/***************************************************************************
 *            vinograd.hpp
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


Tuple<String,DottedRealAssignments,RealVariablesBox,RealVariablesBox,Real,double> VG()
{
    RealVariable x("x"), y("y"), u1("u1"), u2("u2");
    DottedRealAssignments dynamics={dot(x)=pow(y,5)+pow(x,2)*(-x+y)+u1,
                                         dot(y)=pow(y,2)*(-2*x+y)+u2};
    RealVariablesBox inputs={-1/1000_q<=u1<=1/1000_q,-1/1000_q<=u2<=1/1000_q};

    Real e=1/10000000_q;
    RealVariablesBox initial={{1-e<=x<=1+e},{-e<=y<=+e}};

    Real evolution_time=18;
    double step=1.0/16;

    return make_tuple("VG",dynamics,inputs,initial,evolution_time,step);
}
