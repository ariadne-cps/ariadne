/***************************************************************************
 *            lotka-volterra.hpp
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


Tuple<String,DottedRealAssignments,RealVariablesBox,RealVariablesBox,Real,double> LV()
{
    RealVariable x("x"), y("y"), u1("u1"), u2("u2");
    DottedRealAssignments dynamics={dot(x)=u1*x*(1-y),dot(y)=u2*y*(x-1)};
    RealVariablesBox inputs={2.99_dec<=u1<=3.01_dec,0.99_dec<=u2<=1.01_dec};

    RealVariablesBox initial={{x==1.2_dec},{y==1.1_dec}};

    Real evolution_time=10;
    double step=1.0/50;

    return make_tuple("LV",dynamics,inputs,initial,evolution_time,step);
}
