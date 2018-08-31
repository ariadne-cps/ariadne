/***************************************************************************
 *            pi-controller.hpp
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


Tuple<String,DottedRealAssignments,RealVariablesBox,RealVariablesBox,Real,double> PI()
{
    RealVariable v("v"), x("x"), u("u");
    RealExpression dynv = -0.101_dec*(v-20)+1.3203_dec*(x-0.1616_dec)-0.01_dec*pow(v,2);
    DottedRealAssignments dynamics={dot(v)=dynv,dot(x)=-dynv + 3*(20-v) + u};
    RealVariablesBox inputs={-1/10_q<=u<=1/10_q};

    Real e=1/1024_q;
    RealVariablesBox initial={{5<=v<=10},{-e<=x<=+e}};

    Real evolution_time=5;
    double step=1.0/32;

    return make_tuple("PI",dynamics,inputs,initial,evolution_time,step);
}
