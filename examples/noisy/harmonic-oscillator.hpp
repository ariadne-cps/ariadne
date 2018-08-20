/***************************************************************************
 *            harmonic-oscillator.hpp
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


Tuple<String,DottedRealAssignments,RealVariablesBox,RealVariablesBox,Real,double> HO()
{
    RealVariable x("x"), y("y"), u("u");
    DottedRealAssignments dynamics={dot(x)=y+u,dot(y)=-x};
    RealVariablesBox inputs={-4/100_q<=u<=4/100_q};

    Real e=1/10000000_q;
    RealVariablesBox initial={{-e<=x<=e},{-e<=y<=e}};

    Real evolution_time=3.141592_dec;
    double step=1.0/64;

    return make_tuple("HO",dynamics,inputs,initial,evolution_time,step);
}
