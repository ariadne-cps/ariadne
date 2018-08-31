/***************************************************************************
 *            dc-dc.hpp
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


Tuple<String,DottedRealAssignments,RealVariablesBox,RealVariablesBox,Real,double> DC()
{
    Real k0(0.002987);
    Real fp0 = (-11+k0)/600;
    Real fp1 = (k0-1)/15_q;
    Real fq0 = (1-k0)/14;
    Real fq1 = -k0*20/7_q;
    Real gp0 = 1/600_q;
    Real gp1 = 1/15_q;
    Real gq0 = -1/14_q;
    Real gq1 = -20/7_q;

    RealVariable x("x"), y("y"), u1("u1"), u2("u2");
    DottedRealAssignments dynamics={dot(x)=x*fp0+y*fp1+u1*(gp0*x+gp1*y)+u2,
                                    dot(y)=x*fq0+y*fq1+u1*(gq0*x+gq1*y)};
    RealVariablesBox inputs={-2/1000_q<=u1<=2/1000_q,4/15_q<=u2<=6/15_q};

    RealVariablesBox initial={{x==1},{y==5}};

    Real evolution_time=5;
    double step=1.0/10;

    return make_tuple("DC",dynamics,inputs,initial,evolution_time,step);
}
