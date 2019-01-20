/***************************************************************************
 *            laub-loomis.hpp
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


inline Tuple<String,DottedRealAssignments,RealVariablesBox,RealVariablesBox,Real,double> LL()
{
    RealVariable x1("x1"), x2("x2"), x3("x3"), x4("x4"), x5("x5"), x6("x6"), x7("x7"), u("u");
    DottedRealAssignments dynamics={dot(x1)=1.4_dec*x3-0.9_dec*x1,
                                    dot(x2)=2.5_dec*x5-1.5_dec*x2,
                                    dot(x3)=0.6_dec*x7-0.8_dec*x2*x3,
                                    dot(x4)=2-1.3_dec*x3*x4,
                                    dot(x5)=0.7_dec*x1-x4*x5,
                                    dot(x6)=0.3_dec*x1-3.1_dec*x6,
                                    dot(x7)=1.8_dec*x6-1.5_dec*x2*x7,
                                    };
    RealVariablesBox inputs={0.0_dec<=u<=0.0_dec};

    Real W0=2/100_q;
    RealVariablesBox initial={{1.2_dec-W0<=x1<=1.2_dec+W0},{1.05_dec-W0<=x2<=1.05_dec+W0},{1.5_dec-W0<=x3<=1.5_dec+W0},
                              {2.4_dec-W0<=x4<=2.4_dec+W0},{1-W0<=x5<=1+W0},{0.1_dec-W0<=x6<=0.1_dec+W0},{0.45_dec-W0<=x7<=0.45_dec+W0}};

    Real evolution_time=10;
    double step=2.0/100;

    return make_tuple("LL",dynamics,inputs,initial,evolution_time,step);
}
