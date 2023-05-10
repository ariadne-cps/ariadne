/***************************************************************************
 *            chemical-reactor_c.hpp
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

SystemSpecification CHE_c()
{
    RealVariable xA("xA"), xB("xB"), xC("xC"), xD("xD");
    DottedRealAssignments dynamics={dot(xA)=-30*xA*xB-0.4_dec*xA*xC+0.05_dec-0.1_dec*xA,
            dot(xB)=-30*xA*xB+0.05_dec*0.9_dec-0.1_dec*xB,
            dot(xC)=30*xA*xB-0.4_dec*xA*xC-0.1_dec*xC,
            dot(xD)=0.4_dec*xA*xC-0.1_dec*xD};

    Real e=1/1000_q;
    RealExpressionBoundedConstraintSet initial_set={{0<=xA<=e},{0<=xB<=e},{0<=xC<=e},{0<=xD<=e}};

    Real evolution_time=10;

    return {"chemical-reactor",dynamics,initial_set,evolution_time};
}
