/***************************************************************************
 *            vanderpol_c.cpp
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

#include "constrained.hpp"

using namespace std;
using namespace pExplore;

void ariadne_main()
{
    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

    Real x0 = 1.40_dec;
    Real y0 = 2.40_dec;
    Real eps_x0 = 0.15_dec;
    Real eps_y0 = 0.05_dec;
    RealExpressionBoundedConstraintSet initial_set({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0});

    Real evolution_time = 7;

    RealConstant ymax("ymax",2.8_x);
    RealConstant ymin("ymin",-3.0_x);
    RealConstant xmin("xmin",-1.0_x);
    RealConstant xmax("xmax",2.1_x);
    RealConstant rsqr("r^2",2.5_x);
    List<RealExpression> constraints = {y - ymin, x - xmin, ymax - y, xmax - x, sqr(x) + sqr(y) - rsqr};

    constrained_execution("vanderpol",dynamics,constraints,initial_set,evolution_time);
}
