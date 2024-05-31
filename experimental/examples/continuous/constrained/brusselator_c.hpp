/***************************************************************************
 *            brusselator_c.hpp
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

SystemSpecification BRU_c()
{
    RealVariable x("x"), y("y");
    VectorField dynamics({dot(x)=-y-1.5_dec*pow(x,2)-0.5_dec*pow(x,3)-0.5_dec,dot(y)=3*x-y});

    Real e1=5/100_q; Real e2=7/100_q;
    RealExpressionBoundedConstraintSet initial_set({1-e1<=x<=1+e1,1-e2<=y<=1+e2});

    Real evolution_time = 0.75_dec;


    ExactBoxType bx({{-0.816775_x,1.0529_x},{0.00576862_x,1.28336_x}});
    RealBox rbx({{-0.816775_x,1.0529_x},{0.00576862_x,1.28336_x}});
    RealBox geqz({{0.0_x,10000.0_x}});
    RealBox leqz({{-10000.0_x,0.0_x}});
    Figure fig(bx,Projection2d(2,0,1));

    auto x_ = EffectiveScalarMultivariateFunction::coordinate(2,0);
    auto y_ = EffectiveScalarMultivariateFunction::coordinate(2,1);

    BoundedConstraintSet bcs(rbx,EffectiveVectorMultivariateFunction(1,1+sqr((x_-0.0952496158743602606_x)/0.47_x)*-1+sqr((y_-0.667473609477698072_x)/0.33_x)*-1),geqz); // ELLIPSOIDAL FALSE_FOR_SOME
    //BoundedConstraintSet bcs(rbx,EffectiveVectorMultivariateFunction(1,1+sqr((x_-0.095250589076984582_x)/0.474685379404738206_x)*1+sqr((y_-0.66747455080446716_x)/0.575585201101396438_x)*-1),geqz); // HYPERBOIDAL TRUE
    //BoundedConstraintSet bcs(rbx,EffectiveVectorMultivariateFunction(1,1+sqr((x_-0.0589093517952659074_x)/0.57_x)*1+sqr((y_-0.629908768330767366_x)/0.61_x)*-1),geqz); // HYPERBOIDAL TRUE

    //BoundedConstraintSet bcs(rbx,EffectiveVectorMultivariateFunction(1,1+sqr((x_-0.13019831758000422_x)/0.389513781231366818_x)*-1+sqr((y_-0.666815728903910499_x)/0.268921262608703305_x)*-1),geqz); // ELLIPSOIDAL TRUE

    //BoundedConstraintSet bcs(rbx,EffectiveVectorMultivariateFunction(1,1+sqr((x_-0.111673162573729168_x)/0.724985614540928803_x)*1+sqr((y_-0.666513376630953669_x)/0.4_x)*-1),geqz); // HYPERBOIDAL FALSE_FOR_ALL
    //BoundedConstraintSet bcs(rbx,EffectiveVectorMultivariateFunction(1,-1+sqr((x_-0.13140766121263281_x)/0.537646778088166078_x)*1+sqr((y_-0.669168549905571175_x)/0.361273201702065284_x)*1),geqz); // ELLIPSOIDAL FALSE_FOR_ALL
    //BoundedConstraintSet bcs(rbx,EffectiveVectorMultivariateFunction(1,1+sqr((x_-0.0952498254852042204_x)/0.693399134523555616_x)*1+sqr((y_-0.666751732981645251_x)/0.476227212029361858_x)*-1),geqz); // HYPERBOIDAL FALSE_FOR_SOME


    fig.set_line_style(false);
    fig.draw(bcs);
    fig.write("constraints");

    return {"brusselator",dynamics,initial_set,evolution_time};
}
