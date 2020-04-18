/***************************************************************************
 *            test_hybrid_graphics.cpp
 *
 *  Copyright  2011-20  Pieter Collins
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

#include "config.hpp"
#include "../test.hpp"

#include "numeric/logical.hpp"
#include "numeric/real.hpp"
#include "symbolic/variables.hpp"
#include "symbolic/assignment.hpp"

#include "hybrid/hybrid_graphics.hpp"
#include "hybrid/hybrid_set.hpp"
#include "hybrid/hybrid_paving.hpp"

using namespace Ariadne;


Int main(Int argc, char **argv) {
    RealVariable x("x"),y("y"),z("z");
    DiscreteLocation location(1);
    HybridRealBox hbx1(location,{0<=x<=1,2<=y<=3,5<=z<=7});
    HybridExactBox hbx2(location,{x,y,z},ExactBoxType{{1.,2.},{3.,4.},{6.,8.}});

    HybridFigure hfig;
    hfig.set_locations({location});
    hfig.set_bounds(x,-8,8);
    hfig.set_bounds(y,-8,8);
    hfig.set_bounds(z,-8,8);
    hfig.set_variables(x,y);

    hfig << FillColour(1,0,0) << hbx1;
    hfig << FillColour(0,1,0) << hbx2;
    hfig.write("test_hybrid_graphics");
}
