/***************************************************************************
 *            test_hybrid_graphics.cc
 *
 *  Copyright 2011--17  Pieter Collins
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

#include "config.h"
#include "test.h"

#include "numeric/logical.h"
#include "numeric/real.h"
#include "expression/variables.h"
#include "expression/assignment.h"

#include "hybrid/hybrid_graphics.h"
#include "hybrid/hybrid_set.h"

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
