/***************************************************************************
 *            vanderpool.cc
 *
 *  Copyright  2011  Luca Geretti
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

#include "ariadne.h"

using namespace Ariadne;


int main()
{
  	HybridAutomaton system("vanderpool");

  	RealConstant mu("mu",0.5);

    DiscreteState one("one");

    // System variables
    RealVariable x("x");
    RealVariable y("y");
    List<RealVariable> varlist;
    varlist.append(x);
    varlist.append(y);

    RealExpression x_dot = y;
    RealExpression y_dot = mu*y*(1-x*x) - x;
    List<RealExpression> exprlist;
    exprlist.append(x_dot);
    exprlist.append(y_dot);
    VectorFunction dynamic(exprlist, varlist);

    system.new_mode(one,dynamic);

    /// Create a HybridEvolver object
    HybridEvolver evolver;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_cell = Vector<Float>(2,0.4);
    evolver.parameters().hybrid_maximum_step_size[one] = 0.02;
    evolver.verbosity=1;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    Box graphic_box(2, -2.0,2.0, -2.0,2.0);
    HybridTime evolution_time(5.0,1);

    Box initial_box(2, 1.01,1.02, 0.51,0.52);
    HybridEnclosureType initial_enclosure(one,initial_box);

    std::cout << "Computing forward orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(system,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    Figure g;
    g.set_bounding_box(graphic_box);
    array<uint> p(2,0,1);
    g.set_projection_map(ProjectionFunction(p,2));
    g << fill_colour(Colour(0.0,0.5,1.0));
    g << orbit;
    g.write("vanderpool-orbit-forward");

    evolver.parameters().direction = BACKWARD;

    initial_enclosure = HybridEnclosureType(one,orbit.final().bounding_box());

    std::cout << "Computing backward orbit... " << std::flush;
    orbit = evolver.orbit(system,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    Figure g2;
    g2.set_bounding_box(graphic_box);
    g2.set_projection_map(ProjectionFunction(p,2));
    g2 << fill_colour(Colour(0.0,0.5,1.0));
    g2 << orbit;
    g2.write("vanderpool-orbit-backward");

}
