/***************************************************************************
 *            vanderpol.cpp
 *
 *  Copyright  2017  Luca Geretti
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

#include <cstdarg>
#include "ariadne.hpp"

using namespace Ariadne;


int main()
{
    /// Create a HybridAutomton object
    AtomicHybridAutomaton vanderpol("vanderpol");

    /// Create four discrete states
    AtomicDiscreteLocation loc("loc");

    RealConstant mu("mu",1.0_dec);

    RealVariable x("x"), y("y");

    vanderpol.new_mode(loc,{ dot(x) = y, dot(y) = mu * (1 - sqr(x))*y - x});

    /// Finished building the automaton

    cout << "Automaton = " << vanderpol << endl << endl;

    typedef GeneralHybridEvolver::EnclosureType EnclosureType;
    typedef GeneralHybridEvolver::EnclosureListType EnclosureListType;
    typedef GeneralHybridEvolver::OrbitType OrbitType;

    /// Create a HybridEvolver object
    GeneralHybridEvolver evolver(vanderpol);

    /// Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);
    evolver.verbosity = 1;
    std::cout <<  evolver.configuration() << std::endl;

    Real x0(1.40);
    Real y0(2.40);
    Real eps_x0 = 15/100_q;
    Real eps_y0 = 5/100_q;

    HybridSet initial_set(vanderpol|loc,{x0-eps_x0<=x<=x0 +eps_x0,y0-eps_y0<=y<=y0+eps_y0});

    std::cout << "Initial set: " << initial_set << std::endl;
    HybridTime evolution_time(7.0,4);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    plot("vanderpol",Axes2d(-2.5,x,2.5, -3.0,y,3.0), Colour(0.0,0.5,1.0), orbit);
    //plot("vanderpol-tx",Axes2d(0,TimeVariable(),evolution_time.continuous_time(), -2.1,x,2.1), Colour(0.0,0.5,1.0), orbit);
    //plot("vanderpol-ty",Axes2d(0,TimeVariable(),evolution_time.continuous_time(), -3.0,y,3.0), Colour(0.0,0.5,1.0), orbit);
/*
    std::cout << "Discretising orbit" << std::flush;
    HybridGrid grid(vanderpol.state_auxiliary_space());
    HybridGridTreeSet hgts(grid);

    for (ListSet<HybridEnclosure>::ConstIterator it = orbit.reach().begin(); it != orbit.reach().end(); it++)
    {
        std::cout<<"."<<std::flush;
        it->state_auxiliary_set().adjoin_outer_approximation_to(hgts,4);
    }
    std::cout << "done." << std::endl;

    // The following currently fails since auxiliary variables are not tracked
    plot("vanderpol-reach", Axes2d(-2.1,x,2.1, -3.0,y,3.0), Colour(0.0,0.5,1.0), hgts);
    */
}
