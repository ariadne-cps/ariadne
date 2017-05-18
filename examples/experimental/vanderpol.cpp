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
    evolver.configuration().set_maximum_step_size(1e-2);
    evolver.verbosity = 1;
    std::cout <<  evolver.configuration() << std::endl;

    HybridSet initial_set(vanderpol|loc,{x==2,y==0});

    Rational max_time=10.0_q;
    HybridTime evolution_time(max_time,4);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    plot("vanderpol",Axes2d(-2.1,x,2.1, -3.0,y,3.0), Colour(0.0,0.5,1.0), orbit);
    plot("vanderpol-tx",Axes2d(0,TimeVariable(),evolution_time.continuous_time(), -2.1,x,2.1), Colour(0.0,0.5,1.0), orbit);
}
