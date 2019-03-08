/***************************************************************************
 *            bouncingball.cpp
 *
 *  Copyright  2008  Davide Bresolin
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
#include <omp.h>
#include "ariadne.hpp"

using namespace Ariadne;
using std::cout; using std::endl; using std::flush;

inline char activity_symbol(SizeType step) {
    switch (step % 4) {
    case 0: return '\\';
    case 1: return '|';
    case 2: return '/';
    default: return '-';
    }
}

void discretize(HybridGridTreePaving& hgts, GeneralHybridEvolver::OrbitType& orbit, unsigned precision)
{
  int oSize=orbit.reach().size();
  std::cerr<<"\n";
  int index=1;
  for (ListSet<HybridEnclosure>::ConstIterator it = orbit.reach().begin(); it != orbit.reach().end(); it++,index++)
  {
      std::cerr << "\r[" << activity_symbol(static_cast<unsigned>(index)) << "] " << static_cast<int>((index*100)/oSize) << "% " << std::flush;
      it->state_auxiliary_set().adjoin_outer_approximation_to(hgts,precision);
  }
  fprintf(stderr,"\n");
}

Int main(Int argc, const char* argv[])
{
    Nat evolver_verbosity=get_verbosity(argc,argv);

    typedef GeneralHybridEvolver GeneralHybridEvolverType;

    /// Set the system parameters
    Real a = 0.5_dec;  // Coefficient of restitution
    Real g = 9.8_dec;

    /// Set the position and velocity functions.
    RealVariable x("x");
    RealVariable v("v");

    /// Build the Hybrid System

    /// Create a HybridAutomton object
    HybridAutomaton ball;

    /// Create the discrete location
    //DiscreteLocation freefall(StringVariable("ball")|"freefall");
    DiscreteLocation freefall;
    cout << "location = " << freefall << endl << endl;

    /// Create the discrete events
    DiscreteEvent bounce("bounce");
    cout << "event = " << bounce << endl << endl;

    /// Build the automaton
    ball.new_mode(freefall,{dot(x)=v,dot(v)=-g});
    ball.new_guard(freefall,bounce,x<=0,EventKind::IMPACT);
    ball.new_update(freefall,bounce,freefall,{next(x)=x,next(v)=-a*v});
    /// Finished building the automaton

    cout << "Ball = " << ball << endl << endl;
    /// Compute the system evolution

    /// Create a GeneralHybridEvolver object
    GeneralHybridEvolverType evolver(ball);
    evolver.verbosity=evolver_verbosity;

    /// Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(2.0);
    evolver.configuration().set_maximum_step_size(1.0/32);
    std::cout <<  evolver.configuration() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolverType::OrbitType OrbitType;

    Real e(1.0/16);
    HybridSet initial_set(freefall,{2-e<=x<=2+e,-e<=v<=e});
    HybridTime evolution_time(1.5,4);

    std::cout << "Computing evolution... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,Semantics::LOWER);
    std::cout << "done." << std::endl;

    plot("bouncingball-xv",Axes2d(-0.1,x,2.1, -10.1,v,10.1), Colour(0.0,0.5,1.0), orbit);
    plot("bouncingball-tx",Axes2d(0.0,TimeVariable(),1.5,- 0.1,x,2.1), Colour(0.0,0.5,1.0), orbit);

    std::cout << "Discretising orbit" << std::flush;
    HybridGrid grid(ball.state_auxiliary_space());
    HybridGridTreePaving hgts(grid);
    Axes2d x_v_axes(-0.1,x,2.1, -10.1,v,10.1);


    for(unsigned i=2;i<=8;++i)
    {
      auto h = hgts;
      clock_t s_time = clock();
      // run code
      discretize(h,orbit,i);
      // End time
      clock_t e_time = clock();
      float elapsed_time =static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
      std::cout << "instance "<<i<<" in "<<elapsed_time<<" sec" << std::endl;
      char title[32];
      sprintf(title,"%d",i);
      plot(title, x_v_axes, Colour(0.0,0.5,1.0), h);
    }
    std::cout << "done." << std::endl;
}
