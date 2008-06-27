/***************************************************************************
 *            impact_oscillator.cc
 *
 *  Copyright  2008  Pieter Collins
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

#include <iostream>

#include "numeric/float.h"
#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/zonotope.h"
#include "geometry/hybrid_set.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/hybrid_evolver.h"
#include "evaluation/standard_applicator.h"
#include "evaluation/standard_integrator.h"
#include "evaluation/standard_satisfier.h"
#include "evaluation/standard_subdivider.h"
#include "evaluation/orthogonal_reducer.h"
#include "output/epsstream.h"
#include "output/logging.h"
#include "models/henon.h"

using namespace Ariadne;
using namespace std;

/*! \file impact_oscillator.cc
\brief Compute the evolution of an impact oscillator.

Consider the impact oscillator defined by
\f[ \begin{array}{c} \ddot{x} + \zeta \dot{x} + x = \cos(2\pi t/T),\quad x<d; \\[\jot] \dot{x}' = -\lambda \dot{x},\quad x=d.  \f]
This models a ball on a spring with periodic forcing hitting an obstacle.
Here, \f$\zeta\f$ is the damping, \f$d\f$ is the height of the obstacle, \f$T\f$ is the forcing period and \f$\lambda\f$ is the coefficient of restitution.

In order to model this as a hybrid automaton in Ariadne, 
we need to introduce a second discrete transition corresponding to resetting time from \f$1/\omega\f$ to zero,
since we currently only support variables in Euclidean space.

We re-write the system to have with three variables, \f$x\f$, \f$v\f$ and \f$\tau\f$, and evolution
The parameters are \f$(\zeta,T,\lambda,d)\f$.
\f[ \dot{x} = v, \quad \dot{v} = -\zeta v -x +\cos(2\pi t/T) . \f]

See Budd, "Non-Smooth Dynamical Systems and the Grazing Bifurcation", <i>Nonlinear Mathematics and its Applications</i>, Cambridge University Press, 1996.
*/


template<class R,class A, class P>
void dynamic(R& r, const A& x, const P& p) {
  const typename P::value_type twopi=6.2831853071795862;
  r[0]=x[1];
  r[1]=-p[0]*x[1]-x[0]+cos(twopi*x[2]/p[1]);
  r[2]=+1.0;
}

template<class R,class A, class P>
void reset(R& r, const A& x, const P& p) {
  r[0]=x[0];
  r[1]=-p[2]*x[1];
  r[2]=x[2];
}

template<class R,class A, class P>
void clock_reset(R& r, const A& x, const P& p) {
  r[0]=x[0];
  r[1]=x[1];
  r[2]=0.0;
}

template<class R,class A, class P>
void hit(R& r, const A& x, const P& p) {
  return p[3]-r[1];
}

template<class R,class A, class P>
void tick(R& r, const A& x, const P& p) {
  r[0]= r[2]-p[1];
}

ARIADNE_BUILD_FUNCTION(Dynamic,dynamic,3,3,4,255);
ARIADNE_BUILD_FUNCTION(Reset,reset,3,3,4,255);
ARIADNE_BUILD_FUNCTION(Tick,tick,3,3,4,255);
ARIADNE_BUILD_FUNCTION(Hit,hit,1,3,4,255);
ARIADNE_BUILD_FUNCTION(Period,period,1,3,4,255);

template<class R> 
int 
impact_oscillator()
{
  typedef Zonotope<R> ES;

  double maximum_basic_set_radius=0.25;
  double grid_length=0.125;

  EvolutionParameters<R> evolution_parameters;
  evolution_parameters.set_maximum_basic_set_radius(maximum_basic_set_radius);
  evolution_parameters.set_grid_length(grid_length);
  
  StandardApplicator<ES> applicator;
  StandardIntegrator<ES> integrator;
  StandardSatisfier<ES> satisfier;
  StandardSubdivider<ES> subdivider;
  OrthogonalReducer<ES> reducer;
  Evolver<HybridAutomaton<R>,ES> evolver(evolution_parameters,applicator,integrator,satisfier,subdivider,reducer);


  // Define the system
  double parameter_array[4]={1.0,1.0,0.75,0.125};
  Vector<R> parameters(4,parameter_array);
  HybridAutomaton<R> oscillator=HybridAutomaton<R>("Impact oscillator");
  
  DiscreteState location=DiscreteState(0);
  ConstraintSet<R> invariant(Box<R>("[-1,5]x[-1,5]x[-1,4]"));
  oscillator.new_mode(location,VectorField<R>(Dynamic<R>(parameters)),invariant);
  //oscillator.new_forced_transition(DiscreteEvent(1),DiscreteState(0),DiscreteState(0),Reset<R>(parameters),Hit<R>(parameters));
  //oscillator.new_forced_transition(DiscreteEvent(2),DiscreteState(0),DiscreteState(0),Tick<R>(parameters),Period<R>(parameters));

  // Define the initial state
  Box<R> initial_box("[1.499,1.501]x[0.499,0.501]x[0.49,0.51]"); // initial state
  Zonotope<R> initial_zonotope(initial_box);
  HybridBasicSet< Zonotope<R> > initial_set(location,initial_zonotope);

  // Give the evolution time
  Rational time=0.5;

  HybridListSet<ES> reach,evolve;
  make_lpair(reach,evolve)=evolver.reach_evolve(oscillator,initial_set,time,upper_semantics);
  cout << reach <<endl;
  cout << evolve <<endl;
  
  cout << "reach.size()=" << reach.size() << endl;
  cout << "evolve.size()=" << evolve.size() << "  " << flush;


  epsfstream eps;
  Box<R> epsbb=Box<R>("[-4.1,4.1]x[-4.1,4.1]"); // eps bounding box
  eps.open("impact_oscillator.eps",epsbb);
  eps << fill_colour(white);
  eps << line_style(false);
  eps << fill_colour(green);
  eps << fill_colour(blue);
  eps.close();

  return 0;
}

int main() {
  return impact_oscillator<Float64>();
}


