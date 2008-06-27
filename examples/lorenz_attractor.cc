/***************************************************************************
 *            lorenz_attractor.cc
 *
 *  Copyright  2008  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
#include "numeric/approximate_float.h"
#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "geometry/rectangular_set.h"
#include "geometry/zonotope.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/standard_integrator.h"
#include "evaluation/standard_subdivider.h"
#include "evaluation/orthogonal_reducer.h"
#include "evaluation/standard_approximator.h"
#include "evaluation/reachability_analyser.h"
#include "output/epsstream.h"
#include "output/logging.h"
#include "models/lorenz.h"

using namespace Ariadne;
using namespace Ariadne::Models;
using namespace std;

template<class R> int lorenz_attractor();

int main() {
  return lorenz_attractor<Float64>();
}

template<class R> 
int 
lorenz_attractor()
{
  typedef Zonotope<R> ES;

  double maximum_basic_set_radius=0.25;
  double grid_length=0.125;
  int subdivisions=128;

  EvolutionParameters<R> evolution_parameters;
  evolution_parameters.set_maximum_basic_set_radius(maximum_basic_set_radius);
  evolution_parameters.set_grid_length(grid_length);
  
  StandardIntegrator<ES> integrator;
  StandardSubdivider<ES> subdivider;
  OrthogonalReducer<ES> reducer;
  Evolver<VectorField<R>,ES> evolver(evolution_parameters,integrator,subdivider,reducer);

  StandardApproximator<ES> approximator;

  set_applicator_verbosity(0);
  
  Discretiser< VectorField<R>, GridApproximationScheme<R>, ES > discretiser(evolution_parameters,evolver,approximator);
  ReachabilityAnalyser< VectorField<R>, GridApproximationScheme<R> > analyser(evolution_parameters,discretiser);

  double parameter_array[]={8./3.,10,28};
  Point<R> parameters=Point<R>(3,parameter_array);
  R beta=parameters[0];
  R rho=parameters[1];
  R sigma=parameters[2];


  LorenzSystem<R> lorenz=LorenzSystem<R>(parameters);
  Point<R> pt("(3,5)"); 
  smoothness_type s=3;
  cout << "pt="<<pt<<" s="<<s<<endl;
  cout << "lorenz(pt)="<<lorenz(pt)<<endl;
  cout << "lorenz.jacobian(pt)="<<lorenz.jacobian(pt) << endl; 
  cout << "lorenz.derivative(pt,s)="<<lorenz.derivative(pt,s) << endl; 

  Box<R> gbb=Box<R>("[-11.0,5.0]x[-8.0,8.0]") ;
  cout << "gbb=" << gbb << endl;
  FiniteGrid<R> fg=FiniteGrid<R>(gbb,subdivisions); // grid
  cout << "fg=" << fg << endl;
  const Grid<R>& g=fg.grid(); // grid
  Box<R> ir=Box<R>("[1.499,1.501]x[0.499,0.501]"); // initial state
  Box<R> cb=Box<R>("[-4,4]x[-4,4]"); // cutoff box
  Box<R> epsbb=Box<R>("[-4.1,4.1]x[-4.1,4.1]"); // eps bounding box
  
  cb=Box<R>(gbb); // cutoff box
  epsbb=Box<R>(gbb); // eps bounding box
  
  GridMaskSet<R> in=GridMaskSet<R>(fg);
  GridMaskSet<R> bd=GridMaskSet<R>(fg);
  in.adjoin(over_approximation(ir,g));
  bd.adjoin(over_approximation(gbb,g));

  SetInterface< Box<R> >* cr=analyser.chain_reach(lorenz,in);
  GridMaskSet<R>& gmcr=*dynamic_cast<GridMaskSet<R>*>(cr);
  PartitionTreeSet<R> ptcr=PartitionTreeSet<R>(gmcr);

  cout << gmcr <<endl;
  cout << ptcr <<endl;
  
  cout << "gmcr.size()=" << gmcr.size() << endl;
  cout << "ptcr.size()=" << ptcr.size() << "  " << flush;
  cout << "ptcr.capacity()=" << ptcr.capacity() << endl;

  RectangularSet<R> ins(ir);
  RectangularSet<R> cbs(cb);
  SetInterface< Box<R> >* crs=analyser.chain_reach(lorenz,ins);
  cout << "*crs=" << *crs <<endl;
  GridMaskSet<R>& gmcrs=*dynamic_cast<GridMaskSet<R>*>(crs);
  cout << "gmcrs=" << gmcrs <<endl;
  cout << "gmcrs.size()=" << gmcrs.size() << " of " << gmcrs.capacity() << endl;

  epsfstream eps;
  eps.open("henon_chainreach-1.eps",epsbb);
  eps << fill_colour(white) << cb;
  eps << line_style(false);
  eps << fill_colour(green) << gmcr;
  eps << fill_colour(blue) << ir;
  eps << line_style(true);
  eps << ptcr.partition_tree();
  eps.close();

  eps.open("henon_chainreach-2.eps",epsbb);
  eps << line_style(false);
  eps << fill_colour(red) << difference(gmcr.neighbourhood(),gmcr);
  eps << fill_colour(blue) << gmcr.adjoining();
  eps << fill_colour(green) << gmcr;
  eps.close();

  epsbb=Box<R>("[-4.1,4.1]x[-4.1,4.1]"); // eps bounding box
  eps.open("henon_chainreach-3.eps",epsbb);
  eps << line_style(false);
  eps << fill_colour(green) << gmcr;
  eps.close();

  return 0;
}
