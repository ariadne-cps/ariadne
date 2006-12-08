/***************************************************************************
 *            test_hybrid_evolution.cc
 *
 *  Copyright  2006  Alberto Casagrande,  Pieter Collins
 *  Email  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
#include <fstream>
#include <string>

#include "ariadne.h"
#include "real_typedef.h"
#include "geometry/rectangle.h"
#include "geometry/hybrid_set.h"
#include "geometry/set.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/hybrid_automaton.h"
#include "evaluation/applicator.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/hybrid_evolver.h"
#include "output/epsfstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_hybrid_evolution();
  
int main() {
  test_hybrid_evolution<Real>();
  cerr << "INCOMPLETE ";
  return 0;
}

template<class R>
int test_hybrid_evolution() 
{
  
  Rectangle<R> r("[-7,7]x[-7,7]");
  Polyhedron<R> p(r);

  AffineVectorField<R> dynamic(Matrix<R>("[-2,-1;1,-2]"),Vector<R>("[-1,0]"));
  AffineMap<R> reset(Matrix<R>("[5,0;0,5]"),Vector<R>("[0,0]"));
  
  PolyhedralSet<R> invariant(p);
  PolyhedralSet<R> activation(Polyhedron<R>(Rectangle<R>("[-2,2]x[-2,2]")));
  
  HybridAutomaton<R> automaton("Affine automaton");
  DiscreteMode<R>& mode=automaton.new_mode(dynamic,invariant);
  DiscreteTransition<R>& transition=automaton.new_transition(reset,activation,mode,mode);
  
  cout << mode << endl << transition << endl;
  
  Applicator<R> apply;
  LohnerIntegrator<R> lohner(0.1,0.1,0.1);
  HybridEvolver<R> hybrid_evolver(apply,lohner);
  
  Rectangle<R> bounding_box("[-8,8]x[-8,8]");
  FiniteGrid<R> finite_grid(bounding_box,256);
  const Grid<R>& grid=finite_grid.grid();
  
  Rectangle<R> initial_rectangle("[-7,-7]x[-7,-7]");
  HybridGridMaskSet<R> initial_set(1,finite_grid);
  initial_set[0].adjoin(over_approximation(initial_rectangle,grid));
  
  HybridGridMaskSet<R> bounding_set(1,finite_grid);
  bounding_set[0].adjoin(over_approximation(bounding_box,grid));
  
  HybridGridMaskSet<R> chainreach=hybrid_evolver.chainreach(automaton,initial_set,bounding_set);
  
  epsfstream eps("test_hybrid_evolution.eps",bounding_box);
  eps << chainreach[0];
  eps.close();
  
  
  return 0;
}
