/***************************************************************************
 *            test_hybrid_automaton.cc
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
#include "geometry/set.h"
#include "geometry/arbitrary_set.h"
#include "geometry/hybrid_set.h"
#include "geometry/rectangle.h"
#include "geometry/polyhedron.h"
#include "geometry/polyhedral_set.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/hybrid_automaton.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace std;

template<class R> int test_hybrid_automaton();
  
int main() {
  return test_hybrid_automaton<Real>();
}

template<class R>
int test_hybrid_automaton() 
{
  
  Rectangle<R> r("[-1,1]x[-1,1]");
  cout << "r=" << r << endl;

  AffineVectorField<R> dynamic(Matrix<R>("[-0.25,-1.00;1.00,-0.25]"),Vector<R>("[0.00,0.00]"));
  cout << "dynamic=" << dynamic << endl;
  AffineMap<R> reset(Matrix<R>("[-7,0;0,-7]"),Vector<R>("[0,0]"));
  cout << "reset=" << reset << endl;
  
  PolyhedralSet<R> invariant(r);
  cout << "invariant=" << invariant << endl;
  PolyhedralSet<R> activation12(Rectangle<R>("[-0.20,0.00]x[-0.20,0.00]"));
  PolyhedralSet<R> activation21(Rectangle<R>("[0.00,0.20]x[0.00,0.20]"));
  cout << "activation12=" << activation12 << endl;
  cout << "activation21=" << activation21 << endl;
  cout << endl;
  
  HybridAutomaton<R> automaton("Affine test automaton");
  id_type mode1_id=0;
  id_type mode2_id=1;
  const DiscreteMode<R>& mode1=automaton.new_mode(mode1_id,dynamic,invariant);
  const DiscreteMode<R>& mode2=automaton.new_mode(mode2_id,dynamic,invariant);
  id_type event_id=5;
  const DiscreteTransition<R>& transition12=automaton.new_transition(event_id,mode1_id,mode2_id,reset,activation12);
  const DiscreteTransition<R>& transition12=automaton.new_transition(event_id,mode1_id,mode2_id,reset,activation21);
  
  cout << mode1  <<  "\n" << mode2 << "\n" << transition12 << "\n" << transition21 << endl;
  cout << automaton.invariant() << endl;

  return 0;
}
