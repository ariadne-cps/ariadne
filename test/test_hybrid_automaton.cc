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
#include "geometry/polyhedron.h"
#include "geometry/set.h"
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
  
  Rectangle<R> r("[-7,7]x[-7,7]");
  cout << "r=" << r << endl;
  Polyhedron<R> p(r);
  cout << "p=" << p << endl;

  AffineVectorField<R> dynamic(Matrix<R>("[-2,-1;1,-2]"),Vector<R>("[-1,0]"));
  cout << "dynamic=" << dynamic << endl;
  AffineMap<R> reset(Matrix<R>("[5,0;0,5]"),Vector<R>("[0,0]"));
  cout << "reset=" << reset << endl;
  
  PolyhedralSet<R> invariant(p);
  cout << "invariant=" << invariant << endl;
  PolyhedralSet<R> activation(Polyhedron<R>(Rectangle<R>("[-2,2]x[-2,2]")));
  cout << "activation=" << activation << endl;
  
  
  DiscreteMode<R> mode(dynamic,invariant);
  DiscreteTransition<R> transition(reset,activation,mode,mode);
  
  HybridAutomaton<R> automaton("Affine test automaton");
  automaton.add_transition(reset,activation,mode,mode);

  return 0;
}
