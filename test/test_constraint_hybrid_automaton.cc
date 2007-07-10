/***************************************************************************
 *            test_constraint_hybrid_automaton.cc
 *
 *  Copyright  2007  Pieter Collins
 *  Email  Pieter.Collins@cwi.nl
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
#include "test_float.h"
#include "geometry/set_interface.h"
#include "geometry/set_reference.h"
#include "geometry/hybrid_set.h"
#include "geometry/rectangle.h"
#include "geometry/polyhedron.h"
#include "geometry/polyhedral_set.h"
#include "geometry/constraint.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/constraint_hybrid_automaton.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace std;

template<class R> int test_constraint_hybrid_automaton();
  
int main() {
  return test_constraint_hybrid_automaton<Float>();
}

template<class R>
int test_constraint_hybrid_automaton() 
{
  Matrix<R> A("[2,1;1,1]");
  cout << A << "\n\n";


  AffineVectorField<R> dynamic(Matrix<R>("[-0.25,-1.00;1.00,-0.25]"),Vector<R>("[0.00,0.00]"));
  cout << "dynamic=" << dynamic << endl;
  AffineMap<R> reset(Matrix<R>("[-7,0;0,-7]"),Vector<R>("[0,0]"));
  cout << "reset=" << reset << endl;
  
  LinearConstraint<R> constraint1(Vector<R>("[1,0]"),Geometry::greater,R(0));
  LinearConstraint<R> constraint2(Vector<R>("[1,0]"),Geometry::less,R(2));
  LinearConstraint<R> activation(Vector<R>("[1,0]"),Geometry::greater,R(1));
  LinearConstraint<R> guard(Vector<R>("[0,1]"),Geometry::greater,R(1));

  ConstraintHybridAutomaton<R> automaton("Affine test automaton");
  id_type mode1_id=0;
  id_type mode2_id=1;
  const ConstraintDiscreteMode<R>& mode1=automaton.new_mode(mode1_id,dynamic,constraint1);
  const ConstraintDiscreteMode<R>& mode2=automaton.new_mode(mode2_id,dynamic,constraint2);
  id_type event1_id=2;
  id_type event2_id=3;
  const ConstraintDiscreteTransition<R>& transition1=automaton.new_transition(event1_id,mode1_id,mode1_id,reset,activation);
  const ConstraintDiscreteTransition<R>& transition2=automaton.new_forced_transition(event2_id,mode1_id,mode1_id,reset,guard);
  
  cout << mode1  <<  "\n" << mode2 << "\n" << transition1 << "\n" << transition2 << endl;
  cout << automaton << endl;

  return 0;
}
