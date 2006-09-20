/***************************************************************************
 *            test_polyhedron.cc
 *
 *  Copyright  2005-6  Alberto Casagrande,  Pieter Collins
 *  Email casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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

#include <vector>

#include <ppl.hh>

#include "ariadne.h"
#include "base/exception.h"
#include "base/utility.h"
#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"
#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/polyhedron.h"

#include "geometry/polyhedron.tpl"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

using namespace Parma_Polyhedra_Library::IO_Operators;


int
test_ppl_polyhedron()
{
  return 0;
}  


template<typename R>
int 
test_polyhedron() 
{
  cout << "test_polyhedron<" << name<R>() << ">" << endl;
  LinearAlgebra::Matrix<R> A("[1.0,0.875;-1,1.125;0.125,-2.25]");
  LinearAlgebra::Vector<R> b("[1.375,0.5,0.25]");
  Polyhedron<R> ptp;
  
  ptp=Polyhedron<R>(A,b);
  
  return 0;
}



int main() {

  

  test_ppl_polyhedron();
  
  test_polyhedron<Float64>();

  
  cerr << "INCOMPLETE ";
}
