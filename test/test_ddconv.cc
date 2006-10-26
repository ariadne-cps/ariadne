/***************************************************************************
 *            test_ddconv.cc
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

#include "ariadne.h"
#include "numeric/rational.h"
#include "numeric/float64.h"
#include "geometry/ddconv.h"
#include "geometry/ddconv.tpl"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace std;

template<class R> class nt {
  public: typedef R at;
};
  
template<> class nt<Float64> {
  public: typedef Interval<Float64> at;
};
  


template<class R>
int 
test_ddconv() 
{
  cout << "test_ddconv<" << name<R>() << ">" << endl;
  typedef typename nt<R>::at F;
  
  std::vector< LinearAlgebra::Vector<F> > constraints;
  std::vector< LinearAlgebra::Vector<F> > generators;
  std::vector< LinearAlgebra::Vector<F> > new_constraints;
  std::vector< LinearAlgebra::Vector<F> > new_generators;
  constraints.push_back(Vector<F>(Vector<R>("[1,0,-1]")));
  constraints.push_back(Vector<F>(Vector<R>("[2,-1,2]")));
  constraints.push_back(Vector<F>(Vector<R>("[-1,1,0]")));
  constraints.push_back(Vector<F>(Vector<R>("[-1,-1,8]")));

  cout << "constraints=" << constraints << endl;
  ddconv(generators,constraints);
  cout << "generators=" << generators << endl;
  ddconv(new_constraints,generators);
  cout << "new_constraints=" << new_constraints << endl;
  ddconv(new_generators,new_constraints);
  cout << "new_generators=" << generators << endl;
  return 0;
}



int main() {

  test_ddconv<Rational>();
  test_ddconv< Float64 >();

  cerr << "INCOMPLETE ";
}
