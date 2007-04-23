/***************************************************************************
 *            test_grid.cc
 *
 *  9 May 2005
 *  Copyright  2005  Pieter Collins
 *  Email Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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
#include <sstream>
#include <string>

#include "ariadne.h"
#include "test_float.h"
#include "combinatoric/binary_word.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
//#include "geometry/irregular_grid_set.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Combinatoric;
using namespace Ariadne::Geometry;
using namespace std;

template<class R>
int
test_grid_set()
{
  Grid<R> gr1(Rectangle<R>("[0,1]x[0,1]x[0,1]"),LatticeBlock("[1,5]x[-1,3]x[-6,-2]"));
  cout << "gr1=" << gr1 << endl;

  Grid<R> gr2(Rectangle<R>("[0,1]x[0,1]x[0,1]"),LatticeBlock("[1,4]x[-1,2]x[-6,-3]"));
  cout << "gr2=" << gr2 << endl;


  Grid<R> gr(2,0.125);
  cout << "gr=" << gr << endl;

  // Test outer-approximations
  Rectangle<R> r("[-0.125,0.5]x[0.1,0.3]");
  Parallelotope<R> p(r);
  Zonotope<R> z(r);
  cout << "\n" << r << "\n" << z << "\n" << std::endl;

  GridBlock<R> rova=over_approximation(r,gr);

  // Test outer-approximations
  GridBlock<R> roa=outer_approximation(r,gr);
  cout << "roa=" << roa << std::endl;
  GridCellListSet<R> poa=outer_approximation(p,gr);
  cout << "poa=" << poa << std::endl;
  GridCellListSet<R> zoa=outer_approximation(z,gr);
  cout << "zoa=" << zoa << std::endl;

  // Test inner-approximations
  GridBlock<R> rua=inner_approximation(r,gr);
  cout << "rua=" << rua << std::endl;
  GridCellListSet<R> pua=inner_approximation(p,gr);
  cout << "pua=" << pua << std::endl;
  GridCellListSet<R> zua=inner_approximation(z,gr);
  cout << "zua=" << zua << std::endl;

  return 0;
}

template<class R>
int
test_irregular_grid_set()
{
  string input("[0,0.75]x[0,0.625]  [0.625,1]x[1,1.25]  [1.25,1.5]x[1.25,2.5] ");
  stringstream is(input);

  ListSet< Rectangle<R> > ls;
  Rectangle< R > r;
  for (int i=0; i< 3; i++) {
    is >> r;
    cout << "r=" << r << endl;
    ls.push_back(r);
  }
  cout << "ls=" << ls << endl;

  /*  
  IrregularGridMaskSet<R> igms(ls);
  cout << "igms=" << flush << igms << endl;

  cout << ListSet< Rectangle<R> >(igms) << endl;

  ListSet< Rectangle<R> > ls1,ls2;
  ls1.push_back(ls[0]);
  ls1.push_back(ls[2]);
  ls2.push_back(ls[1]);
  cout << "ls1=" << ls1 << "\n" << "ls2=" << ls2 << std::endl;
  
  IrregularGrid<R> igr1(ls1);
  cout << "Finished constructing irregular grid" << std::endl;
  cout << "igr1=" << igr1 << std::endl;

  IrregularGrid<R> igr2(ls2);
  cout << "igr2=" << igr2 << std::endl;

  IrregularGrid<R> igrj(igr1,igr2);
  cout << "igrj=" << igrj << "\n";

  IrregularGrid<R> gr1(ls1);
  cout << "gr1=" << gr1 << std::endl;

  IrregularGrid<R> gr2(ls2);
  cout << "gr2=" << gr2 << std::endl;
  */

  return 0;
}


int main() {
  test_grid_set<Float>();
  test_irregular_grid_set<Float>();
  cerr << "INCOMPLETE ";
  return 0;
}
