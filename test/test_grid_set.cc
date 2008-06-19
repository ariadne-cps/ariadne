/***************************************************************************
 *            test_grid.cc
 *
 *
 *  Copyright  2005  Pieter Collins
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
#include <fstream>
#include <sstream>
#include <string>

#include "ariadne.h"
#include "test_float.h"
#include "combinatoric/binary_word.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
//#include "geometry/irregular_grid_set.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

template<class R>
int
test_grid_set()
{
  Grid<R> gr1(Box<R>("[0,1]x[0,1]x[0,1]"),LatticeBlock("[1,5]x[-1,3]x[-6,-2]"));
  cout << "gr1=" << gr1 << endl;

  Grid<R> gr2(Box<R>("[0,1]x[0,1]x[0,1]"),LatticeBlock("[1,4]x[-1,2]x[-6,-3]"));
  cout << "gr2=" << gr2 << endl;


  Grid<R> gr(2,0.125);
  cout << "gr=" << gr << endl;

  // Test outer-approximations
  Box<R> bx("[-0.125,0.5]x[0.1,0.3]");
  Box<R> r(bx);
  Zonotope<R> z(bx);
  cout << "\nr=" << r << "\nz=" << z << "\n" << std::endl;

  GridBlock<R> rova=over_approximation(bx,gr);

  // Test outer-approximations
  GridBlock<R> roa=outer_approximation(r,gr);
  cout << "roa=" << roa << std::endl;
  GridCellListSet<R> zoa=outer_approximation(z,gr);
  cout << "zoa=" << zoa << std::endl;

  // Test inner-approximations
  GridBlock<R> rua=inner_approximation(r,gr);
  cout << "rua=" << rua << std::endl;
  GridCellListSet<R> zua=inner_approximation(z,gr);
  cout << "zua=" << zua << std::endl;

  // test GridMaskSet clone
  {
    Grid<R> g(Vector<R>("[0.25,0.25]"));
    Box<R> bb("[-2,2]x[-2,2]");
    Box<R> r("[-1.375,0.625]x[0.5,1.375]");
    
    GridMaskSet<R> gms(g,bb);
    gms.adjoin_outer_approximation(r);
    GridMaskSet<R>* gms_clone=gms.clone();
    ARIADNE_TEST_ASSERT(gms_clone->size()==gms.size());
    ARIADNE_TEST_EVALUATE(*gms_clone);
    delete gms_clone;
  }

  return 0;
}

template<class R>
int
test_irregular_grid_set()
{
  string input("[0,0.75]x[0,0.625]  [0.625,1]x[1,1.25]  [1.25,1.5]x[1.25,2.5] ");
  stringstream is(input);

  ListSet< Box<R> > ls;
  Box< R > r;
  for (int i=0; i< 3; i++) {
    is >> r;
    cout << "r=" << r << endl;
    ls.adjoin(r);
  }
  cout << "ls=" << ls << endl;

  /*  
  IrregularGridMaskSet<R> igms(ls);
  cout << "igms=" << flush << igms << endl;

  cout << ListSet< Box<R> >(igms) << endl;

  ListSet< Box<R> > ls1,ls2;
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
  test_grid_set<Flt>();
  test_irregular_grid_set<Flt>();
  cerr << "INCOMPLETE ";
  return 0;
}
