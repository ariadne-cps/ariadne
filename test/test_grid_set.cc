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
#include "real_typedef.h"
#include "combinatoric/binary_word.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/grid_set.h"
#include "geometry/list_set.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

template<class R>
int
test_grid_set()
{

  ListSet<R,Rectangle> ls;
  Rectangle< R > r;
  
  string input("[0,0.75]x[0,0.625]  [0.625,1]x[1,1.25]  [1.25,1.5]x[1.25,2.5] ");
  stringstream is(input);
  for (int i=0; i< 3; i++) {
    is >> r;
    cout << "r=" << r << endl;
    ls.push_back(r);
  }
  cout << "ls=" << ls << endl;

  RegularGrid<R> gr(2,0.125);

  GridCellListSet<R> gcls(ls);
  cout << "gcls=" << gcls << endl;
  cout << "ListSet(gcls)=" << ListSet<R,Rectangle>(gcls) << endl;
  cout << "gcls.lattice_set().bounding_block()=" << gcls.lattice_set().bounding_block() << endl;


  
  GridMaskSet<R> gms(ls);
  cout << "gms=" << flush << gms << endl;

  GridMaskSet<R> gcms(gcls);
  cout << gcms << endl;
  GridCellListSet<R> gclms(gms);
  cout << gclms << endl;

  cout << ListSet<R,Rectangle>(gms) << endl;

  ListSet<R,Rectangle> ls1,ls2;
  ls1.push_back(ls[0]);
  ls1.push_back(ls[2]);
  ls2.push_back(ls[1]);
  
  IrregularGrid<R> igr1(ls1);
  IrregularGrid<R> igr2(ls2);

  IrregularGrid<R> igrj(igr1,igr2);

  cout << igr1 << "\n" << igr2 << "\n" << igrj << "\n";

  IrregularGrid<R> gr1(ls1);

  IrregularGrid<R> gr2(ls2);

  GridCellListSet<R> gcls1(ls1);
  cout << gcls1 << "\n";

  FiniteGrid<R> fgr=FiniteGrid<R>(igrj,igrj.lattice_block());
  //GridMaskSet<R> gms1(fgr,grlsj1);
  //GridMaskSet<R> gms2(fgr,grlsj2);
  //cout << regular_intersection(gms1,gms2);
  //cout << join(gms1,gms2);

  // Test over-approximations
  r=Rectangle<R>("[-0.125,0.5]x[0.1,0.3]");
  Parallelotope<R> p(r);
  Zonotope<R> z(r);
  GridBlock<R> roa=over_approximation(r,gr);
  GridCellListSet<R> poa=over_approximation(p,gr);
  GridCellListSet<R> zoa=over_approximation(z,gr);

  return 0;
}

int main() {
  test_grid_set<Real>();
  cerr << "INCOMPLETE ";
  return 0;
}
