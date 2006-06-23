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
#include <string>

#include "ariadne.h"
#include "real_typedef.h"
#include "base/exception.h"
#include "base/utility.h"
#include "base/binary_word.h"
#include "numeric/numerical_types.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/grid_set.h"
#include "geometry/list_set.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;


int main() {

  cout << "test_grid: " << flush;
  ofstream clog("test_grid.log");

  ListSet<Real,Rectangle> ls;
  Rectangle< Real > r;
  
  string input("[0,0.75]x[0,0.625]  [0.625,1]x[1,1.25]  [1.25,1.5]x[1.25,2.5] ");
  stringstream is(input);
  for (int i=0; i< 3; i++) {
    is >> r;
    clog << "r=" << r << endl;
    ls.push_back(r);
  }
  clog << ls << endl;

  RegularGrid<Real> gr(2,0.125);
  GridRectangleListSet<Real> grls(gr,ls);
  clog << grls << endl;
  clog << ListSet<Real,Rectangle>(grls) << endl;

  GridCellListSet<Real> gcls(grls);
  clog << gcls << endl;
  clog << ListSet<Real,Rectangle>(gcls) << endl;

  GridMaskSet<Real> gms(grls);
  clog << gms << endl;

  GridMaskSet<Real> gcms(gcls);
  clog << gcms << endl;
  GridCellListSet<Real> gclms(gms);
  clog << gclms << endl;

  clog << ListSet<Real,Rectangle>(gms) << endl;

  ListSet<Real,Rectangle> ls1,ls2;
  ls1.push_back(ls[0]);
  ls1.push_back(ls[2]);
  ls2.push_back(ls[1]);

  IrregularGrid<Real> igr1(ls1);
  IrregularGrid<Real> igr2(ls2);

  IrregularGrid<Real> igrj(igr1,igr2);

  clog << igr1 << "\n" << igr2 << "\n" << igrj << "\n";

  IrregularGrid<Real> gr1(ls1);
  GridRectangleListSet<Real> grls1(gr1,ls1);
  GridRectangleListSet<Real> grlsj1(igrj,grls1);
  clog << grlsj1 << "\n";

  IrregularGrid<Real> gr2(ls2);
  GridRectangleListSet<Real> grls2(gr2,ls2);
  GridRectangleListSet<Real> grlsj2(igrj,grls2);
  clog << grlsj2 << "\n";

  GridCellListSet<Real> gcls1(grls1);
  clog << gcls1 << "\n";
  GridRectangleListSet<Real> grlsc1(igrj,gcls1);
  clog << grlsc1 << "\n";

  FiniteGrid<Real> fgr=FiniteGrid<Real>(igrj,igrj.bounding_box());
  //GridMaskSet<Real> gms1(fgr,grlsj1);
  //GridMaskSet<Real> gms2(fgr,grlsj2);
  //clog << regular_intersection(gms1,gms2);
  //clog << join(gms1,gms2);

  clog.close();
  cout << "INCOMPLETE\n";

  return 0;
}
