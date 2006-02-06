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
#include "exception.h"
#include "utility.h"
#include "numerical_type.h"
#include "point.h"
#include "rectangle.h"
#include "binary_word.h"
#include "grid_set.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

template class Rectangle<Rational>;
template class Grid<Rational>;
template class FiniteGrid<Rational>;
template class InfiniteGrid<Rational>;

template class GridCell<Rational>;
template class GridRectangle<Rational>;

template class GridRectangleListSet<Rational>;
template class GridCellListSet<Rational>;
template class GridMaskSet<Rational>;

template class Rectangle<Dyadic>;
template class Grid<Dyadic>;
template class FiniteGrid<Dyadic>;
template class InfiniteGrid<Dyadic>;

template class GridCell<Dyadic>;
template class GridRectangle<Dyadic>;

template class GridRectangleListSet<Dyadic>;
template class GridCellListSet<Dyadic>;
template class GridMaskSet<Dyadic>;

int main() {

  cout << "test_grid: " << flush;
  ofstream clog("test_grid.log");

  ListSet<Rational,Rectangle> ls;
  Rectangle< Rational > r;
  
  string input("[[0,5/6],[0,4/3]], [[2/3,1],[1,4/3]], [[4/3,3/2],[4/3,5/2]] ");
  stringstream is(input);
  for (int i=0; i< 4; i++) {
	  is >> r;
  	ls.push_back(r);
  }
  clog << ls << endl;

  FiniteGrid<Rational> gr(ls);
  GridRectangleListSet<Rational> grls(gr,ls);
  clog << grls << endl;
  clog << ListSet<Rational,Rectangle>(grls) << endl;

  GridCellListSet<Rational> gcls(grls);
  clog << gcls << endl;
  clog << ListSet<Rational,Rectangle>(gcls) << endl;

  GridMaskSet<Rational> gms(grls);
  clog << gms << endl;

  GridMaskSet<Rational> gcms(gcls);
  clog << gcms << endl;
  GridCellListSet<Rational> gclms(gms);
  clog << gclms << endl;

  clog << ListSet<Rational,Rectangle>(gms) << endl;

  ListSet<Rational,Rectangle> ls1,ls2;
  ls1.push_back(ls[0]);
  ls1.push_back(ls[2]);
  ls2.push_back(ls[1]);

  FiniteGrid<Rational> fg1(ls1);
  FiniteGrid<Rational> fg2(ls2);

  FiniteGrid<Rational> fgj(fg1,fg2);

  clog << fg1 << "\n" << fg2 << "\n" << fgj << "\n";
  clog << FiniteGrid<Rational>::index_translation(fg1,fgj) << "\n";
  clog << FiniteGrid<Rational>::index_translation(fg2,fgj) << "\n";

  FiniteGrid<Rational> gr1(ls1);
  GridRectangleListSet<Rational> grls1(gr1,ls1);
  GridRectangleListSet<Rational> grlsj1(fgj,grls1);
  clog << grlsj1 << "\n";

  FiniteGrid<Rational> gr2(ls2);
  GridRectangleListSet<Rational> grls2(gr2,ls2);
  GridRectangleListSet<Rational> grlsj2(fgj,grls2);
  clog << grlsj2 << "\n";

  GridCellListSet<Rational> gcls1(grls1);
  clog << gcls1 << "\n";
  GridRectangleListSet<Rational> grlsc1(fgj,gcls1);
  clog << grlsc1 << "\n";

  GridMaskSet<Rational> gms1(fgj,grlsj1);
  GridMaskSet<Rational> gms2(fgj,grlsj2);
  clog << regular_intersection(gms1,gms2);
  clog << join(gms1,gms2);

  clog.close();
  cout << "PASS\n";

  return 0;
}
