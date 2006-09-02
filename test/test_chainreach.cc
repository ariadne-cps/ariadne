/***************************************************************************
 *            test_chainreach.cc
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#include <fstream>

#include "real_typedef.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "system/henon_map.h"
#include "evaluation/apply.h"
#include "output/epsfstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Postscript;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace std;

template class Rectangle< Real >;
template class Point< Real >;

int main() {
  cout << "test_chainreach: " << flush;
  ofstream clog("test_chainreach.log");

  Point<Real> params=Point<Real>("(1.5,0.875)");
  Real a=params[0];
  Real b=params[1];

  HenonMap<Real> h=HenonMap<Real>(a,b);
  Rectangle<Real> gbb=Rectangle<Real>("[-11,5]x[-8,8]") ;
  FiniteGrid<Real> g=FiniteGrid<Real>(gbb,128); // grid
  Rectangle<Real> ir=Rectangle<Real>("[1.499,1.501]x[0.499,0.501]"); // initial state
  Rectangle<Real> cb=Rectangle<Real>("[-4,4]x[-4,4]"); // cutoff box
  Rectangle<Real> epsbb=Rectangle<Real>("[-4.1,4.1]x[-4.1,4.1]"); // eps bounding box
  
  cb=Rectangle<Real>(gbb); // cutoff box
  epsbb=Rectangle<Real>(gbb); // eps bounding box
  
  GridMaskSet<Real> in=GridMaskSet<Real>(g);
  GridMaskSet<Real> bd=GridMaskSet<Real>(g);
  in.adjoin(over_approximation(ir,g));
  bd.adjoin(over_approximation(gbb,g));


  GridMaskSet<Real> cr=chainreach(h,in,bd);
  PartitionTreeSet<Real> ptcr=PartitionTreeSet<Real>(cr);

  clog << cr <<endl;
  clog << ptcr <<endl;
  
  clog << "cr.size()=" << cr.size() << endl;
  clog << "ptcr.size()=" << ptcr.size() << "  " << flush;
  clog << "ptcr.capacity()=" << ptcr.capacity() << endl;

  epsfstream eps=epsfstream("test_chainreach-1.eps",epsbb);
  eps.set_pen_colour("black");
  eps.set_fill_colour("white");
  eps << cb;
  eps.set_line_style(false);
  eps.set_fill_colour("green");
  eps << cr;
  eps.set_line_style(true);
  eps.set_fill_style(false);
  eps << ptcr.partition_tree();
  eps.set_line_style(true);
  eps.set_fill_colour("blue");
  eps << ir;
  eps.close();

  GridMaskSet<Real> gmcr=GridMaskSet<Real>(cr);
  eps.open("test_chainreach-2.eps",epsbb);
  eps.set_line_style(false);
  eps.set_fill_colour("red");
  eps << difference(gmcr.neighbourhood(),gmcr);
  eps.set_fill_colour("blue");
  eps << gmcr.adjoining();
  eps.set_fill_colour("green");
  eps << gmcr;
  eps.close();
  
  clog.close();
  cout << "PASS\n";

  return 0;

}
