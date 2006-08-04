/***************************************************************************
 *            test_partition_tree.cc
 *
 *  1 July 2005
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
#include "linear_algebra/matrix.h"
#include "linear_algebra/vector.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/list_set.h"
#include "geometry/partition_tree_set.h"
#include "utility/epsfstream.h"

#include "test.h"

using namespace std;
using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Postscript;

int main() {

  cout << "test_partition_tree: " << flush;
  ofstream clog("test_partition_tree.log");

  std::vector<Ariadne::dimension_type> seqa;

  Rectangle<Real> bb;
  SubdivisionSequence seq(2);
  BinaryWord bnw;
  BooleanArray bna;
  BinaryTree bnt;
  BooleanArray bla;

  string input("[0,1]x[0,3]" "[0,1]" "[0,0,1,0,0,1,1,1,0,1,1]" "[1,0,1,1,0,1]" 
               "[0.125,0.25]" "[2,1;0.5,1] " 
               " [-4,4]x[-4,4]  [1,1,0,1,0,0,1,0,1,0,1,1,0,0,0] [0,1,0,1,0,1,0,1]");
  stringstream is(input);
  is >> bb >> seqa;
  is >> bnw;
  bna=BooleanArray(bnw.begin(),bnw.begin()+bnw.size());
  bnt=BinaryTree(bna);
  is >> bnw;
  bla=BooleanArray(bnw.begin(),bnw.begin()+bnw.size());
  seq=SubdivisionSequence(seqa.begin(),seqa.begin(),seqa.end());

  clog << bb << "  " << seq << "  " << bna << "  " << bla << " " << endl;

  PartitionScheme<Real> pg(bb,seq);
  PartitionTree<Real> pt(pg,bnt);
  PartitionTreeSet<Real> pts(pg,bnt,bla);

  PartitionTree<Real>::const_iterator ptree_iter=pt.begin();
  PartitionTree<Real>::const_iterator ptree_end=pt.end();

  PartitionTreeSet<Real>::const_iterator ptreeset_iter=pts.begin();
  PartitionTreeSet<Real>::const_iterator ptreeset_end=pts.end();

  IntervalVector<Real> iv(2);
  iv[0]=Interval<Real>(-0.5,0.5);
  iv[1]=Interval<Real>(-0.5,0.5);
  Vector<Real> c;
  Matrix<Real> A;
  is >> c >> A >> bb;
  Parallelotope<Real> pltp(c,A);
  seq=SubdivisionSequence(2);
  pg=PartitionScheme<Real>(bb,seq);
  uint dpth=12;
  PartitionTreeSet<Real> ptsoa=over_approximation< Real, Parallelotope<Real>  >(pltp,pg,dpth);
  PartitionTreeSet<Real> ptsua=under_approximation< Real, Parallelotope<Real>  >(pltp,pg,dpth);
  PartitionTreeSet<Real> ptsouta=outer_approximation< Real, Parallelotope<Real>  >(pltp,pg,dpth);
  PartitionTreeSet<Real> ptsina=inner_approximation< Real, Parallelotope<Real>  >(pltp,pg,dpth);
  RegularGrid<Real> rg(2,Real(0.125));
  Grid<Real>& g=rg;
  LatticeRectangle lr(2);
  is >> lr;
  FiniteGrid<Real> fg(g,bb);
  ListSet<Real,Rectangle> lsua=ptsua;
  GridMaskSet<Real> gmsouta=under_approximation(ptsouta,fg);
  GridMaskSet<Real> gmsina=over_approximation(ptsina,fg);
  
  
  clog << pt << endl;
  clog << pts << endl;
  clog << ListSet<Real,Rectangle>(pts) << endl;
  clog << pltp << endl;
  clog << ptsua << endl;
  clog << ptsoa << endl;
  clog << ptsina << endl;
  clog << ptsouta << endl;
  clog << gmsina << endl;
  clog << gmsouta << endl;

  clog.close();

  epsfstream<Real> eps;
  eps.open("test_partition_tree-1.eps",bb);
  eps.set_fill_colour("red");
  eps << ptsouta;
  eps.set_fill_colour("blue");
  eps << ptsina;
  eps.set_fill_style(false);
  eps << pltp;
  eps.close();
  
  eps.open("test_partition_tree-2.eps",bb);
  eps.set_fill_colour("red");
  eps << ptsoa;
  eps.set_fill_colour("blue");
  eps << ptsua;
  eps.set_fill_style(false);
  eps << pltp;
  eps.close();
  
  eps.open("test_partition_tree-3.eps",bb);
  eps.set_fill_colour("red");
  eps << gmsouta;
  eps.set_fill_colour("blue");
  eps << gmsina;
  eps.set_fill_style(false);
  eps << pltp;
  eps.close();
  
  cout << "PASS\n";

  return 0;
}
