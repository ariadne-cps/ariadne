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
#include "base/stlio.h"
#include "base/utility.h"
#include "combinatoric/binary_word.h"
#include "numeric/numerical_types.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/vector.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/list_set.h"
#include "geometry/partition_tree_set.h"
#include "output/epsfstream.h"

#include "test.h"

using namespace std;
using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Combinatoric;
using namespace Ariadne::Geometry;
using namespace Ariadne::Postscript;

int main() {


  

  std::vector<Ariadne::dimension_type> seqa;

  Rectangle<Real> bb("[0,1]x[0,3]");
  SubdivisionSequence seq(2);
  BinaryWord bnw;
  BooleanArray bna;
  BinaryTree bnt;
  BooleanArray bla;

  string input(" [1,1,0,1,0,0,1,0,1,0,1,1,0,0,0] [0,1,0,1,0,1,0,1]");
  stringstream is(input);

  seqa.resize(2); seqa[0]=0; seqa[1]=1;
  seq=SubdivisionSequence(seqa.begin(),seqa.begin(),seqa.end());
  bnw=BinaryWord("[0,0,1,0,0,1,1,1,0,1,1]");
  bna=BooleanArray(bnw.begin(),bnw.begin()+bnw.size());
  bnt=BinaryTree(bna);
  bnw=BinaryWord("[1,0,1,1,0,1]");
  bla=BooleanArray(bnw.begin(),bnw.begin()+bnw.size());

  cout << "bb=" << bb << "  seq=" << seq << "  bna=" << bna << "  bla=" << bla << endl;

  PartitionScheme<Real> pg(bb,seq);
  PartitionTree<Real> pt(pg,bnt);
  PartitionTreeSet<Real> pts(pg,bnt,bla);
  ListSet<Real,Rectangle> rls(pts);

  PartitionTree<Real>::const_iterator ptree_iter=pt.begin();
  PartitionTree<Real>::const_iterator ptree_end=pt.end();

  PartitionTreeSet<Real>::const_iterator ptreeset_iter=pts.begin();
  PartitionTreeSet<Real>::const_iterator ptreeset_end=pts.end();

  Vector< Interval<Real> > iv(2);
  iv[0]=Interval<Real>(-0.5,0.5);
  iv[1]=Interval<Real>(-0.5,0.5);
  Vector<Real> c("[0.125,0.25]");
  Matrix<Real> A("[2,1;0.5,1]");
  bb=Rectangle<Real>("[-4,4]x[-4,4]");
  Parallelotope<Real> pltp(c,A);
  cout << "pltp=" << pltp << endl << "bb=" << bb << endl;


  seq=SubdivisionSequence(2);
  cout << "seq=" << seq << " " << seq.body_size() << " " << seq.tail_size() << " " << seq.dimension() << endl;

  
  cout << "pt=" << pt << endl;
  cout << "pts=" << pts << endl;
  cout << "rls=" << rls << endl;
  cout << "pltp=" << pltp << endl;
  pg=PartitionScheme<Real>(bb,seq);
  uint dpth=8;
  PartitionTreeSet<Real> ptsua=under_approximation< Real, Parallelotope<Real>  >(pltp,pg,dpth);
  cout << "ptsua=" << ptsua << endl;
  PartitionTreeSet<Real> ptsoa=over_approximation< Real, Parallelotope<Real>  >(pltp,pg,dpth);
  cout << "ptsoa=" << ptsoa << endl;
  PartitionTreeSet<Real> ptsina=inner_approximation< Real, Parallelotope<Real>  >(pltp,pg,dpth);
  cout << "ptsina=" << ptsina << endl;
  PartitionTreeSet<Real> ptsouta=outer_approximation< Real, Parallelotope<Real>  >(pltp,pg,dpth);
  cout << "ptsouta=" << ptsouta << endl;

  RegularGrid<Real> rg(2,Real(0.125));
  Grid<Real>& g=rg;
  FiniteGrid<Real> fg(g,bb);
  ListSet<Real,Rectangle> lsua=ptsua;
  cout << "lsua=" << lsua << endl;
  GridMaskSet<Real> gmsina=under_approximation(ptsina,fg);
  cout << "gmsina=" << gmsina << endl;
  GridMaskSet<Real> gmsouta=over_approximation(ptsouta,fg);
  cout << "gmsouta=" << gmsouta << endl;

  

  epsfstream eps;
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
  
  return 0;
}
