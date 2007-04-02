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
#include <sstream>
#include <fstream>
#include <string>

#include "ariadne.h"
#include "test_float.h"
#include "base/stlio.h"
#include "combinatoric/binary_word.h"
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
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Combinatoric;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;

int main() {


  

  std::vector<Ariadne::dimension_type> seqa;

  Rectangle<Float> bb("[0,1]x[0,3]");
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

  PartitionScheme<Float> pg(bb,seq);
  PartitionTree<Float> pt(pg,bnt);
  PartitionTreeSet<Float> pts(pg,bnt,bla);
  ListSet< Rectangle<Float> > rls(pts);

  PartitionTree<Float>::const_iterator ptree_iter=pt.begin();
  PartitionTree<Float>::const_iterator ptree_end=pt.end();

  PartitionTreeSet<Float>::const_iterator ptreeset_iter=pts.begin();
  PartitionTreeSet<Float>::const_iterator ptreeset_end=pts.end();

  Vector< Interval<Float> > iv(2);
  iv[0]=Interval<Float>(-0.5,0.5);
  iv[1]=Interval<Float>(-0.5,0.5);
  Vector<Float> c("[0.125,0.25]");
  Matrix<Float> A("[2,1;0.5,1]");
  bb=Rectangle<Float>("[-4,4]x[-4,4]");
  Rectangle<Float> rect("[-1.5,2.5]x[0.875,2.25]");
  Parallelotope<Float> pltp(c,A);
  cout << "pltp=" << pltp << endl << "bb=" << bb << endl;


  seq=SubdivisionSequence(2);
  cout << "seq=" << seq << " " << seq.body_size() << " " << seq.tail_size() << " " << seq.dimension() << endl;

  assert(subset(pltp,bb));
  assert(!subset(bb,pltp));
  assert(!disjoint(pltp,bb));
  
  cout << "pt=" << pt << endl;
  cout << "pts=" << pts << endl;
  cout << "rect=" << rect << endl;
  cout << "pltp=" << pltp << endl;
  cout << "rls=" << rls << endl;
  pg=PartitionScheme<Float>(bb,seq);
  uint dpth=8;
  PartitionTreeSet<Float> ptsrua=under_approximation< Float, Rectangle<Float>  >(rect,pg,dpth);
  cout << "under_approximation(rect,pg,dpth)=" << ptsrua << endl;
  PartitionTreeSet<Float> ptsroa=over_approximation< Float, Rectangle<Float>  >(rect,pg,dpth);
  cout << "over_approximation(rect,pg,dpth)=" << ptsrua << endl;

  PartitionTreeSet<Float> ptsua=under_approximation< Float, Parallelotope<Float>  >(pltp,pg,dpth);
  cout << "ptsua=" << ptsua << endl;
  PartitionTreeSet<Float> ptsoa=over_approximation< Float, Parallelotope<Float>  >(pltp,pg,dpth);
  cout << "ptsoa=" << ptsoa << endl;
  PartitionTreeSet<Float> ptsina=inner_approximation< Float, Parallelotope<Float>  >(pltp,pg,dpth);
  cout << "ptsina=" << ptsina << endl;
  PartitionTreeSet<Float> ptsouta=outer_approximation< Float, Parallelotope<Float>  >(pltp,pg,dpth);
  cout << "ptsouta=" << ptsouta << endl;

  Grid<Float> g(2,Float(0.125));
  FiniteGrid<Float> fg(g,bb);
  ListSet< Rectangle<Float> > lsua=ptsua;
  cout << "lsua=" << lsua << endl;
  cout << "ptsua.size()=" << ptsua.size() << ", lsua.size()=" << lsua.size() << endl;
  // Under approximation of partition tree sets not currently implemented
  //GridMaskSet<Float> gmsina=under_approximation(ptsina,fg);
  //cout << "gmsina=" << gmsina << endl;
  GridMaskSet<Float> gmsouta=over_approximation(ptsouta,fg);
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
  //eps << gmsina;
  eps.set_fill_style(false);
  eps << pltp;
  eps.close();
  
  return 0;
}
