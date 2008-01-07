/***************************************************************************
 *            test_partition_tree.cc
 *
 *  Copyright  2005-7  Pieter Collins
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
#include "geometry/box.h"
#include "geometry/box_list_set.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "geometry/partition_tree_set.h"
#include "output/epsstream.h"

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

  Box<Flt> bb("[0,1]x[0,3]");
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

  PartitionScheme<Flt> pg(bb,seq);
  PartitionTree<Flt> pt(pg,bnt);
  PartitionTreeSet<Flt> pts(pg,bnt,bla);
  BoxListSet<Flt> rls(pts);

  PartitionTree<Flt>::const_iterator ptree_iter=pt.begin();
  PartitionTree<Flt>::const_iterator ptree_end=pt.end();

  PartitionTreeSet<Flt>::const_iterator ptreeset_iter=pts.begin();
  PartitionTreeSet<Flt>::const_iterator ptreeset_end=pts.end();

  Vector< Interval<Flt> > iv(2);
  iv[0]=Interval<Flt>(-0.5,0.5);
  iv[1]=Interval<Flt>(-0.5,0.5);
  Point<Flt> c("(0.125,0.25)");
  Matrix<Flt> A("[2,1;0.5,1]");
  bb=Box<Flt>("[-4,4]x[-4,4]");
  Box<Flt> bx("[-1.5,2.5]x[0.875,2.25]");
  Zonotope<Flt> z(c,A);
  cout << "z=" << z << endl << "bb=" << bb << endl;


  seq=SubdivisionSequence(2);
  cout << "seq=" << seq << " " << seq.body_size() << " " << seq.tail_size() << " " << seq.dimension() << endl;

  assert((bool)(subset(z,bb)));
  assert((bool)(!subset(bb,z)));
  assert((bool)(!disjoint(z,bb)));
  
  cout << "pt=" << pt << endl;
  cout << "pts=" << pts << endl;
  cout << "bx=" << bx << endl;
  cout << "z=" << z << endl;
  cout << "rls=" << rls << endl;
  pg=PartitionScheme<Flt>(bb,seq);
  uint dpth=8;
  PartitionTreeCell<Flt> ptcrova=over_approximation(bx,pg);
  cout << "over_approximation(rect,pg,dpth)=" << ptcrova << endl; 

  PartitionTreeSet<Flt> ptsina=inner_approximation(z,pg,dpth);
  cout << "ptsina=" << ptsina << endl;
  PartitionTreeSet<Flt> ptsouta=outer_approximation(z,pg,dpth);
  cout << "ptsouta=" << ptsouta << endl;

  Grid<Flt> g(2,Flt(0.125));
  FiniteGrid<Flt> fg(g,bb);
  BoxListSet<Flt> lsina=ptsina;
  cout << "lsina=" << lsina << endl;
  cout << "ptsina.size()=" << ptsina.size() << ", lsina.size()=" << lsina.size() << endl;
  // Under approximation of partition tree sets not currently implemented
  //GridMaskSet<Flt> gmsina=under_approximation(ptsina,fg);
  //cout << "gmsina=" << gmsina << endl;
  GridMaskSet<Flt> gmsouta=outer_approximation(ptsouta,fg);
  cout << "gmsouta=" << gmsouta << endl;

  

  epsfstream eps;
  eps.open("test_partition_tree-1.eps",bb);
  eps << fill_colour(red) << ptsouta;
  eps << fill_colour(blue) << ptsina;
  eps << fill_colour(transparant) << z;
  eps.close();
  
  eps.open("test_partition_tree-2.eps",bb);
  eps << fill_colour(red) << gmsouta;
  //eps << fill_colour(blue) << gmsina;
  eps << fill_colour(transparant) << z;
  eps.close();
  
  return 0;
}
