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
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/list_set.h"
#include "geometry/partition_tree_set.h"

#include "test.h"

using namespace std;
using namespace Ariadne;
using namespace Ariadne::Geometry;

int main() {

  cout << "test_partition_grid: " << flush;
  ofstream clog("test_partition_grid.log");

  std::vector<Ariadne::dimension_type> seqa;

  Rectangle<Real> r;
  SubdivisionSequence seq(2);
  BinaryWord bnw;
  BooleanArray bna;
  BinaryTree bnt;
  BooleanArray bla;

  string input("[0,1]x[0,3]" "[0,1]" "[0,0,1,0,0,1,1,1,0,1,1]" "[1,0,1,1,0,1]");
  stringstream is(input);
  is >> r >> seqa >> bnw; //>> bla;
  bna=BooleanArray(bnw.begin(),bnw.begin()+bnw.size());
  bnt=BinaryTree(bna);
  seq=SubdivisionSequence(seqa.begin(),seqa.begin(),seqa.end());

  clog << r << "  " << seq << "  " << bna << "  " << bla << " " << endl;

  PartitionScheme<Real> pg(r,seq);
  PartitionTree<Real> pt(pg,bnt);
  PartitionTreeSet<Real> pts(pg,bnt,bla);

  PartitionTree<Real>::const_iterator ptree_iter=pt.begin();
  PartitionTree<Real>::const_iterator ptree_end=pt.end();

  PartitionTreeSet<Real>::const_iterator ptreeset_iter=pts.begin();
  PartitionTreeSet<Real>::const_iterator ptreeset_end=pts.end();

  clog << pt << endl;
  clog << pts << endl;
  clog << ListSet<Real,Rectangle>(pts) << endl;

  clog.close();
  cout << "PASS\n";

  return 0;
}
