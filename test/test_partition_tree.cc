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
#include "exception.h"
#include "utility.h"
#include "numerical_type.h"
#include "state.h"
#include "rectangle.h"
#include "binary_word.h"
#include "partition_tree.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

template class Rectangle<Dyadic>;
template class PartitionScheme<Dyadic>;
template class PartitionTreeCell<Dyadic>;
template class PartitionTree<Dyadic>;
template class PartitionTreeSet<Dyadic>;


int main() {

  cout << "test_partition_grid: " << flush;
  ofstream clog("test_partition_grid.log");
//  ostream& clog=cerr;

  std::vector<dimension_type> seqa;

  Rectangle<Dyadic> r;
  SubdivisionSequence seq;
  BinaryArray bna;
  BooleanArray bla;

  string input("[[0,1],[0,3]]" "[0,1]" "[0,0,1,0,0,1,1,1,0,1,1]" "[1,0,1,1,0,1]");
  stringstream is(input);
  is >> r >> seqa >> bna >> bla;
  seq=SubdivisionSequence(seqa.begin(),seqa.begin(),seqa.end());

  clog << r << "  " << seq << "  " << bna << "  " << bla << " " << endl;

  PartitionScheme<Dyadic> pg(r,seq);
  PartitionTree<Dyadic> pt(pg,bna);
  PartitionTreeSet<Dyadic> pts(pg,bna,bla);

  PartitionTree<Dyadic>::const_iterator ptree_iter=pt.begin();
  PartitionTree<Dyadic>::const_iterator ptree_end=pt.end();

  PartitionTreeSet<Dyadic>::const_iterator ptreeset_iter=pts.begin();
  PartitionTreeSet<Dyadic>::const_iterator ptreeset_end=pts.end();

  clog << pt << endl;
  clog << pts << endl;
  clog << ListSet<Dyadic,Rectangle>(pts) << endl;

  clog.close();
  cout << "PASS\n";

  return 0;
}
