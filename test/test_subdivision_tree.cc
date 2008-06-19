/***************************************************************************
 *            test_subdivision_tree.cc
 *
 *  Copyright  2005-6  Pieter Collins
 *  Email casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
#include <string>

#include "ariadne.h"
#include "base/stlio.h"
#include "combinatoric/subdivision_tree_set.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

int test_subdivision_tree();

int main() {
  test_subdivision_tree();
  cerr << "INCOMPLETE ";
  return 0;
}

int 
test_subdivision_tree()
{
  vector<dimension_type> seqa;

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

  //cout << "  seq=" << seq << "  bna=" << bna << "  bla=" << bla << endl;
  
  SubdivisionTree st(seq,bnt);
  
  SubdivisionTreeSet sts(seq,bnt,bla);
    
  return 0;
}
