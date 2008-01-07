/***************************************************************************
 *            test_binary_word.cc
 *
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
#include <string>

#include "ariadne.h"

#include "base/stlio.h"
#include "combinatoric/binary_word.h"
#include "combinatoric/binary_tree.h"

#include "base/utility.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Combinatoric;
using namespace std;


int test_binary_word();


int main() {
  test_binary_word();
  return 0;
}

  
int test_binary_word() {
  cout << __PRETTY_FUNCTION__ << endl;
  
  string istr = "[0,1,1,0,1,0] [0,1,1,0,1,0] 011010 011 1101 111 1010 ";
  stringstream iss(istr);

  BinaryWord bw0,bw1,bw2,bw3,bw4,bw5,bw6;
  BinaryWordList bwl;
  BinaryWordFixedSizeList bwfsl(12);
  BinaryTree bwt;
  
  vector<bool> v;

  iss >> v;
  bw0=BinaryWord(v);
  
  iss >> bw1 >> bw2 >> bw3 >> bw4 >> bw5 >> bw6;
  
  ARIADNE_TEST_ASSERT(bw1.size()==6);
  
  cout << bw0 << " "  << bw1 << " " << bw2 << endl;
  ARIADNE_TEST_ASSERT(bw0==bw1);
  ARIADNE_TEST_ASSERT(bw1==bw2);
  ARIADNE_TEST_ASSERT(bw1.is_prefix(bw1));
  ARIADNE_TEST_ASSERT(bw1.is_subword(bw1));
  ARIADNE_TEST_ASSERT(bw3!=bw1);
  ARIADNE_TEST_ASSERT(bw3.is_prefix(bw1));
  ARIADNE_TEST_ASSERT(bw3.is_subword(bw1));
  ARIADNE_TEST_ASSERT(bw4!=bw1);
  ARIADNE_TEST_ASSERT(!bw4.is_prefix(bw1));
  ARIADNE_TEST_ASSERT(bw4.is_subword(bw1));
  ARIADNE_TEST_ASSERT(!bw5.is_subword(bw1));
  ARIADNE_TEST_ASSERT(!bw1.is_subword(bw5));
  ARIADNE_TEST_ASSERT(bw6.is_subword(bw1));
  
  bw5=bw4;
  ARIADNE_TEST_ASSERT(bw5==bw4);
  bw5=BinaryWord("011010");
  ARIADNE_TEST_ASSERT(bw5==bw1);

  return ARIADNE_TEST_FAILURES;
}
