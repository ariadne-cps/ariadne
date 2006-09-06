/***************************************************************************
 *            test_binary_word.cc
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
#include <string>

#include "ariadne.h"

#include "base/stlio.h"
#include "base/exception.h"
#include "base/binary_word.h"
#include "base/binary_tree.h"
#include "numeric/numerical_types.h"

#include "base/utility.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

template class Rectangle< Rational >;

int main() {
  cout << "test_binary_word: " << flush;

  string istr = "[0,1,1,0,1,0] ";
  stringstream iss(istr);

  BinaryWord bw;
  BinaryWordList bwl;
  BinaryWordFixedSizeList bwfsl(12);
  BinaryTree bwt;
  
  vector<bool> v;

  iss >> v;
  bw=BinaryWord(v);
  
  cout << "INCOMPLETE\n";

    return 0;
}
