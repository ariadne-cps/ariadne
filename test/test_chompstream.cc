/***************************************************************************
 *            test_chompfstream.cc
 *
 *  24 June 2005
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
#include "combinatoric/lattice_set.h"
#include "combinatoric/lattice_map.h"
#include "output/chompstream.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

int main() {

  LatticeBlock lbb("[-1,3]x[-1,3]x[-1,3]");
  cout << lbb << endl;
  LatticeMaskSet lms(lbb);
  cout << lms << endl;
  
  LatticeBlock lb("[-1,2]x[1,3]x[0,1]");
  lms.adjoin(lb);
  LatticeCell lc=*lb.begin();
  
  LatticeMultiMap lmm(3,2);
  lmm.adjoin_to_image(LatticeCell("(-1,1,0)"),LatticeCell("(0,0)"));
  lmm.adjoin_to_image(LatticeCell("(-1,1,0)"),LatticeCell("(0,1)"));
  lmm.adjoin_to_image(LatticeCell("(-1,2,0)"),LatticeCell("(0,1)"));
  
  ofstream ofs;
  chompstream chomps(ofs);
  ofs.open("test_chompfstream.hom");
  chomps << lc << "\n\n";
  chomps << lms << "\n";
  chomps << lmm << "\n";
  ofs.close();

  chomps.redirect(cout);
  chomps << lc << "\n\n";
  chomps << lms << "\n";
  chomps << lmm << "\n";

  return 0;
}
