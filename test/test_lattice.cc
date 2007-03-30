/***************************************************************************
 *            test_lattice.cc
 *
 *  Copyright  2006  Pieter Collins, Alberto Casagrande
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

#include "declarations.h"
#include "test_float.h"

#include "combinatoric/lattice_set.h"
#include "combinatoric/lattice_map.h"

using namespace Ariadne;
using namespace Combinatoric;
using namespace System;
using namespace std;

int main() {

  LatticeCell lc("(0,2,-1)");
  cout << "lc=" << lc << endl;
  LatticeBlock lb("[0,5]x[-2,3]x[-1,1]");
  cout << "lb=" << lb << endl;
  
  
  IndexArray argary=IndexArray(3);
  argary[0]=-1;
  argary[1]=2;
  argary[2]=3;
  LatticeCell arglc(argary);
  cout << "arglc=" << arglc << endl;
  
 
  LatticeCellListSet imglcls=LatticeCellListSet(2);
  cout << "imglcls.dimension()=" << imglcls.dimension() << endl;
  cout << "imglcls=" << imglcls << endl;
  IndexArray imgary=IndexArray(2);
  imgary[0]=4;
  imgary[1]=2;
  imglcls.adjoin(LatticeCell(imgary));
  imgary[0]=5;
  imglcls.adjoin(LatticeCell(imgary));
  cout << "imglcls=" << imglcls << endl;
  
  LatticeMultiMap lmap=LatticeMultiMap(3,2);
  cout << "lmap=" << lmap << endl;
  lmap.adjoin_to_image(arglc,imglcls);
  cout << "lmap=" << lmap << endl;
  argary[0]=0;
  arglc=LatticeCell(argary);
  lmap.adjoin_to_image(arglc,imglcls);
  cout << "lmap=" << lmap << endl;
  imgary[1]=3;
  lmap.adjoin_to_image(arglc,LatticeCell(imgary));
  cout << "lmap=" << lmap << endl << endl;
  
  argary[0]=0;
  argary[1]=0;
  argary[2]=0;
  LatticeCell arglc2=LatticeCell(argary);
  
  cout << arglc2 << " " << lmap(arglc2) << endl;
  cout << arglc2 << " " << lmap.image(arglc2) << endl;
  cout << arglc << " " << lmap(arglc) << endl;
  cout << arglc << " " << lmap.image(arglc) << endl;
  cout << endl;
  
  cout << LatticeBlock(arglc) << " " << lmap(LatticeBlock(arglc)) << endl;
  LatticeCellListSet lcls(2);
  cout << lmap(arglc2) << endl;
  lcls=lmap.image(arglc2);
  cout << arglc2 << " " << lcls << endl;
  lcls=lmap(arglc);
  cout << lcls << endl;
  lcls.adjoin(LatticeCell(imgary));
  cout << lcls << endl;
  lcls=lmap(LatticeCellListSet(arglc));
  cout << lmap << endl;

  //assert(lmap.inverse().inverse() == lmap);

  LatticeMaskSet lms(lcls);
  cout << lmap.weak_preimage(lms);
  cout << lmap.inverse().image(lms);
  //assert(lmap.weak_preimage(lms) == lmap.inverse().image(lms));
  cout << lmap.inverse().image(lms);
  cout << lmap.strong_preimage(lms);

  
  return 0;
}
