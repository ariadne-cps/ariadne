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
#include "real_typedef.h"

#include "geometry/lattice_set.h"
#include "system/lattice_map.h"

using namespace Ariadne;
using namespace Geometry;
using namespace System;
using namespace std;

int main() {
  cout << "test_lattice: " << flush;
  ofstream clog("test_lattice.log");
  
  IndexArray argary=IndexArray(3);
  argary[0]=-1;
  argary[1]=2;
  argary[2]=3;
  LatticeCell arglc(argary);
  
 
  LatticeCellListSet imglcls=LatticeCellListSet(2);
  IndexArray imgary=IndexArray(2);
  imgary[0]=4;
  imgary[1]=2;
  imglcls.adjoin(LatticeCell(imgary));
  imgary[0]=5;
  imglcls.adjoin(LatticeCell(imgary));
  
  LatticeMap lm=LatticeMap(3,2);
  lm.adjoin_to_image(arglc,imglcls);
  argary[0]=0;
  arglc=LatticeCell(argary);
  lm.adjoin_to_image(arglc,imglcls);
  imgary[1]=3;
  lm.adjoin_to_image(arglc,LatticeCell(imgary));
  clog << lm << endl << endl;
  
  argary[0]=0;
  argary[1]=0;
  argary[2]=0;
  LatticeCell arglc2=LatticeCell(argary);
  
  clog << arglc2 << " " << lm(arglc2) << endl;
  clog << arglc2 << " " << lm.apply(arglc2) << endl;
  clog << arglc << " " << lm(arglc) << endl;
  clog << arglc << " " << lm.apply(arglc) << endl;
  clog << endl;
  
  clog << LatticeRectangle(arglc) << " " << lm(LatticeRectangle(arglc)) << endl;
  LatticeCellListSet lcls(0);
  clog << lm(arglc2) << endl;
  lcls=lm.apply(arglc2);
  clog << arglc2 << " " << lcls << endl;
  lcls=lm(arglc);
  clog << lcls << endl;
  lcls.adjoin(LatticeCell(imgary));
  clog << lcls << endl;
  lcls=lm(LatticeCellListSet(arglc));
  clog << lm << endl;

  
  clog.close();
  cout << "PASS" << endl;
  
  return 0;
}
