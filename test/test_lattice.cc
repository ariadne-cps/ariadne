
#include "declarations.h"
#include "real_typedef.h"

#include "geometry/lattice_set.h"
#include "system/lattice_map.h"

using namespace Ariadne;
using namespace Geometry;
using namespace Evaluation;
using namespace std;

int main() {
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
  cout << lm << endl << endl;
  
  argary[0]=0;
  argary[1]=0;
  argary[2]=0;
  LatticeCell arglc2=LatticeCell(argary);
  
  cout << arglc2 << " " << lm(arglc2) << endl;
  cout << arglc2 << " " << lm.apply(arglc2) << endl;
  cout << arglc << " " << lm(arglc) << endl;
  cout << arglc << " " << lm.apply(arglc) << endl;
  cout << endl;
  
  cout << LatticeRectangle(arglc) << " " << lm(LatticeRectangle(arglc)) << endl;
  LatticeCellListSet lcls(0);
  cout << lm(arglc2) << endl;
  lcls=lm.apply(arglc2);
  cout << arglc2 << " " << lcls << endl;
  lcls=lm(arglc);
  cout << lcls << endl;
  lcls.adjoin(LatticeCell(imgary));
  cout << lcls << endl;
  lcls=lm(LatticeCellListSet(arglc));
  cout << lm << endl;
  cout << "Passed" << endl;
  
  return 0;
}
