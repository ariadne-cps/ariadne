#include <iostream>
#include <fstream>
#include <cstdio>


#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "system/henon_map.h"
#include "evaluation/newton.h"

#include "test.h"
#include "real_typedef.h"

using namespace std;
using namespace Ariadne;

int main() {
  cout << "test_interval_newton: " << flush;
  ofstream clog("test_partition_tree.log");

  System::HenonMap<Real> h(1.5,0.375);
  Geometry::Rectangle<Real> r("[-2.125,-2]x[-2.125,-2]");
  Real e=1e-10;
  
  clog << h << r << e << endl;
  
  Geometry::Rectangle<Real> fr=Evaluation::interval_newton(System::DifferenceMap<Real>(h),r,e);
  test_assert(fr.radius()<e,"interval_newton radius");
  clog << fr << endl;

  Geometry::Point<Real> fp("(-2.09201281589026466534,-2.09201281589026466534)");
  test_assert(fr.contains(fp),"interval_newton result");

  fp=h(fp);
  
  cout << "Pass\n";

  return 0;
}
