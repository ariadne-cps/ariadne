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

  

  System::HenonMap<Real> h(1.5,0.375);
  Geometry::Rectangle<Real> r("[-2.125,-2]x[-2.125,-2]");
  Real e=1e-10;
  
  cout << h << r << e << endl;
  
  Geometry::Rectangle<Real> fr;
  try {
    fr=Evaluation::interval_newton(System::DifferenceMap<Real>(h),r,e);
  }
  catch(Evaluation::EvaluationException e) {
    cout << "No solution found" << endl;
    throw e;
  }
  assert(fr.radius()<e);
  cout << fr << endl;

  Geometry::Point<Real> fp("(-2.09201281589026466534,-2.09201281589026466534)");
  assert(fr.contains(fp));

  fp=h(fp);
  


  return 0;
}
