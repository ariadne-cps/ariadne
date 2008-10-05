#include <iostream>
using namespace std;

#include "numeric.h"
#include "differential_variable.h"

int main()
{
  DifferentialVariable x=DifferentialVariable(3,2,1.0,1);
  cout << x << x*x << endl;
}

