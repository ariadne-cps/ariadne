#include "blas.hpp"
#include "blas_output.hpp"

using namespace BLAS;
using std::cout;

int main() {
  complex<double> X[3] = {complex<double>(1,1),complex<double>(0,2),3};
  complex<double> Y[3] = {4,5,6};
  complex<double> A[12] = {1,2,3,4,5,6,7,8,9,10,11,12};
  int m=3;
  int n=3;

  //axpy(3,0.1,X,1,Y,1);

  cout << vector(3,X) << "\n\n";
  cout << vector(3,Y) << "\n\n";
  cout << matrix(3,3,A) << "\n";
  cout << vector(3,X) << "\n";
  cout << vector(3,Y) << "\n";

  gemv(RowMajor, NoTrans, 3,3,complex<double>(0.01),A,3,X,1,complex<double>(2.0),Y,1);
  cout << vector(3,Y) << "\n\n";

  cout << "Exiting\n";

  return 0;
}
  
