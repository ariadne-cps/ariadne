#include <numeric.h>
#include "vector.h"
#include "matrix.h"

template class boost::numeric::ublas::matrix<Ariadne::Float>;
template class boost::numeric::ublas::matrix<Ariadne::Interval>;

namespace Ariadne {

Matrix<Float> midpoint(const Matrix<Interval>& A) {
  Matrix<Float> R(A.row_size(),A.column_size());
  for(uint i=0; i!=A.row_size(); ++i) {
    for(uint j=0; j!=A.row_size(); ++j) {
      R[i][j]=A[i][j].midpoint();
    }
  }
  return R;
}


template class Matrix<Float>;
template class Matrix<Interval>;

} // namespace Ariadne

