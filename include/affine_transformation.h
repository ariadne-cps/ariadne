#ifndef ARIADNE_AFFINE_TRANSFORMATION_H
#define ARIADNE_AFFINE_TRANSFORMATION_H

#include "vector.h"
#include "matrix.h"

namespace Ariadne {

template<class X> class AffineTransformation;

//! The affine transformation \f$x\mapsto A(x-c)+b \f$.
template<class X>
class AffineTransformation
{
  Vector<Float> _centre;
  Vector<X> _b;
  Matrix<X> _A;
 public:
  AffineTransformation(const Vector<Float>& c, const Vector<X>& b, const Matrix<X>& A)
    : _centre(c), _b(b), _A(A) { assert(A.column_size()==c.size()); assert(A.row_size()==b.size()); }
  uint result_size() const { return _b.size(); }
  uint argument_size() const { return _centre.size(); }
  const Vector<Float>& centre() const { return _centre; }
  const Vector<Float>& c() const { return _centre; }
  const Vector<X>& b() const { return _b; }
  const Matrix<X>& A() const { return _A; }
};

} // namespace Ariadne

#endif /* ARIADNE_AFFINE_TRANSFORMATION_H */

