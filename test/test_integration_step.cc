#include "ariadne.h"
#include "numeric/numerical_types.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "system/affine_vector_field.h"
#include "evaluation/lohner_integrator.h"

#include "test.h"

using namespace Ariadne;

typedef MPFloat Real;

int main() {
  Evaluation::C1LohnerIntegrator<Real> lohner=Evaluation::C1LohnerIntegrator<Real>(0.125,0.5,0.0625);
  Geometry::Rectangle<Real> r=Geometry::Rectangle<Real>("[0.96,1.04]x[0.46,0.54]");
  Geometry::Parallelotope<Real> p=Geometry::Parallelotope<Real>(r);
  LinearAlgebra::Matrix<Real> A=LinearAlgebra::Matrix<Real>("[-0.25,-1;+1,-0.25]");
  LinearAlgebra::Vector<Real> b=LinearAlgebra::Vector<Real>("[0,0]");
  System::AffineVectorField<Real> avf=System::AffineVectorField<Real>(A,b);

  Real step_size=Real(0.125);
  Geometry::Parallelotope<Real> next_paral=lohner.integration_step(avf,p,step_size);
  std::cerr << next_paral << "\n";
    
  return 0;
}
