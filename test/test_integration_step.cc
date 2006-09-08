/***************************************************************************
 *            test_integration_step.cc
 *
 *  Copyright  2006  Pieter Collins
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

#include <fstream>

#include "ariadne.h"
#include "real_typedef.h"
#include "numeric/numerical_types.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "system/affine_vector_field.h"
#include "evaluation/lohner_integrator.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

int main() {

  
  
  Evaluation::C1LohnerIntegrator<Real> lohner=Evaluation::C1LohnerIntegrator<Real>(0.125,0.5,0.0625);
  Geometry::Rectangle<Real> r=Geometry::Rectangle<Real>("[0.96,1.04]x[0.46,0.54]");
  cout << "r=" << r << endl;
  Geometry::Parallelotope<Real> p=Geometry::Parallelotope<Real>(r);
  cout << "p=" << p << endl;
  LinearAlgebra::Matrix<Real> A=LinearAlgebra::Matrix<Real>("[-0.25,-1;+1,-0.25]");
  cout << "A=" << A << endl;
  LinearAlgebra::Vector<Real> b=LinearAlgebra::Vector<Real>("[0,0]");
  cout << "b=" << b << endl;
  System::AffineVectorField<Real> avf=System::AffineVectorField<Real>(A,b);
  cout << "avf=" << avf << endl;

  Real h=Real(0.125);
  cout << "h=" << h << endl;
  Geometry::Parallelotope<Real> next_paral=lohner.integration_step(avf,p,h);
  cout << next_paral << "\n";
  
  

  
  return 0;
}
