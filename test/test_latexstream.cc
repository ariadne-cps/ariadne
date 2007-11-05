/***************************************************************************
 *            test_latexstream.cc
 *
 *  Copyright  2007  Pieter Collins
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
#include <sstream>
#include <string>

#include "ariadne.h"
#include "test_float.h"

#include "numeric/integer.h"
#include "numeric/rational.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix.h"
#include "geometry/rectangle.h"
#include "output/latexstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace std;

int main() {

  Integer z(42);
  Rational q(23,5);
  Flt x(1.25);

  Vector<Flt> v("[-1,2,0.5]");
  Matrix<Flt> A("[-1,2,0.5;0.5,-1.0,2.0]");
  Rectangle<Flt> r("[-1,1]x[0.5,1.5]");

  latexfstream tex;
  tex.open("test_latexstream-1.tex");
  tex << "\\[ " << z << "\\ " << q << "\\ " << x << "\\]\n";
  tex << "\\[ " << A << " \\quad " << v << "\\]\n";
  tex.close();
  
  tex.open("test_latexstream-2.tex","\\usepackage{amssymb}");
  tex << "\\[ " << r << " \\subset\\mathbb{R}^2" << "\\]\n";
  tex.close();
  
  return 0;
}
