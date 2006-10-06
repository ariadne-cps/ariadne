/***************************************************************************
 *            test_tensor.cc
 *
 *  Copyright  2006  Pieter Collins, Alberto Casagrande
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

#define NO_CBLAS

#include <iostream>
#include <fstream>

#include "declarations.h"
#include "real_typedef.h"
#include "numeric/numerical_types.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/tensor.h"

using namespace std;
using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;

template<class R>
int test_tensor();

int main() {
  test_tensor<Rational>();
  cout << "Incomplete";
  return 0;
}

template<class R>
int
test_tensor()
{
  return 0;
}
