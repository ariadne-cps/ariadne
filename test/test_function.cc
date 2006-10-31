/***************************************************************************
 *            test_function.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#include <cassert>
#include <fstream>
#include <string>
#include <iomanip>

#include "numeric/function.h"
#include "numeric/float64.h"
#include "numeric/mpfloat.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

typedef Float64(*unary_func)(const Float64&);


template<class R>
int
test_inverse_pair(
  std::string name, 
  R(*fnl)(const R&),
  R(*fnu)(const R&),
  R(*ifnl)(const R&),
  R(*ifnu)(const R&) )
{
  cout << name << endl;
  R o=1;
  R iml=fnl(o);
  R imu=fnu(o);
  R ol=ifnl(iml);
  R ou=ifnu(imu);
  cout << iml << " <= " << name << "(1) <= " << imu << endl;
  cout << ol << " <=    1   <= " << ou << endl;
  assert(iml<imu);
  assert(ol<=o);
  assert(o<=ou);
  return 0;
}
  
  
template<class R>
void
test_function()
{
  cout << __PRETTY_FUNCTION__ << endl;
  
  cout << setprecision(20);
  mpf_set_default_prec (128);
  
  test_inverse_pair("exp",&exp_down<R>,&exp_up<R>,&log_down<R>,&log_up<R>);
  test_inverse_pair("sin",&sin_down<R>,&sin_up<R>,&asin_down<R>,&asin_up<R>);
  test_inverse_pair("cos",&cos_down<R>,&cos_up<R>,&acos_down<R>,&acos_up<R>);
  test_inverse_pair("tan",&tan_down<R>,&tan_up<R>,&atan_down<R>,&atan_up<R>);
  test_inverse_pair("sinh",&sinh_down<R>,&sinh_up<R>,&asinh_down<R>,&asinh_up<R>);
  test_inverse_pair("cosh",&cosh_down<R>,&cosh_up<R>,&acosh_down<R>,&acosh_up<R>);
  test_inverse_pair("tanh",&tanh_down<R>,&tanh_up<R>,&atanh_down<R>,&atanh_up<R>);
  return;
}

int main() {

  test_function<Float64>();
  test_function<MPFloat>();

  return 0;
}
