/***************************************************************************
 *            test_arithmetic.cc
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

#include <fstream>
#include <iomanip>

#include <gmpxx.h>
#include <mpfr.h>
#include <boost/numeric/interval.hpp>

#include "numeric/arithmetic.h"
#include "numeric/interval.h"
#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"

#include "test.h"

using namespace std;
using namespace Ariadne::Numeric;
using Ariadne::name;
using Ariadne::Rational;
using Ariadne::MPFloat;
using Ariadne::Float64;

void test_boost_rounding(); 
template<class R> void test_arithmetic();
template<> void test_arithmetic<Rational>();

int main() {

  cout << setprecision(20);
  mpf_set_default_prec (8);

  test_boost_rounding();
  
  test_arithmetic<Rational>();
  test_arithmetic<MPFloat>();
  test_arithmetic<Float64>();
 
  cerr << "INCOMPLETE ";

  return 0;
}

void
test_boost_rounding() 
{
  double x=1;
  double y=3;
  double zl,zu;
  
  { 
    boost::numeric::interval_lib::rounded_arith_std<double> rnd;
    zl=rnd.div_down(x,y);
    zu=rnd.div_up(x,y);
  }
  cout << zl << " <= " << x << "/" << y << " <= " << zu << endl;
}

  
template<class R>
void
test_arithmetic()
{
  cout << "test_arithmetic<" << name<R>() << ">" << endl;
  
  R f1(1.25);
  R f2(2.25);
  R f3;
  R f4;
  //cout << "prec(f1)=" << f1.get_prec() << endl;
  
  f3=add_down(f1,f2);
  f4=add_up(f1,f2);
  cout << f3 << " <= " << f1 << " + " << f2 << " <= " << f4 << endl;
  f3=sub_down(f1,f2);
  f4=sub_up(f1,f2);
  cout << f3 << " <= " << f1 << " - " << f2 << " <= " << f4 << endl;
  f3=mul_down(f1,f2);
  f4=mul_up(f1,f2);
  cout << f3 << " <= " << f1 << " * " << f2 << " <= " << f4 << endl;
  f3=div_down(f1,f2);
  f4=div_up(f1,f2);
  cout << f3 << " <= " << f1 << " / " << f2 << " <= " << f4 << endl;
  assert(f3<f4);
  assert(Rational(f3)<=Rational(5,9));
  assert(Rational(f4)>=Rational(5,9));
  cout << mul_down(f3,f2) << " <= (" << f1 << "/" << f2 << ")*" << f2 << " <= " << mul_up(f4,f2) << endl;
  assert(mul_down(f3,f2)<f1);
  assert(mul_up(f4,f2)>f1);
  
  f3=mul_down(f3,f2);
  f4=mul_up(f4,f2);
  cout << f3 << " <= " << f1 << " <= " << f4 << endl;
  
  R h=0.5;
  R e=0.5;
  R o=1.0;
  R l=add_down(o,e);
  R u=add_up(o,e);
  int n=0;
  while(l==u && n!=256) {
    e=mul_approx(e,h);
    ++n;
    l=add_down(o,e);
    u=add_up(o,e);
  }
  cout << l << " <= " << o << " + " << e << " <= " << u << endl;
  
  R z=0.0;
  R t=3.0;
  R tl=div_down(o,t);
  R tu=div_up(o,t);
  cout << tl << " <= " << tu << endl;
  R ol=mul_down(tl,t);
  R ou=mul_up(tu,t);
  cout << ol << " <= " << o << " <= " << ou << endl;
  R zl=sub_down(ol,o);
  R zu=sub_up(ou,o);
  cout << "zl=" << zl << "  zu=" << zu << endl;
  assert(tl<tu);
  assert(ol<o);
  assert(o<ou);
  assert(zl<z);
  assert(z<zu);
  //assert(zl<=z && z <=zu);
  
  Interval<R> io(1);
  Interval<R> it(3);
  Interval<R> iao=(io/it)*it;
  cout << iao << endl;
  assert(contains_value(iao,o));
  Interval<R> iaz=iao-io;
  cout << iaz << endl;
  assert(contains_value(iaz,z)); 
  cout << endl;

  return;
}


template<>
void
test_arithmetic<Rational>()
{
  typedef Rational R;
  
  R f1(1.25);
  R f2(2.25);
  R f3;

  f3=add(f1,f2);
  cout << f1 << " + " << f2 << " = " << f3 << endl;
  assert(f3==R(7,2));
  f3=sub(f1,f2);
  cout << f1 << " - " << f2 << " = " << f3 << endl;
  assert(f3==R(-1,1));
  f3=mul(f1,f2);
  cout << f1 << " * " << f2 << " = " << f3 << endl;
  assert(f3==R(45,16));
  f3=div(f1,f2);
  cout << f1 << " / " << f2 << " = " << f3 << endl;
  assert(f3==R(5,9));

}
