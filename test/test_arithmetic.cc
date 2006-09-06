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

template<typename R> void test_arithmetic(ofstream&);

int main() {

  cout << "test_arithmetic: " << flush;
  ofstream clog("test_arithmetic.log");
  clog << setprecision(20);
  mpf_set_default_prec (8);

  test_arithmetic<Float64>(clog);
  test_arithmetic<MPFloat>(clog);
  test_arithmetic<Rational>(clog);
  
  clog.close();
  cout << "INCOMPLETE\n";

  return 0;
}

template<typename Real>
void
test_arithmetic(ofstream& clog)
{
  clog << "test_arithmetic<" << name<Real>() << ">" << endl;
  
  Real f1(1.25);
  Real f2(2.25);
  Real f3;
  Real f4;
  //clog << "prec(f1)=" << f1.get_prec() << endl;
  
  f3=add_down(f1,f2);
  f4=add_up(f1,f2);
  clog << f3 << " <= " << f1 << " + " << f2 << " <= " << f4 << endl;
  f3=sub_down(f1,f2);
  f4=sub_up(f1,f2);
  clog << f3 << " <= " << f1 << " - " << f2 << " <= " << f4 << endl;
  f3=mul_down(f1,f2);
  f4=mul_up(f1,f2);
  clog << f3 << " <= " << f1 << " * " << f2 << " <= " << f4 << endl;
  f3=div_down(f1,f2);
  f4=div_up(f1,f2);
  clog << f3 << " <= " << f1 << " / " << f2 << " <= " << f4 << endl;
  f3=mul_down(f3,f2);
  f4=mul_up(f4,f2);
  clog << f3 << " <= " << f1 << " <= " << f4 << endl;
  
  Real h=0.5;
  Real e=0.5;
  Real o=1.0;
  Real l=add_down(o,e);
  Real u=add_up(o,e);
  int n=0;
  while(l==u && n!=256) {
    e=e*h;
    ++n;
    l=add_down(o,e);
    u=add_up(o,e);
  }
  clog << l << " <= " << o << " + " << e << " <= " << u << endl;
  
  Real z=0.0;
  Real t=3.0;
  Real tl=div_down(o,t);
  Real tu=div_up(o,t);
  clog << tl << " <= " << tu << endl;
  Real ol=mul_down(tl,t);
  Real ou=mul_up(tu,t);
  clog << ol << " <= " << o << " <= " << ou << endl;
  Real zl=sub_down(ol,o);
  Real zu=sub_up(ou,o);
  clog << "zl=" << zl << "  zu=" << zu << endl;
  assert(tl<=tu);
  assert(ol<=o);
  assert(o<=ou);
  assert(zl<=z);
  assert(z<=zu);
  //assert(zl<=z && z <=zu);
  
  Interval<Real> io(1);
  Interval<Real> it(3);
  Interval<Real> iao=(io/it)*it;
  clog << iao << endl;
  assert(in(o,iao));
  Interval<Real> iaz=iao-io;
  clog << iaz << endl;
  assert(in(z,iaz)); 
  clog << endl;
  
/*
  boost::numeric::interval<Real> io(1);
  boost::numeric::interval<Real> it(3);
  boost::numeric::interval<Real> iao=(io/it)*it;
  clog << iao << endl;
  assert(in(o,iao));
  boost::numeric::interval<Real> iaz=iao-io;
  clog << iaz << endl;
  assert(in(z,iaz));
*/
}
