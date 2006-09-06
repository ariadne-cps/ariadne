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

#include <fstream>
#include <iomanip>

#include "numeric/function.h"
#include "numeric/float64.h"
#include "numeric/mpfloat.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

template<typename R>
void
test_function(ostream& clog) 
{
  clog << "test_function<" << name<R>() << ">" << endl;
  
  clog << setprecision(20);
  mpf_set_default_prec (8);

  { 
    R z(0);
    R o(1);
    R t(2);
    R fl,fu,zl,zu,ol,ou,tl,tu;
    
    clog << "sqrt" << endl;
    fl=sqrt_down(t);
    fu=sqrt_up(t);
    clog << fl << " <= " << fu << endl;
    assert(fl <= fu);
    tl=mul_down(fl,fl);
    tu=mul_up(fu,fu);
    clog << tl << " <= " << t << " <= " << tu << endl;
    assert(tl <= tu);
    assert(tl <= t);
    assert(t <= tu);

    clog << "exp" << endl;
    fl=exp_down(o);
    fu=exp_up(o);
    clog << fl << " <= " << fu << endl;
    assert(fl <= fu);
    ol=log_down(fl);
    ou=log_up(fu);
    clog << ol << " <= " << o << " <= " << ou << endl;
    assert(ol <= o && o <= ou);
    zl=log_down(ol);
    zu=log_up(ou);
    clog << zl << " <= " << z << " <= " << zu << endl;
    assert(zl <= z && z <= zu);
    
    clog << "sin" << endl;
    fl=sin_down(o);
    fu=sin_up(o);
    clog << fl << " <= " << fu << endl;
    assert(fl <= fu);
    ol=asin_down(fl);
    ou=asin_up(fl);
    clog << ol << " <= " << ou << endl;
    assert(ol <= ou);
    assert(ol <= o);
    assert(o <= ou);
  }
  
  return;
}

int main() {
  cout << "test_function: " << flush;
  ofstream clog("test_function.log");

  test_function<Float64>(clog);
  test_function<MPFloat>(clog);
  
  clog.close();
  cout << "INCOMPLETE\n";

  return 0;
}
