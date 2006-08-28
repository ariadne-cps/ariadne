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

#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/rounded_arith.hpp>
#include <boost/numeric/interval/io.hpp>

#include "real_typedef.h"
#include "numeric/function.h"

#include "test.h"

using namespace std;
using namespace Ariadne;
using boost::numeric::interval;

template<typename R>
void
test_function(ostream& clog) 
{
  clog << setprecision(20);
  mpf_set_default_prec (8);

  { 
    Numeric::rounding<R> rnd;
    R z(0);
    R o(1);
      
    R eol=rnd.exp_down(o);
    R eou=rnd.exp_up(o);
    clog << eol << " < " << eou << endl;
    assert(eol < eou);
    R ol=rnd.log_down(eol);
    R ou=rnd.log_up(eou);
    clog << ol << " < " << o << " < " << ou << endl;
    assert(ol < o && o < ou);
    R zl=rnd.log_down(ol);
    R zu=rnd.log_up(ou);
    clog << zl << " < " << z << " < " << zu << endl;
    assert(zl < z && z < zu);
  }
  
  return;
}

int main() {
  cout << "test_function: " << flush;
  ofstream clog("test_function.log");

  test_function<Float64>(clog);
  //test_function<MPFloat>(clog);
  
  clog.close();
  cout << "INCOMPLETE\n";

  return 0;
}
