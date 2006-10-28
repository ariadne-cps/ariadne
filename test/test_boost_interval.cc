/***************************************************************************
 *            test_boost_interval.cc
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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/io.hpp>

using namespace std;

int test_boost_interval();

int main() {
  cout << setprecision(20);
  test_boost_interval();
  
  return 0;
}


int
test_boost_interval()
{
  cout << "test_boost_interval<double>" << endl;
  using namespace boost::numeric;
  
  interval<double> o(1.0);
  interval<double> t(3.0);
  interval<double> odt=o/t;
  interval<double> oa=odt*t;
  cout << o << " / " << t << " = " << odt << endl;
  cout << o << " in " << oa << endl;
  assert(odt.lower()<odt.upper());
  assert(oa.lower()<o.lower());
  assert(oa.upper()>o.upper());
  
  return 0;
}
