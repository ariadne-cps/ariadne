/***************************************************************************
 *            test_boost.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
#include <cassert>

#include <boost/numeric/interval.hpp>

#include "test/test.h"

using namespace std;


class TestBoost {
 public:
	void test_rounding() 
	{
		double x=1;
		double y=3;
		double zl,zu;
		
		{ 
			boost::numeric::interval_lib::rounded_arith_std<double> rnd;
			zl=rnd.div_down(x,y);
			zu=rnd.div_up(x,y);
		}
		if(!(zl<zu)) {
			cerr << "Warning: boost::numeric::interval_lib::rounded_arith_std<double> does not round correctly\n";
		}
		ARIADNE_TEST_ASSERT(zl<zu);
	}

  void test()
  {
    ARIADNE_TEST_CALL(test_rounding());
  }
};

int main() 
{
  TestBoost().test();
  return ARIADNE_TEST_FAILURES;
}




