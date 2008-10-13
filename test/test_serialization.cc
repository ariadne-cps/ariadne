/***************************************************************************
 *            test_serialization.cc
 *
 *  Copyright  2008  Pieter Collins
 *
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

#include "serialization.h"
#include "array.h"
#include "stlio.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

class TestSerialization
{
 public:
  int test() {
    test_array();
    return 0;
  }

  int test_array() {
    double data[]={5,23,42,111,242};
    const array<double> oary1(data,data+5);
    const array<double> oary2(data,data+3);
    ofstream ofs("test_serialization-array.txt");
    text_oarchive txtoa(ofs);
    txtoa << oary1 << oary2;
    ofs.close();
    
    array<double> iary1(100),iary2(1);
    ifstream ifs("test_serialization-array.txt");
    text_iarchive txtia(ifs);
    txtia >> iary1 >> iary2;
    ofs.close();
    
    ARIADNE_TEST_EQUAL(oary1,iary1);
    ARIADNE_TEST_EQUAL(oary2,iary2);

    return 0;
  }

};

  
int main() {
  return TestSerialization().test();
}
