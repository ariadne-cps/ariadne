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

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/string.hpp>

#include "config.h"
#include "macros.h"
#include "serialization.h"
#include "array.h"
#include "numeric.h"
#include "float.h"
#include "interval.h"
#include "stlio.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

namespace Ariadne {
void serialize(boost::archive::text_oarchive& a, Float& flt, const unsigned int v);
void serialize(boost::archive::text_iarchive& a, Float& flt, const unsigned int v);
}

class TestSerialization
{
  public:
    void test() {
        test_array();
        test_numeric();
    }

    void test_array() {
        double data[]={5,23,42,111,242};
        const Array<double> oary1(data,data+5);
        const Array<double> oary2(data,data+3);
        std::ofstream ofs("test_serialization-Array.txt");
        boost::archive::text_oarchive txtoa(ofs);
        txtoa << oary1 << oary2;
        ofs.close();

        Array<double> iary1(100),iary2(1);
        std::ifstream ifs("test_serialization-Array.txt");
        boost::archive::text_iarchive txtia(ifs);
        txtia >> iary1 >> iary2;
        ifs.close();

        ARIADNE_TEST_EQUAL(oary1,iary1);
        ARIADNE_TEST_EQUAL(oary2,iary2);
    }

    void test_numeric() {

        const Float nan = 0.0/0.0;
        const Float inf = 1.0/0.0;

        std::ofstream ofs("test_serialization-numeric.txt");
        boost::archive::text_oarchive txtoa(ofs);

        double xary[] = { 0.0, 1.0, 4.2, 1e-72, 1.2e+72 };
        const Array<Float> oxary(xary,xary+5);
        for(uint i=0; i!=oxary.size(); ++i) { txtoa << oxary[i]; }

        // Test output of special values
        ARIADNE_TEST_TRY( txtoa << inf );
        txtoa << static_cast<const double&>(1.0);
        ARIADNE_TEST_TRY( txtoa << nan );
        ARIADNE_TEST_TRY( txtoa << static_cast<const Float&>(-inf) );

        ofs.close();

        Array<Float> ixary(oxary.size());
        std::ifstream ifs("test_serialization-numeric.txt");
        boost::archive::text_iarchive txtia(ifs);

        for(uint i=0; i!=ixary.size(); ++i) {
            txtia >> ixary[i];
            ARIADNE_TEST_EQUALS(ixary[i],oxary[i]);
        }

        Float input_nan, input_inf, input_one;
        // Test input of constant "inf"
        try {
            txtia >> input_inf;
            if(input_inf==0.0) { ARIADNE_TEST_WARN("Inputing 'inf' floating-point value from archive yields 0.0"); }
            else { ARIADNE_TEST_EQUALS( input_inf, inf ); }
        } catch(...) {
            ARIADNE_TEST_WARN("Cannot input 'inf' floating-point value from archive.");
        }
        try {
            txtia >> input_one;
            ARIADNE_TEST_EQUALS( input_one, 1.0 );
        } catch(...) {
            ARIADNE_TEST_WARN("Error in archive read leaves archive in an invalid state.");
        }

        try {
            txtia >> input_nan;
        } catch(...) {
        }
        ARIADNE_TEST_ASSERT( std::isnan(input_nan.get_d()) );
        ifs.close();

    }

};

#include <sstream>

int main() {
    TestSerialization().test();
    return ARIADNE_TEST_FAILURES;
}
