/***************************************************************************
 *            test_multi_index.cpp
 *
 *  Copyright  2007-9  Pieter Collins
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

#include "config.h"
#include "algebra/multi_index.hpp"

#include "test.hpp"

using namespace std;
using namespace Ariadne;


class TestMultiIndex
{
  public:
    TestMultiIndex() { }

    Void test() {
        ARIADNE_TEST_CALL(test_constructor());
        ARIADNE_TEST_CALL(test_comparison());
        ARIADNE_TEST_CALL(test_addition());
        ARIADNE_TEST_CALL(test_increment());
    }

    Void test_constructor() {
        DegreeType a1[3]={2,0,3};
        MultiIndex i1(3u,a1);
        ARIADNE_TEST_EQUAL((Int)i1.number_of_variables(),3);
        ARIADNE_TEST_EQUAL((Int)i1.degree(),5);
        ARIADNE_TEST_EQUAL((Int)i1[0],2);
        ARIADNE_TEST_EQUAL((Int)i1[1],0);
        ARIADNE_TEST_EQUAL((Int)i1[2],3);

        DegreeType a2[5]={2,0,3,0,5};
        MultiIndex i2(5,a2);
        ARIADNE_TEST_EQUAL((Int)i2.number_of_variables(),5);
        ARIADNE_TEST_EQUAL((Int)i2.degree(),10);
        ARIADNE_TEST_EQUAL((Int)i2[3],0);
        ARIADNE_TEST_EQUAL((Int)i2[4],5);

        MultiIndex i3=MultiIndex::unit(9u,2u);
        ARIADNE_TEST_EQUAL((Int)i3.number_of_variables(),9);
        ARIADNE_TEST_EQUAL((Int)i3.degree(),1);
        ARIADNE_TEST_EQUAL((Int)i3[2],1);
        ARIADNE_TEST_EQUAL((Int)i3[8],0);
    }

    Void test_comparison() {
        DegreeType a1[3]={2,0,3};
        DegreeType a2[3]={1,0,4};
        DegreeType a3[3]={3,0,1};
        ARIADNE_TEST_CONSTRUCT(MultiIndex,i1,(3,a1));
        ARIADNE_TEST_CONSTRUCT(MultiIndex,i2,(3,a2));
        ARIADNE_TEST_CONSTRUCT(MultiIndex,i3,(3,a3));

        ARIADNE_TEST_ASSERT(graded_less(i1,i2));
        ARIADNE_TEST_ASSERT(graded_less(i3,i1));
        ARIADNE_TEST_ASSERT(graded_less(i3,i2));
        ARIADNE_TEST_ASSERT(!graded_less(i2,i1));

        ARIADNE_TEST_ASSERT(lexicographic_less(i2,i1));
        ARIADNE_TEST_ASSERT(lexicographic_less(i1,i3));
        ARIADNE_TEST_ASSERT(lexicographic_less(i2,i3));

        ARIADNE_TEST_ASSERT(reverse_lexicographic_less(i2,i1));
        ARIADNE_TEST_ASSERT(reverse_lexicographic_less(i1,i3));
        ARIADNE_TEST_ASSERT(reverse_lexicographic_less(i2,i3));

    }

    Void test_addition() {
        DegreeType a1[3]={2,0,3}, a2[3]={1,1,0}, a3[3]={3,1,3};
        ARIADNE_TEST_EQUAL(MultiIndex(3,a1)+MultiIndex(3,a2),MultiIndex(3,a3));
    }

    Void test_increment() {
        MultiIndex a(4);
        SizeType n=0;
        while(a.degree()<=5) {
            MultiIndex b=a; ++a; ++n;
            ARIADNE_TEST_BINARY_PREDICATE(graded_less,b,a);
            MultiIndex::IndexType d=0; for(MultiIndex::SizeType i=0; i!=a.size(); ++i) { d+=a[i]; }
            ARIADNE_TEST_EQUAL(a.degree(),d);
        }
        ARIADNE_ASSERT_EQUAL(n,126);
    }

};

Int main() {
  TestMultiIndex().test();
  return ARIADNE_TEST_FAILURES;
}

