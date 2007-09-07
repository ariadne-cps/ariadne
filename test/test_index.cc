/***************************************************************************
 *            test_polynomial.cc
 *
 *  Copyright  2007  Pieter Collins
 *  Email Pieter.Collins@cwi.nl
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

#include "function/sorted_index.h"
#include "function/multi_index.h"
#include "function/position_index.h"

#include "test.h"

using namespace std;
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Function;

class TestSortedIndex
{
 public:
  int test_construct() {
    uint a1[4]={1,1,3,3};
    uint a2[4]={1,1,3,3};
    SortedIndex s1(4,4,a1);
    SortedIndex s2(5,4,a2);
    
    ARIADNE_EVALUATE(s1);
    ARIADNE_EVALUATE(SortedIndex(4,4,a2));
    ARIADNE_EVALUATE(s2);
    
    return 0;
  }

  int test_compare() {
    uint a1[3]={1,1,3};
    uint a2[3]={1,2,2};
    uint a3[4]={2,3,1,2};
    uint a4[4]={1,2,2,3};
    
    ARIADNE_CONSTRUCT(SortedIndex,s1,(4,3,a1));
    ARIADNE_CONSTRUCT(SortedIndex,s2,(4,3,a2));
    ARIADNE_CONSTRUCT(SortedIndex,s3,(4,4,a3));
    ARIADNE_CONSTRUCT(SortedIndex,s4,(4,4,a4));
    SortedIndex s5(5,4,a4);
    
    ARIADNE_TEST_ASSERT(!(s1<s1));
    ARIADNE_TEST_ASSERT(s1<s2);
    ARIADNE_TEST_ASSERT(s3==s4);
    //ARIADNE_FAIL(s4==s5);

    return 0;
  }

  int test_addition() {
    uint a1[6]={1,1,3,4,6,6};
    uint a2[5]={1,2,2,3,5};
    uint a3[11]={1,1,1,2,2,3,3,4,5,6,6};
    
    SortedIndex s1(7,6,a1);
    SortedIndex s2(7,5,a2);
    SortedIndex s3(7,11,a3);
    
    ARIADNE_TEST_ASSERT(s1+s2==s3);
    ARIADNE_TEST_ASSERT((s1+=s2)==s3);
    return 0;
  }
  
  int test_increment() {
    uint a1[5]={0,0,2,2,2};
    uint a2[5]={0,0,2,2,3};
    uint a3[5]={0,0,2,3,3};
    uint a4[5]={0,0,3,3,3};
    uint a5[5]={0,1,1,1,1};
    uint a6[5]={3,3,3,3,3};
    uint a7[6]={0,0,0,0,0,0};
    ARIADNE_CONSTRUCT(SortedIndex,s0,(4));
    ARIADNE_CONSTRUCT(SortedIndex,s1,(4,5,a1));
    ARIADNE_CONSTRUCT(SortedIndex,s2,(4,5,a2));
    ARIADNE_CONSTRUCT(SortedIndex,s3,(4,5,a3));
    ARIADNE_CONSTRUCT(SortedIndex,s4,(4,5,a4));
    ARIADNE_CONSTRUCT(SortedIndex,s5,(4,5,a5));
    ARIADNE_CONSTRUCT(SortedIndex,s6,(4,5,a6));
    ARIADNE_CONSTRUCT(SortedIndex,s7,(4,6,a7));
    ARIADNE_EVALUATE(++s0);
    ARIADNE_TEST_ASSERT(++s1==s2);
    ARIADNE_TEST_ASSERT(++s2==s3);
    ARIADNE_TEST_ASSERT(++s3==s4);
    ARIADNE_TEST_ASSERT(++s4==s5);
    ARIADNE_TEST_ASSERT(++s5!=s6);
    ARIADNE_TEST_ASSERT(++s6==s7);
    
    return 0;
  }

  int test_position() {
    SortedIndex s(4);
    uint n=256;
    for(uint i=0; i!=n; ++i) {
      if(i!=s.position()) {
        ARIADNE_EVALUATE(s);
        ARIADNE_CHECK(s.position(),i);
      }
      ++s;
    }
    ARIADNE_EVALUATE(s);
    ARIADNE_CHECK(s.position(),n);
    return 0;
  }

  int test() {
    test_construct();
    test_compare();
    test_addition();
    test_increment();
    test_position();
    return 0;
  }

};





class TestMultiIndex
{
 public:
  TestMultiIndex() {
  }
  
  int test() {
    test_constructor();
    test_comparison();
    test_addition();
    test_position();
    test_conversion();
    return 0;
  }

  int test_constructor() {
    size_type a1[3]={2,0,3};
    MultiIndex i1(3,a1);
    ARIADNE_EVALUATE(i1);
    ARIADNE_TEST_ASSERT(i1.number_of_variables()==3);
    ARIADNE_TEST_ASSERT(i1.degree()==5);
    ARIADNE_TEST_ASSERT(i1.position()==44);
    return 0;
  }

  int test_comparison() {
    size_type a1[3]={2,0,3};
    size_type a2[3]={1,0,4};
    size_type a3[3]={3,0,1};
    MultiIndex i1(3,a1);
    MultiIndex i2(3,a2);
    MultiIndex i3(3,a3);
    ARIADNE_TEST_ASSERT(a2<a1);
    ARIADNE_TEST_ASSERT(a3<a1);
    ARIADNE_TEST_ASSERT(a3<a2);
    ARIADNE_TEST_ASSERT(!(a1<a2));
    return 0;
  }

  int test_addition() {
    uint a1[3]={2,0,3}, a2[3]={1,1,0}, a3[3]={3,1,3};
    ARIADNE_TEST_ASSERT(MultiIndex(3,a1)+MultiIndex(3,a2)==MultiIndex(3,a3));
    return 0;
  }
  
  int test_position() {
    MultiIndex s(4);
    uint n=256;
    for(uint i=0; i!=n; ++i) {
      if(i!=s.position()) {
        ARIADNE_EVALUATE(s);
        ARIADNE_CHECK(s.position(),i);
      }
      ++s;
    }
    ARIADNE_EVALUATE(s);
    ARIADNE_CHECK(s.position(),n);
    return 0;
  }

  int test_conversion() {
    size_type a1[6]={3,0,2,0,0,2};
    size_type a2[4]={3,0,2,1};
    
    ARIADNE_CONSTRUCT(SortedIndex,s1,(4,6,a1));
    ARIADNE_CONSTRUCT(MultiIndex,i1,(4,a2));
    ARIADNE_CONSTRUCT(PositionIndex,p1,(4,143));
    ARIADNE_CONSTRUCT(PositionIndex,p2,(4,14));

    ARIADNE_TEST_ASSERT(MultiIndex(s1)==i1);
    ARIADNE_TEST_ASSERT(SortedIndex(i1)==s1);
    ARIADNE_TEST_ASSERT(PositionIndex(s1)==p1);
    ARIADNE_TEST_ASSERT(PositionIndex(i1)==p1);
    ARIADNE_TEST_ASSERT(MultiIndex(p1)==i1);
    ARIADNE_EVALUATE(SortedIndex(p2));
    ARIADNE_EVALUATE(MultiIndex(PositionIndex(4,1234)))
    ARIADNE_TEST_ASSERT(SortedIndex(p1)==s1);

    return 0;
  }

  int test_compose() {
    return 0;
  }

};

int main() {
  TestSortedIndex().test();
  TestMultiIndex().test();
  cout << ARIADNE_TEST_FAILURES << " failures\n";
  return 0;
}

