/***************************************************************************
 *            test_hybrid_set.cc
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

#include "config.h"
#include "test.h"

#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "function/function.h"

#include "geometry/function_set.h"
#include "hybrid/hybrid_set.h"
#include "output/graphics.h"

using namespace std;
using namespace Ariadne;


class TestHybridSet {
  public:
    Void test();
  private:
    Void test_hybrid_list_set();
    Void test_hybrid_grid_set();
};

Void
TestHybridSet::test()
{
    ARIADNE_TEST_CALL(test_hybrid_list_set());
    //ARIADNE_TEST_CALL(test_hybrid_grid_set());
}


Void
TestHybridSet::test_hybrid_grid_set()
{
    // Test compilation without correctness
    Grid g;
    GridTreeSet gts;
    MonolithicHybridSpace hspc;
    DiscreteLocation loc;
    HybridGrid hg(hspc);
    HybridGridTreeSet hgts(hg);
    Figure fig;
    hg.has_location(loc);
    g=hg[loc];

    hgts.has_location(loc);
    GridTreeSet const& gtscr = hgts[loc];
    GridTreeSet& gtsr = hgts[loc];
    draw(fig,hgts[loc]);
}


Void
TestHybridSet::test_hybrid_list_set()
{
    HybridListSet<ExactBoxType> hls;
    DiscreteLocation loc1(123);
    DiscreteLocation loc2(105);
    RealVariable x("x"), y("y")
    ;
    RealSpace spc1({x});
    RealSpace spc2({x,y});

    ExactBoxType bx1={{0,1}};
    ExactBoxType bx2={{2,3}};
    ExactBoxType bx3={{1,2},{2,3}};
    ExactBoxType bx4={{4,5},{5,6}};
    ExactBoxType bx5={{6,7},{8,9}};
    ARIADNE_TEST_EXECUTE(hls.adjoin(loc1,spc1,bx1));
    ARIADNE_TEST_FAIL(hls.adjoin(loc1,bx3)); // Should fail due to incompatible dimensions
    ARIADNE_TEST_FAIL(hls.adjoin(loc2,bx3)); // Should fail due to unspecified space
    ARIADNE_TEST_EXECUTE(hls.adjoin(loc2,spc2,bx3));
    ARIADNE_TEST_EXECUTE(hls.adjoin(loc2,bx4));
    ARIADNE_TEST_PRINT(hls);
    ARIADNE_TEST_EXECUTE(hls.adjoin(HybridBoxType(loc1,spc1,bx2)));
    ARIADNE_TEST_EXECUTE(hls.adjoin(loc2,bx5));

    ARIADNE_TEST_PRINT(hls);

    HybridListSet<ExactBoxType>::ConstIterator iter=hls.begin();
    ARIADNE_TEST_EQUAL(*iter,HybridBoxType(loc2,spc2,bx3));
    ++iter;
    ARIADNE_TEST_EQUAL(*iter,HybridBoxType(loc2,spc2,bx4));
    ++iter;
    ARIADNE_TEST_EQUAL(*iter,HybridBoxType(loc2,spc2,bx5));
    ++iter;
    ARIADNE_TEST_EQUAL(*iter,HybridBoxType(loc1,spc1,bx1));
    ++iter;
    ARIADNE_TEST_EQUAL(*iter,HybridBoxType(loc1,spc1,bx2));
    ++iter;
    ARIADNE_TEST_ASSERT(iter==hls.end());

    ListSet<ExactBoxType> ls1; ls1.adjoin(bx1); ls1.adjoin(bx2);
    ListSet<ExactBoxType> ls2; ls2.adjoin(bx3); ls2.adjoin(bx4); ls2.adjoin(bx5);

    HybridListSet<ExactBoxType>::LocationsConstIterator loc_iter=hls.locations_begin();
    ARIADNE_TEST_EQUAL(loc_iter->first,loc2);
    ARIADNE_TEST_EQUAL(loc_iter->second.second,ls2);
    ++loc_iter;
    ARIADNE_TEST_EQUAL(loc_iter->first,loc1);
    ARIADNE_TEST_EQUAL(loc_iter->second.second,ls1);
    ++loc_iter;
    ARIADNE_TEST_ASSERT(loc_iter==hls.locations_end());

}


Int main() {
    TestHybridSet().test();
    return ARIADNE_TEST_FAILURES;
}

