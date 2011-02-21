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

#include "test.h"

#include "vector.h"
#include "matrix.h"
#include "function.h"

#include "function_set.h"
#include "hybrid_set.h"
#include "graphics.h"

using namespace std;
using namespace Ariadne;


class TestHybridSet {
  public:
    void test();
  private:
    void test_hybrid_image_set();
    void test_hybrid_list_set();
    void test_hybrid_grid_set();
};

void
TestHybridSet::test()
{
    ARIADNE_TEST_CALL(test_hybrid_image_set());
    ARIADNE_TEST_CALL(test_hybrid_list_set());
    //ARIADNE_TEST_CALL(test_hybrid_grid_set());
}


void
TestHybridSet::test_hybrid_image_set()
{
    HybridImageSet his;
    Box bx=Box::unit_box(2);
    Matrix<Float> A=Matrix<Float>::identity(2);
    Vector<Float> b1=Vector<Float>::unit(2,0);
    Vector<Float> b2=Vector<Float>::unit(2,1);
    DiscreteLocation loc1(123);
    DiscreteLocation loc2(105);
    ImageSet ims1(bx,VectorAffineFunction(A,b1));
    ImageSet ims2(bx,VectorAffineFunction(A,b2));
    his.insert(make_pair(loc1,ims1));
    his[loc2]=ims2;

    HybridImageSet::locations_const_iterator iter=his.locations_begin();
    ARIADNE_ASSERT_EQUAL(iter->second,ims2);
    ++iter;
    ARIADNE_ASSERT_EQUAL(iter->second,ims1);
}

void
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
    draw(fig,hgts);
}


void
TestHybridSet::test_hybrid_list_set()
{
    HybridListSet<Box> hls;
    DiscreteLocation loc1(123);
    DiscreteLocation loc2(105);
    Box bx1=make_box("[0,1]");
    Box bx2=make_box("[2,3]");
    Box bx3=make_box("[1,2]x[2,3]");
    Box bx4=make_box("[4,5]x[5,6]");
    Box bx5=make_box("[6,7]x[8,9]");
    ARIADNE_TEST_EXECUTE(hls[loc1].adjoin(bx1));
    ARIADNE_TEST_FAIL(hls[loc1].adjoin(bx3)); // Should fail due to incompatible dimensions
    ARIADNE_TEST_EXECUTE(hls.insert(make_pair(loc2,ListSet<Box>(bx3))));
    ARIADNE_TEST_EXECUTE(hls[loc2].adjoin(bx4));
    ARIADNE_TEST_EXECUTE(hls.adjoin(HybridBox(loc1,bx2)));
    ARIADNE_TEST_EXECUTE(hls.adjoin(loc2,bx5));

    ARIADNE_TEST_PRINT(hls);

    HybridListSet<Box>::const_iterator iter=hls.begin();
    ARIADNE_TEST_EQUAL(*iter,HybridBox(loc2,bx3));
    ++iter;
    ARIADNE_TEST_EQUAL(*iter,HybridBox(loc2,bx4));
    ++iter;
    ARIADNE_TEST_EQUAL(*iter,HybridBox(loc2,bx5));
    ++iter;
    ARIADNE_TEST_EQUAL(*iter,HybridBox(loc1,bx1));
    ++iter;
    ARIADNE_TEST_EQUAL(*iter,HybridBox(loc1,bx2));
    ++iter;
    ARIADNE_TEST_ASSERT(iter==hls.end());

    ListSet<Box> ls1; ls1.adjoin(bx1); ls1.adjoin(bx2);
    ListSet<Box> ls2; ls2.adjoin(bx3); ls2.adjoin(bx4); ls2.adjoin(bx5);

    HybridListSet<Box>::locations_const_iterator loc_iter=hls.locations_begin();
    ARIADNE_TEST_EQUAL(loc_iter->first,loc2);
    ARIADNE_TEST_EQUAL(loc_iter->second,ls2);
    ++loc_iter;
    ARIADNE_TEST_EQUAL(loc_iter->first,loc1);
    ARIADNE_TEST_EQUAL(loc_iter->second,ls1);
    ++loc_iter;
    ARIADNE_TEST_ASSERT(loc_iter==hls.locations_end());

}


int main() {
    TestHybridSet().test();
    return ARIADNE_TEST_FAILURES;
}

