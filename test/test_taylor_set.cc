/***************************************************************************
 *            test_taylor_set.cc
 *
 *  Copyright 2008  Pieter Collins
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
#include <iomanip>
#include "multi_index.h"
#include "taylor_model.h"
#include "zonotope.h"
#include "taylor_set.h"
#include "grid_set.h"
#include "graphics.h"

#include "test.h"
using namespace std;
using namespace Ariadne;

class TestTaylorSet {
  public:
    void test();
  private:
    void test_linearise();
    void test_discretise();
    void test_split();
    void test_recondition();
    void test_subdivide();
    void test_subsume();
    void test_draw();
};

GridTreeSet compute_outer_approximation(const TaylorSet& ts, const Grid& g, uint subd, uint depth);

void
TestTaylorSet::test()
{
    ARIADNE_TEST_CALL(test_linearise());
    ARIADNE_TEST_CALL(test_discretise());
    ARIADNE_TEST_CALL(test_split());
    ARIADNE_TEST_CALL(test_recondition());
    ARIADNE_TEST_CALL(test_subsume());
    ARIADNE_TEST_CALL(test_subdivide());
    ARIADNE_TEST_CALL(test_draw());
}


void
TestTaylorSet::test_subsume()
{
    TaylorSet ts=TaylorSet(2,2,2, 0.0,1.0,0.5,0.0,0.0,0.0, 0.25, 0.0,0.5,1.0,1.0,0.0,0.0, 0.375);
    TaylorSet cts=TaylorSet(2,4,2, 0.0, 1.0,0.5,0.25,0.0,  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,  0.0,
                                   0.0, 0.5,1.0,0.0,0.375, 1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,  0.0);
    ARIADNE_TEST_EQUAL(ts.subsume(),cts);
}


void
TestTaylorSet::test_linearise()
{
    TaylorSet ts(2,2,2, 0.0,1.0,0.25,0.0,0.0,0.0, 0.0, 0.0,0.5,1.0,1.0,0.0,0.0, 0.0);

    Zonotope z=zonotope(ts);
    Box b=ts.bounding_box();

    ARIADNE_TEST_EQUAL(ts.linearise(),TaylorSet(2,2,1, 0.0,1.0,0.25, 0.0, 0.0,0.5,1.0, 1.0));

    Box bounding_box=ts.bounding_box()+Vector<Interval>(2,Interval(-1,1));
    plot("test_taylor_set-linearise",PlanarProjectionMap(2,0,1),bounding_box,
         Colour(1,0,0),b,Colour(1,0,1),z,Colour(0,0,1),ts);
}


void
TestTaylorSet::test_discretise()
{
    TaylorSet ts(2,2,2, 0.0,1.0,0.25,0.0,0.0,0.0, 0.0, 0.0,0.5,1.0,1.0,0.0,0.0, 0.0);

    Grid grid(2);
    uint height(5);
    uint depth(4);
    GridTreeSet gts=ts.discretise(grid,depth);


    BooleanArray box_tree=make_binary_word("1010111110101011011100100010010110110010010110010101100100111011111100100011100100001111100100011100100001011111100100001000110110110110010010110010010110110010010110010011110010001110011011001000000");
    BooleanArray box_mask=make_binary_word("0000000101101000101001011010001011010111010110101110010110100111101011010111010110101010010110100000");
    GridTreeSet expected_box_discretisation=GridTreeSet(grid,height,box_tree,box_mask);

    BooleanArray affine_tree=make_binary_word("101011111010101101110010001001011011010110001011001010001110111100101100100010101110010000110111101100100001011110010001001111111000100011101000001110011011001000000");
    BooleanArray affine_mask=make_binary_word("00000001011010001101001010000100101101010110111010101101001001010010110010110100000");
    GridTreeSet expected_affine_discretisation=GridTreeSet(grid,height,affine_tree,affine_mask);

    BooleanArray rough_tree=make_binary_word("101011110101101101110001100010110010100011101101101100100101100100001110111101000101000111101100100110100001111100001011000101101000101111100100011011001000000");
    BooleanArray rough_mask=make_binary_word("00000001110100101000000101001010011011010111010110001010110100100010100110100000");
    GridTreeSet expected_rough_discretisation=GridTreeSet(grid,height,rough_tree,rough_mask);

    BooleanArray nonlinear_tree=make_binary_word("101011110101101101110001100010110010100011101101101100100101100100001110111101000101100010110100011110110010011100100001111100100010110000101111110001000111010011001000000");
    BooleanArray nonlinear_mask=make_binary_word("00000001110100101000000101001010011011001011011110101010001010011010001010010110100000");
    GridTreeSet expected_nonlinear_discretisation=GridTreeSet(grid,height,nonlinear_tree,nonlinear_mask);

    ARIADNE_TEST_BINARY_PREDICATE(subset,expected_nonlinear_discretisation,gts);
    ARIADNE_TEST_BINARY_PREDICATE(subset,gts,expected_box_discretisation);

    if(gts==expected_box_discretisation) {
        ARIADNE_TEST_WARN("Using box over approximation for discretisation of TaylorSet; large errors expected.");
    } else if(gts==expected_affine_discretisation) {
        ARIADNE_TEST_WARN("Using affine over approximation for discretisation of TaylorSet; moderate errors expected.");
    } else if(!subset(gts,expected_rough_discretisation)) {
        ARIADNE_TEST_WARN("Discretisation of TaylorSet may not be sufficiently accurate.");
    }

    Box bounding_box=ts.bounding_box()+Vector<Interval>(2,Interval(-1,1));
    plot("test_taylor_set-discretise",PlanarProjectionMap(2,0,1),bounding_box,
         Colour(1,0,0),gts,Colour(1,0,1),expected_nonlinear_discretisation,Colour(0,0,1),ts);

    plot("test_taylor_set-expected",bounding_box,
         Colour(1.00,0,0),expected_affine_discretisation,
         Colour(0.75,0,0),expected_rough_discretisation,
         Colour(0.55,0,0),expected_nonlinear_discretisation);
}


void
TestTaylorSet::test_split()
{
    TaylorSet ts(2,2,2, 0.0,1.0,0.25,0.0,0.0,0.0, 0.0, 0.0,0.5,1.0,1.0,0.0,0.0, 0.0);

    TaylorSet ts1,ts2,ts3,ts4,ts5,ts6;
    make_lpair(ts1,ts2)=ts.split();
    make_lpair(ts3,ts4)=ts2.split();
    make_lpair(ts5,ts6)=ts4.split();

    ARIADNE_TEST_EQUAL(ts.split().first,
        TaylorSet(2,2,2, -0.5,0.5,0.25,0.0,0.0,0.0, 0.0, +0.0,-0.25,1.0,0.25,0.0,0.0, 0.0));
    ARIADNE_TEST_EQUAL(ts.split().second,
        TaylorSet(2,2,2, +0.5,0.5,0.25,0.0,0.0,0.0, 0.0, +0.5,+0.75,1.0,0.25,0.0,0.0, 0.0));

    Box bounding_box=ts.bounding_box()+Vector<Interval>(2,Interval(-1,1));
    plot("test_taylor_set-split",bounding_box,
         Colour(0,0.0,1),ts1,
         Colour(0,0.4,1),ts3,
         Colour(0,0.8,1),ts5,
         Colour(0,0.9,1),ts6);

}

void
TestTaylorSet::test_recondition()
{
    TaylorSet ts(2,2,2, 0.,4.,0.,0.5,0.,0., 0.0,  0.,3.,1.,0.,0.,0., 0.0);
    TaylorSet rts=ts.recondition();
    ARIADNE_TEST_PRINT(ts);
    ARIADNE_TEST_PRINT(rts);
    Box bounding_box(2, -10.,+10., -10.,+10.);
    plot("test_taylor_set-recondition-1",bounding_box,Colour(0,1,1),rts,Colour(1,0,1),ts);

    //ts=TaylorSet(2,3,2, 0.,4.,0.,0.,0.5,0.,0.,0.,0., 0.0,
    //                    0.,3.,1.,0.,0.0,0.,0.,0.,0., 0.0);
    //rts=ts.recondition();
    //plot("test_taylor_set-recondition-2",bounding_box,Colour(0,1,1),rts,Colour(1,0,1),ts);

    TaylorSet sts1,sts2,srts1,srts2;;
    make_lpair(sts1,sts2)=ts.split();
    make_lpair(srts1,srts2)=rts.split();

    TaylorSet sts11,sts12,sts21,sts22,srts11,srts12,srts21,srts22;
    make_lpair(sts11,sts12)=sts1.split();
    make_lpair(sts21,sts22)=sts2.split();
    make_lpair(srts11,srts12)=srts1.split();
    make_lpair(srts21,srts22)=srts2.split();

   plot("test_taylor_set-split-norecondition",bounding_box,
         Colour(0,1,1),sts11,Colour(0,1,1),sts12,Colour(1,0,1),sts21,Colour(1,0,1),sts22);
    plot("test_taylor_set-split-recondition",bounding_box,
         Colour(0,1,1),srts11,Colour(0,1,1),srts12,Colour(1,0,1),srts21,Colour(1,0,1),srts22);
}

void
TestTaylorSet::test_subdivide()
{
}

void
TestTaylorSet::test_draw()
{
    TaylorSet ts(2,2,2, 0.0,1.0,0.25,1.0,0.0,0.0, 0.0, 0.0,0.5,1.0,0.0,0.0,0.0, 0.0);
    GridTreeSet gts=compute_outer_approximation(ts,Grid(2),6,6);

    plot("test_taylor_set-draw-1",gts.bounding_box(),Colour(1,0,1),gts,Colour(0,0,1),ts);
    plot("test_taylor_set-draw-2.png",ts);

}


GridTreeSet
compute_outer_approximation(const TaylorSet& set, const Grid& grid, uint subd, uint depth)
{
    assert(set.dimension()==2);
    assert(set.generators_size()==2);

    GridTreeSet gts(grid);

    Float rad=1.0/(1<<subd);
    Vector<Interval> v(2);
    for(Float x=-1.0; x!=1.0; x+=rad) {
        for(Float y=-1.0; y!=1.0; y+=rad) {
            v[0]=Interval(x,x+rad); v[1]=Interval(y,y+rad);
            gts.adjoin_outer_approximation(Box(evaluate(set.models(),v)),depth);
        }
    }
    gts.recombine();
    return gts;
}


int main() {
    TestTaylorSet().test();
    return ARIADNE_TEST_FAILURES;
}
