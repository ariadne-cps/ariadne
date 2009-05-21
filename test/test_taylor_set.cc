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
    TaylorSet ts1=TaylorSet(2,2,2, 0.0,1.0,0.5,0.0,0.0,0.0, 0.25, 0.0,0.5,1.0,1.0,0.0,0.0, 0.375);
    TaylorSet cts1=TaylorSet(2,4,2, 0.0, 1.0,0.5,0.25,0.0,  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,  0.0,
                                    0.0, 0.5,1.0,0.0,0.375, 1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,  0.0);
    ARIADNE_TEST_EQUAL(ts1.subsume(),cts1);

    TaylorSet ts2=TaylorSet(2,2,2, 0.0,1.0,0.5,0.0,0.0,0.0, 0.0, 0.0,0.5,1.0,1.0,0.0,0.0, 0.375);
    TaylorSet cts2=TaylorSet(2,3,2, 0.0, 1.0,0.5,0.0,  0.0,0.0,0.0,0.0,0.0,0.0,  0.0,
                                    0.0, 0.5,1.0,0.375, 1.0,0.0,0.0,0.0,0.0,0.0,  0.0);
    ARIADNE_TEST_EQUAL(ts2.subsume(),cts2);
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

void plot(const char* filename, const TaylorSet& set) {
    Figure fig;
    fig.set_bounding_box(set.bounding_box());
    fig.set_line_width(0.0);
    draw(fig,set);
    fig.write(filename);
}


void
TestTaylorSet::test_discretise()
{
    TaylorSet ts(2,2,2, 0.0,1.0,0.25,0.0,0.0,0.0, 0.03125, 0.0,0.5,1.0,1.0,0.0,0.0, 0.0);

    Grid grid(2);
    uint height(5);
    uint depth(6);
    GridTreeSet gts=ts.discretise(grid,depth);

    //ListSet<Box> bxfine=ts.discretise(1.0/(1<<(depth/2+3)));
    //GridTreeSet gtsfine=outer_approximation(bxfine,grid,depth+4);
    //ARIADNE_TEST_PRINT(gtsfine);

    ListSet<Box> bxls=ts.discretise(1.0/(1<<depth/2));
    GridTreeSet bxgts=outer_approximation(bxls,grid,depth);

    // The discretisation which would be obtained using the bounding box.
    BooleanArray box_tree=make_binary_word("1010111110101011011100100010010110110010010110010101100100111011111100100011100100001111100100011100100001011111100100001000110110110110010010110010010110110010010110010011110010001110011011001000000");
    BooleanArray box_mask=make_binary_word("0000000101101000101001011010001011010111010110101110010110100111101011010111010110101010010110100000");
    GridTreeSet expected_box_discretisation=GridTreeSet(grid,height,box_tree,box_mask);

    // The discretisation which would be obtained using an affine over-approximation.
    BooleanArray affine_tree=make_binary_word("101011111010101101110010001001011011010110001011001010001110111100101100100010101110010000110111101100100001011110010001001111111000100011101000001110011011001000000");
    BooleanArray affine_mask=make_binary_word("00000001011010001101001010000100101101010110111010101101001001010010110010110100000");
    GridTreeSet expected_affine_discretisation=GridTreeSet(grid,height,affine_tree,affine_mask);

    // A rough discretisation, which should be a superset of the computed discretisation (test for accuracy).
    BooleanArray rough_tree=make_binary_word("101011110101101101110001100010110010100011101101101100100101100100001110111101000101000111101100100110100001111100001011000101101000101111100100011011001000000");
    BooleanArray rough_mask=make_binary_word("00000001110100101000000101001010011011010111010110001010110100100010100110100000");
    GridTreeSet expected_rough_discretisation=GridTreeSet(grid,height,rough_tree,rough_mask);

    // A very good discretisation, which should be a subset of the computed discretisation (test for correctness).
    BooleanArray fine_tree=make_binary_word("101011110101101101110101010101100011010110110001110000111101010011000001110110111000100101101000010111111000100011101000001011011110111000100101101000010111110001000111010000101101101110001001011010000011101101101101010101011001001011011001001011110010001110010001011011011110010001110010001011001001011011011001000000111111111111011000110100000000110110111110000101010000011011011111100010001110100000010101010101101100100101100100111111101110100010101000101111000101100001011011010100001011110111000001011100010110001011011010000101101111011000011101000101100000110110110110110010010110000011111101101000011101100001110100001011011000010110111001010000101111010110001011001010001011101100001110010100011111111101110001001011010000101111100010001110100001011010101001110101001101000111100000111011011100100000111111100100000001111101101010011010001110101001011000001111111001000000000111111011011100100000111110101001011000010111100001011000001110111111000010110001011010001011010100011011010110001011001010001110110111000000111011100010110000111011010000111011000011101000101100010110110110110100000010111111011110110000111001010000011111010110001011001010001011101100001110010100000111111110110000111010000110101100010110010100111011000010111001010001011110110000111001010000111011011100100000111111100100000000000");
    BooleanArray fine_mask=make_binary_word("00000000000011000011011100010111100101010010001101011010011000101010010010101011010010001010100100000000000001010001010010110101100001011010110010100001010000001100111111111010101011111110101101001111111111010110100101101010110100101000101001101011010100101110111000010010110010001111101011011111011100100101101110110010101001100101101010100100101010010101001001010101101001111110111011001000000101000010100000011110110011101100001010000000000101000010101101011010110100110101101001001101011100101101011001010001010110110010001101101001101000001000000000100101010001001011010101001001010100011011010011001011010100100101010100010010101000010100001010000000000");
    GridTreeSet expected_fine_discretisation=GridTreeSet(grid,height,fine_tree,fine_mask);

    ARIADNE_ASSERT(subset(expected_fine_discretisation,expected_rough_discretisation));
    ARIADNE_ASSERT(subset(expected_rough_discretisation,expected_affine_discretisation));

    //ARIADNE_TEST_PRINT(gts.measure());

    ARIADNE_TEST_BINARY_PREDICATE(subset,expected_fine_discretisation,gts);
    ARIADNE_TEST_BINARY_PREDICATE(subset,gts,expected_rough_discretisation);

    if(gts==expected_box_discretisation) {
        ARIADNE_TEST_WARN("Using box over approximation for discretisation of TaylorSet; large errors expected.");
    } else if(gts==expected_affine_discretisation) {
        ARIADNE_TEST_WARN("Using affine over approximation for discretisation of TaylorSet; moderate errors expected.");
    } else if(!subset(gts,expected_rough_discretisation)) {
        ARIADNE_TEST_WARN("Discretisation of TaylorSet may not be sufficiently accurate.");
    }

    Box bounding_box=ts.bounding_box()+Vector<Interval>(2,Interval(-1,1));

    plot("test_taylor_set-discretise-box",PlanarProjectionMap(2,0,1),bounding_box,
         Colour(1,0,1),bxls);
    plot("test_taylor_set-discretise-box_grid",PlanarProjectionMap(2,0,1),bounding_box,
         Colour(1,0,1),bxgts);
    plot("test_taylor_set-discretise-grid",PlanarProjectionMap(2,0,1),bounding_box,
         //Colour(1,0,0),gts,Colour(1,0,1),expected_nonlinear_discretisation,Colour(0,0,1),ts);
         Colour(1,0,0),expected_rough_discretisation,Colour(1,0,1),gts,Colour(0.5,0,0),expected_fine_discretisation);

    plot("test_taylor_set-discretise-expected",bounding_box,
         Colour(1.00,0,0),expected_box_discretisation,
         Colour(0.85,0,0),expected_affine_discretisation,
         Colour(0.65,0,0),expected_rough_discretisation,
         Colour(0.50,0,0),expected_fine_discretisation);
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
    std::cerr<<"INCOMPLETE "<<std::flush;
    return ARIADNE_TEST_FAILURES;
}
