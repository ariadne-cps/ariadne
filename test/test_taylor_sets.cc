/***************************************************************************
 *            test_taylor_sets.cc
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
#include "function.h"
#include "constraint.h"
#include "taylor_set.h"
#include "affine_set.h"
#include "grid_set.h"
#include "graphics.h"

#include "test.h"
using namespace std;
using namespace Ariadne;

class TestTaylorImageSet {
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

GridTreeSet compute_outer_approximation(const TaylorImageSet& ts, const Grid& g, uint subd, uint depth);

void
TestTaylorImageSet::test()
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
TestTaylorImageSet::test_subsume()
{
    TaylorImageSet ts1=TaylorImageSet(2,2,2, 0.0,1.0,0.5,0.0,0.0,0.0, 0.25, 0.0,0.5,1.0,1.0,0.0,0.0, 0.375);
    TaylorImageSet cts1=TaylorImageSet(2,4,2, 0.0, 1.0,0.5,0.25,0.0,  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,  0.0,
                                    0.0, 0.5,1.0,0.0,0.375, 1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,  0.0);
    ARIADNE_TEST_EQUAL(ts1.subsume(),cts1);

    TaylorImageSet ts2=TaylorImageSet(2,2,2, 0.0,1.0,0.5,0.0,0.0,0.0, 0.0, 0.0,0.5,1.0,1.0,0.0,0.0, 0.375);
    TaylorImageSet cts2=TaylorImageSet(2,3,2, 0.0, 1.0,0.5,0.0,  0.0,0.0,0.0,0.0,0.0,0.0,  0.0,
                                    0.0, 0.5,1.0,0.375, 1.0,0.0,0.0,0.0,0.0,0.0,  0.0);
    ARIADNE_TEST_EQUAL(ts2.subsume(),cts2);
}


void
TestTaylorImageSet::test_linearise()
{
/*
    TaylorImageSet ts(2,2,2, 0.0,1.0,0.25,0.0,0.0,0.0, 0.0, 0.0,0.5,1.0,1.0,0.0,0.0, 0.0);

    Zonotope z=zonotope(ts);
    Box b=ts.bounding_box();

    ARIADNE_TEST_EQUAL(ts.linearise(),TaylorImageSet(2,2,1, 0.0,1.0,0.25, 0.0, 0.0,0.5,1.0, 1.0));

    Box bounding_box=ts.bounding_box()+Vector<Interval>(2,Interval(-1,1));
    plot("test_taylor_set-linearise",PlanarProjectionMap(2,0,1),bounding_box,
         Colour(1,0,0),b,Colour(1,0,1),z,Colour(0,0,1),ts);
*/
}

void plot(const char* filename, const TaylorImageSet& set) {
    Figure fig;
    fig.set_bounding_box(set.bounding_box());
    fig.set_line_width(0.0);
    draw(fig,set);
    fig.write(filename);
}


void
TestTaylorImageSet::test_discretise()
{
    Grid grid(2);
    uint height(5);
    uint depth(4);

    TaylorImageSet ts(2,2,2, 0.0,1.0,0.25,0.0,0.0,0.0, 0.015625, 0.0,0.5,1.0,1.0,0.0,0.0, 0.015625);
    GridTreeSet grid_discretisation=ts.discretise(grid,depth);

    TaylorImageSet tsz(2,2,2, 0.0,1.0,0.25,0.0,0.0,0.0, 0.031251, 0.0,0.5,1.0,1.0,0.0,0.0, 0.00001);
    ListSet<Box> rough_boxes_approximation=ts.discretise(1.0/(1<<(depth-2)));
    ListSet<Box> boxes_approximation=ts.discretise(1.0/(1<<(depth-2)));
    ListSet<Box> fine_boxes_approximation=ts.discretise(1.0/(1<<(depth+3)));
    GridTreeSet boxes_discretisation=outer_approximation(boxes_approximation,grid,depth-1);

    //GridTreeSet computed_box_discretisation=outer_approximation(ts.bounding_box(),grid,depth-3);
    //ARIADNE_TEST_PRINT(computed_box_discretisation);
    //GridTreeSet computed_affine_discretisation=outer_approximation(zonotope_over_approximation,grid,depth-2);
    //ARIADNE_TEST_PRINT(computed_affine_discretisation);
    //GridTreeSet computed_rough_discretisation=outer_approximation(rough_boxes_approximation,grid,depth-2);
    //ARIADNE_TEST_PRINT(computed_rough_discretisation);
    //GridTreeSet computed_fine_discretisation=outer_approximation(fine_boxes_approximation,grid,depth+2);
    //ARIADNE_TEST_PRINT(computed_fine_discretisation);

    // The discretisation which would be obtained using the bounding box.
    BooleanArray box_tree=make_binary_word("101011111010101010101100100101101101110010011001001011100100110010010110111001001100100101110010011001001110110110110010010110010010110110010010110010010110110110010000111111110010011001001110010011001001111001001100100111001001100100111110010011001001110010011001001111001001100100111001001100100111111001001100100111001001100100111100100110010011100100110010011111001001100100111001001100100111100100110010011100100110010011111110010011001000111100100110010000111111001001100100011110010011001000000");
    BooleanArray box_mask=make_binary_word("00000000010100011111111011111111001111111101111101000001010010100010100101000010100111111111111111111111111111111111111111111111111111111111111111111111111111111111111101011111010111111111111111111111010111110101111111101111111100111111110111110100000");
    GridTreeSet expected_box_discretisation=GridTreeSet(grid,height,box_tree,box_mask);

    // The discretisation which would be obtained using an affine over-approximation.
    BooleanArray affine_tree=make_binary_word("101011111010101101110010001001011011010110001011001010001110111100101100100010101110010000110111101100100001011110010001001111111000100011101000001110011011001000000");
    BooleanArray affine_mask=make_binary_word("00000001011010001101001010000100101101010110111010101101001001010010110010110100000");
    GridTreeSet expected_affine_discretisation=GridTreeSet(grid,height,affine_tree,affine_mask);

    // A rough discretisation, which should be a superset of the computed discretisation (test for accuracy).
    BooleanArray rough_tree=
        make_binary_word("1010111110101010101010100101101101011000101100101000111011011011001001011001001011011000001110111101000101000111101100100110100001111100001011000101101000101111100100011011001000000");
    BooleanArray rough_mask=
        make_binary_word("0000000000100011010010100000010100101000100011011010111010110001010110100100010100110100000");
    GridTreeSet expected_rough_discretisation=GridTreeSet(grid,height,rough_tree,rough_mask);

    // A very good discretisation, which should be a subset of the computed discretisation (test for correctness).
    BooleanArray fine_tree=make_binary_word("101011111010101010101010101010101010110010010110110111101010101010110010010110110010010110010011110101010011011011001001011001001100110110010001110110010011001100000111011011001001011001010001011110110000111001010000101101111010110001011001010001011101100001110010100010110110101100010110010100001110110110110110110110110010010110010011010011011001000110101010110010011110010001110000101111011011001001011001001011011001001000111011001001011001000101101111011001001011001000111011001001011001000101111011001001011000101100101001011011001001011001001011011011011011001001011000000011110110110111101100001110010100000011111111011001001011010001011110010001110001001011011010110010000101111110110010010110010001110110010001001011001010010110110010010110001011011011001001011010110010000101111011110110010010110010001110110010001000111011001001011010001011110010000010111101101110001001011010000111011110110010010110010000011101110110010001001011010110010000101101111101100100010001110101100100001011011101100100010010110101100100011111111101011000101100101000101110110000111001010001010101110101100100100111011001001011001000111011010110010000111110110010010110010001111001100000011110101110101100100100111011001001011001000111110100100111001000111101100100000111110000001111101010010110010011110110010000011101100100101100100000011111111101100100101100100011101100001110010100101100100011101101100100101100100101101100001011110010101100100101101100100101100100001110111101101100100001110010000011111110100011100100001111100010001110010101100100101101100100101100100011101101101100100000111110010001110010100101100100101101101100100001011011011110010101100100010101101100100101100100001011111111111101100100010001110101100100000111011101100100010010110101100100000111011111011001000100011101011001000010110111011001000100101101011001000001111011111100010001110100000111111101100100010001110010000111110010001110010000111110111001000011111001000111001000001111111000000000000");
    BooleanArray fine_mask=make_binary_word("00000000000000000101000000000010100010100101000010001010010101001011001010101111000101001010011101101010110001101001010101101101010100011010010100000000000101001010010010110000010101011011100001010010100010101100101001011000010100101100101001011000101001011010100010100101000000101001000001111101101010111111010110110010110101011100101110001010010110010101011010100010100101111010110010111000010100101100101010011010110110010110011101010110111000101001011001001010101100101101100101010010010110110010101011001011011010010101011011010101111111010101101011010000110100011010110100101000001111101010110101101001101010100110100001000001110110101101000011010110100000110101101001101101010110100111010110101110110010110101110101101000111101011010111101001010001010110101101011101011010011110101110101101010110101111010110000101101010111010110100000010101001001011001001010101100101100100101010010010110110010101011001011001101011010011001010100101000101001010000101000101001010000100000000000");
    GridTreeSet expected_fine_discretisation=GridTreeSet(grid,height,fine_tree,fine_mask);

    Box bounding_box=ts.bounding_box()+Vector<Interval>(2,Interval(-1,1));

    plot("test_taylor_set-discretise-box",PlanarProjectionMap(2,0,1),bounding_box,
         Colour(1,0,1),boxes_approximation);
    plot("test_taylor_set-discretise-box_grid",PlanarProjectionMap(2,0,1),bounding_box,
         Colour(1,0,1),boxes_discretisation);
    plot("test_taylor_set-discretise-grid",PlanarProjectionMap(2,0,1),bounding_box,
         Colour(1,0,0),expected_rough_discretisation,
         Colour(1,0,1),grid_discretisation,
         Colour(0.5,0,0),expected_fine_discretisation);

    plot("test_taylor_set-discretise-expected",bounding_box,
         Colour(1.00,0,0),expected_box_discretisation,
         //Colour(0.85,0,0),expected_affine_discretisation,
         Colour(0.65,0,0),expected_rough_discretisation,
         Colour(0.50,0,0),expected_fine_discretisation);

    //plot("test_taylor_set-discretise-computed",bounding_box,
    //     Colour(1.00,0,0),computed_box_discretisation,
    //     Colour(0.85,0,0),computed_affine_discretisation);
    //     Colour(0.65,0,0),computed_rough_discretisation,
    //     Colour(0.50,0,0),computed_fine_discretisation);

    ARIADNE_ASSERT(subset(expected_fine_discretisation,expected_rough_discretisation));
    //ARIADNE_ASSERT(subset(expected_rough_discretisation,expected_affine_discretisation));

    ARIADNE_TEST_BINARY_PREDICATE(subset,expected_fine_discretisation,grid_discretisation);
    ARIADNE_TEST_BINARY_PREDICATE(subset,grid_discretisation,expected_rough_discretisation);

    //ARIADNE_TEST_BINARY_PREDICATE(subset,computed_fine_discretisation,gts);
    //ARIADNE_TEST_BINARY_PREDICATE(subset,gts,computed_rough_discretisation);

    if(grid_discretisation==expected_box_discretisation) {
        ARIADNE_TEST_WARN("Using box over approximation for discretisation of TaylorImageSet; large errors expected.");
    } else if(grid_discretisation==expected_affine_discretisation) {
        ARIADNE_TEST_WARN("Using affine over approximation for discretisation of TaylorImageSet; moderate errors expected.");
    } else if(!subset(grid_discretisation,expected_rough_discretisation)) {
        ARIADNE_TEST_WARN("Discretisation of TaylorImageSet may not be sufficiently accurate.");
    }



}


void
TestTaylorImageSet::test_split()
{
    TaylorImageSet ts(2,2,2, 0.0,1.0,0.25,0.0,0.0,0.0, 0.0, 0.0,0.5,1.0,1.0,0.0,0.0, 0.0);

    TaylorImageSet ts1,ts2,ts3,ts4,ts5,ts6;
    make_lpair(ts1,ts2)=ts.split();
    make_lpair(ts3,ts4)=ts2.split();
    make_lpair(ts5,ts6)=ts4.split();

    ARIADNE_TEST_EQUAL(ts.split().first,
        TaylorImageSet(2,2,2, -0.5,0.5,0.25,0.0,0.0,0.0, 0.0, +0.0,-0.25,1.0,0.25,0.0,0.0, 0.0));
    ARIADNE_TEST_EQUAL(ts.split().second,
        TaylorImageSet(2,2,2, +0.5,0.5,0.25,0.0,0.0,0.0, 0.0, +0.5,+0.75,1.0,0.25,0.0,0.0, 0.0));

    Box bounding_box=ts.bounding_box()+Vector<Interval>(2,Interval(-1,1));
    plot("test_taylor_set-split",bounding_box,
         Colour(0,0.0,1),ts1,
         Colour(0,0.4,1),ts3,
         Colour(0,0.8,1),ts5,
         Colour(0,0.9,1),ts6);

    // Test split with an error term
    ts=TaylorImageSet(2,2,2, 0.5,1.0,0.25,0.0,0.0,0.0, 2.5, 0.0,0.5,1.0,1.0,0.0,0.0, 1.0);
    ARIADNE_TEST_EQUAL(ts.split().first,
        TaylorImageSet(2,2,2, -0.75,1.0,0.25,0.0,0.0,0.0, 1.25, 0.0,0.5,1.0,1.0,0.0,0.0, 1.0));
    ARIADNE_TEST_EQUAL(ts.split().second,
        TaylorImageSet(2,2,2, 1.75,1.0,0.25,0.0,0.0,0.0, 1.25, 0.0,0.5,1.0,1.0,0.0,0.0, 1.0));


}

void
TestTaylorImageSet::test_recondition()
{
    TaylorImageSet ts(2,2,2, 0.,4.,0.,0.5,0.,0., 0.0,  0.,3.,1.,0.,0.,0., 0.0);
    TaylorImageSet rts=ts.recondition();
    ARIADNE_TEST_PRINT(ts);
    ARIADNE_TEST_PRINT(rts);
    Box bounding_box(2, -10.,+10., -10.,+10.);
    plot("test_taylor_set-recondition-1",bounding_box,Colour(0,1,1),rts,Colour(1,0,1),ts);

    //ts=TaylorImageSet(2,3,2, 0.,4.,0.,0.,0.5,0.,0.,0.,0., 0.0,
    //                    0.,3.,1.,0.,0.0,0.,0.,0.,0., 0.0);
    //rts=ts.recondition();
    //plot("test_taylor_set-recondition-2",bounding_box,Colour(0,1,1),rts,Colour(1,0,1),ts);

    TaylorImageSet sts1,sts2,srts1,srts2;;
    make_lpair(sts1,sts2)=ts.split();
    make_lpair(srts1,srts2)=rts.split();

    TaylorImageSet sts11,sts12,sts21,sts22,srts11,srts12,srts21,srts22;
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
TestTaylorImageSet::test_subdivide()
{
}

void
TestTaylorImageSet::test_draw()
{
    TaylorImageSet ts(2,2,2, 0.0,1.0,0.25,1.0,0.0,0.0, 0.0, 0.0,0.5,1.0,0.0,0.0,0.0, 0.0);
    GridTreeSet gts=compute_outer_approximation(ts,Grid(2),6,6);

    plot("test_taylor_set-draw-1",gts.bounding_box(),Colour(1,0,1),gts,Colour(0,0,1),ts);

}


GridTreeSet
compute_outer_approximation(const TaylorImageSet& set, const Grid& grid, uint subd, uint depth)
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

/*
typedef Vector<Float> FloatVector;
List<FloatVector> scatter(const Box& D, const VectorFunction& f, const List<ScalarFunction>& g) {
    ARIADNE_ASSERT(D.size()==2);
    List<FloatVector> res;
    const uint N=16;
    FloatVector s(2);
    for(uint i=0; i<=N; ++i) {
        s[0]=((N-i)*D[0].lower()+i*D[0].upper())/N;
        for(uint j=0; j<=N; ++j) {
            s[1]=((N-j)*D[1].lower()+j*D[1].upper())/N;
            bool valid=true;
            for(uint k=0; k!=g.size(); ++k) {
                const ScalarFunction& gk=g[k];
                if(gk(s)>0.0) { valid=false; break; }
            }
            if(valid) { res.append(s); }
        }
    }
    return res;
}
*/

class TestTaylorConstrainedImageSet
{
  private:
    Figure figure;
  public:
    void test() {
        figure.set_bounding_box(Box(2, -4.0,+4.0, -4.0,+4.0));
        ARIADNE_TEST_CALL(test_affine_approximation());
        ARIADNE_TEST_CALL(test_disjoint());
        ARIADNE_TEST_CALL(test_outer_approximation());
        ARIADNE_TEST_CALL(test_draw());
    }

    void test_affine_approximation() {
        ScalarFunction s=ScalarFunction::coordinate(2,0);
        ScalarFunction t=ScalarFunction::coordinate(2,1);
        ScalarFunction x=ScalarFunction::coordinate(2,0);
        ScalarFunction y=ScalarFunction::coordinate(2,1);
        //Box domain(2, -0.5,1.5, -0.5,0.5);
        Box domain(2, -1.0,+1.0, -1.0,+1.0);

        TaylorConstrainedImageSet set(domain,(s+0.5*t,t+sqr(s)));
        set.new_negative_constraint(s+t-0.25);
        GridTreeSet paving(set.dimension());
        set.subdivision_adjoin_outer_approximation_to(paving,5u);
        paving.recombine();

        figure.set_fill_colour(1.0,0.0,0.0);
        figure.draw(set.affine_over_approximation());
        figure.set_fill_colour(0.25,0.75,0.75);
        figure.draw(paving);
        figure.set_fill_colour(0.0,1.0,1.0);
        figure.draw(set);
        figure.set_fill_colour(1.0,1.0,0.0);
        figure.draw(set.affine_approximation());
        figure.write("test_taylor_set-affine_approximation");
        figure.clear();
    }


    void test_disjoint() {
        ScalarFunction s=ScalarFunction::coordinate(2,0);
        ScalarFunction t=ScalarFunction::coordinate(2,1);
        ScalarFunction x=ScalarFunction::coordinate(2,0);
        ScalarFunction y=ScalarFunction::coordinate(2,1);
        Box domain(2, -0.5,1.5, -0.5,0.5);

        TaylorConstrainedImageSet set(domain,(s,t+sqr(s)));
        set.new_negative_constraint(s+t);
        ARIADNE_TEST_PRINT(set);
        Box bx1(2, 1.0,2.0, 0.0,1.0);
        ARIADNE_TEST_PRINT(bx1);
        ARIADNE_TEST_ASSERT(set.disjoint(bx1));
    }

    void test_outer_approximation() {
        ScalarFunction s=ScalarFunction::coordinate(2,0);
        ScalarFunction t=ScalarFunction::coordinate(2,1);
        ScalarFunction x=ScalarFunction::coordinate(2,0);
        ScalarFunction y=ScalarFunction::coordinate(2,1);
        //Box domain(2, -0.5,1.5, -0.5,0.5);
        Box domain(2, -1.0,+1.0, -1.0,+1.0);
        uint accuracy = 4u;
        uint high_accuracy = accuracy+1u;

        TaylorConstrainedImageSet set(domain,(s+0.5*t,t+sqr(s)));
        //TaylorConstrainedImageSet set(domain,(s,t));
        set.new_negative_constraint(s+t-0.25);
        ARIADNE_TEST_PRINT(set);

        GridTreeSet subdivision_paving(set.dimension());
        set.subdivision_adjoin_outer_approximation_to(subdivision_paving,high_accuracy);
        subdivision_paving.recombine();
        ARIADNE_TEST_EQUALS(subdivision_paving.measure(),3.375);

        figure.set_bounding_box(Box(2, -2.0,+2.0, -2.0,+2.0));
        figure.set_fill_colour(0.0,0.5,1.0);
        figure.draw(subdivision_paving);
        figure.write("test_taylor_set-subdivision_paving");
        figure.clear();


        GridTreeSet affine_paving(set.dimension());
        set.affine_adjoin_outer_approximation_to(affine_paving,accuracy);
        affine_paving.recombine();

        ARIADNE_TEST_ASSERT(subset(subdivision_paving,affine_paving));
        if(affine_paving.measure()/subdivision_paving.measure() > 1.10) {
            ARIADNE_TEST_WARN("TaylorConstrainedImageSet::affine_adjoin_outer_approximation_to(...) yields poor approximation");
        }

        GridTreeSet constraint_paving(set.dimension());
        set.constraint_adjoin_outer_approximation_to(constraint_paving,accuracy);
        constraint_paving.recombine();
        figure.set_fill_colour(1.0,0.0,0.0);
        figure.draw(constraint_paving);
        figure.set_fill_colour(0.0,0.5,1.0);
        figure.draw(subdivision_paving);
        figure.write("test_taylor_set-constraint_paving");
        figure.clear();

        ARIADNE_TEST_ASSERT(subset(subdivision_paving,constraint_paving));
        ARIADNE_TEST_LESS(constraint_paving.measure()/subdivision_paving.measure() , 1.10);
    }

    void test_draw(const std::string& str, const TaylorConstrainedImageSet& set, uint acc) {
        ARIADNE_TEST_PRINT(set);
        figure.clear();
        GridTreeSet paving(set.dimension());
        set.adjoin_outer_approximation_to(paving,acc);
        figure.set_fill_opacity(1.0);
        figure.set_fill_colour(red);
        figure.draw(paving);
        figure.set_fill_colour(green);
        figure.set_fill_opacity(0.5);
        figure.draw(set);
        figure.write((std::string("test_taylor_set-draw-")+str).c_str());
        figure.clear();
    }

    void test_draw() {
        ScalarFunction s=ScalarFunction::coordinate(2,0);
        ScalarFunction t=ScalarFunction::coordinate(2,1);
        ScalarFunction x=ScalarFunction::coordinate(2,0);
        ScalarFunction y=ScalarFunction::coordinate(2,1);
        uint accuracy = 3u;
        uint high_accuracy = 5u;

        // Test drawing with outer approximation computed using domain subdivision
        Box domain(2, -1.0,+1.0, -1.0,+1.0);
        TaylorConstrainedImageSet set(domain,(s+0.5*t,t+sqr(s)));
        set.new_negative_constraint(s+t-0.25);
        GridTreeSet paving(set.dimension());
        set.subdivision_adjoin_outer_approximation_to(paving,high_accuracy);
        paving.recombine();

        figure.set_bounding_box(Box(2, -2.0,+2.0, -2.0,+2.0));
        figure.set_fill_opacity(1.0);
        figure.set_fill_colour(0.0,0.5,1.0);
        figure.draw(paving);
        figure.set_fill_colour(green);
        figure.set_fill_opacity(0.5);
        figure.draw(set);
        figure.write("test_taylor_set-draw");
        figure.clear();


        // Draw a variety of shapes
        figure.set_bounding_box(Box(2, -4.0,+4.0, -4.0,+4.0));
        test_draw("box",TaylorConstrainedImageSet(Box(2,-1.01,1.01,-1.01,1.01),(s,t)),accuracy);
        test_draw("polytope",TaylorConstrainedImageSet(Box(2,-2.05,2.05,-1.05,1.05),(s,t),(s+t<=1.5)),accuracy);
        test_draw("dome",TaylorConstrainedImageSet(Box(2,-1.0,1.0,-1.0,1.0),(s,t),(s*s+t<=0.251)),accuracy);
        test_draw("disc",TaylorConstrainedImageSet(Box(2,-1.0,1.0,-1.0,1.0),(s,t),(s*s+t*t<=0.751)),accuracy);
        test_draw("parallelotope",TaylorConstrainedImageSet(Box(2,-1.0,1.0,-1.0,1.0),(2*s+t,s+t)),accuracy);
        test_draw("ellipse",TaylorConstrainedImageSet(Box(2,-1.0,1.0,-1.0,1.0),(2*s+t,s+t),(s*s+t*t<=0.75)),accuracy);
        test_draw("concave",TaylorConstrainedImageSet(Box(2,-1.01,1.01,-1.01,1.01),(s,0.25*s*s+t),(2*s+0.25*s*s+t-0.5<=0)),accuracy);

        test_draw("empty",TaylorConstrainedImageSet(Box(2,-1.01,1.01,-1.01,1.01),(s,0.25*s*s+t),(1.25+s+t*t<=0)),accuracy);
        test_draw("curve",TaylorConstrainedImageSet(Box(2,-1.01,1.01,-1.01,1.01),(s,t),(s*s-t==0)),accuracy);


    }
};



int main() {
    //TestTaylorImageSet().test();
    TestTaylorConstrainedImageSet().test();
    std::cerr<<"INCOMPLETE "<<std::flush;
    return ARIADNE_TEST_FAILURES;
}
