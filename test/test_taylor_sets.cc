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



class TestTaylorConstrainedImageSet
{
  private:
    Figure figure;
  public:
    void test() {
        figure.set_bounding_box(Box(2, -4.0,+4.0, -4.0,+4.0));
        ARIADNE_TEST_CALL(test_outer_approximation()); return;
        ARIADNE_TEST_CALL(test_disjoint());
        //ARIADNE_TEST_CALL(test_discretise());
        ARIADNE_TEST_CALL(test_affine_approximation());
        ARIADNE_TEST_CALL(test_outer_approximation());
        ARIADNE_TEST_CALL(test_split());
        ARIADNE_TEST_CALL(test_recondition());
        ARIADNE_TEST_CALL(test_subsume());
        ARIADNE_TEST_CALL(test_draw());
    }

    void test_subsume()
    {
        Interval e(-1.0,+1.0);
        IntervalVector d1=IntervalVector::unit_box(2);
        VectorTaylorFunction x=VectorTaylorFunction::identity(d1);
        VectorTaylorFunction f1( ( x[0]+0.5*x[1]+0.25*e , 0.5*x[0]+x[1]+x[0]*x[0]+0.375*e ) );
        IntervalVector d2=IntervalVector::unit_box(4);
        VectorTaylorFunction y=VectorTaylorFunction::identity(d2);
        VectorTaylorFunction f2( ( y[0]+0.5*y[1]+0.25*y[2] , 0.5*y[0]+y[1]+y[0]*y[0]+0.375*y[3] ) );
        TaylorConstrainedImageSet ts1=TaylorConstrainedImageSet(f1);
        TaylorConstrainedImageSet cts1=TaylorConstrainedImageSet(f2);
        //ARIADNE_TEST_EQUAL(ts1.subsume(),cts1);

        f1=VectorTaylorFunction( ( x[0]+0.5*x[1] , 0.5*x[0]+x[1]+x[0]*x[0]+0.375*e ) );
        d2=IntervalVector::unit_box(3);
        y=VectorTaylorFunction::identity(d2);
        f2=VectorTaylorFunction( ( y[0]+0.5*y[1] , 0.5*y[0]+y[1]+y[0]*y[0]+0.375*y[2] ) );
        TaylorConstrainedImageSet ts2=TaylorConstrainedImageSet(f1);
        TaylorConstrainedImageSet cts2=TaylorConstrainedImageSet(f2);
        //ARIADNE_TEST_EQUAL(ts2.subsume(),cts2);
    }




    void test_discretise()
    {
        Grid grid(2);
        uint depth(4);

        Interval e(-1.0,+1.0);
        VectorTaylorFunction x=VectorTaylorFunction::identity(IntervalVector::unit_box(2));

        TaylorConstrainedImageSet ts( VectorTaylorFunction(( x[0]+0.25*x[1]+0.015625*e , 0.5*x[0]+x[1]+x[0]*x[0]+0.015625*e )) );
        TaylorConstrainedImageSet tsz( VectorTaylorFunction(( x[0]+0.25*x[1]+0.031251*e, 0.5*x[0]+x[1]+x[0]*x[0]+0.00001*e )) );

        GridTreeSet box_discretisation=outer_approximation(ts.bounding_box(),grid,depth-3);
        ARIADNE_TEST_PRINT(box_discretisation);
        AffineSet as=ts.affine_over_approximation();
        GridTreeSet affine_discretisation(grid);
        as.adjoin_outer_approximation_to(affine_discretisation,depth-2);
        //GridTreeSet affine_discretisation=outer_approximation(as,grid,depth-2);
        ARIADNE_TEST_PRINT(affine_discretisation);
        GridTreeSet rough_discretisation=outer_approximation(ts,grid,depth-2);
        ARIADNE_TEST_PRINT(rough_discretisation);
        GridTreeSet fine_discretisation=outer_approximation(ts,grid,depth+2);
        ARIADNE_TEST_PRINT(fine_discretisation);

        Box bounding_box=ts.bounding_box();
        plot("test_taylor_set-discretise",bounding_box,
             Colour(1.00,0,0),box_discretisation,
             Colour(0.85,0,0),affine_discretisation,
             Colour(0.65,0,0),rough_discretisation,
             Colour(0.50,0,0),fine_discretisation);

        ARIADNE_TEST_BINARY_PREDICATE(subset,fine_discretisation,rough_discretisation);
        ARIADNE_TEST_BINARY_PREDICATE(subset,rough_discretisation,box_discretisation);

    }


    void test_split()
    {
        Interval e(-1.0,+1.0);
        Vector<IntervalTaylorModel> s=IntervalTaylorModel::variables(2);
/*

        IntervalVector dom=IntervalVector(2, -1.0,+1.0, -1.0,1.0);
        TaylorConstrainedImageSet ts(VectorTaylorFunction( dom,(s[0]+0.25*s[1],0.5*s[0]+s[1]+s[0]*s[0])) );

        TaylorConstrainedImageSet ts1,ts2,ts3,ts4,ts5,ts6;
        make_lpair(ts1,ts2)=ts.split(0);
        make_lpair(ts3,ts4)=ts2.split(1);
        make_lpair(ts5,ts6)=ts4.split(0);

        IntervalVector dom1=IntervalVector(2, -1.0,0.0, -1.0,1.0);
        IntervalVector dom2=IntervalVector(2, 0.0,+1.0, -1.0,1.0);
        ARIADNE_TEST_EQUAL(ts.split().first,
            TaylorConstrainedImageSet(VectorTaylorFunction( dom1, (-0.5+0.5*s[0]+0.25*s[1], -0.25*s[0]+s[1]+0.25*s[0]*s[0]) ) ) );
        ARIADNE_TEST_EQUAL(ts.split().second,
            TaylorConstrainedImageSet(VectorTaylorFunction( dom2, (+0.5+0.5*s[0]+0.25*s[1], +0.5+0.75*s[0]+s[1]+0.25*s[0]*s[0]) ) ) );

        Box bounding_box=ts.bounding_box()+Vector<Interval>(2,Interval(-1,1));
        plot("test_taylor_set-split",bounding_box,
            Colour(0,0.0,1),ts1,
            Colour(0,0.4,1),ts3,
            Colour(0,0.8,1),ts5,
            Colour(0,0.9,1),ts6);

        // Test split with an error term
        ts=TaylorConstrainedImageSet(2,2,2, 0.5,1.0,0.25,0.0,0.0,0.0, 2.5, 0.0,0.5,1.0,1.0,0.0,0.0, 1.0);
        ARIADNE_TEST_EQUAL(ts.split().first,
            TaylorConstrainedImageSet(2,2,2, -0.75,1.0,0.25,0.0,0.0,0.0, 1.25, 0.0,0.5,1.0,1.0,0.0,0.0, 1.0));
        ARIADNE_TEST_EQUAL(ts.split().second,
            TaylorConstrainedImageSet(2,2,2, 1.75,1.0,0.25,0.0,0.0,0.0, 1.25, 0.0,0.5,1.0,1.0,0.0,0.0, 1.0));

*/
    }

    void test_recondition() {
/*
        Interval e(-1.0,+1.0);
        VectorTaylorFunction x=VectorTaylorFunction::identity(IntervalVector::unit_box(2));
        TaylorConstrainedImageSet ts( ( 4.0*x[0]+0.5*x[0]*x[0] , 3.0*x[0]+x[1] ) );
        TaylorConstrainedImageSet rts=ts.recondition();
        ARIADNE_TEST_PRINT(ts);
        ARIADNE_TEST_PRINT(rts);
        Box bounding_box(2, -10.,+10., -10.,+10.);
        plot("test_taylor_set-recondition-1",bounding_box,Colour(0,1,1),rts,Colour(1,0,1),ts);
*/
    }

    void test_affine_approximation() {
        RealScalarFunction s=RealScalarFunction::coordinate(3,0);
        RealScalarFunction t=RealScalarFunction::coordinate(3,1);
        RealScalarFunction u=RealScalarFunction::coordinate(3,2);
        RealScalarFunction x=RealScalarFunction::coordinate(2,0);
        RealScalarFunction y=RealScalarFunction::coordinate(2,1);
        Real e(Interval(-1,+1));
        //Box domain(2, -0.5,1.5, -0.5,0.5);
        Box domain(3, -1.0,+1.0, -1.0,+1.0, -1.0,+1.0);

        //TaylorConstrainedImageSet set(domain,(s+0.5*t+(1.0/256)*s*t+t*t,t+0.375*sqr(s)+0.25*e),(0.625*s*s+t-1<=0, 0.25*s+0.5*t==0, -0.25*s+0.5*t+0.125*e==0));

        //TaylorConstrainedImageSet set(domain,(s+0.5*t,t),(s+t<=1.0,t*t-s<=1.0,t+u+0.0625*e==0)); //,(t-4<=0)); //, -0.25*s+0.5*t+0.125*e==0));//, 0.25*s+0.5*t==0
        TaylorConstrainedImageSet set(domain,(s+t+0.25*s*s-0.25*s*t+0.125*t*t,-s+t),(u<=2,-s+t+u-1+0*e/16==0)); //,(t-4<=0)); //, -0.25*s+0.5*t+0.125*e==0));//, 0.25*s+0.5*t==0

        ARIADNE_TEST_PRINT(set);
        AffineSet affine_set=set.affine_over_approximation();
        ARIADNE_TEST_PRINT(affine_set);

        Vector< Affine<Float> > a=Affine<Float>::variables(4);
        AffineSet expected_affine_set(IntervalVector::unit_box(4),(a[0]+a[1]+0.625*a[3],-a[0]+a[1]));
        expected_affine_set.new_inequality_constraint((-2.0)+a[2]);
        expected_affine_set.new_equality_constraint(-a[0]+a[1]+a[2]-1.0);
        ARIADNE_TEST_PRINT(expected_affine_set);

        ARIADNE_TEST_EQUALS(affine_set,expected_affine_set);

        DRAWING_ACCURACY=1;
        figure.set_fill_colour(0.5,1.0,1.0);
        figure.draw(affine_set);
        figure.set_fill_colour(0.0,0.5,0.5);
        figure.set_fill_opacity(0.5);
        figure.draw(set);
        figure.write("test_taylor_set-affine_over_approximation");
        figure.clear();

        return;

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
        RealScalarFunction s=RealScalarFunction::coordinate(2,0);
        RealScalarFunction t=RealScalarFunction::coordinate(2,1);
        RealScalarFunction x=RealScalarFunction::coordinate(2,0);
        RealScalarFunction y=RealScalarFunction::coordinate(2,1);
        Box domain(2, -0.5,1.5, -0.5,0.5);

        TaylorConstrainedImageSet set(domain,(s,t+sqr(s)));
        set.new_negative_constraint(s+t);
        ARIADNE_TEST_PRINT(set);
        Box bx1(2, 1.0,2.0, 0.0,1.0);
        ARIADNE_TEST_PRINT(bx1);
        ARIADNE_TEST_ASSERT(set.disjoint(bx1));
    }

    void test_outer_approximation() {
        RealScalarFunction s=RealScalarFunction::coordinate(2,0);
        RealScalarFunction t=RealScalarFunction::coordinate(2,1);
        RealScalarFunction x=RealScalarFunction::coordinate(2,0);
        RealScalarFunction y=RealScalarFunction::coordinate(2,1);
        //Box domain(2, -0.5,1.5, -0.5,0.5);
        Box domain(2, -1.0,+1.0, -1.0,+1.0);
        uint accuracy = 4u;

        TaylorConstrainedImageSet set(domain,(s+0.5*t,t+sqr(s)));
        set.new_negative_constraint(s+t-0.25);
        ARIADNE_TEST_PRINT(set);

        Grid grid(2);
        int height(3);
        BinaryWord tree(
            "11010110111010101010101010101010011011010101010101010011011010101100011101010011000011101101100011100001110110001110000001110110110101010011010100111000011101011000111010011000001110110110001110000111011000111000000111111110100000000011101101101101101000000011111111101100001110100010110000111011010000111011000011101000101100001110110110100000111110110000111010001011000011101101000011101100001110100010110000111011011011010000001111111011000011101000101100001110110100001110110000111010001011000011101101101000001111101100001110100010110000111011010000111011000011101000101100010110110110110110110100000000111011011010101010101010101010010110110110110110000000001111111111011011011000111000011101100011100000111011010100100110101010011110010000001111111011001001000111110010001110000000001111011111101000111001000011100101100100000111111101100001011100100010110010110000101111101000101010001011101000100001010101101010110110011000110100010111010110100100111100000111011100001110100100001111111011011011000000111110101010001011101000101010001011111000101100001011010100001110110110111000000011111110110000111010001011000010111010100010111000101000101111110000101100001011110000101101000101101101101101100000001111101101010101100011101010100110101001110000111100000111011010101001101010011100001101011001100011110000011111110000000011111101011000111100000111101000111001000001111111001000111000001110101000101101000001111111101000111100000000001111111111010001011000110101001101000101101100001110110000111101001100010110001110100011010011000101101101100000111011011010000011111110000000000111111111111011000011101000101100001110110100001110110000111010001011000011101101101000001111101100001110100010110000111011010000111011000011101000101100001110110110110100000011111110110000111010001011000011101101000011101100001110100010110000111011011010000011111011000011101000101100001110110100001110110000111010001011000110111101010110100100101101110000111010010001111111000000001111010101001111000001110110011000100001110110101010100110101100011110100110000011101101010100110101001110000111010101001101010011100001111100000011111100000001111111100100011110000000001110110110110110100000001111111111100000000000000");

        BinaryWord leaves(
            "000000000000001000000000100000110001011100011011100110111110000001000101110001100101111000110111001101111100111111110000010000001101101001101100100011011010011011000100001101101001101100100011011010011011000010000011011010011011001000110110100110110001000011011010011011001000110110100110100000001000000000000000000001000000100000000000110111001101111000010100001010111110010101101011011111111110110101110100101111001001010110010010101011010110101101111111111101001100111110101000011000110100000001000010101101011010110101001010010100111010111100100101100100101010101010101101011010110101101100000010000011111100111101110100010000111110111010001110100100001000000011100100001100101000010100101001010100100001100100000000010110010101010110001000010010110100100100101101000010000001000010000000001101101001101100100011011010011011000100001101101001101100100011011010011011000010000011011010011011001000110110100110110001000011011010011011001000110110100110111111101011100011010010000000111101000011010010001111110111001101000011111011101000111101110100010000010000001010010000000000000100000010000000000000");

        GridTreeSet expected_very_high_accuracy_paving(grid,height,tree,leaves);

        DRAWING_ACCURACY=3;
        figure.set_bounding_box(Box(2, -2.0,+2.0, -1.5,+2.5));

        figure.set_fill_opacity(0.5);
        figure.set_fill_colour(1.0,0.5,0.5);
        figure.draw(expected_very_high_accuracy_paving);
        figure.set_fill_colour(0.0,1.0,1.0);
        figure.draw(set);
        figure.write("test_taylor_set-paving-expected");
        figure.clear();


        GridTreeSet subdivision_paving(set.dimension());
        set.subdivision_adjoin_outer_approximation_to(subdivision_paving,accuracy);
        subdivision_paving.recombine();
        ARIADNE_TEST_PRINT(subdivision_paving);

        ARIADNE_TEST_ASSERT(subset(expected_very_high_accuracy_paving,subdivision_paving));
        if(subdivision_paving.measure()/expected_very_high_accuracy_paving.measure() > 1.15) {
            ARIADNE_TEST_WARN("TaylorConstrainedImageSet::subdivision_outer_approximation(...) may yield poor approximation "<<
                              "(factor "<<subdivision_paving.measure()/expected_very_high_accuracy_paving.measure()<<")");
        }

        DRAWING_ACCURACY=2;
        figure.set_fill_colour(1.0,0.5,0.5);
        figure.draw(subdivision_paving);
        figure.set_fill_colour(0.0,1.0,1.0);
        figure.draw(set);
        figure.write("test_taylor_set-paving-subdivision");
        figure.clear();


        GridTreeSet affine_paving(set.dimension());
        set.affine_adjoin_outer_approximation_to(affine_paving,accuracy);
        affine_paving.recombine();
        ARIADNE_TEST_PRINT(affine_paving);

        ARIADNE_TEST_ASSERT(subset(expected_very_high_accuracy_paving,affine_paving));
        if(affine_paving.measure()/expected_very_high_accuracy_paving.measure() > 1.20) {
            ARIADNE_TEST_WARN("TaylorConstrainedImageSet::affine_outer_approximation(...) may yield poor approximation "<<
                              "(factor "<<affine_paving.measure()/expected_very_high_accuracy_paving.measure()<<")");
        }

        figure.set_fill_colour(1.0,0.0,0.0);
        figure.draw(difference(expected_very_high_accuracy_paving,affine_paving));
        figure.set_fill_colour(1.0,0.5,0.5);
        figure.draw(affine_paving);
        figure.set_fill_colour(0.0,1.0,1.0);
        figure.draw(set);
        figure.set_fill_colour(0.0,1.0,1.0);
        figure.write("test_taylor_set-paving-affine");
        figure.clear();


        GridTreeSet constraint_paving(set.dimension());
        set.constraint_adjoin_outer_approximation_to(constraint_paving,accuracy);
        constraint_paving.recombine();
        ARIADNE_TEST_PRINT(constraint_paving);

        ARIADNE_TEST_ASSERT(subset(expected_very_high_accuracy_paving,constraint_paving));
        if(constraint_paving.measure()/expected_very_high_accuracy_paving.measure() > 1.10) {
            ARIADNE_TEST_WARN("TaylorConstrainedImageSet::constraint_outer_approximation(...) yields poor approximation "<<
                              "(factor "<<constraint_paving.measure()/expected_very_high_accuracy_paving.measure()<<")");
        }

        figure.set_fill_colour(1.0,0.5,0.5);
        figure.draw(constraint_paving);
        figure.set_fill_colour(0.0,1.0,1.0);
        figure.draw(set);
        figure.write("test_taylor_set-paving-constraint");
        figure.clear();

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

        RealScalarFunction s=RealScalarFunction::coordinate(2,0);
        RealScalarFunction t=RealScalarFunction::coordinate(2,1);
        RealScalarFunction x0=RealScalarFunction::coordinate(3,0);
        RealScalarFunction x1=RealScalarFunction::coordinate(3,1);
        RealScalarFunction x2=RealScalarFunction::coordinate(3,2);
        uint accuracy = 3u;

        // Draw a variety of shapes
        figure.set_bounding_box(Box(2, -4.0,+4.0, -4.0,+4.0));

        ARIADNE_TEST_TRY( test_draw("box",TaylorConstrainedImageSet(Box(2,-1.01,1.01,-1.01,1.01),(s,t)),accuracy) );
        ARIADNE_TEST_TRY( test_draw("polytope",TaylorConstrainedImageSet(Box(2,-2.05,2.05,-1.05,1.05),(s,t),(s+t<=1.5)),accuracy ));
        ARIADNE_TEST_TRY( test_draw("dome",TaylorConstrainedImageSet(Box(2,-1.0,1.0,-1.0,1.0),(s,t),(s*s+t<=0.251)),accuracy) );
        ARIADNE_TEST_TRY( test_draw("disc",TaylorConstrainedImageSet(Box(2,-1.0,1.0,-1.0,1.0),(s,t),(s*s+t*t<=0.751)),accuracy) );
        ARIADNE_TEST_TRY( test_draw("parallelotope",TaylorConstrainedImageSet(Box(2,-1.0,1.0,-1.0,1.0),(2*s+t,s+t)),accuracy) );
        ARIADNE_TEST_TRY( test_draw("ellipse",TaylorConstrainedImageSet(Box(2,-1.0,1.0,-1.0,1.0),(2*s+t,s+t),(s*s+t*t<=0.75)),accuracy) );
        ARIADNE_TEST_TRY( test_draw("concave",TaylorConstrainedImageSet(Box(2,-1.01,1.01,-1.01,1.01),(s,0.25*s*s+t),(2*s+0.25*s*s+t-0.5<=0)),accuracy) );

        ARIADNE_TEST_TRY( test_draw("empty",TaylorConstrainedImageSet(Box(2,-1.01,1.01,-1.01,1.01),(s,0.25*s*s+t),(1.25+s+t*t<=0)),accuracy) );
        ARIADNE_TEST_TRY( test_draw("curve",TaylorConstrainedImageSet(Box(2,-1.01,1.01,-1.01,1.01),(s,t),(s*s-t==0)),accuracy) );

        // The following example is an enclosure from a hybrid evolution
        TaylorConstrainedImageSet enclosure(
            Box(3, -0.125,0.25, -0.125,0.125, 0.0,2.0),
            ( 8*(x2+0.0645161*x1-0.0322581*x0+0.0645162), -8+8*(0.5*x2+1.03226*x1-0.516129*x0+1.03226)),
            (0.0645161*x1-1.03226*x0+0.0645161<=0, 0.0645161*x1-1.03226*x0+0.0645161<=0, x2+0.0645161*x1-1.03226*x0+0.0645161==0));

        ARIADNE_TEST_TRY( test_draw("enclosure",enclosure,accuracy) );

    }
};



int main() {
    //TestTaylorConstrainedImageSet().test();
    TestTaylorConstrainedImageSet().test();
    std::cerr<<"INCOMPLETE "<<std::flush;
    return ARIADNE_TEST_FAILURES;
}
