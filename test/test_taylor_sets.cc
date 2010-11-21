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
        ARIADNE_TEST_CALL(test_disjoint());
        ARIADNE_TEST_CALL(test_discretise());
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
        uint height(5);
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
        RealScalarFunction s=RealScalarFunction::coordinate(2,0);
        RealScalarFunction t=RealScalarFunction::coordinate(2,1);
        RealScalarFunction x=RealScalarFunction::coordinate(2,0);
        RealScalarFunction y=RealScalarFunction::coordinate(2,1);
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
        RealScalarFunction s=RealScalarFunction::coordinate(2,0);
        RealScalarFunction t=RealScalarFunction::coordinate(2,1);
        RealScalarFunction x=RealScalarFunction::coordinate(2,0);
        RealScalarFunction y=RealScalarFunction::coordinate(2,1);
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
    //TestTaylorConstrainedImageSet().test();
    TestTaylorConstrainedImageSet().test();
    std::cerr<<"INCOMPLETE "<<std::flush;
    return ARIADNE_TEST_FAILURES;
}
