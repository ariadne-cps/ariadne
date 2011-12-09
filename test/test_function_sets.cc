/***************************************************************************
 *            test_function_sets.cc
 *
 *  Copyright  2009  Pieter Collins
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

#include "config.h"
#include "function.h"
#include "box.h"
#include "grid_set.h"
#include "affine_set.h"
#include "function_set.h"
#include "graphics.h"

#include "test.h"

using namespace Ariadne;
using namespace std;


class TestConstrainedImageSet
{
  private:
    Figure figure;
  public:
    void test() {
        figure.set_bounding_box(Box(2, -4.0,+4.0, -4.0,+4.0));
        ARIADNE_TEST_CALL(test_draw());
    }

    void test_constructor() {
        List<RealScalarFunction> s=RealScalarFunction::coordinates(3);
        List<RealScalarFunction> x=RealScalarFunction::coordinates(2);

        RealBoxSet d(3,RealIntervalSet(-1,+2));
        RealConstrainedImageSet set(d,(s[0],0.25*s[0]*s[0]+s[1]+0.5*s[2]));
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);
        set.new_space_constraint(x[0]+x[1]<=2.0);
        set.apply((x[0]+x[1],x[0]-x[1]*x[1]));
    }

    void test_geometry() {
        List<RealScalarFunction> p=RealScalarFunction::coordinates(1);
        List<RealScalarFunction> s=RealScalarFunction::coordinates(3);
        List<RealScalarFunction> x=RealScalarFunction::coordinates(2);

        // Test the polytope
        RealConstrainedImageSet polytope((RealIntervalSet(-2,+2),RealIntervalSet(-2,+2)),(x[0],x[1]));
        polytope.new_parameter_constraint(x[0]+1.5*+x[1]<=1);
        ARIADNE_TEST_ASSERT(polytope.separated( (Interval(1.0,2.0),Interval(0.5,1.0)) ));
        ARIADNE_TEST_ASSERT(polytope.overlaps( (Interval(0.0,1.0),Interval(0.5,1.0)) ));

        // Test the unit disc
        RealConstrainedImageSet disc((RealIntervalSet(-2,+2),RealIntervalSet(-2,+2)),(x[0],x[1]));
        disc.new_parameter_constraint(x[0]*x[0]+x[1]*x[1]<=1);
        ARIADNE_TEST_ASSERT(disc.overlaps( (Interval(-0.5,0.5),Interval(0.25,0.5)) ));
        ARIADNE_TEST_ASSERT(disc.separated( (Interval(1,2),Interval(0.5,1)) ));
        ARIADNE_TEST_ASSERT(disc.overlaps( (Interval(0.75,2),Interval(0.5,1)) ));

        // Test a one-dimensional parabolic set
        RealConstrainedImageSet parabola(RealBoxSet(1u,RealIntervalSet(-1,+1)),(p[0],p[0]*p[0]));
        ARIADNE_TEST_PRINT(parabola);
        ARIADNE_TEST_ASSERT(parabola.separated( (Interval(0,0.5),Interval(0.5,1)) ));
        ARIADNE_TEST_ASSERT(parabola.overlaps( (Interval(0.75,2),Interval(0.5,1)) ));

        // Test whether the second iterate of the Henon map intersects a box
        RealBoxSet d(2,RealIntervalSet(-0.5,+0.5));
        RealVectorFunction h((1.5-x[0]*x[0]-0.375*x[1],x[0]));
        RealVectorFunction f=compose(h,h);
        RealConstrainedImageSet set(d,f);
        //set.new_parameter_constraint(0<=x[0]+x[1]<=1);

        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_ASSERT(set.separated(Box(2, 1.0,1.25, 1.0,1.25)));
        ARIADNE_TEST_ASSERT(set.overlaps(Box(2, -1.0,-0.875, 1.375,1.625)));
        ARIADNE_TEST_PRINT(f(Point(2,0.375,-0.375)));

        RealConstrainedImageSet idisc(RealBoxSet((RealIntervalSet(-2,+2),RealIntervalSet(-2,+2))),RealVectorFunction((x[0],x[1])));

    }

    void test_separated() {
        List<RealScalarFunction> s=RealScalarFunction::coordinates(3);
        List<RealScalarFunction> x=RealScalarFunction::coordinates(2);

        RealBoxSet d(3,RealIntervalSet(-1.1,+2.1));
        RealConstrainedImageSet set(d,(s[0],0.25*s[0]*s[0]+s[1]+0.5*s[2]));
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);

        Box b(set.bounding_box());
        List<Box> stack(1u,b);
        while(!stack.empty()) {
            Box bx=stack.back();
            stack.pop_back();
            if(bx.radius()>=0.25) {
                Pair<Box,Box> sbx=bx.split();
                stack.append(sbx.first); stack.append(sbx.second);
            } else {
                tribool overlaps = set.overlaps(bx);
                if(overlaps) {
                    figure.set_fill_colour(0.0,0.5,0.5);
                    figure.draw(bx);
                } else if(possibly(overlaps)) {
                    figure.set_fill_colour(0.0,1.0,1.0);
                    figure.draw(bx);
                }
            }
        }
        figure.write("test_constrained_image_set-separated");
        figure.clear();
    }

    void test_approximation() {
        List<RealScalarFunction> s=RealScalarFunction::coordinates(3);
        List<RealScalarFunction> x=RealScalarFunction::coordinates(2);

        RealBoxSet d(3,RealIntervalSet(-1,+2));
        RealConstrainedImageSet set(d,(s[0],0.25*s[0]*s[0]+s[1]+0.5*s[2]));
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);
        set.new_space_constraint(x[0]+x[1]<=2.0);
        ARIADNE_TEST_PRINT(set);
        GridTreeSet paving(2);
        uint depth=2;
        paving.adjoin_outer_approximation(set,depth);
        set.adjoin_outer_approximation_to(paving,depth);
        figure.draw(paving);
    }


    void test_split() {
        RealScalarFunction o=RealScalarFunction::constant(3,1.0);
        RealScalarFunction s0=RealScalarFunction::coordinate(3,0);
        RealScalarFunction s1=RealScalarFunction::coordinate(3,1);
        RealScalarFunction s2=RealScalarFunction::coordinate(3,2);
        RealScalarFunction x0=RealScalarFunction::coordinate(2,0);
        RealScalarFunction x1=RealScalarFunction::coordinate(2,1);
        RealBoxSet d(3,RealIntervalSet(-1,+1));
        RealConstrainedImageSet set(d,(s0,s1+0.5*s2*s2));
        set.new_parameter_constraint(s0+0.75*s1+s2<=0.0);

        RealConstrainedImageSet subset1,subset2;
        make_lpair(subset1,subset2)=set.split(0);
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_PRINT(subset1);
        ARIADNE_TEST_PRINT(subset2);
        ARIADNE_TEST_PRINT(set.split(1));

        RealConstrainedImageSet subset11,subset12,subset21,subset22;
        make_lpair(subset11,subset12)=subset1.split(0);
        make_lpair(subset21,subset22)=subset2.split(0);
        ARIADNE_TEST_PRINT(subset11);
        subset11.apply(RealVectorFunction((x0+2.5,x1)));
        ARIADNE_TEST_PRINT(subset11);

        set.apply(RealVectorFunction((x0-2.5,x1)));
        Figure figure;
        figure.set_bounding_box(Box(2, -4.0,+4.0, -4.0,+4.0));
        figure.set_fill_colour(1.0,1.0,1.0);
        figure.draw(set.bounding_box());
        figure.set_fill_colour(0.75,0.75,0.75);
        figure.draw(set.affine_approximation());
        figure.set_fill_colour(0.5,0.5,0.5);
        figure.draw(subset1.affine_approximation());
        figure.draw(subset2.affine_approximation());
        figure.set_fill_colour(0.25,0.25,0.25);
        figure.draw(subset11.affine_approximation());
        figure.draw(subset12.affine_approximation());
        figure.draw(subset21.affine_approximation());
        figure.draw(subset22.affine_approximation());
        figure.write("test_constrained_image_set-split");
    }

    void test_affine_approximation() {
        // Test conversionn is exact for the affine set -2<x<1; 0<y<2 3x+y<1
        List<RealScalarFunction> s=RealScalarFunction::coordinates(2);
        RealBoxSet d( (RealIntervalSet(-2.0,1.0),RealIntervalSet(0.0,2.0)) );
        RealConstrainedImageSet set(d,(s[0],s[1]));
        set.new_parameter_constraint(3*s[0]+s[1]<=1);
        IntervalAffineConstrainedImageSet affine_set=set.affine_approximation();
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_PRINT(affine_set);
        //ARIADNE_TEST_PRINT(set.affine_approximation());
    }

    void test_draw(const std::string& str, const RealConstrainedImageSet& set, uint acc) {
        figure.clear();
        figure.set_bounding_box(Box(2, -1.75,+1.75,-1.5,+2.0));
        GridTreeSet paving(set.dimension());
        set.adjoin_outer_approximation_to(paving,acc+1);
        figure.set_fill_opacity(1.0);
        figure.set_fill_colour(red);
        paving.recombine();
        figure.draw(paving);
        figure.set_fill_colour(green);
        figure.set_fill_opacity(0.5);
        typedef RealConstrainedImageSet CIS;
        CIS s1,s2,s3,s4,s5,s6,s7,s8, s9,s10,s11,s12,s13,s14,s15,s16;;
        make_lpair(s1,s2)=set.split();
        make_lpair(s3,s4)=s2.split(); make_lpair(s1,s2)=s1.split();
        make_lpair(s7,s8)=s4.split(); make_lpair(s5,s6)=s3.split(); make_lpair(s3,s4)=s2.split(); make_lpair(s1,s2)=s1.split();
        make_lpair(s15,s16)=s8.split(); make_lpair(s13,s14)=s7.split(); make_lpair(s11,s12)=s6.split(); make_lpair(s9,s10)=s5.split();
        make_lpair(s7,s8)=s4.split(); make_lpair(s5,s6)=s3.split(); make_lpair(s3,s4)=s2.split(); make_lpair(s1,s2)=s1.split();
        figure.draw(s1); figure.draw(s2); figure.draw(s3); figure.draw(s4);
        figure.draw(s5); figure.draw(s6); figure.draw(s7); figure.draw(s8);
        figure.draw(s9); figure.draw(s10); figure.draw(s11); figure.draw(s12);
        figure.draw(s13); figure.draw(s14); figure.draw(s15); figure.draw(s16);
        figure.write((std::string("test_function_set-draw-")+str).c_str());
        figure.clear();
    }

    void test_draw2(const std::string& str, const RealConstrainedImageSet& set, uint acc) {
        figure.clear();
        figure.set_bounding_box(Box(2, -1.75,+1.75,-1.5,+2.0));
        GridTreeSet paving(set.dimension());
        set.adjoin_outer_approximation_to(paving,acc+3);
        figure.set_fill_opacity(1.0);
        figure.set_fill_colour(red);
        paving.recombine();
        figure.draw(paving);
        figure.set_fill_colour(green);
        figure.set_fill_opacity(0.5);
        typedef RealConstrainedImageSet CIS;
        CIS s1,s2,s3,s4,s5,s6,s7,s8, s9,s10,s11,s12,s13,s14,s15,s16;;
        make_lpair(s1,s2)=set.split();
        make_lpair(s3,s4)=s2.split(); make_lpair(s1,s2)=s1.split();
        make_lpair(s7,s8)=s4.split(); make_lpair(s5,s6)=s3.split(); make_lpair(s3,s4)=s2.split(); make_lpair(s1,s2)=s1.split();
        make_lpair(s15,s16)=s8.split(); make_lpair(s13,s14)=s7.split(); make_lpair(s11,s12)=s6.split(); make_lpair(s9,s10)=s5.split();
        make_lpair(s7,s8)=s4.split(); make_lpair(s5,s6)=s3.split(); make_lpair(s3,s4)=s2.split(); make_lpair(s1,s2)=s1.split();
        figure.draw(s1); figure.draw(s2); figure.draw(s3); figure.draw(s4);
        figure.draw(s5); figure.draw(s6); figure.draw(s7); figure.draw(s8);
        figure.draw(s9); figure.draw(s10); figure.draw(s11); figure.draw(s12);
        figure.draw(s13); figure.draw(s14); figure.draw(s15); figure.draw(s16);
        figure.write((std::string("test_function_set-draw-")+str).c_str());
        figure.clear();
    }

    void test_draw() {
        RealScalarFunction s=RealScalarFunction::coordinate(2,0);
        RealScalarFunction t=RealScalarFunction::coordinate(2,1);
        RealScalarFunction x=RealScalarFunction::coordinate(2,0);
        RealScalarFunction y=RealScalarFunction::coordinate(2,1);
        uint acc = 2u;

        test_draw("ellipse",RealConstrainedImageSet(RealBoxSet(2,RealIntervalSet(-1.0,1.0)),(2*s+t,s+t),(s*s+t*t<=0.75)),acc+1u);
        //test_draw("concave",RealConstrainedImageSet(Box(2,-1.01,1.01,-1.01,1.01),(s,1.0*s*s+t),(2*s+0.25*s*s+t-2.0<=0)),acc);
    }
};




int main(int argc, const char* argv[])
{
    TestConstrainedImageSet().test();
    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

