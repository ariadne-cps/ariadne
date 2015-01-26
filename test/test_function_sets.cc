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
#include "function/function.h"
#include "geometry/box.h"
#include "geometry/grid_set.h"
#include "geometry/affine_set.h"
#include "geometry/function_set.h"
#include "output/graphics.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

class TestConstrainedImageSet
{
  private:
    Figure figure;
  public:
    Void test() {
        figure.set_bounding_box(ExactBox{{-4.0,+4.0},{-4.0,+4.0}});
        ARIADNE_TEST_CALL(test_constructor());
        ARIADNE_TEST_CALL(test_geometry());
        ARIADNE_TEST_CALL(test_separated());
        ARIADNE_TEST_CALL(test_split());
        ARIADNE_TEST_CALL(test_affine_approximation());
        ARIADNE_TEST_CALL(test_draw());
    }

    Void test_constructor() {
        List<EffectiveScalarFunction> s=EffectiveScalarFunction::coordinates(3);
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);

        BoxSet d(3,IntervalSet(-1,+2));
        EffectiveConstrainedImageSet set(d,(s[0],s[0]*s[0]/4+s[1]+s[2]/2));
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);
        set.new_space_constraint(x[0]+x[1]<=2);
        set.apply((x[0]+x[1],x[0]-x[1]*x[1]));
    }

    Void test_geometry() {
        List<EffectiveScalarFunction> p=EffectiveScalarFunction::coordinates(1);
        List<EffectiveScalarFunction> s=EffectiveScalarFunction::coordinates(3);
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        ExactBox box1(2);
        ExactBox box2(2);
        ExactBox box3(2);
        Figure fig;
        Colour set_colour(0,0,1);
        Colour box_colour(1,0,1);

        // Test the polytope
        EffectiveConstrainedImageSet polytope((IntervalSet(-2,+2),IntervalSet(-2,+2)),(x[0],x[1]));
        polytope.new_parameter_constraint(x[0]+Real(1.5)*+x[1]<=1);
        box1=ExactBox( (ExactInterval(1.0,2.0),ExactInterval(0.5,1.0)) );
        box2=ExactBox( (ExactInterval(0.0,1.0),ExactInterval(0.5,1.0)) );
        ARIADNE_TEST_ASSERT(polytope.separated(box1));
        //ARIADNE_TEST_ASSERT(polytope.overlaps(box2));

        plot("test_function_sets-geometry-polytope",widen(polytope.bounding_box(),0.5),set_colour,polytope,box_colour,box1,box_colour,box2);

        // Test the unit disc
        EffectiveConstrainedImageSet disc((IntervalSet(-2,+2),IntervalSet(-2,+2)),(x[0],x[1]));
        disc.new_parameter_constraint(x[0]*x[0]+x[1]*x[1]<=1);
        box1=ExactBox( (ExactInterval(-0.5,0.5),ExactInterval(0.25,0.5)) );
        box2=ExactBox( (ExactInterval(1,2),ExactInterval(0.5,1)) );
        box3=ExactBox( (ExactInterval(0.75,2),ExactInterval(-1.0,-0.5)) );
        //ARIADNE_TEST_ASSERT(disc.overlaps(box1));
        ARIADNE_TEST_ASSERT(disc.separated(box2));
        //ARIADNE_TEST_ASSERT(disc.overlaps(box3));

        plot("test_function_sets-geometry-disc",widen(disc.bounding_box(),0.5),set_colour,disc,box_colour,box1,box_colour,box2,box_colour,box3);

        // Test a one-dimensional parabolic set
        EffectiveConstrainedImageSet parabola(BoxSet(1u,IntervalSet(-1,+1)),(p[0],p[0]*p[0]));
        box1=ExactBox( (ExactInterval(0,0.5),ExactInterval(0.5,1)) );
        box2=ExactBox( (ExactInterval(0.75,2),ExactInterval(0.5,1)) );
        ARIADNE_TEST_PRINT(parabola);
        ARIADNE_TEST_ASSERT(parabola.separated(box1));
        //ARIADNE_TEST_ASSERT(parabola.overlaps(box2));

        plot("test_function_sets-geometry-parabola",widen(parabola.bounding_box(),0.5),set_colour,parabola,box_colour,box1,box_colour,box2);

        // Test whether the second iterate of the Henon map intersects a box
        Real half(0.5);
        BoxSet d(2,IntervalSet(-half,+half));
        Real a(1.5); Real b(0.375);
        EffectiveVectorFunction h((a-x[0]*x[0]-b*x[1],x[0]));
        EffectiveVectorFunction f=compose(h,h);
        EffectiveConstrainedImageSet set(d,f);
        set.new_parameter_constraint(0<=x[0]+x[1]<=1);

        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_ASSERT(set.separated(ExactBox{{1.0,1.25},{1.0,1.25}}));
        //ARIADNE_TEST_ASSERT(set.overlaps(ExactBox(2, -1.0,-0.875, 1.375,1.625)));
        ARIADNE_TEST_PRINT(f(ExactPoint{0.375,-0.375}));


        ValidatedConstrainedImageSet idisc(ExactBox({{-2.0,+2.0},{-2.0,+2.0}}),(x[0],x[1]));
        idisc.new_parameter_constraint(x[0]*x[0]+x[1]*x[1]<=1);
        box1=ExactBox( (ExactInterval(-0.5,0.5),ExactInterval(0.25,0.75)) );
        box1=ExactBox( (ExactInterval(-0.5,0.75),ExactInterval(0.25,0.75)) );
        box2=ExactBox( (ExactInterval(1,2),ExactInterval(0.5,1)) );
        box3=ExactBox( (ExactInterval(0.75,2),ExactInterval(-1.0,-0.5)) );
        //ARIADNE_TEST_ASSERT(idisc.overlaps(box1));
        ARIADNE_TEST_ASSERT(idisc.separated(box2));
        //ARIADNE_TEST_ASSERT(idisc.overlaps(box3));
        plot("test_function_sets-geometry-idisc",widen(idisc.bounding_box(),0.5),set_colour,idisc,box_colour,box1,box_colour,box2,box_colour,box3);
    }

    Void test_separated() {
        List<EffectiveScalarFunction> s=EffectiveScalarFunction::coordinates(3);
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);

        BoxSet d(3,IntervalSet(Decimal(-1.1),Decimal(+2.1)));
        EffectiveConstrainedImageSet set(d,(s[0],s[0]*s[0]/4+s[1]+s[2]/2));
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);

        Figure figure;
        figure.set_bounding_box(widen(set.bounding_box(),0.5));
        ExactBox b(make_exact_box(set.bounding_box()));
        List<ExactBox> stack(1u,b);
        while(!stack.empty()) {
            ExactBox bx=stack.back();
            stack.pop_back();
            if(bx.radius()>=0.125) {
                Pair<ExactBox,ExactBox> sbx=bx.split();
                stack.append(sbx.first); stack.append(sbx.second);
            } else {
                Tribool overlaps = set.overlaps(bx);
                Tribool separated = set.separated(bx);
            if(definitely(overlaps)) {
                    figure.set_fill_colour(0.0,0.5,0.0);
                } else if(definitely(separated)) {
                    figure.set_fill_colour(0.75,1.0,0.75);
                } else {
                    figure.set_fill_colour(0.55,0.75,0.5);
            }
                figure.draw(bx);
            }
        }
        figure.set_fill_opacity(0.5);
        figure.set_fill_colour(0.0,0.5,0.5);
        figure.draw(set);
        figure.write("test_function_set-separated");
        figure.clear();
    }

    Void test_approximation() {
        List<EffectiveScalarFunction> s=EffectiveScalarFunction::coordinates(3);
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);

        BoxSet d(3,IntervalSet(-1,+2));
        EffectiveConstrainedImageSet set(d,(s[0],s[0]*s[0]/4+s[1]+s[2]/2));
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);
        set.new_space_constraint(x[0]+x[1]<=2.0);
        ARIADNE_TEST_PRINT(set);
        GridTreeSet paving(2);
        Nat depth=2;
        paving.adjoin_outer_approximation(set,depth);
        set.adjoin_outer_approximation_to(paving,depth);
        figure.draw(paving);
    }


    Void test_split() {
        EffectiveScalarFunction o=EffectiveScalarFunction::constant(3,1);
        EffectiveScalarFunction s0=EffectiveScalarFunction::coordinate(3,0);
        EffectiveScalarFunction s1=EffectiveScalarFunction::coordinate(3,1);
        EffectiveScalarFunction s2=EffectiveScalarFunction::coordinate(3,2);
        EffectiveScalarFunction x0=EffectiveScalarFunction::coordinate(2,0);
        EffectiveScalarFunction x1=EffectiveScalarFunction::coordinate(2,1);
        BoxSet d(3,IntervalSet(-1,+1));
        EffectiveConstrainedImageSet set(d,(s0,s1+s2*s2/2));
        set.new_parameter_constraint(s0+Real(0.75)*s1+s2<=Real(0));

        EffectiveConstrainedImageSet subset1,subset2;
        make_lpair(subset1,subset2)=set.split(0);
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_PRINT(subset1);
        ARIADNE_TEST_PRINT(subset2);
        ARIADNE_TEST_PRINT(set.split(1));

        EffectiveConstrainedImageSet subset11,subset12,subset21,subset22;
        make_lpair(subset11,subset12)=subset1.split(0);
        make_lpair(subset21,subset22)=subset2.split(0);
        ARIADNE_TEST_PRINT(subset11);
        subset11.apply(EffectiveVectorFunction((x0+Real(2.5),x1)));
        ARIADNE_TEST_PRINT(subset11);

        set.apply(EffectiveVectorFunction((x0-Real(2.5),x1)));
        Figure figure;
        figure.set_bounding_box(ExactBox({{-4.0,+4.0},{-4.0,+4.0}}));
        figure.set_fill_colour(1.0,1.0,1.0);
        figure.draw(make_exact_box(set.bounding_box()));
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
        figure.write("test_function_set-split");
    }

    Void test_affine_approximation() {
        // Test conversionn is exact for the affine set -2<x<1; 0<y<2 3x+y<1
        List<EffectiveScalarFunction> s=EffectiveScalarFunction::coordinates(2);
        BoxSet d( (IntervalSet(-2,1),IntervalSet(0,2)) );
        EffectiveConstrainedImageSet set(d,(s[0],s[1]));
        set.new_parameter_constraint(3*s[0]+s[1]<=1);
        ValidatedAffineConstrainedImageSet affine_set=set.affine_approximation();
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_PRINT(affine_set);
        //ARIADNE_TEST_PRINT(set.affine_approximation());
    }

    Void test_draw(const StringType& str, const EffectiveConstrainedImageSet& set, Nat acc) {
        figure.clear();
        figure.set_bounding_box(ExactBox({{-2.75,+2.75},{-1.5,+2.0}}));
        GridTreeSet paving(set.dimension());
        set.adjoin_outer_approximation_to(paving,acc+1);
        figure.set_fill_opacity(1.0);
        figure.set_fill_colour(red);
        paving.recombine();
        figure.draw(paving);
        figure.set_fill_colour(green);
        figure.set_fill_opacity(0.5);
        typedef EffectiveConstrainedImageSet CIS;
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
        figure.write((StringType("test_function_set-draw-")+str).c_str());
        figure.clear();
    }

    Void test_draw2(const StringType& str, const EffectiveConstrainedImageSet& set, Nat acc) {
        figure.clear();
        figure.set_bounding_box(ExactBox{{-1.75,+1.75},{-1.5,+2.0}});
        GridTreeSet paving(set.dimension());
        set.adjoin_outer_approximation_to(paving,acc+3);
        figure.set_fill_opacity(1.0);
        figure.set_fill_colour(red);
        paving.recombine();
        figure.draw(paving);
        figure.set_fill_colour(green);
        figure.set_fill_opacity(0.5);
        typedef EffectiveConstrainedImageSet CIS;
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
        figure.write((StringType("test_function_set-draw-")+str).c_str());
        figure.clear();
    }

    Void test_draw() {
        EffectiveScalarFunction s=EffectiveScalarFunction::coordinate(2,0);
        EffectiveScalarFunction t=EffectiveScalarFunction::coordinate(2,1);
        EffectiveScalarFunction x=EffectiveScalarFunction::coordinate(2,0);
        EffectiveScalarFunction y=EffectiveScalarFunction::coordinate(2,1);
        Nat acc = 2u;

        test_draw("ellipse",EffectiveConstrainedImageSet(BoxSet(2,IntervalSet(-1,1)),(2*s+t,s+t),(s*s+t*t<=0.75)),acc+1u);
        //test_draw("concave",EffectiveConstrainedImageSet(ExactBox(2,-1.01,1.01,-1.01,1.01),(s,1.0*s*s+t),(2*s+0.25*s*s+t-2.0<=0)),acc);
    }
};




Int main(Int argc, const char* argv[])
{
    TestConstrainedImageSet().test();
    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

