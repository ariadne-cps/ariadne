/***************************************************************************
 *            test_function_sets.cpp
 *
 *  Copyright  2009-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>

#include "config.hpp"
#include "function/function.hpp"
#include "function/taylor_model.hpp"
#include "algebra/algebra.hpp"
#include "geometry/box.hpp"
#include "geometry/grid_paving.hpp"
#include "geometry/affine_set.hpp"
#include "geometry/function_set.hpp"
#include "function/formula.hpp"
#include "output/graphics.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestConstrainedImageSet
{
  private:
    Figure figure;
  public:
    Void test() {
        figure.set_bounding_box(ExactBoxType{{-4.0,+4.0},{-4.0,+4.0}});
        ARIADNE_TEST_CALL(test_separated()); return;
        ARIADNE_TEST_CALL(test_constructor());
        ARIADNE_TEST_CALL(test_domain());
        ARIADNE_TEST_CALL(test_geometry());
        ARIADNE_TEST_CALL(test_separated());
        ARIADNE_TEST_CALL(test_split());
        ARIADNE_TEST_CALL(test_affine_approximation());
        ARIADNE_TEST_CALL(test_draw());
    }

    Void test_constructor() {
        List<EffectiveScalarMultivariateFunction> s=EffectiveScalarMultivariateFunction::coordinates(3);
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);

        RealBox d(3,RealInterval(-1,+2));
        EffectiveConstrainedImageSet set(d,{s[0],s[0]*s[0]/4+s[1]+s[2]/2});
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);
        set.new_space_constraint(x[0]+x[1]<=2);
        set.apply({x[0]+x[1],x[0]-x[1]*x[1]});
    }

    Void test_domain() {
        List<EffectiveScalarMultivariateFunction> s=EffectiveScalarMultivariateFunction::coordinates(3);
        RealBox d(3,RealInterval(Decimal(-1.1),Decimal(+2.1)));
        EffectiveConstrainedImageSet set(d,{s[0],s[0]*s[0]/4+s[1]+s[2]/2});
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);

        ExactBoxType bx(cast_exact_box(set.bounding_box()));
        LowerKleenean overlaps = set.overlaps(bx);
        ValidatedLowerKleenean check_overlaps=overlaps.check(Effort::get_default());
        ARIADNE_TEST_ASSERT(definitely(not check_overlaps));
        LowerKleenean separated = set.separated(bx);
        ValidatedLowerKleenean check_separated=separated.check(Effort::get_default());
        ARIADNE_TEST_ASSERT(definitely(check_separated));
    }

    Void test_geometry() {
        List<EffectiveScalarMultivariateFunction> p=EffectiveScalarMultivariateFunction::coordinates(1);
        List<EffectiveScalarMultivariateFunction> s=EffectiveScalarMultivariateFunction::coordinates(3);
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        ExactBoxType box1(2);
        ExactBoxType box2(2);
        ExactBoxType box3(2);
        Figure fig;
        Colour set_colour(0,0,1);
        Colour box_colour(1,0,1);

        // Test the polytope
        EffectiveConstrainedImageSet polytope({{-2,+2},{-2,+2}},{x[0],x[1]});
        polytope.new_parameter_constraint(x[0]+Real(1.5)*+x[1]<=1);
        box1=ExactBoxType( {ExactIntervalType(1.0,2.0),ExactIntervalType(0.5,1.0)} );
        box2=ExactBoxType( {ExactIntervalType(0.0,1.0),ExactIntervalType(0.5,1.0)} );
        ARIADNE_TEST_ASSERT(definitely(polytope.separated(box1),Effort::get_default()));
        //ARIADNE_TEST_ASSERT(polytope.overlaps(box2));

        plot("test_function_sets-geometry-polytope",widen(polytope.bounding_box(),0.5_exact),set_colour,polytope,box_colour,box1,box_colour,box2);

        // Test the unit disc
        EffectiveConstrainedImageSet disc({RealInterval(-2,+2),RealInterval(-2,+2)},{x[0],x[1]});
        disc.new_parameter_constraint(x[0]*x[0]+x[1]*x[1]<=1);
        box1=ExactBoxType( {ExactIntervalType(-0.5,0.5),ExactIntervalType(0.25,0.5)} );
        box2=ExactBoxType( {ExactIntervalType(1,2),ExactIntervalType(0.5,1)} );
        box3=ExactBoxType( {ExactIntervalType(0.75,2),ExactIntervalType(-1.0,-0.5)} );
        //ARIADNE_TEST_ASSERT(disc.overlaps(box1));
        ARIADNE_TEST_ASSERT(definitely(disc.separated(box2),Effort::get_default()));
        //ARIADNE_TEST_ASSERT(disc.overlaps(box3));

        plot("test_function_sets-geometry-disc",widen(disc.bounding_box(),0.5_exact),set_colour,disc,box_colour,box1,box_colour,box2,box_colour,box3);

        // Test a one-dimensional parabolic set
        EffectiveConstrainedImageSet parabola({RealInterval(-1,+1)},{p[0],p[0]*p[0]});
        box1=ExactBoxType( {ExactIntervalType(0,0.5),ExactIntervalType(0.5,1)} );
        box2=ExactBoxType( {ExactIntervalType(0.75,2),ExactIntervalType(0.5,1)} );
        ARIADNE_TEST_PRINT(parabola);
        ARIADNE_TEST_ASSERT(parabola.separated(box1));
        //ARIADNE_TEST_ASSERT(parabola.overlaps(box2));

        plot("test_function_sets-geometry-parabola",widen(parabola.bounding_box(),0.5_exact),set_colour,parabola,box_colour,box1,box_colour,box2);

        // Test whether the second iterate of the Henon map intersects a box
        Real half(0.5);
        RealBox d(2,RealInterval(-half,+half));
        Real a(1.5); Real b(0.375);
        EffectiveVectorMultivariateFunction h={a-x[0]*x[0]-b*x[1],x[0]};
        EffectiveVectorMultivariateFunction f=compose(h,h);
        EffectiveConstrainedImageSet set(d,f);
        set.new_parameter_constraint(0<=x[0]+x[1]<=1);

        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_ASSERT(set.separated(ExactBoxType{{1.0,1.25},{1.0,1.25}}));
        //ARIADNE_TEST_ASSERT(set.overlaps(ExactBoxType(2, -1.0,-0.875, 1.375,1.625)));

        ARIADNE_TEST_PRINT((ExactPointType{0.375,-0.375}));
        ARIADNE_TEST_PRINT(f(ExactPointType{0.375,-0.375}));


        ValidatedConstrainedImageSet idisc(ExactBoxType({{-2.0,+2.0},{-2.0,+2.0}}),{x[0],x[1]});
        idisc.new_parameter_constraint(x[0]*x[0]+x[1]*x[1]<=1);
        box1=ExactBoxType( {ExactIntervalType(-0.5,0.5),ExactIntervalType(0.25,0.75)} );
        box1=ExactBoxType( {ExactIntervalType(-0.5,0.75),ExactIntervalType(0.25,0.75)} );
        box2=ExactBoxType( {ExactIntervalType(1,2),ExactIntervalType(0.5,1)} );
        box3=ExactBoxType( {ExactIntervalType(0.75,2),ExactIntervalType(-1.0,-0.5)} );
        //ARIADNE_TEST_ASSERT(idisc.overlaps(box1));
        ARIADNE_TEST_ASSERT(idisc.separated(box2));
        //ARIADNE_TEST_ASSERT(idisc.overlaps(box3));
        plot("test_function_sets-geometry-idisc",widen(idisc.bounding_box(),0.5_exact),set_colour,idisc,box_colour,box1,box_colour,box2,box_colour,box3);
    }

    Void test_separated() {
        List<EffectiveScalarMultivariateFunction> s=EffectiveScalarMultivariateFunction::coordinates(3);
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);

        RealBox d(3,RealInterval(Decimal(-1.1),Decimal(+2.1)));
//        RealBox d(3,RealInterval(Decimal(-0.1015625),Decimal(+2.1015625)));
        EffectiveConstrainedImageSet set(d,{s[0],s[0]*s[0]/4+s[1]+s[2]/2});
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);

        figure.clear();
        figure.set_bounding_box(widen(set.bounding_box(),0.5_exact));
        ExactBoxType b(cast_exact_box(set.bounding_box()));
        List<ExactBoxType> stack(1u,b);
        while(!stack.empty()) {
            ExactBoxType bx=stack.back();
            stack.pop_back();
            if(possibly(bx.radius()>=0.125_exact)) {
                Pair<ExactBoxType,ExactBoxType> sbx=bx.split();
                stack.append(sbx.first); stack.append(sbx.second);
            } else {
                ValidatedLowerKleenean overlaps = check(set.overlaps(bx),Effort::get_default());
                LowerKleenean separated = set.separated(bx);
            if(definitely(overlaps)) {
                    figure.set_fill_colour(0.0,0.5,0.0);
                } else if(definitely(separated,Effort::get_default())) {
                    figure.set_fill_colour(0.75,1.0,0.75);
                } else {
                    figure.set_fill_colour(0.55,0.75,0.5);
                }
                figure << bx;
            }
        }
        figure.set_fill_opacity(0.5);
        figure.set_fill_colour(0.0,0.5,0.5);
        figure.draw(set);
        figure.write("test_function_set-separated");
        figure.clear();


        // Regression test
        EffectiveConstrainedImageSet quadratic_set(RealBox({{-1.0,1.0},{-1.0,1.0}}));
        quadratic_set.apply( {2*x[0]+x[1]+x[0]*x[0]/4,x[0]+x[1]} );
        ExactBoxType box{{0.750000,1.00000},{0.00000,0.250000}};
        ARIADNE_TEST_ASSERT(not definitely(set.separated(box).check(Effort::get_default())));
    }

    Void test_approximation() {
        List<EffectiveScalarMultivariateFunction> s=EffectiveScalarMultivariateFunction::coordinates(3);
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);

        RealBox d(3,RealInterval(-1,+2));
        EffectiveConstrainedImageSet set(d,{s[0],s[0]*s[0]/4+s[1]+s[2]/2});
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);
        set.new_space_constraint(x[0]+x[1]<=2.0);
        ARIADNE_TEST_PRINT(set);
        GridTreePaving paving(2);
        Nat depth=2;
        paving.adjoin_outer_approximation(set,depth);
        set.adjoin_outer_approximation_to(paving,depth);
        figure.draw(paving);
    }


    Void test_split() {
        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(3,1);
        EffectiveScalarMultivariateFunction s0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction s1=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction s2=EffectiveScalarMultivariateFunction::coordinate(3,2);
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction x1=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveVectorMultivariateFunction translation;
        RealBox d(3,RealInterval(-1,+1));
        EffectiveConstrainedImageSet set(d,{s0,s1+s2*s2/2});
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
        translation={x0-Real(2.5),x1};
        subset11.apply(translation);
        ARIADNE_TEST_PRINT(subset11);

        translation={x0-Real(2.5),x1};
        set.apply(translation);
        figure.clear();
        figure.set_bounding_box(ExactBoxType({{-4.0,+4.0},{-4.0,+4.0}}));
        figure.set_fill_colour(1.0,1.0,1.0);
        figure.draw(cast_exact_box(set.bounding_box()));
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
        List<EffectiveScalarMultivariateFunction> s=EffectiveScalarMultivariateFunction::coordinates(2);
        RealBox d={RealInterval(-2,1),RealInterval(0,2)};
        EffectiveConstrainedImageSet set(d,{s[0],s[1]});
        set.new_parameter_constraint(3*s[0]+s[1]<=1);
        ValidatedAffineConstrainedImageSet affine_set=set.affine_approximation();
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_PRINT(affine_set);
        //ARIADNE_TEST_PRINT(set.affine_approximation());
    }

    Void test_draw(const StringType& str, const EffectiveConstrainedImageSet& set, Nat acc) {
        figure.clear();
        figure.set_bounding_box(ExactBoxType({{-2.75,+2.75},{-1.5,+2.0}}));
        GridTreePaving paving(set.dimension());
        set.adjoin_outer_approximation_to(paving,acc+1);
        figure.set_fill_opacity(1.0);
        figure.set_fill_colour(red);
        paving.recombine();
        figure.draw(paving);
        figure.set_fill_colour(green);
        figure.set_fill_opacity(0.5);
        typedef EffectiveConstrainedImageSet CIS;
        CIS s1,s2,s3,s4,s5,s6,s7,s8, s9,s10,s11,s12,s13,s14,s15,s16;
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
        figure.set_bounding_box(ExactBoxType{{-1.75,+1.75},{-1.5,+2.0}});
        GridTreePaving paving(set.dimension());
        set.adjoin_outer_approximation_to(paving,acc+3);
        figure.set_fill_opacity(1.0);
        figure.set_fill_colour(red);
        paving.recombine();
        figure.draw(paving);
        figure.set_fill_colour(green);
        figure.set_fill_opacity(0.5);
        typedef EffectiveConstrainedImageSet CIS;
        CIS s1,s2,s3,s4,s5,s6,s7,s8, s9,s10,s11,s12,s13,s14,s15,s16;
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
        EffectiveScalarMultivariateFunction s=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(2,1);
        Nat acc = 2u;

        test_draw("ellipse",EffectiveConstrainedImageSet(RealBox(2,RealInterval(-1,1)),{2*s+t,s+t},{s*s+t*t<=0.75}),acc+1u);
        //test_draw("concave",EffectiveConstrainedImageSet(ExactBoxType(2,-1.01,1.01,-1.01,1.01),{s,1.0*s*s+t),(2*s+0.25*s*s+t-2.0<=0}),acc);
    }
};




Int main(Int argc, const char* argv[])
{
    TestConstrainedImageSet().test();
    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

