/***************************************************************************
 *            test_affine_sets.cpp
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
#include "geometry/box.hpp"
#include "geometry/affine_set.hpp"
#include "geometry/function_set.hpp"
#include "geometry/grid_paving.hpp"
#include "output/graphics.hpp"
#include "output/geometry2d.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

static const Colour colour(0.5,1.0,1.0);
static const Colour expected_colour(1.0,0.25,0.25);

namespace Ariadne {
inline decltype(auto) operator<=(Affine<FloatDPBounds>const& af, Dyadic w) { return af <= FloatDPValue(w,dp); }
inline FloatDPValue operator"" _ex (long double x) { return FloatDPValue((double)x); }
inline FloatDPBounds operator/(Int n1, FloatDPValue x2) { return FloatDPValue(n1)/x2; }
}

struct FloatDPValueVector2d : FloatDPValueVector, Vector2d {
    FloatDPValueVector2d(double vx, double vy) : FloatDPValueVector{FloatDPValue(vx),FloatDPValue(vy)}, Vector2d(vx,vy) { }
};

class TestAffineSet
{
  private:
    Figure figure;
    ExactBoxType D;
    Vector< Affine<FloatDPBounds> > x;
    ValidatedAffineConstrainedImageSet set;
  public:
    TestAffineSet() : set(Matrix<FloatDPValue>(2,2),Vector<FloatDPValue>(2)) { }

    Void test_pure_constraint() {
        figure.clear();

        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2);

        set=ValidatedAffineConstrainedImageSet(D, {x[0]+2,x[0]+x[1]} );

        set.new_parameter_constraint(x[0]+0.5_ex*x[1]+0.75_ex <= 0.0_ex);
        set.new_parameter_constraint(-0.5_ex*x[0]+1.5_ex*x[1]+0.5_ex <= 0.0_ex);

        ARIADNE_TEST_PRINT(set);
        figure.draw(set);

        set=ValidatedAffineConstrainedImageSet(D, {x[0]+x[1]+2,x[0]+x[1]+2});
        set.new_parameter_constraint(-0.5_ex*x[0]-1.5_ex*x[1]+0.5_ex <= 0.0_ex);

        ARIADNE_TEST_PRINT(set);
        figure.draw(set);
    }

    Void test_constrained_image() {
        D = ExactBoxType::unit_box(3);
        x = Affine<FloatDPBounds>::variables(3);

        set=ValidatedAffineConstrainedImageSet(D, {2*x[0]+3*x[1]+x[2],x[0]+x[1]-2} );

        set.new_parameter_constraint(1.0_ex*x[0]+0.5_ex*x[1]+0.25_ex*x[2]<=0.75_ex);
        set.new_parameter_constraint(-2.0_ex*x[0]-3.0_ex*x[1]-1.0_ex*x[2]<=3.00_ex);

        ARIADNE_TEST_PRINT(set);
        Polytope2d expected_set({{-3.0,-3.0}, {-3.0,-3.66667}, {2.0,-2.0}, {4.0,-1.0},
                                     {3.0,-0.5}, {1.0,-1.0}, {-2.0,-2.0}, {-3.0,-2.5}});

        figure.clear(); figure.set_bounding_box(ExactBoxType{{-7.,+7.},{-5.,+1.}}); figure.set_fill_opacity(0.5);
        figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        figure.write("test_affine_set-constrained_image");
    }


    Void test_separated() {
        D=ExactBoxType::unit_box(3);
        x=Affine<FloatDPBounds>::variables(3);
        set=ValidatedAffineConstrainedImageSet(D, {-0.9375_ex+0.0625_ex*x[0]+0.5_ex*x[1],0.87890625_ex-0.1171875_ex*x[0]+x[1]+0.00390625_ex*x[2]});
        set.new_parameter_constraint(-1.1875_ex+0.0625_ex*x[0]+x[1]<=0.0_ex);
        ARIADNE_TEST_PRINT(set);
        ExactBoxType cell1({{-1.0,-0.9375},{0.8125, 0.875}}); // subset
        ExactBoxType cell2({{-0.5625,-0.50},{1.4375,1.5}}); // overlaps
        ExactBoxType cell3({{-0.875,-0.75},{0.625,0.7578125}}); // touches at (-0.875,0.7578125);
        ExactBoxType cell4({{-1.1850,-1.125},{0.0625,0.125}}); // almost touches
        ExactBoxType cell5({{-0.9375,-0.875},{0.4375,0.5}}); // disjoint
        ExactBoxType cell6({{-1.5,-1.375},{0.5,0.625}}); // disjoint; regression test

        ARIADNE_TEST_ASSERT(!definitely(set.separated(cell1)));
        ARIADNE_TEST_ASSERT(!definitely(set.separated(cell2)));
        ARIADNE_TEST_ASSERT(!definitely(set.separated(cell3)));
        //ARIADNE_TEST_ASSERT(!definitely(affine_set.overlaps(cell3)));
        ARIADNE_TEST_ASSERT(possibly(set.separated(cell4)));
        ARIADNE_TEST_ASSERT(definitely(set.separated(cell5)));
        ARIADNE_TEST_ASSERT(definitely(set.separated(cell6)));

        figure.clear();
        figure.set_bounding_box(widen(set.bounding_box(),0.125_x));
        figure << set
               << fill_colour(0,0,1) << cell1
               << fill_colour(0,1,0) << cell2
               << fill_colour(1,1,0) << cell3
               << fill_colour(1,0,0) << cell4
               << fill_colour(1,0,0) << cell5
               << fill_colour(1,0,0) << cell6
               ;
        figure.write("test_affine_set-separated");

    }

    Void test_outer_approximation() {

        {
            D=ExactBoxType::unit_box(2);
            x=Affine<FloatDPBounds>::variables(2);
            set=ValidatedAffineConstrainedImageSet(D,{x[0],x[1]});
            set.new_parameter_constraint(x[0]+x[1]<=-0.5_ex);
            ARIADNE_TEST_PRINT(set);

            GridTreePaving paving = set.outer_approximation(Grid(2),4);
            paving.recombine();
            std::cout<<std::setprecision(17);
            ARIADNE_TEST_PRINT(paving);

            figure.clear();
            figure.set_bounding_box(widen(set.bounding_box(),0.125_x));
            figure.set_fill_opacity(0.5);
            figure.set_fill_colour(1.0,0.0,0.0);
            figure.draw(paving);
            figure.set_fill_colour(0.0,0.5,1.0);
            figure.draw(set);
            figure.write("test_affine_set-outer_approximation-4");
        }

        {
            // The following set has difficulties
            set=ValidatedAffineConstrainedImageSet(D,{x[0],x[1]});
            set.new_parameter_constraint(2.0_ex*x[0]+1.0_ex*x[1]<=0.5_ex);
            set.new_parameter_constraint(-0.5_ex*x[0]+1.0_ex*x[1]<=0.75_ex);
            set.new_parameter_constraint(-0.5_ex*x[0]-1.0_ex*x[1]<=0.875_ex);

            figure.clear();
            figure.set_bounding_box(widen(ExactBoxType{{-1.0,+1.0},{-1.0,+1.0}},+0.125_x));
            figure.set_fill_opacity(0.5);
            figure.set_fill_colour(1.0,0.0,0.0);
            figure.draw(set.outer_approximation(Grid(2),3));
            figure.set_fill_colour(0.0,0.5,1.0);
            figure.draw(set);

            figure.write("test_affine_set-outer_approximation-2");
            figure.clear();
        }

        {
            D=ExactBoxType::unit_box(3);
            x=Affine<FloatDPBounds>::variables(3);

            set=ValidatedAffineConstrainedImageSet( D, {2.0_ex*x[0]+3.0_ex*x[1]+1.0_ex*x[2]+0.05_ex, x[0]+x[1]+2.051_ex} );

            set.new_parameter_constraint(+1.01_ex*x[0]+0.51_ex*x[1]+0.251_ex*x[2]<=0.751_ex);
            set.new_parameter_constraint(-0.51_ex*x[0]+1.51_ex*x[1]+0.01_ex*x[2]<=0.51_ex);
            set.new_parameter_constraint(-2.01_ex*x[0]-3.01_ex*x[1]-1.02_ex*x[2]<=3.01_ex);

            Grid grid(2);
            GridTreePaving paving(grid);
            set.adjoin_outer_approximation_to(paving,3);
            paving.recombine();

            figure.clear();
            figure.set_bounding_box(ExactBoxType{{-7.0,+7.0},{-1.0,+5.0}});
            figure.set_fill_colour(1.0,0.0,0.0);
            figure.set_fill_opacity(0.5);
            figure.draw(paving);
            figure.set_fill_colour(0.0,0.5,1.0);
            figure.draw(set);
            figure.write("test_affine_set-outer_approximation-1");
            figure.clear();
        }

        {
            set=ValidatedAffineConstrainedImageSet(ExactBoxType::unit_box(3),{-0.9375_ex+0.0625_ex*x[0]+0.5_ex*x[1],0.87890625_ex-0.1171875_ex*x[0]+x[1]+0.00390625_ex*x[2]});
            set.new_parameter_constraint(0.0625_ex*x[0]+x[1]<=1.1875_ex);
            ARIADNE_TEST_PRINT(set);

            GridTreePaving paving = set.outer_approximation(Grid(2),4);
            paving.recombine();
            std::cout<<std::setprecision(17);
            ARIADNE_TEST_PRINT(paving);

            figure.clear();
            figure.set_bounding_box(widen(set.bounding_box(),+0.125_x));
            figure.set_fill_opacity(0.5);
            figure.set_fill_colour(1.0,0.0,0.0);
            figure.draw(paving);
            figure.set_fill_colour(0.0,0.5,1.0);
            figure.draw(set);
            figure.write("test_affine_set-outer_approximation-3");
        }


    }

    Void test_empty() {
        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2);

        set=ValidatedAffineConstrainedImageSet(D, {0.125_ex*x[0]+1.0_ex, 0.125_ex*x[1]+1.0_ex} );
        set.new_parameter_constraint(-1.0_ex*x[1] <= -2.0_ex);

        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_ASSERT(set.is_empty());

        figure.clear();
        figure.draw(set);
        figure.write("test_affine_set-empty");
    }

    Void test_invalid() {
        {
            D = ExactBoxType::unit_box(1);
            x = Affine<FloatDPBounds>::variables(2);
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D, {0.125_ex*x[0]+1.0_ex, 0.125_ex*x[1]+1.0_ex} ));
        }
        {
            D = ExactBoxType::unit_box(2);
            x = Affine<FloatDPBounds>::variables(2);
            auto a = Affine<FloatDPBounds>::variables(1);
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D, {0.125_ex*x[0]+1.0_ex, 0.125_ex*x[1]+1.0_ex}, {a[0]<=-0.25_ex} ));
        }
    }

    Void test_uniform_error() {
        FloatDPValue e=0.25_ex;
        FloatDPBounds E = FloatDPBounds(-e,+e);
        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2);
        ValidatedAffineConstrainedImageSet set1=ValidatedAffineConstrainedImageSet(D, {0.5_ex*x[0]+0.25_ex*x[1]+E/2, 0.25_ex*x[0]-0.25_ex*x[1]+E} );
        D = ExactBoxType::unit_box(4);
        x = Affine<FloatDPBounds>::variables(4);
        ValidatedAffineConstrainedImageSet set2=ValidatedAffineConstrainedImageSet(D, {0.5_ex*x[0]+0.25_ex*x[1]+e*x[2]/2, 0.25_ex*x[0]-0.25_ex*x[1]+e*x[3]} );

        figure.clear();
        figure.set_bounding_box(ExactBoxType{{-1.0,+1.0},{-1.0,+1.0}});
        figure.set_fill_opacity(0.5);
        figure.set_fill_colour(0,1,0);
        figure.draw(set2);
        figure.set_fill_colour(0,0,1);
        figure.draw(set1);
        figure.write("test_affine_set-uniform_error");

        Grid g(2);
        GridTreePaving paving1=outer_approximation(set1,g,4);
        GridTreePaving paving2=outer_approximation(set2,g,4);

        figure.clear();
        figure.set_fill_colour(1,0,0);
        figure.draw(paving1);
        figure.set_fill_colour(0,0,1);
        figure.draw(set1);
        figure.write("test_affine_set-uniform_error-paving");

        ARIADNE_TEST_ASSERT(paving1==paving2);

    }


    Void test_draw() {
        FloatDPValueVector2d offsets{0.0,0.0}; // Offsets
        FloatDPValueVector& o=offsets;
        Polytope2d expected_set({{0.0,0.0}});
        Vector< Affine<FloatDPBounds> > a;
        ExactBoxType dom;
        figure.clear();
        figure.set_bounding_box(ExactBoxType{{-1.0,+15.0},{-1.0,+15.0}});
        figure.set_fill_colour(0.5,1.0,1.0);
        figure.set_fill_opacity(0.5);

        {
            // Draw overlapping sets to check colours
            figure << fill_colour(expected_colour) << Polytope2d({{0.,-0.5},{2.,-0.5},{2.,0.},{0.,0.}})
                   << fill_colour(colour) << Polytope2d({{1.,-0.5},{3.,-0.5},{3.,0.},{1.,0.}});
        }

        {
            // Test draw of nondegenerate two-dimensional image set
            a=Affine<FloatDPBounds>::variables(2);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPValueVector2d{1.0,13.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.5_ex*a[0],o[1]+0.25_ex*a[1]});
            expected_set=Polytope2d({{-0.5,-0.25}, {+0.5,-0.25}, {+0.5,+0.25},{-0.5,+0.25}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of nondegenerate two-dimensional constrained image set
            a=Affine<FloatDPBounds>::variables(2);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPValueVector2d(4.0,13.0);
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.5_ex*a[0],o[1]+0.25_ex*a[1]},{a[0]+a[1]<=0.5_ex});
            expected_set=Polytope2d({{-0.5,-0.25}, {+0.5,-0.25}, {+0.5,-0.125}, {-0.25,+0.25},{-0.5,+0.25}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of two-dimensional image set with errors
            a=Affine<FloatDPBounds>::variables(2);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPValueVector2d{7.0,13.0};
            FloatDPBounds e=FloatDPBounds(-1.0,+1.0);

            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.25_ex*a[0]+0.25_ex*a[1]+0.5_ex*e,o[1]+0.5_ex*a[0]-0.5_ex*a[1]});
            expected_set=Polytope2d({{-1.0,0.0}, {-0.5,-1.0}, {+0.5,-1.0}, {+1.0,0.0}, {+0.5,+1.0}, {-0.5,+1.0}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of nondegenerate one-dimensional image set
            a=Affine<FloatDPBounds>::variables(1);
            dom=ExactBoxType::unit_box(1);
            offsets=FloatDPValueVector2d{1.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.5_ex*a[0],o[1]+0.25_ex*a[0]},{0.0_ex*a[0]<=1.0_ex,a[0]==0.75_ex});
            expected_set=Polytope2d({{0.375, 0.1875}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of one-dimensional image set
            a=Affine<FloatDPBounds>::variables(1);
            dom=ExactBoxType::unit_box(1);
            offsets=FloatDPValueVector2d{4.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.5_ex*a[0],o[1]+0.25_ex*a[0]});
            expected_set=Polytope2d({{-0.5,-0.25}, {0.5,0.25}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of one-dimensional constraint set
            a=Affine<FloatDPBounds>::variables(2);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPValueVector2d{7.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{0.0_ex*a[0]<=1.0_ex,a[0]+a[1]==0.5_ex});
            expected_set=Polytope2d({{-0.5,+1.0}, {+1.0,-0.5}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of one-dimensional constraint set with one proper constraint on boundary
            a=Affine<FloatDPBounds>::variables(2);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPValueVector2d{10.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{1.0_ex<=a[0], a[1]<=0.5_ex});
            expected_set=Polytope2d({{+1.0,-1.0}, {+1.0,+0.5}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set with constraint a<=0 and a>=0
            a=Affine<FloatDPBounds>::variables(2);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPValueVector2d{13.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+0.5_ex*a[1]<=0.75_ex,0.75_ex<=a[0]+0.5_ex*a[1]});
            expected_set=Polytope2d({{+0.25,+1.0}, {+1.0,-0.5}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set with constraint 0<=a<=0
            a=Affine<FloatDPBounds>::variables(2);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPValueVector2d{13.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{0.75_ex<=a[0]+0.5_ex*a[1]<=0.75_ex});
            expected_set=Polytope2d({{+0.25,+1.0}, {+1.0,-0.5}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set with degenerate constraints at corner of box
            a=Affine<FloatDPBounds>::variables(2);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPValueVector2d{1.0,4.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+2*a[1]<=3,1.5_ex*(a[0]+a[1])<=3,2*a[0]+a[1]<=3});
            expected_set=Polytope2d({{-1.0,-1.0}, {+1.0,-1.0}, {+1.0,+1.0}, {-1.0,+1.0}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set with degenerate constraints near corner of box
            a=Affine<FloatDPBounds>::variables(2);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPValueVector2d{4.0,4.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+2*a[1]<=2,1.5_ex*(a[0]+a[1])<=2,2*a[0]+a[1]<=2});
            expected_set=Polytope2d({{-1.0,-1.0}, {+1.0,-1.0}, {+1.0,0.0}, {+0.667,+0.667}, {0.0,+1.0}, {-1.0,+1.0}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set repeated constraints
            a=Affine<FloatDPBounds>::variables(2);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPValueVector2d{7.0,4.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+2.0_ex*a[1]<=2.0_ex,a[0]*(1/3.0_ex)+a[1]*(2/3.0_ex)<=(2/3.0_ex)});
            expected_set=Polytope2d({{-1.0,-1.0}, {+1.0,-1.0}, {+1.0,0.5}, {0.0,+1.0}, {-1.0,+1.0}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of projection of three-dimensional box with single equality constraint
            a=Affine<FloatDPBounds>::variables(3);
            dom=ExactBoxType::unit_box(3);
            offsets=FloatDPValueVector2d{10.0,4.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+2.0_ex*a[1]+a[2]==1.5_ex});
            expected_set=Polytope2d({{+1.0,-0.25}, {+1.0,+0.75}, {+0.5,+1.0}, {-1.0,+1.0}, {-1.0,+0.75}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of two-dimensional set with nondegenerate inequality and equality constraints
            a=Affine<FloatDPBounds>::variables(3);
            dom=ExactBoxType::unit_box(3);
            offsets=FloatDPValueVector2d{1.0,7.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.3_ex*a[0]+0.20_ex*a[1]+0.05_ex*a[2],o[1]-0.10_ex*a[0]+0.1_ex*a[1]+0.05_ex*a[2]},{a[1]-a[2]<=-0.25_ex,a[0]+a[1]+a[2]==0.5_ex});
            expected_set=Polytope2d({{+0.0,0.0}}) + offsets; // Unknown
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        figure.write("test_affine_set-draw");
    }

    Void test() {
        figure.set_bounding_box(ExactBoxType{{-4.0_ex,+4.0_ex},{-4.0_ex,+4.0_ex}});
        ARIADNE_TEST_CALL(test_empty());
        ARIADNE_TEST_CALL(test_invalid());
        ARIADNE_TEST_CALL(test_pure_constraint());
        ARIADNE_TEST_CALL(test_constrained_image());
        ARIADNE_TEST_CALL(test_separated());
        ARIADNE_TEST_CALL(test_outer_approximation());
        ARIADNE_TEST_CALL(test_uniform_error());
        ARIADNE_TEST_CALL(test_draw());
    }

};



Int main(Int argc, const char* argv[])
{
    TestAffineSet().test();
    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

