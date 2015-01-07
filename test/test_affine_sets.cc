/***************************************************************************
 *            test_affine_sets.cc
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
#include "output/geometry2d.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

struct ExactFloatVector2d : ExactFloatVector, Vector2d {
    ExactFloatVector2d(double x, double y) : ExactFloatVector{ExactFloat(x),ExactFloat(y)}, Vector2d(x,y) { }
};

struct Polytope2d
    : public DrawableInterface
{
    List<Point2d> points;
  public:
    Polytope2d(uint n, ...) {
        assert(n>=1); va_list args; va_start(args,n);
        for(uint i=0; i!=n; ++i) {
            // NOTE: Need to store values first since order of evaluating function arguments is undefined
            double x=va_arg(args,double); double y=va_arg(args,double);
            this->points.push_back(Point2d(x,y));
        } va_end(args);
    }

    virtual Polytope2d* clone() const { return new Polytope2d(*this); }
    virtual uint dimension() const { return 2u; }

    virtual void draw(CanvasInterface& canvas, const Projection2d& p) const {
        if(points.size()==1) { canvas.dot(points[0].x,points[0].y); return; }
        canvas.move_to(points[0].x,points[0].y);
        for(uint i=1; i!=points.size(); ++i) {
            canvas.line_to(points[i].x,points[i].y);
        }
        canvas.line_to(points[0].x,points[0].y);
        canvas.fill();
    }

    Polytope2d operator+(const Vector2d& v) {
        Polytope2d r(*this); for(uint i=0; i!=r.points.size(); ++i) { r.points[i]+=v; } return r;
    }
};

static const Colour colour(0.5,1.0,1.0);
static const Colour expected_colour(1.0,0.25,0.25);

namespace Ariadne {
ExactFloat operator"" _ex (long double x) { return ExactFloat((double)x); }
ValidatedFloat operator/(int n1, ExactFloat x2) { return ExactFloat(n1)/x2; }
}

class TestAffineSet
{
  private:
    Figure figure;
    ExactBox D;
    Vector< ValidatedAffineFunction > x;
    ValidatedAffineConstrainedImageSet set;
  public:
    TestAffineSet() : set(Matrix<ExactFloat>(2,2),Vector<ExactFloat>(2)) { }

    void test_pure_constraint() {
        figure.clear();

        D = ExactBox::unit_box(2);
        x = ValidatedAffineFunction::variables(2);

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

    void test_constrained_image() {
        D = ExactBox::unit_box(3);
        x = Affine<ValidatedNumber>::variables(3);

        set=ValidatedAffineConstrainedImageSet(D, {2*x[0]+3*x[1]+x[2],x[0]+x[1]-2} );

        set.new_parameter_constraint(1.0_ex*x[0]+0.5_ex*x[1]+0.25_ex*x[2]<=0.75_ex);
        set.new_parameter_constraint(-2.0_ex*x[0]-3.0_ex*x[1]-1.0_ex*x[2]<=3.00_ex);

        ARIADNE_TEST_PRINT(set);
        Polytope2d expected_set(8, -3.0,-3.0, -3.0,-3.66667, 2.0,-2.0, 4.0,-1.0,
                                    3.0,-0.5, 1.0,-1.0, -2.0,-2.0, -3.0,-2.5);

        figure.clear(); figure.set_bounding_box(ExactBox{{-7.,+7.},{-5.,+1.}}); figure.set_fill_opacity(0.5);
        figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        figure.write("test_affine_set-constrained_image");
    }


    void test_separated() {
        D=ExactBox::unit_box(3);
        x=Affine<ValidatedNumber>::variables(3);
        ValidatedAffineConstrainedImageSet affine_set(D, {-0.9375_ex+0.0625_ex*x[0]+0.5_ex*x[1],0.87890625_ex-0.1171875_ex*x[0]+x[1]+0.00390625_ex*x[2]});
        affine_set.new_parameter_constraint(-1.1875_ex+0.0625_ex*x[0]+x[1]<=0.0_ex);
        ARIADNE_TEST_PRINT(affine_set);
        ExactBox cell1({{-1.0,-0.9375},{0.8125, 0.875}}); // subset
        ExactBox cell2({{-0.5625,-0.50},{1.4375,1.5}}); // overlaps
        ExactBox cell3({{-0.875,-0.75},{0.625,0.7578125}}); // touches at (-0.875,0.7578125);
        ExactBox cell4({{-1.1850,-1.125},{0.0625,0.125}}); // almost touches
        ExactBox cell5({{-0.9375,-0.875},{0.4375,0.5}}); // disjoint
        ExactBox cell6({{-1.5,-1.375},{0.5,0.625}}); // disjoint; regression test

        ARIADNE_TEST_ASSERT(definitely(!affine_set.separated(cell1)));
        ARIADNE_TEST_ASSERT(definitely(!affine_set.separated(cell2)));
        ARIADNE_TEST_ASSERT(!definitely(affine_set.separated(cell3)));
        ARIADNE_TEST_ASSERT(possibly(affine_set.separated(cell4)));
        ARIADNE_TEST_ASSERT(definitely(affine_set.separated(cell5)));
        ARIADNE_TEST_ASSERT(definitely(affine_set.separated(cell6)));

        figure.clear();
        figure.set_bounding_box(widen(affine_set.bounding_box(),0.125));
        figure << affine_set
               << fill_colour(0,0,1) << cell1
               << fill_colour(0,1,0) << cell2
               << fill_colour(1,1,0) << cell3
               << fill_colour(1,0,0) << cell4
               << fill_colour(1,0,0) << cell5
               << fill_colour(1,0,0) << cell6
               ;
        figure.write("test_affine_set-separated");

    }

    void test_outer_approximation() {

        {
            D=ExactBox::unit_box(2);
            x=Affine<ValidatedNumber>::variables(2);
            ValidatedAffineConstrainedImageSet set(D,{x[0],x[1]});
            set.new_parameter_constraint(x[0]+x[1]<=-0.5_ex);
            ARIADNE_TEST_PRINT(set);

            GridTreeSet paving = set.outer_approximation(Grid(2),4);
            paving.recombine();
            std::cout<<std::setprecision(17);
            ARIADNE_TEST_PRINT(paving);

            Figure figure;
            figure.set_bounding_box(widen(set.bounding_box(),+0.125));
            figure.set_fill_opacity(0.5);
            figure.set_fill_colour(1.0,0.0,0.0);
            figure.draw(paving);
            figure.set_fill_colour(0.0,0.5,1.0);
            figure.draw(set);
            figure.write("test_affine_set-outer_approximation-4");
        }

        {
            // The following set has difficulties
            ValidatedAffineConstrainedImageSet set(D,{x[0],x[1]});
            set.new_parameter_constraint(2.0_ex*x[0]+1.0_ex*x[1]<=0.5_ex);
            set.new_parameter_constraint(-0.5_ex*x[0]+1.0_ex*x[1]<=0.75_ex);
            set.new_parameter_constraint(-0.5_ex*x[0]-1.0_ex*x[1]<=0.875_ex);

            Figure figure;
            figure.set_bounding_box(widen(ExactBox{{-1.0,+1.0},{-1.0,+1.0}},+0.125));
            figure.set_fill_opacity(0.5);
            figure.set_fill_colour(1.0,0.0,0.0);
            figure.draw(set.outer_approximation(Grid(2),3));
            figure.set_fill_colour(0.0,0.5,1.0);
            figure.draw(set);

            figure.write("test_affine_set-outer_approximation-2");
            figure.clear();
        }

        {
            D=ExactBox::unit_box(3);
            x=Affine<ValidatedNumber>::variables(3);

            set=ValidatedAffineConstrainedImageSet( D, {2.0_ex*x[0]+3.0_ex*x[1]+1.0_ex*x[2]+0.05_ex, x[0]+x[1]+2.051_ex} );

            set.new_parameter_constraint(+1.01_ex*x[0]+0.51_ex*x[1]+0.251_ex*x[2]<=0.751_ex);
            set.new_parameter_constraint(-0.51_ex*x[0]+1.51_ex*x[1]+0.01_ex*x[2]<=0.51_ex);
            set.new_parameter_constraint(-2.01_ex*x[0]-3.01_ex*x[1]-1.02_ex*x[2]<=3.01_ex);

            Grid grid(2);
            GridTreeSet paving(grid);
            set.adjoin_outer_approximation_to(paving,3);
            paving.recombine();

            Figure figure;
            figure.set_bounding_box(ExactBox{{-7.0,+7.0},{-1.0,+5.0}});
            figure.set_fill_colour(1.0,0.0,0.0);
            figure.set_fill_opacity(0.5);
            figure.draw(paving);
            figure.set_fill_colour(0.0,0.5,1.0);
            figure.draw(set);
            figure.write("test_affine_set-outer_approximation-1");
            figure.clear();
        }

        {
            ValidatedAffineConstrainedImageSet set(ExactBox::unit_box(3),{-0.9375_ex+0.0625_ex*x[0]+0.5_ex*x[1],0.87890625_ex-0.1171875_ex*x[0]+x[1]+0.00390625_ex*x[2]});
            set.new_parameter_constraint(0.0625_ex*x[0]+x[1]<=1.1875_ex);
            ARIADNE_TEST_PRINT(set);

            GridTreeSet paving = set.outer_approximation(Grid(2),4);
            paving.recombine();
            std::cout<<std::setprecision(17);
            ARIADNE_TEST_PRINT(paving);

            Figure figure;
            figure.set_bounding_box(widen(set.bounding_box(),+0.125));
            figure.set_fill_opacity(0.5);
            figure.set_fill_colour(1.0,0.0,0.0);
            figure.draw(paving);
            figure.set_fill_colour(0.0,0.5,1.0);
            figure.draw(set);
            figure.write("test_affine_set-outer_approximation-3");
        }


    }


    void test_empty() {
        D = ExactBox::unit_box(2);
        x = Affine<ValidatedNumber>::variables(2);

        set=ValidatedAffineConstrainedImageSet(D, {0.125_ex*x[0]+1.0_ex, 0.125_ex*x[1]+1.0_ex} );
        set.new_parameter_constraint(-1.0_ex*x[1] <= -2.0_ex);

        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_ASSERT(set.empty());

        figure.clear();
        figure.draw(set);
        figure.write("test_affine_set-empty");
    }

    void test_uniform_error() {
        double e=0.25;
        ExactInterval I = ExactInterval(-1,+1);
        ValidatedNumber E = ValidatedNumber(-e,+e);
        D = ExactBox::unit_box(2);
        x = Affine<ValidatedNumber>::variables(2);
        ValidatedAffineConstrainedImageSet set1=ValidatedAffineConstrainedImageSet(D, {0.5_ex*x[0]+0.25_ex*x[1]+E/2, 0.25_ex*x[0]-0.25_ex*x[1]+E} );
        D = ExactBox::unit_box(4);
        x = Affine<ValidatedNumber>::variables(4);
        ValidatedAffineConstrainedImageSet set2=ValidatedAffineConstrainedImageSet(D, {0.5_ex*x[0]+0.25_ex*x[1]+e*x[2]/2, 0.25_ex*x[0]-0.25_ex*x[1]+e*x[3]} );

        figure.clear();
        figure.set_bounding_box(ExactBox{{-1.0,+1.0},{-1.0,+1.0}});
        figure.set_fill_opacity(0.5);
        figure.set_fill_colour(0,1,0);
        figure.draw(set2);
        figure.set_fill_colour(0,0,1);
        figure.draw(set1);
        figure.write("test_affine_set-uniform_error");

        Grid g(2);
        GridTreeSet paving1=outer_approximation(set1,g,4);
        GridTreeSet paving2=outer_approximation(set2,g,4);

        figure.clear();
        figure.set_fill_colour(1,0,0);
        figure.draw(paving1);
        figure.set_fill_colour(0,0,1);
        figure.draw(set1);
        figure.write("test_affine_set-uniform_error-paving");

        ARIADNE_TEST_ASSERT(paving1==paving2);

    }


    void test_draw() {
        ExactFloatVector2d offsets{0.0,0.0}; // Offsets
        ExactFloatVector& o=offsets;
        double& ox=(double&)offsets[0]; double& oy=(double&)offsets[1];
        Polytope2d expected_set(1,0.0,0.0);
        Vector< Affine<ValidatedNumber> > a;
        Vector<ExactInterval> dom;
        figure.clear();
        figure.set_bounding_box(ExactBox{{-1.0,+15.0},{-1.0,+15.0}});
        figure.set_fill_colour(0.5,1.0,1.0);
        figure.set_fill_opacity(0.5);

        {
            // Draw overlapping sets to check colours
            figure << fill_colour(expected_colour) << Polytope2d(4,0.,-0.5,2.,-0.5,2.,0.,0.,0.)
                   << fill_colour(colour) << Polytope2d(4,1.,-0.5,3.,-0.5,3.,0.,1.,0.);
        }

        {
            // Test draw of nondegenerate two-dimensional image set
            a=Affine<ValidatedNumber>::variables(2);
            dom=ExactBox::unit_box(2);
            offsets=ExactFloatVector2d{1.0,13.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.5_ex*a[0],o[1]+0.25_ex*a[1]});
            expected_set=Polytope2d(4, -0.5,-0.25, +0.5,-0.25, +0.5,+0.25,-0.5,+0.25) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of nondegenerate two-dimensional constrained image set
            a=Affine<ValidatedNumber>::variables(2);
            dom=ExactBox::unit_box(2);
            offsets=ExactFloatVector2d(4.0,13.0);
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.5_ex*a[0],o[1]+0.25_ex*a[1]},{a[0]+a[1]<=0.5_ex});
            expected_set=Polytope2d(5, -0.5,-0.25, +0.5,-0.25, +0.5,-0.125, -0.25,+0.25,-0.5,+0.25) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of two-dimensional image set with errors
            a=Affine<ValidatedNumber>::variables(2);
            dom=ExactBox::unit_box(2);
            offsets=ExactFloatVector2d{7.0,13.0};
            ValidatedNumber e=ValidatedNumber(-1.0,+1.0);

            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.25_ex*a[0]+0.25_ex*a[1]+0.5_ex*e,o[1]+0.5_ex*a[0]-0.5_ex*a[1]});
            expected_set=Polytope2d(6, -1.0,0.0, -0.5,-1.0, +0.5,-1.0, +1.0,0.0, +0.5,+1.0, -0.5,+1.0) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of nondegenerate one-dimensional image set
            a=Affine<ValidatedNumber>::variables(1);
            dom=ExactBox::unit_box(1);
            offsets=ExactFloatVector2d{1.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.5_ex*a[0],o[1]+0.25_ex*a[0]},{0*a[0]<=1.0_ex,a[0]==0.75_ex});
            expected_set=Polytope2d(1, 0.375, 0.1875) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of one-dimensional image set
            offsets=ExactFloatVector2d{4.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.5_ex*a[0],o[1]+0.25_ex*a[0]});
            expected_set=Polytope2d(2, -0.5,-0.25, 0.5,0.25) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of one-dimensional constraint set
            a=Affine<ValidatedNumber>::variables(2);
            dom=ExactBox::unit_box(2);
            offsets=ExactFloatVector2d{7.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{0*a[0]<=1.0_ex,a[0]+a[1]==0.5_ex});
            expected_set=Polytope2d(2, -0.5,+1.0, +1.0,-0.5) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of one-dimensional constraint set with one proper constraint on boundary
            offsets=ExactFloatVector2d{10.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{1.0_ex<=a[0], a[1]<=0.5_ex});
            expected_set=Polytope2d(2, +1.0,-1.0, +1.0,+0.5) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set with constraint a<=0 and a>=0
            offsets=ExactFloatVector2d{13.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+0.5_ex*a[1]<=0.75_ex,0.75_ex<=a[0]+0.5_ex*a[1]});
            expected_set=Polytope2d(2, +0.25,+1.0, +1.0,-0.5) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set with constraint 0<=a<=0
            offsets=ExactFloatVector2d{13.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{0.75_ex<=a[0]+0.5_ex*a[1]<=0.75_ex});
            expected_set=Polytope2d(2, +0.25,+1.0, +1.0,-0.5) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set with degenerate constraints at corner of box
            offsets=ExactFloatVector2d{1.0,4.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+2*a[1]<=3,1.5_ex*(a[0]+a[1])<=3,2*a[0]+a[1]<=3});
            expected_set=Polytope2d(4, -1.0,-1.0, +1.0,-1.0, +1.0,+1.0, -1.0,+1.0) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set with degenerate constraints near corner of box
            offsets=ExactFloatVector2d{4.0,4.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+2*a[1]<=2,1.5_ex*(a[0]+a[1])<=2,2*a[0]+a[1]<=2});
            expected_set=Polytope2d(6, -1.0,-1.0, +1.0,-1.0, +1.0,0.0, +0.667,+0.667, 0.0,+1.0, -1.0,+1.0) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set repeated constraints
            offsets=ExactFloatVector2d{7.0,4.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+2*a[1]<=2,a[0]*(1/3.0_ex)+a[1]*(2/3.0_ex)<=(2/3.0_ex)});
            expected_set=Polytope2d(5, -1.0,-1.0, +1.0,-1.0, +1.0,0.5, 0.0,+1.0, -1.0,+1.0) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of projection of three-dimensional box with single equality constraint
            a=Affine<ValidatedNumber>::variables(3);
            dom=ExactBox::unit_box(3);
            offsets=ExactFloatVector2d{10.0,4.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+2*a[1]+a[2]==1.5_ex});
            expected_set=Polytope2d(5, +1.0,-0.25, +1.0,+0.75, +0.5,+1.0, -1.0,+1.0, -1.0,+0.75) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }
        {
            // Test draw of two-dimensional set with nondegenerate inequality and equality constraints
            offsets=ExactFloatVector2d{1.0,7.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.3_ex*a[0]+0.20_ex*a[1]+0.05_ex*a[2],o[1]-0.10_ex*a[0]+0.1_ex*a[1]+0.05_ex*a[2]},{a[1]-a[2]<=-0.25_ex,a[0]+a[1]+a[2]==0.5_ex});
            expected_set=Polytope2d(1, +0.0,0.0) + offsets; // Unknown
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        figure.write("test_affine_set-draw");
    }

    void test() {
        figure.set_bounding_box(ExactBox{{-4.0_ex,+4.0_ex},{-4.0_ex,+4.0_ex}});
        ARIADNE_TEST_CALL(test_empty());
        ARIADNE_TEST_CALL(test_pure_constraint());
        ARIADNE_TEST_CALL(test_constrained_image());
        ARIADNE_TEST_CALL(test_separated());
        ARIADNE_TEST_CALL(test_outer_approximation());
        ARIADNE_TEST_CALL(test_uniform_error());
        ARIADNE_TEST_CALL(test_draw());
    }

};



int main(int argc, const char* argv[])
{
    TestAffineSet().test();
    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

