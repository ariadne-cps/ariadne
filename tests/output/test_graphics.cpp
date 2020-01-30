/***************************************************************************
 *            test_graphics.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include "config.hpp"
#include "../test.hpp"

#include "function/function.hpp"
#include "output/graphics.hpp"
#include "geometry/point.hpp"
#include "geometry/box.hpp"
#include "geometry/curve.hpp"
#include "geometry/function_set.hpp"
#include "geometry/grid_paving.hpp"

using namespace Ariadne;


inline EffectiveScalarMultivariateFunction operator+(EffectiveScalarMultivariateFunction f, double c) { return f+Real(c); }
inline EffectiveScalarMultivariateFunction operator*(double c, EffectiveScalarMultivariateFunction f) { return Real(c)*f; }


Int main(Int argc, char **argv)
{

    ExactBoxType bx1(2); bx1[0]=ExactIntervalType(-0.2,0.2); bx1[1]=ExactIntervalType(-0.1,0.10);
    ExactBoxType bx2(2); bx2[0]=ExactIntervalType(0.1,0.3); bx2[1]=ExactIntervalType(0.05,0.15);
    ExactBoxType bx3(2); bx3[0]=ExactIntervalType(0.2,0.4); bx3[1]=ExactIntervalType(0.10,0.25);
    ExactBoxType bx4(2); bx4[0]=ExactIntervalType(0.25,0.5); bx4[1]=ExactIntervalType(0.20,0.50);
    ExactBoxType bx5(2); bx5[0]=ExactIntervalType(0.4,0.8); bx5[1]=ExactIntervalType(0.40,1.1);

    //Zonotope z1(z1c,z1g);
    //Polytope p1=polytope(z1);
    Real p(0.5);
    EffectiveVectorMultivariateFunction x=EffectiveVectorMultivariateFunction::identity(3);
    EffectiveVectorMultivariateFunction afn1={0.05*x[0]+0.05*x[2]+0.15,0.05*x[1]+0.05*x[2]+0.6};
    ValidatedConstrainedImageSet s1(ExactBoxType::unit_box(3),afn1);
    ApproximateBoxType bbx1=widen(s1.bounding_box(),0.25_x);

    EffectiveVectorMultivariateFunction rf(1u, sqr(x[0])+sqr(x[1])-sqr(p));
    ConstraintSet cs1(rf,RealBox(1u,RealInterval(-1,0)));

    {
        double h=10000;
        Figure g;
        g.set_bounding_box(ExactBoxType{{-1.,+1.},{-1*h,+1*h}});
        g.set_fill_colour(0.5,1.0,1.0);
        g.set_line_width(10);
        g << ExactBoxType({{-0.5,+0.0},{-0.5*h, +0.5*h}});
        g.set_line_width(1);
        g << ExactBoxType({{0.25,+0.75},{-0.5*h, +0.5*h}});
        g.write("test_graphics-canvas");
    }

    Figure g;
    g << fill_colour(0.5,1.0,1.0)
      << line_colour(0.0,0.0,0.0)
      << bx1
      << bx2
      << bx3
      << bx4
      << bx5;
    g.write("test_graphics-bx1");
    //g.display();
    g.clear();

    g.set_fill_colour(1.0,0.5,1.0);
    g.draw(bx1);
    g.draw(bx2);
    g.set_fill_colour(magenta);
    g.draw(bx5);
    g.write("test_graphics-bx2");
    g.clear();


    ExactBoxType bx2d(2); bx2d[0]=ExactIntervalType(0.2,0.4); bx2d[1]=ExactIntervalType(0.2,0.5);
    ExactBoxType bx3d(3); bx3d[0]=ExactIntervalType(0.2,0.4); bx3d[1]=ExactIntervalType(0.2,0.5); bx3d[2]=ExactIntervalType(0.2,0.7);
    g.set_projection(3,0,1);
    g.set_bounding_box(bx3d.bounding_box());
    g.draw(bx3d);
    g.write("test_graphics-bx3");
    g.clear();

    g.set_projection_map(PlanarProjectionMap(2,0,1));
    g.set_bounding_box(bbx1);

    ARIADNE_TEST_PRINT(s1);
    g << fill_colour(0.0,0.5,0.5)
      << s1;
    g.write("test_graphics-set");
    g.clear();

    InterpolatedCurve cv(0,Point<FloatDPValue>(2,FloatDPValue(0.0)));
    for(Int i=1; i<=10; ++i) {
        Point<FloatDPValue> pt(2); pt[0]=FloatDPValue(i/10.); pt[1]=FloatDPValue(i*i/100.);
        cv.insert(i,pt);
    }
    g.set_bounding_box(cv.bounding_box());
    g.set_line_colour(1,0,0);
    g.draw(cv);
    g.set_line_colour(0,0,0);
    g.write("test_graphics-curve");
    g.clear();

    GridTreePaving gts(2);
    gts.adjoin_outer_approximation(bx1, 4);
    gts.adjoin_outer_approximation(bx2, 5);
    gts.adjoin_outer_approximation(bx3, 6);
    gts.adjoin_outer_approximation(bx4, 7);
    gts.adjoin_outer_approximation(bx5, 8);
    gts.recombine();

    ExactBoxType bbox(2); bbox[0]=ExactIntervalType(-2,2); bbox[1]=ExactIntervalType(-2,2);
    g.set_bounding_box(bbox);
    g << gts;
    g.write("test_graphics-paving");
    g.clear();

}

