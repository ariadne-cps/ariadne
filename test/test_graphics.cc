/***************************************************************************
 *            test_graphics.cc
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
 
#include "function.h"
#include "graphics.h"
#include "point.h"
#include "box.h"
#include "zonotope.h"
#include "polytope.h"
#include "curve.h"
#include "approximate_taylor_model.h"
#include "function_set.h"
#include "grid_set.h"

using namespace Ariadne;


struct RadiusSquare : FunctionData<1,2,1> {
    template<class R, class A, class P>
    void compute(R& r, const A& x, const P& p) const {
        r[0]=sqr(x[0])+sqr(x[1])-sqr(p[0]);
    }
};
                   


int main(int argc, char **argv) 
{

    Box bx1(2); bx1[0]=Interval(-0.2,0.2); bx1[1]=Interval(-0.1,0.10);
    Box bx2(2); bx2[0]=Interval(0.1,0.3); bx2[1]=Interval(0.05,0.15);
    Box bx3(2); bx3[0]=Interval(0.2,0.4); bx3[1]=Interval(0.10,0.25);
    Box bx4(2); bx4[0]=Interval(0.25,0.5); bx4[1]=Interval(0.20,0.50);
    Box bx5(2); bx5[0]=Interval(0.4,0.8); bx5[1]=Interval(0.40,1.1);
    double z1cdata[]={0.15,0.6}; double z1gdata[]={0.05,0.0,0.05, 0.0,0.05,0.05};
    Vector<Float> z1c(2,z1cdata);
    Matrix<Float> z1g(2,3,z1gdata);
    Zonotope z1(z1c,z1g);
    Vector<Float> ts1c=z1c-Vector<Float>(2,Float(0.25));
    Matrix<Float> ts1g=z1g;
    AffineFunction afn1(ts1g,ts1c);
    ApproximateTaylorModel ts1(Box(3,Interval(-1,1)),afn1,1,0);

    Function<RadiusSquare> radius(Vector<Float>(1u,0.5));
    ConstraintSet cs1(Box(1u,Interval(-1,0)),radius);
    

    Figure g;
    g << fill_colour(0.5,1.0,1.0)
      << bx1
      << bx2
      << bx3
      << bx4
      << bx5;
    g << fill_colour(0.0,0.5,0.5)
      << polytope(z1)
      << polytope(ts1);
    g.write("test_graphics-bx1");
    
    g << fill_colour(1.0,1.0,0.5);
    //g.display();

    g.clear();

    g.set_fill_colour(1.0,0.5,1.0);
    g.draw(bx1);
    g.draw(bx2);
    g.set_fill_colour(magenta);
    g.draw(bx5);
    g.write("test_graphics-bx2");


    Box bx2d(2); bx2d[0]=Interval(0.2,0.4); bx2d[1]=Interval(0.2,0.5);
    Box bx3d(3); bx3d[0]=Interval(0.2,0.4); bx3d[1]=Interval(0.2,0.5); bx3d[2]=Interval(0.2,0.7);
    g.clear();
    g.set_projection_map(ProjectionFunction(2,3,1));
    g.draw(bx3d);
    g.write("test_graphics-bx3");
   

    InterpolatedCurve cv(Point(2,0.0));
    for(int i=1; i<=10; ++i) {
        Point pt(2); pt[0]=i/10.; pt[1]=sqr(pt[0]);
        cv.insert(i,pt);
    }

    g.clear();
    g.set_line_colour(1,0,0);
    g.draw(cv);
    g.set_line_colour(0,0,0);
    g.write("test_graphics-cv");

    GridTreeSet gts(2);
    gts.adjoin_outer_approximation(ImageSet(bx1), 6);
    gts.adjoin_outer_approximation(ImageSet(bx2), 7);
    gts.adjoin_outer_approximation(ImageSet(bx3), 8);
    gts.adjoin_outer_approximation(ImageSet(bx4), 9);
    gts.adjoin_outer_approximation(ImageSet(bx5),10);
    gts.recombine();

    g.clear();
    Box bbox(2); bbox[0]=Interval(-1,1); bbox[1]=Interval(-1,1);
    g.set_bounding_box(bbox);
    g << gts;
    g.write("test_graphics-gts1");

    g.clear();
    g.set_bounding_box(Box());
    draw(g,gts);
    g.write("test_graphics-gts2");
    

}

