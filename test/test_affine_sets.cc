/***************************************************************************
 *      test_affine_sets.cc
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

#include "function.h"
#include "box.h"
#include "grid_set.h"
#include "affine_set.h"
#include "function_set.h"
#include "graphics.h"

#include "test.h"

using namespace Ariadne;
using namespace std;


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

    virtual void draw(CanvasInterface& canvas) const {
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

class TestAffineSet
{
  private:
    Figure figure;
    Matrix<Float> G;
    Vector<Float> h;
    Vector<Float> a;
    Float b;
    AffineSet set;
  public:
    TestAffineSet() : set(Matrix<Float>(2,2),Vector<Float>(2)) { }

    void test_pure_constraint() {
     figure.clear();

     G=Matrix<Float>(2,2, 1.0,1.0, 0.0,1.0);
     h=Vector<Float>(2, 2.0,0.0);
     set=AffineSet(G,h);

     a.resize(2);
     a[0]=1.0; a[1]=0.5; b=0.75;
     set.new_inequality_constraint(a,b);
     a[0]=-0.5; a[1]=1.5; b=0.5;
     set.new_inequality_constraint(a,b);

     ARIADNE_TEST_PRINT(set);
     figure.draw(set);

     G[0][1]=1.0;
     h[0]=2.0; h[1]=2.0;
     set=AffineSet(G,h);
     //a[0]=1.0; a[1]=0.5; b=0.75;
     //set.new_inequality_constraint(a,b);
     a[0]=-0.5; a[1]=-1.5; b=0.5;
     set.new_inequality_constraint(a,b);

     ARIADNE_TEST_PRINT(set);
     figure.draw(set);
    }

    void test_constrained_image() {
     G=Matrix<Float>(2,3, 2.0,3.0,1.0, 1.0,1.0,0.0);
     //G=Matrix<Float>(2,3, 2.0,3.25,1.0, 1.0,1.0,0.0);
     h=Vector<Float>(2, 0.0,-2.0);
     set=AffineSet(G,h);

     a.resize(3);
     a[0]=1.0; a[1]=0.5; a[2]=0.25; b=0.75;
     set.new_inequality_constraint(a,b);
     a[0]=-0.5; a[1]=1.5; a[2]=0.0; b=0.5;
     //set.new_inequality_constraint(a,b);
     a[0]=-2.0; a[1]=-3.5; a[2]=-1.0; b=3.0;
     a[0]=-2.0; a[1]=-3.0; a[2]=-1.0; b=3.0;
     set.new_inequality_constraint(a,b);

     ARIADNE_TEST_PRINT(set);
     Polytope2d expected_set(8, -3.0,-3.0, -3.0,-3.66667, 2.0,-2.0, 4.0,-1.0,
      3.0,-0.5, 1.0,-1.0, -2.0,-2.0, -3.0,-2.5);

     figure.clear(); figure.set_bounding_box(Box(2,-7.,+7.,-5.,+1.)); figure.set_fill_opacity(0.5);
     figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
     figure.write("test_affine_set-constrained_image");
    }


    void test_outer_approximation() {
     G=Matrix<Float>(2,3, 2.0,3.0,1.0, 1.0,1.0,0.0);
     h=Vector<Float>(2, 0.05,2.051);
     set=AffineSet(G,h);
     a.resize(3);

     a[0]=1.01; a[1]=0.51; a[2]=0.251; b=0.751;
     set.new_inequality_constraint(a,b);
     a[0]=-0.51; a[1]=1.51; a[2]=0.01; b=0.51;
     set.new_inequality_constraint(a,b);
     a[0]=-2.01; a[1]=-3.5; a[2]=-1.02; b=3.01;
     a[0]=-2.01; a[1]=-3.01; a[2]=-1.02; b=3.01;
     set.new_inequality_constraint(a,b);

     Grid grid(2);
     GridTreeSet paving(grid);
     set.adjoin_outer_approximation_to(paving,3);
     paving.recombine();

     Figure figure;
     figure.set_bounding_box(Box(2, -4.0,+4.0, -4.0,+4.0));
     figure.set_fill_colour(1.0,0.0,0.0);
     figure.set_fill_opacity(0.5);
     figure.draw(paving);
     figure.set_fill_colour(0.0,0.5,1.0);
     figure.draw(set);


     // The following set has difficulties
     G=Matrix<Float>::identity(2);
     h=Vector<Float>::zero(2);
     set=AffineSet(G,h);

     a.resize(2);
     a[0]=2.0; a[1]=1.0; b=0.5;
     set.new_inequality_constraint(a,b);
     a[0]=-0.5; a[1]=1.0; b=0.75;
     set.new_inequality_constraint(a,b);
     a[0]=-0.5; a[1]=-1.0; b=0.875;
     set.new_inequality_constraint(a,b);

     figure.set_fill_colour(1.0,0.0,0.0);
     figure.draw(set.outer_approximation(Grid(2),3));
     figure.set_fill_colour(0.0,0.5,1.0);
     figure.draw(set);

     figure.write("test_affine_set-outer_approximation");

    }


    void test_empty() {
     G=Matrix<Float>(2,2, 0.125,0.0, 0.0,0.125);
     h=Vector<Float>(2, 1.0,1.0);
     set=AffineSet(G,h);
     a.resize(2);

     a[0]=0.0; a[1]=-1.0; b=-2.0;
     set.new_inequality_constraint(a,b);

     ARIADNE_TEST_PRINT(set);
     ARIADNE_TEST_ASSERT(set.empty());

     figure.clear();
     figure.draw(set);
     figure.write("test_affine_set-empty");
    }



    void test_draw() {
     Vector2d offsets(0.0,0.0); // Offsets
     double& ox=offsets.x; double& oy=offsets.y;
     Polytope2d expected_set(1,0.0,0.0);
     Vector< Affine<Float> > a;
     Vector<Interval> dom;
     figure.clear();
     figure.set_bounding_box(Box(2, -0.25,+14.25, -0.25,+11.25));
     figure.set_fill_colour(0.5,1.0,1.0);
     figure.set_fill_opacity(0.5);

     {
      // Draw overlapping sets to check colours
      figure << fill_colour(expected_colour) << Polytope2d(4,0.,-0.25,2.,-0.25,2.,0.,0.,0.)
    << fill_colour(colour) << Polytope2d(4,1.,-0.25,3.,-0.25,3.,0.,1.,0.);
     }

     {
      // Test draw of nondegenerate zero-dimensional image set
      a=Affine<Float>::variables(1);
      dom=Vector<Interval>::unit_box(1);
      offsets=Vector2d(1.0,1.0);
      set=AffineSet(dom, (offsets.x+0.5*a[0],offsets.y+0.25*a[0]),(0*a[0]-1),(a[0]-0.75));
      expected_set=Polytope2d(1, 0.375, 0.1875) + offsets;
      figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
     }

     {
      // Test draw of one-dimensional image set
      a=Affine<Float>::variables(1);
      dom=Vector<Interval>::unit_box(1);
      offsets=Vector2d(4.0,1.0);
      set=AffineSet(dom, (offsets.x+0.5*a[0],offsets.y+0.25*a[0]));
      expected_set=Polytope2d(2, -0.5,-0.25, 0.5,0.25) + offsets;
      figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
     }

     {
      // Test draw of one-dimensional constraint set
      a=Affine<Float>::variables(2);
      dom=Vector<Interval>::unit_box(2);
      offsets=Vector2d(7.0,1.0);
      set=AffineSet(dom, (offsets.x+a[0],offsets.y+a[1]),(0*a[0]-1.0),(a[0]+a[1]-0.5));
      expected_set=Polytope2d(2, -0.5,+1.0, +1.0,-0.5) + offsets;
      figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
     }

     {
      // Test draw of one-dimensional constraint set with one proper constraint on boundary
      a=Affine<Float>::variables(2);
      dom=Vector<Interval>::unit_box(2);
      offsets=Vector2d(10.0,1.0);
      set=AffineSet(dom, (offsets.x+a[0],offsets.y+a[1]),(-a[0]+1.0, a[1]-0.5));
      expected_set=Polytope2d(2, +1.0,-1.0, +1.0,+0.5) + offsets;
      figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
     }

     {
      // Test draw of set with constraint a<=0 and a>=0
      offsets=Vector2d(13.0,1.0);
      set=AffineSet(dom, (ox+a[0],oy+a[1]),(a[0]+0.5*a[1]-0.75,-a[0]-0.5*a[1]+0.75));
      expected_set=Polytope2d(2, +0.25,+1.0, +1.0,-0.5) + offsets;
      figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
     }

     {
      // Test draw of set with degenerate constraints at corner of box
      offsets=Vector2d(1.0,4.0);
      set=AffineSet(dom, (ox+a[0],oy+a[1]),(a[0]+2*a[1]-3,1.5*(a[0]+a[1])-3,2*a[0]+a[1]-3));
      expected_set=Polytope2d(4, -1.0,-1.0, +1.0,-1.0, +1.0,+1.0, -1.0,+1.0) + offsets;
      figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
     }

     {
      // Test draw of set with degenerate constraints near corner of box
      offsets=Vector2d(4.0,4.0);
      set=AffineSet(dom, (ox+a[0],oy+a[1]),(a[0]+2*a[1]-2,1.5*(a[0]+a[1])-2,2*a[0]+a[1]-2));
      expected_set=Polytope2d(6, -1.0,-1.0, +1.0,-1.0, +1.0,0.0, +0.667,+0.667, 0.0,+1.0, -1.0,+1.0) + offsets;
      figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
     }

     {
      // Test draw of set repeated constraints
      offsets=Vector2d(7.0,4.0);
      set=AffineSet(dom, (ox+a[0],oy+a[1]),(a[0]+2*a[1]-2,a[0]*(1/3.0)+a[1]*(2/3.0)-(2/3.0)));
      expected_set=Polytope2d(5, -1.0,-1.0, +1.0,-1.0, +1.0,0.5, 0.0,+1.0, -1.0,+1.0) + offsets;
      figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
     }

     {
      // Test draw of projection of three-dimensional box with single equality constraint
      a=Affine<Float>::variables(3);
      dom=Vector<Interval>::unit_box(3);
      offsets=Vector2d(10.0,4.0);
      set=AffineSet(dom, (ox+a[0],oy+a[1]),(0*a[0]-1),(a[0]+2*a[1]+a[2]-1.5));
      expected_set=Polytope2d(5, +1.0,-0.25, +1.0,+0.75, +0.5,+1.0, -1.0,+1.0, -1.0,+0.75) + offsets;
      figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
     }
     {
      // Test draw of two-dimensional set with nondegenerate inequality and equality constraints
      offsets=Vector2d(1.0,7.0);
      set=AffineSet(dom, (ox+0.3*a[0]+0.20*a[1]+0.05*a[2],oy-0.10*a[0]+0.1*a[1]+0.05*a[2]),(a[1]-a[2]+0.25),(a[0]+a[1]+a[2]-0.5));
      expected_set=Polytope2d(1, +0.0,0.0) + offsets; // Unknown
      figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
     }
     figure.write("test_affine_set-draw");
    }

    void test() {
     figure.set_bounding_box(Box(2, -4.0,+4.0, -4.0,+4.0));
     ARIADNE_TEST_CALL(test_empty());
     ARIADNE_TEST_CALL(test_pure_constraint());
     ARIADNE_TEST_CALL(test_constrained_image());
     ARIADNE_TEST_CALL(test_outer_approximation());
     ARIADNE_TEST_CALL(test_draw());
    }

};



int main(int argc, const char* argv[])
{
    TestAffineSet().test();
    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

