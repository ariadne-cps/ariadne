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

#include "function.h"
#include "box.h"
#include "grid_set.h"
#include "affine_set.h"
#include "function_set.h"
#include "graphics.h"

#include "test.h"

using namespace Ariadne;
using namespace std;




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

        figure.set_fill_colour(1.0,0.0,0.0);
        figure.draw(set);
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

        Figure figure;
        figure.set_bounding_box(Box(2, -4.0,+4.0, -4.0,+4.0));
        figure.set_fill_colour(1.0,0.0,0.0);
        figure.draw(paving);
        figure.set_fill_colour(0.0,0.0,1.0);
        figure.draw(set);


        // The following set has difficulties
        G=Matrix<Float>::identity(2);
        h=Vector<Float>::zero(2);
        set=AffineSet(G,h);
        figure.set_fill_colour(0.75,0.75,0.75);
        figure.draw(set);

        a.resize(2);
        a[0]=2.0; a[1]=1.0; b=0.5;
        set.new_inequality_constraint(a,b);
        a[0]=-0.5; a[1]=1.0; b=0.75;
        set.new_inequality_constraint(a,b);
        a[0]=-0.5; a[1]=-1.0; b=0.875;
        set.new_inequality_constraint(a,b);

        figure.set_fill_colour(1.0,0.25,0.25);
        figure.draw(set.outer_approximation(Grid(2),3));

        figure.set_fill_colour(0.5,0.5,0.5);
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

    void test() {
        figure.set_bounding_box(Box(2, -4.0,+4.0, -4.0,+4.0));
        ARIADNE_TEST_CALL(test_empty());
        ARIADNE_TEST_CALL(test_pure_constraint());
        ARIADNE_TEST_CALL(test_constrained_image());
        figure.write("test_affine_set-draw");
        figure.clear();
        ARIADNE_TEST_CALL(test_outer_approximation());
    }

};



int main(int argc, const char* argv[])
{
    TestAffineSet().test();
    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

