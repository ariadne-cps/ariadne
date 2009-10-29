/***************************************************************************
 *            test_textplot.cc
 *
 *  Copyright 2009  Davide Bresolin
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
#include "textplot.h"
#include "point.h"
#include "box.h"
#include "zonotope.h"
#include "polytope.h"
#include "curve.h"
#include "taylor_set.h"
#include "function_set.h"
#include "grid_set.h"
#include "hybrid_set.h"

#include "user_function.h"

using namespace Ariadne;


struct RadiusSquare : VectorFunctionData<1,2,1> {
    template<class R, class A, class P>
    static void compute(R& r, const A& x, const P& p) {
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
    // 3d Boxes
    Box bx6(3, 0.2,0.4, 0.10,0.25, 0.3, 0.5);
    Box bx7(3, 0.25,0.5, 0.20,0.50, 0.1, 0.2);
    Box bx8(3, 0.4,0.8, 0.40,1.1, 0.4,1.25);
    double z1cdata[]={0.15,0.6}; double z1gdata[]={0.05,0.0,0.05, 0.0,0.05,0.05};
    Vector<Float> z1c(2,z1cdata);
    Matrix<Float> z1g(2,3,z1gdata);
    Zonotope z1(z1c,z1g);
    Vector<Float> ts1c=z1c-Vector<Float>(2,Float(0.25));
    Matrix<Float> ts1g=z1g;
    VectorAffineFunction afn1(ts1g,ts1c);
    TaylorSet ts1(afn1,Box::unit_box(3));

    VectorUserFunction<RadiusSquare> radius(Vector<Float>(1u,0.5));
    ConstraintSet cs1(Box(1u,Interval(-1,0)),radius);
    
    std::cout << "Testing boxes.." << std::endl;
    TextPlot g("test_textplot-bx1.txt");
    g << bx1
      << bx2
      << bx3
      << bx4
      << bx5;
    g.close();

    g.open("test_textplot-bx2.txt");
    g.draw(bx6);
    g.draw(bx7);
    g.draw(bx8);
    g.close();

/*
 *  Zonotopes and TaylorSets are not supported by textplot, skipping test
 *   
    std::cout << "Testing zonotopes and TaylorSets.." << std::endl;    
    g.open("test_textplot-zts.txt");
    g << z1
      << ts1;
    g.close();
 *
 */
    
    std::cout << "Testing interpolated curves.." << std::endl;    
    InterpolatedCurve cv(Point(2,0.0));
    for(int i=1; i<=10; ++i) {
        Point pt(2); pt[0]=i/10.; pt[1]=sqr(pt[0]);
        cv.insert(i,pt);
    }

    g.open("test_textplot-cv.txt");
    g << cv;
    g.close();
    
    std::cout << "Testing list sets.." << std::endl;
    ListSet<Box> ls;
    ls.adjoin(bx1);
    ls.adjoin(bx2);
    ls.adjoin(bx3);
    
    g.open("test_textplot-ls.txt");
    g << ls;
    g.close();

    std::cout << "Testing grid sets.." << std::endl;    
    GridTreeSet gts(2);
    gts.adjoin_outer_approximation(ImageSet(bx1), 2);
    gts.adjoin_outer_approximation(ImageSet(bx2), 3);
    gts.adjoin_outer_approximation(ImageSet(bx3), 4);
    gts.adjoin_outer_approximation(ImageSet(bx4), 5);
    gts.adjoin_outer_approximation(ImageSet(bx5), 6);
    gts.recombine();

    g.open("test_textplot-gts.txt");
    g << gts;
    g.close();
    
    std::cout << "Testing hybrid basic sets.." << std::endl;
    DiscreteState q1(1);
    DiscreteState q2(2);
    DiscreteState q3(3);
    HybridBasicSet<Box> hbx1(q1,bx1);
    HybridBasicSet<Box> hbx2(q2,bx2);
    HybridBasicSet<Box> hbx3(q3,bx3);
    g.open("test_textplot-hbx.txt");
    g << hbx1 << hbx2 << hbx3;
    g.close();

    std::cout << "Testing hybrid list sets.." << std::endl;
    HybridListSet<Box> hls;
    hls.adjoin(q1,bx1);
    hls.adjoin(q2,bx2);
    hls.adjoin(q3,bx3);
    g.open("test_textplot-hls.txt");
    g << hls;
    g.close();
    
     std::cout << "Testing hybrid grid sets.." << std::endl;
    HybridGridTreeSet hgts;
    hgts.insert(make_pair(q1,GridTreeSet(2)));
    hgts[q1].adjoin_outer_approximation(ImageSet(bx1), 2);
    hgts.insert(make_pair(q2,GridTreeSet(2)));
    hgts[q2].adjoin_outer_approximation(ImageSet(bx2), 3);
    hgts.insert(make_pair(q3,GridTreeSet(2)));
    hgts[q3].adjoin_outer_approximation(ImageSet(bx3), 4);
    hgts.recombine();
    g.open("test_textplot-hgts.txt");
    g << hgts;
    g.close();
    
}

