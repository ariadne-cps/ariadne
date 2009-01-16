/***************************************************************************
 *            test_taylor_set.cc
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
 
#include <iostream>
#include "multi_index.h"
#include "taylor_variable.h"
#include "zonotope.h"
#include "taylor_set.h"
#include "grid_set.h"
#include "graphics.h"

#include "test.h"
using namespace std;
using namespace Ariadne;

class TestTaylorSet {
  public:
    void test();
  private:
    void test_outer_approximation();
};

void
TestTaylorSet::test()
{
    ARIADNE_TEST_CALL(test_outer_approximation());
}

void 
TestTaylorSet::test_outer_approximation()
{
    MultiIndex e0=MultiIndex::unit(2,0);
    MultiIndex e1=MultiIndex::unit(2,1);
    Vector<TaylorVariable> tv(2,TaylorVariable(2));
    tv[0][e0]=1.0;
    tv[0][e1]=0.25;
    tv[1][e0]=0.5;
    tv[1][e1]=1.0;
    tv[1][2u*e0]=1.0;

    TaylorSet ts=TaylorSet(tv);
    Zonotope z=zonotope(ts);

    Grid grid(2);
    uint depth(4);
    GridTreeSet gts=outer_approximation(ts,grid,depth);
    gts.recombine();

    Box bounding_box=ts.bounding_box()+Vector<Interval>(2,Interval(-1,1));
    plot("test_taylor_set-ts",bounding_box,Colour(0,1,1),ts);
    plot("test_taylor_set-z",bounding_box,Colour(0,1,0),z);
    plot("test_taylor_set-gts",bounding_box,Colour(0,1,1),gts);
    plot("test_taylor_set",PlanarProjectionMap(2,0,1),bounding_box,Colour(0,1,0),gts,Colour(1,0,1),z,Colour(0,1,1),ts);

    uint height(5);
    BooleanArray tree=make_binary_word("10101111010110110101010010110010010110110010010110100011101101010101001011011000001111111100000011011010000111101010011010001111000001111100001011000101101000111110110000111001000011111001000111100000000");
    BooleanArray word=make_binary_word("000000000100101000101001000000001000100001111111011111101100100001010110100100001001010001010010000000");
    GridTreeSet expected_outer_approximation=GridTreeSet(grid,height,tree,word);
    outer_approximation(ts,Grid(2),4);
    ARIADNE_TEST_CHECK(gts,expected_outer_approximation);
   
}


    
int main() {
    TestTaylorSet().test();
    return ARIADNE_TEST_FAILURES;
}
