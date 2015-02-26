/***************************************************************************
 *            test_box.cc
 *
 *  Copyright  2005-14  Alberto Casagrande, Pieter Collins
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

#include "numeric/module.h"
#include "geometry/box.h"

#include "test.h"

using namespace Ariadne;
using std::cout; using std::endl;

template<class I> Box<I> make_box(InitializerList<Interval<Flt>> lst) {
    return Box<I>(Array<I>(lst));
}

template<class BX> class TestBox {
    typedef typename BX::IntervalType I;
    typedef typename I::LowerBoundType L;
    typedef typename I::UpperBoundType U;
    typedef typename I::MidpointType X;
    typedef I IntervalType;
    typedef Point<X> PointType;
    typedef BX BoxType;
  public:
    void test() {
        ARIADNE_TEST_CALL(test_geometric_predicates());
    }

    void test_geometric_predicates() {

    PointType s1({1.0_x,1.0_x});
    PointType s2({1.5_x,1.5_x});
    PointType s3({1.375_x,1.375_x});
    PointType s4({2.0_x,2.0_x});
    PointType s5({0.75_x,0.125_x});

    ARIADNE_TEST_CONSTRUCT( BoxType, r0, (2) );
    ARIADNE_TEST_CONSTRUCT( BoxType, r1, ({{0.0_x,1.0_x},{0.0_x,1.0_x}}) );
    ARIADNE_TEST_CONSTRUCT( BoxType, r2, ({{-0.5_x,1.5_x},{-0.375_x,0.5_x}}) );
    ARIADNE_TEST_CONSTRUCT( BoxType, r3, ({{-0.25_x,0.625_x},{0.375_x,1.5_x}}) );
    ARIADNE_TEST_CONSTRUCT( BoxType, r4, ({{0.5625_x,1.125_x},{0.4375_x,1.375_x}}) );
    ARIADNE_TEST_CONSTRUCT( BoxType, r5, ({{0.0_x,1.1875_x},{0.4375_x,1.0_x}}) );
    ARIADNE_TEST_CONSTRUCT( BoxType, r6, ({{0.0_x,1.0_x},{0.0_x,0.5_x}}) );

    ARIADNE_TEST_EQUALS(r2.measure(), 1.75_x);

    ARIADNE_TEST_ASSERT(r2.centre()==PointType({0.5_x,0.0625_x}));

    ARIADNE_TEST_ASSERT(decide(!disjoint(r1,r1)));
    ARIADNE_TEST_ASSERT(definitely(subset(r1,r1)));

    ARIADNE_TEST_ASSERT(definitely(contains(r1,s1)));
    ARIADNE_TEST_ASSERT(decide(!contains(r1,s3)));
    ARIADNE_TEST_ASSERT(decide(contains(r1,s5)));
    //ARIADNE_TEST_ASSERT(!subset_of_interior(r1,r1));
    //ARIADNE_TEST_ASSERT(subset_of_open_cover(r1,cover1));
    //ARIADNE_TEST_ASSERT(!subset_of_open_cover(r1,cover2));

    ARIADNE_TEST_EXECUTE( r4=intersection(r1,r2) );
    ARIADNE_TEST_ASSERT(r4==r6);

    ARIADNE_TEST_EXECUTE( r1[1]=IntervalType(0.0_x,0.5_x) );
    ARIADNE_TEST_ASSERT(r1==r6);


    r1=BoxType({{0,1},{-1,1},{-1,0}});

    }
};


int main() {
    TestBox<ExactFloatBox>().test();

    return ARIADNE_TEST_FAILURES;
}
