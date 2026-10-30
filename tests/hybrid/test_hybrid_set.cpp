/***************************************************************************
 *            test_hybrid_set.cpp
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

#include <iostream>
#include <fstream>

#include "config.hpp"
#include "macro/testing.hpp"

#include "numeric/rational.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "function/function.hpp"

#include "set/function_set.hpp"
#include "hybrid/hybrid_set.hpp"
#include "hybrid/hybrid_paving.hpp"
#include "labelling/figure.hpp"

using namespace std;
using namespace Ariadne;


class TestHybridSet {
public:
    Void test() {
        ARIADNE_TEST_CALL(test_hybrid_list_set());
        ARIADNE_TEST_CALL(test_hybrid_paving_separation());
        ARIADNE_TEST_CALL(test_hybrid_initial_outer_approximation());
        //ARIADNE_TEST_CALL(test_compile_incorrect_hybrid_paving());
    }

    Void test_hybrid_list_set() {
        HybridListSet<ExactBoxType> hls;
        DiscreteLocation loc1(123);
        DiscreteLocation loc2(105);
        RealVariable x("x"), y("y");
        RealSpace spc1({x});
        RealSpace spc2({x, y});

        ExactBoxType bx1 = {{0, 1}};
        ExactBoxType bx2 = {{2, 3}};
        ExactBoxType bx3 = {{1, 2},
                            {2, 3}};
        ExactBoxType bx4 = {{4, 5},
                            {5, 6}};
        ExactBoxType bx5 = {{6, 7},
                            {8, 9}};
        ARIADNE_TEST_EXECUTE(hls.adjoin(loc1, spc1, bx1));
        ARIADNE_TEST_FAIL(hls.adjoin(loc1, bx3)); // Should fail due to incompatible dimensions
        ARIADNE_TEST_FAIL(hls.adjoin(loc2, bx3)); // Should fail due to unspecified space
        ARIADNE_TEST_EXECUTE(hls.adjoin(loc2, spc2, bx3));
        ARIADNE_TEST_EXECUTE(hls.adjoin(loc2, bx4));
        ARIADNE_TEST_PRINT(hls);
        ARIADNE_TEST_EXECUTE(hls.adjoin(HybridExactBox(loc1, spc1, bx2)));
        ARIADNE_TEST_EXECUTE(hls.adjoin(loc2, bx5));

        ARIADNE_TEST_PRINT(hls);

        HybridListSet<ExactBoxType>::ConstIterator iter = hls.begin();
        ARIADNE_TEST_EQUAL(*iter, HybridExactBox(loc2, spc2, bx3));
        ++iter;
        ARIADNE_TEST_EQUAL(*iter, HybridExactBox(loc2, spc2, bx4));
        ++iter;
        ARIADNE_TEST_EQUAL(*iter, HybridExactBox(loc2, spc2, bx5));
        ++iter;
        ARIADNE_TEST_EQUAL(*iter, HybridExactBox(loc1, spc1, bx1));
        ++iter;
        ARIADNE_TEST_EQUAL(*iter, HybridExactBox(loc1, spc1, bx2));
        ++iter;
        ARIADNE_TEST_ASSERT(iter == hls.end());

        ListSet<ExactBoxType> ls1;
        ls1.adjoin(bx1);
        ls1.adjoin(bx2);
        ListSet<ExactBoxType> ls2;
        ls2.adjoin(bx3);
        ls2.adjoin(bx4);
        ls2.adjoin(bx5);

        HybridListSet<ExactBoxType>::LocationsConstIterator loc_iter = hls.locations_begin();
        ARIADNE_TEST_EQUAL(loc_iter->first, loc2);
        ARIADNE_TEST_EQUAL(loc_iter->second.second, ls2);
        ++loc_iter;
        ARIADNE_TEST_EQUAL(loc_iter->first, loc1);
        ARIADNE_TEST_EQUAL(loc_iter->second.second, ls1);
        ++loc_iter;
        ARIADNE_TEST_ASSERT(loc_iter == hls.locations_end());
    }

    // ValidatedLowerKleenean(false) maps to UNLIKELY rather than FALSE:
    // such lower predicates cannot prove the negation with not possibly().
    // A location absent from the paving is disjoint from every box there.
    // In particular, checking separation must never dereference map::end().
    Void test_hybrid_paving_separation() {
        DiscreteLocation loc(123);
        RealVariable x("x");
        RealSpace state_space({x});
        MonolithicHybridSpace hybrid_space;
        hybrid_space.new_location(loc,state_space);
        HybridGridTreePaving paving{HybridGrid(hybrid_space)};
        HybridExactBox box(loc,state_space,ExactBoxType({{0,1}}));
        ARIADNE_TEST_ASSERT(definitely(paving.separated(box)));
        ARIADNE_TEST_EQUAL(paving.overlaps(box).repr(), LogicalValue::UNLIKELY);
        ARIADNE_TEST_EQUAL(paving.covers(box).repr(), LogicalValue::UNLIKELY);

        // A present but empty Euclidean paving is separated as well.
        paving[loc];
        ARIADNE_TEST_ASSERT(definitely(paving.separated(box)));
        ARIADNE_TEST_EQUAL(paving.overlaps(box).repr(), LogicalValue::UNLIKELY);
        ARIADNE_TEST_EQUAL(paving.covers(box).repr(), LogicalValue::UNLIKELY);
    }

    // Isolate the bounding-box paving from hybrid evolution.
    Void test_hybrid_initial_outer_approximation() {
        DiscreteLocation loc;
        RealVariable x("x"), v("v");
        RealSpace state_space({x,v});
        MonolithicHybridSpace hybrid_space;
        hybrid_space.new_location(loc,state_space);
        Dyadic height=2;
        Dyadic radius=1/pow(two,6);
        HybridBoundedConstraintSet initial(loc,
            {height-radius<=x<=height+radius,-radius<=v<=radius});
        HybridGridTreePaving paving{HybridGrid(hybrid_space)};
        paving.adjoin_outer_approximation(initial,4u);
        ARIADNE_TEST_ASSERT(paving.nontrivial_location(loc));
        ARIADNE_TEST_ASSERT(paving.size()>0u);
        // The domain is a RealBox; use an outward enclosure before
        // converting to an exact box for the paving containment check.
        const HybridUpperBoxes bounds=initial.bounding_box();
        ExactBoxType initial_box=cast_exact_box(bounds.euclidean_set(loc,state_space));
        ARIADNE_TEST_ASSERT(paving.euclidean_set(loc).superset(initial_box));
    }

    Void test_compile_incorrect_hybrid_paving() {
        Grid g;
        GridTreePaving gts;
        MonolithicHybridSpace hspc;
        DiscreteLocation loc;
        HybridGrid hg(hspc);
        HybridGridTreePaving hgts(hg);
        Figure fig(hgts[loc].bounding_box(), Projection2d(2, 0, 1));
        hg.has_location(loc);
        g = hg[loc];

        hgts.has_location(loc);
        GridTreePaving const &gtscr = hgts[loc];
        draw(fig, gtscr);
    }
};

Int main() {
    TestHybridSet().test();
    return ARIADNE_TEST_FAILURES;
}

