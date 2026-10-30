/***************************************************************************
 *            test_verify.cpp
 *
 *  Copyright  2006-20  Pieter Collins
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
#include <string>

#include "config.hpp"
#include "function/taylor_function.hpp"
#include "function/formula.hpp"
#include "algebra/algebra.hpp"
#include "set/function_set.hpp"
#include "paving/grid_paving.hpp"
#include "dynamics/enclosure.hpp"
#include "dynamics/orbit.hpp"
#include "dynamics/vector_field_evolver.hpp"
#include "dynamics/reachability_analyser.hpp"
#include "expression-set/expression_set.hpp"
#include "integration/integrator.hpp"
#include "labelling/figure.hpp"
#include "drawing/drawing_command_line_interface.hpp"

#include "macro/testing.hpp"

using namespace std;
using namespace Ariadne;

class TestVerifySafety
{

    typedef VectorFieldEvolver EvolverType;
    typedef EvolverType::SystemType SystemType;
    typedef SystemType::TimeType TimeType;
    typedef SystemType::StateSpaceType StateSpaceType;
    typedef EvolverType::EnclosureType EnclosureType;
    typedef RealExpressionBoundedConstraintSet SymbolicSetType;
    typedef BoundedConstraintSet SetType;
    typedef ReachabilityAnalyser<SystemType> AnalyserType;
    typedef AnalyserType::BoundingDomainType BoundingDomainType;
    typedef AnalyserType::StorageType StorageType;

    SystemType system;
    AnalyserType analyser;
    RealVariable x;
    SymbolicSetType symbolic_initial_set;
    SetType initial_set;
    SymbolicSetType symbolic_safe_set;
    SetType safe_set;
  public:

    static SystemType build_system()
    {
        std::cout<<std::setprecision(20);
        std::cerr<<std::setprecision(20);
        std::clog<<std::setprecision(20);
        RealVariable x("x");
        // Contract toward the centre of a grid cell.  Its outer
        // approximation reaches a fixed point without accumulating
        // spurious neighbouring cells as a stationary flow can.
        SystemType sys(List<DottedRealAssignment>{dot(x)=0.375_dec-x});
        cout << "Done building system\n";
        return sys;
    }

    static AnalyserType build_analyser(const SystemType& system)
    {
        GradedTaylorSeriesIntegrator integrator(StepMaximumError(1e-2_pr));

        EvolverType evolver(system,integrator);

        AnalyserType analyser(evolver);
        analyser.configuration().set_maximum_grid_fineness(2);
        cout << "Done building analyser\n";
        return analyser;
    }

    TestVerifySafety()
        : system(build_system()),
          analyser(build_analyser(system)),
          x("x"),
          symbolic_initial_set({0.3125_dec<=x<=0.34375_dec}),
          initial_set(symbolic_initial_set.euclidean_set(system.state_space())),
          symbolic_safe_set({-2_dec<=x<=2_dec}),
          safe_set(symbolic_safe_set.euclidean_set(system.state_space()))
    {
        cout << "Done creating initial and safe sets\n" << endl;

        cout << "system=" << system << endl;
        cout << "initial_set=" << initial_set << endl;
        cout << "safe_set=" << safe_set << endl;
    }

    Void test_verify_safety_no_bounding() {
        analyser.configuration().set_bounding_domain_ptr(nullptr);
        auto safety_certificate=analyser.verify_safety(initial_set,safe_set);
        ARIADNE_TEST_PRINT(safety_certificate.is_safe);
        ARIADNE_TEST_PRINT(safety_certificate.chain_reach_set.size());
        ARIADNE_TEST_PRINT(safety_certificate.safe_set.size());
        ARIADNE_TEST_ASSERT(definitely(safety_certificate.is_safe));
    }

    Void test_verify_safety_crossing_bounding() {
        // Restriction must reject an initial set entirely outside the bounding domain.
        analyser.configuration().set_bounding_domain(BoundingDomainType(analyser.system().dimension(),{1,2}));
        ARIADNE_TEST_FAIL(analyser.verify_safety(initial_set,safe_set));
    }

    Void test_verify_safety_crossing_safe() {
        analyser.configuration().set_bounding_domain_ptr(nullptr);
        // Initial cells are outside this safe set, so no proof is possible.
        SymbolicSetType smaller_symbolic_safe_set({1_dec<=x<=2_dec});
        auto smaller_safe_set = smaller_symbolic_safe_set.euclidean_set(system.state_space());

        auto safety_certificate=analyser.verify_safety(initial_set,smaller_safe_set);
        ARIADNE_TEST_ASSERT(is_indeterminate(safety_certificate.is_safe));
    }


    Void test() {
        ARIADNE_TEST_CALL(test_verify_safety_no_bounding());
        ARIADNE_TEST_CALL(test_verify_safety_crossing_bounding());
        ARIADNE_TEST_CALL(test_verify_safety_crossing_safe());
    }

};


Int main(Int argc, const char* argv[])
{
    if (not DrawingCommandLineInterface().acquire(argc,argv)) return -1;

    TestVerifySafety().test();
    return ARIADNE_TEST_FAILURES;
}

