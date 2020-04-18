/***************************************************************************
 *            test_hybrid_automaton.cpp
 *
 *  Copyright  2009-20  Pieter Collins
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
#include "../test.hpp"

#include "symbolic/expression.hpp"
#include "symbolic/valuation.hpp"
#include "hybrid/hybrid_automata.hpp"

using namespace Ariadne;

class TestHybridSystem {
  public:
    Void test();
  private:
    Void test_build_hybrid_system();
    Void test_static_analysis();
  private:
    RestrictiveHybridAutomaton _system;
};

Void
TestHybridSystem::test()
{
    std::clog<<std::boolalpha;
    std::cerr<<std::boolalpha;
    ARIADNE_TEST_CALL(test_build_hybrid_system());
    ARIADNE_TEST_CALL(test_static_analysis());
}


Void
TestHybridSystem::test_build_hybrid_system()
{
    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0_dec);
    RealConstant hmin("hmin",5.5_dec);
    RealConstant hmax("hmax",8.0_dec);
    RealConstant delta("delta",0.05_dec);
    RealConstant lambda("lambda",0.02_dec);
    RealConstant rate("rate",0.3_dec);

    // Declare the system variables
    RealVariable height("height");
    RealVariable alpha("alpha");

    // Create the tank object
    RestrictiveHybridAutomaton tank_component;
    // The water level is always given by the same dynamic
    // The inflow is controlled by the valve alpha, the outflow depends on the
    // pressure, which is proportional to the water height.
    tank_component.new_dynamic(DiscretePredicate(true),{dot(height)=-lambda*height+rate*alpha});
    tank_component.disable_events(DiscretePredicate(true),complement(EventSet()));
    ARIADNE_TEST_PRINT(tank_component.mode(DiscreteLocation()));
    ARIADNE_TEST_PRINT(tank_component);

    // Describe the valve model

    // Declare the events we use in the valve automaton
    DiscreteEvent start_opening("start_opening");
    DiscreteEvent start_closing("start_closing");
    DiscreteEvent finished_opening("finished_opening");
    DiscreteEvent finished_closing("finished_closing");

    // Declare the locations we use
    StringConstant open("open");
    StringConstant opening("opening");
    StringConstant closed("closed");
    StringConstant closing("closing");

    RealVariable beta("beta");
    StringVariable valve("valve");

    RestrictiveHybridAutomaton valve_component;

    valve_component.disable_events(valve==open || valve==opening,complement(EventSet{start_closing,finished_opening}));
    valve_component.disable_events((valve==closed || valve==closing),complement(EventSet{start_opening,finished_closing}));
    valve_component.disable_events(valve==open,finished_opening);
    valve_component.disable_events(valve==closed,finished_closing);

    // Since alpha is a known constant when the valve is open or closed,
    // specify alpha by an algebraic equation.
    valve_component.new_auxiliary(valve==open,{let(alpha)=1});
    valve_component.new_auxiliary(valve==closed,{let(alpha)=-1});
    // Specify the differential equation for how the valve opens/closes.
    valve_component.new_dynamic(valve==opening,{dot(alpha)=1/T});
    valve_component.new_dynamic(valve==closing,{dot(alpha)=-1/T});

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    //valve_automaton.new_invariant(open,start_closing,height<=hmax || (height>=hmin && !(height<=hmin+1)));
    valve_component.new_invariant(valve==open,start_closing,height<=hmax);
    valve_component.new_invariant(valve==opening,start_closing,height<=hmax);
    valve_component.new_invariant(valve==opening,finished_opening,alpha<=1);
    valve_component.new_invariant(valve==closed,start_opening,height>=hmin);
    valve_component.new_invariant(valve==closing,start_opening,height>=hmin);
    valve_component.new_invariant(valve==closing,finished_closing,alpha>=0);

    valve_component.new_transition(valve==closed,start_opening,(next(valve)=opening),{next(alpha)=alpha},height<=hmin);
    valve_component.new_transition(valve==closing,start_opening,next(valve)=opening,{next(alpha)=alpha},height<=hmin);
    valve_component.new_transition(valve==open,start_closing,(next(valve)=closing),{next(alpha)=alpha},height>=hmax);
    valve_component.new_transition(valve==opening,start_closing,(next(valve)=closing),{next(alpha)=alpha},height>=hmax);

    // Set the transitions for when the valve finished opening.
    // Since alpha is defined by an algebraic equation in the new mode,
    // it may not be specified in the reset.
    valve_component.new_transition(valve==opening,finished_opening,next(valve)=open,alpha>=1);
    valve_component.new_transition(valve==closing,finished_closing,next(valve)=closed,alpha<=0);

    ARIADNE_TEST_PRINT(valve_component.mode({valve|opening}));
    ARIADNE_TEST_PRINT(valve_component.mode({valve|closing}));
    ARIADNE_TEST_PRINT(valve_component.mode({valve|open}));
    ARIADNE_TEST_PRINT(valve_component.mode({valve|closed}));
    ARIADNE_TEST_PRINT(valve_component);

    RestrictiveHybridAutomaton watertank_system=compose({tank_component,valve_component});
    std::cout << "watertank_system:\n" << watertank_system << "\n";

    _system=watertank_system;
    _system.verbosity=0;

    Set<DiscreteLocation> reachable = valve_component.reachable_locations(valve|open);
    std::cerr << reachable;
    for(Set<DiscreteLocation>::ConstIterator iter=reachable.begin(); iter!=reachable.end(); ++iter) {
        std::cerr<<valve_component.mode(*iter)<<"\n";
        valve_component.check_mode(*iter);
    }

}





Void TestHybridSystem::test_static_analysis()
{
    DiscreteLocation valve_opening(StringVariable("valve")|"opening");
    ARIADNE_TEST_EXECUTE(this->_system.check_mode(valve_opening));
    //ARIADNE_TEST_EXECUTE(_system.check_mode((tank|draining,valve|open)));
    //ARIADNE_TEST_EXECUTE(_system.check_reachable_modes((tank|draining,valve|opening)));
    //ARIADNE_TEST_EXECUTE(_system.discrete_reachability((tank|draining,valve|opening)));

    //DiscreteLocation invalid_location=(tank|"nonexistent",valve|"opening");
    //ARIADNE_TEST_FAIL(_system.discrete_reachability(invalid_location));
}



class TestHybridAutomaton {
  public:
    Void test();

    Void test_distinct_modes();
    Void test_overspecified_location();
    Void test_overspecified_dynamic();
    Void test_algebraic_loops();
    Void test_nonexistent_mode();
    Void test_multiple_guard();
    Void test_multiple_transition();
    Void test_overspecified_reset();

    Void test_underspecified_mode();
    Void test_underspecified_reset();


    Void test_build_hybrid_system();
    Void test_build_intensional_hybrid_automaton();
    Void test_static_analysis();
  private:
    CompositeHybridAutomaton _system;
};

Void
TestHybridAutomaton::test()
{
    std::clog<<std::boolalpha;
    std::cerr<<std::boolalpha;
    ARIADNE_TEST_CALL(test_distinct_modes());
    ARIADNE_TEST_CALL(test_overspecified_location());
    ARIADNE_TEST_CALL(test_overspecified_dynamic());
    ARIADNE_TEST_CALL(test_algebraic_loops());
    ARIADNE_TEST_CALL(test_nonexistent_mode());
    ARIADNE_TEST_CALL(test_multiple_guard());
    ARIADNE_TEST_CALL(test_multiple_transition());
    ARIADNE_TEST_CALL(test_overspecified_reset());
    ARIADNE_TEST_CALL(test_underspecified_mode());
    ARIADNE_TEST_CALL(test_underspecified_reset());
    return;
    ARIADNE_TEST_CALL(test_build_hybrid_system());
    ARIADNE_TEST_CALL(test_static_analysis());
}


Void
TestHybridAutomaton::test_distinct_modes()
{
    HybridAutomaton system;
    StringVariable q1("q1");
    StringVariable q2("q2");
    ARIADNE_TEST_EXECUTE(system.new_mode(q1|"s1"));
    ARIADNE_TEST_EXECUTE(system.new_mode(q1|"s2"));
    ARIADNE_TEST_THROWS(system.new_mode(q2|"s3"),IndistinguishableModeError);
}

Void
TestHybridAutomaton::test_overspecified_location()
{
    StringVariable q1("q1");
    StringVariable q2("q2");
    StringVariable q3("q3");

    HybridAutomaton system1;
    ARIADNE_TEST_EXECUTE(system1.new_mode(q1|"s1"));
    DiscreteLocation loc1({q1|"s1"});
    ARIADNE_TEST_ASSERT(system1.has_mode(loc1));
    ARIADNE_TEST_ASSERT(system1.has_partial_mode(loc1));
    ARIADNE_TEST_PRINT(system1.mode(loc1));
    DiscreteLocation loc1e({q1|"s1",q3|"s3"});
    ARIADNE_TEST_ASSERT(not system1.has_mode(loc1e));
    ARIADNE_TEST_ASSERT(system1.has_partial_mode(loc1e));
    ARIADNE_TEST_PRINT(system1.mode(loc1e));

    HybridAutomaton system2;
    ARIADNE_TEST_EXECUTE(system2.new_mode(q2|"s2"));
    CompositeHybridAutomaton system({system1,system2});
    DiscreteLocation loc({q1|"s1",q2|"s2"});
    ARIADNE_TEST_ASSERT(system.has_mode(loc));
    ARIADNE_TEST_ASSERT(system.has_partial_mode(loc));
    ARIADNE_TEST_PRINT(system.mode(loc));
    DiscreteLocation loce({q1|"s1",q2|"s2",q3|"s3"});
    ARIADNE_TEST_ASSERT(not system.has_mode(loce));
    ARIADNE_TEST_ASSERT(system.has_partial_mode(loce));
    ARIADNE_TEST_PRINT(system.mode(loce));
}

Void
TestHybridAutomaton::test_overspecified_dynamic()
{
    HybridAutomaton system;
    RealVariable x("x");
    RealVariable y("y");
    ARIADNE_TEST_EXECUTE( system.new_mode( {let(x)=y},{dot(y)=x} ) );
    ARIADNE_TEST_THROWS( system.new_mode( {let(x)=1,let(x)=x+2} ),SystemSpecificationError);
    ARIADNE_TEST_THROWS( system.new_mode( {dot(x)=1,dot(x)=0} ),SystemSpecificationError);
    ARIADNE_TEST_THROWS( system.new_mode( {let(x)=1},{dot(x)=0} ),SystemSpecificationError);
}

Void
TestHybridAutomaton::test_algebraic_loops()
{
    HybridAutomaton system;
    RealVariable x("x");
    RealVariable y("y");
    RealVariable z("z");
    ARIADNE_TEST_THROWS( HybridAutomaton().new_mode( {let(x)=y+1,let(y)=x} ),AlgebraicLoopError);
    ARIADNE_TEST_THROWS( HybridAutomaton().new_mode( {let(x)=2*x+1} ),AlgebraicLoopError);
    ARIADNE_TEST_EXECUTE( HybridAutomaton().new_mode( {let(x)=y,let(y)=z} ) );
    ARIADNE_TEST_EXECUTE( HybridAutomaton().new_mode( {let(x)=y,let(y)=z},{dot(z)=x+y+z} ) );
}

Void
TestHybridAutomaton::test_nonexistent_mode()
{
    HybridAutomaton system;
    StringVariable q("q");
    StringVariable r("r");
    DiscreteEvent e("e");
    RealVariable x("x");
    ARIADNE_TEST_EXECUTE( system.new_mode( {q|"1"},{dot(x)=0} ) );
    ARIADNE_TEST_EXECUTE( system.new_mode( {q|"2",r|"1"} ) );
    ARIADNE_TEST_EXECUTE( system.new_mode( {q|"2",r|"2"} ) );
    ARIADNE_TEST_PRINT( system );
    ARIADNE_TEST_THROWS( system.new_action( {q|"2"}, e, x>=0 ), NonExistentModeError );
    ARIADNE_TEST_THROWS( system.new_transition( {q|"1"}, e, {q|"3"} ), NonExistentModeError );
    ARIADNE_TEST_THROWS( system.new_transition( {q|"3"}, e, {q|"1"} ), NonExistentModeError );
    ARIADNE_TEST_THROWS( system.new_transition( {q|"1"}, e, {q|"1",r|"2"} ), NonExistentModeError );
    ARIADNE_TEST_EXECUTE( system.new_transition( {q|"1"}, e, {q|"2",r|"1"} ) );
}

Void
TestHybridAutomaton::test_multiple_guard()
{
    HybridAutomaton system;
    DiscreteLocation q1(1);
    DiscreteLocation q2(2);
    DiscreteEvent e1(1);
    DiscreteEvent e2(2);
    RealVariable x("x");
    RealVariable y("y");
    Dyadic c(0.125);
    ARIADNE_TEST_EXECUTE( system.new_mode( q1,(dot(x)=1,dot(y)=1) ) );
    ARIADNE_TEST_EXECUTE( system.new_mode( q2,(dot(x)=2) ) );
    ARIADNE_TEST_EXECUTE( system.new_action( q1, x<=c, e1, x>=0 ) );
    ARIADNE_TEST_EXECUTE( system.new_action( q2, x<=c, e1, x>=0 ) );
    ARIADNE_TEST_EXECUTE( system.new_action( q1, y<=c, e2, y>=0 ) );
    ARIADNE_TEST_THROWS( system.new_action( q1, x+y<=c, e1, x+y>=0 ), MultipleGuardError );
}

Void
TestHybridAutomaton::test_multiple_transition()
{
    HybridAutomaton system;
    DiscreteLocation q1(1);
    DiscreteLocation q2(2);
    DiscreteEvent e1(1);
    DiscreteEvent e2(2);
    RealVariable x("x");
    ARIADNE_TEST_EXECUTE( system.new_mode( q1,(dot(x)=1) ) );
    ARIADNE_TEST_EXECUTE( system.new_mode( q2,(dot(x)=2) ) );
    ARIADNE_TEST_EXECUTE( system.new_transition( q1, e1, q1, (prime(x)=x) ) );
    ARIADNE_TEST_EXECUTE( system.new_transition( q2, e1, q1, (prime(x)=x) ) );
    ARIADNE_TEST_EXECUTE( system.new_transition( q1, e2, q1, (prime(x)=x) ) );
    ARIADNE_TEST_THROWS( system.new_transition( q1, e1, q2, (prime(x)=x) ), MultipleTransitionError );
}

Void
TestHybridAutomaton::test_overspecified_reset()
{
    HybridAutomaton system;
    DiscreteLocation q1(1);
    DiscreteEvent e1(1);
    DiscreteEvent e2(2);
    DiscreteEvent e3(3);
    DiscreteEvent e4(4);
    DiscreteEvent e5(5);
    RealVariable x("x");
    RealVariable y("y");
    RealVariable z("z");
    ARIADNE_TEST_EXECUTE( system.new_mode( q1,{let(x)=1}, {dot(y)=1} ) );
    ARIADNE_TEST_EXECUTE( system.new_transition( q1, e1, q1 ) ); // OK, underspecified
    ARIADNE_TEST_EXECUTE( system.new_transition( q1, e2, q1, {prime(y)=y+1} ) );
    ARIADNE_TEST_EXECUTE( system.new_transition( q1, e3, q1, {prime(y)=y+1,prime(z)=x} ) ); // OK, z in another component
    ARIADNE_TEST_THROWS( system.new_transition( q1, e4, q1, {prime(x)=x+1} ), OverspecifiedResetError );
    ARIADNE_TEST_THROWS( system.new_transition( q1, e5, q1, {prime(y)=2*x+4,prime(y)=2*(x+2)} ), OverspecifiedResetError );
}

Void
TestHybridAutomaton::test_underspecified_mode()
{
    HybridAutomaton system;
    DiscreteLocation qs(0);
    DiscreteLocation q1(1);
    DiscreteLocation q2(2);
    DiscreteLocation q3(3);
    DiscreteLocation q4(4);
    DiscreteLocation q5(5);
    DiscreteEvent e("e");
    RealVariable x("x");
    RealVariable y("y");
    RealVariable z("z");

    // The following mode should be valid
    ARIADNE_TEST_EXECUTE( system.new_mode( qs,{let(x)=y+1},{dot(y)=x+y} ) );
    ARIADNE_TEST_EXECUTE( system.new_action( qs, x+y<=1, e, x+y>=0 ) );
    ARIADNE_TEST_EXECUTE( system.new_update( qs, e, qs, {prime(y)=x+y} ) );
    ARIADNE_TEST_EXECUTE( system.check_mode( qs ) );
    // Test invalidity of the following modes
    ARIADNE_TEST_EXECUTE( system.new_mode( q1,{let(x)=y} ) );
    ARIADNE_TEST_THROWS( system.check_mode( q1 ), UnderspecifiedDynamicError );
    ARIADNE_TEST_EXECUTE( system.new_mode( q2,{dot(x)=x+y} ) );
    ARIADNE_TEST_THROWS( system.check_mode( q2 ), UnderspecifiedDynamicError );
    ARIADNE_TEST_EXECUTE( system.new_mode( q3,{dot(x)=x} ) );
    ARIADNE_TEST_EXECUTE( system.new_action( q3, y<=1, e, x>=0 ) );
    ARIADNE_TEST_THROWS( system.check_mode( q3 ), UnderspecifiedConstraintError );
    ARIADNE_TEST_EXECUTE( system.new_mode( q4,{dot(x)=x} ) );
    ARIADNE_TEST_EXECUTE( system.new_action( q4, x<=1, e, y>=0 ) );
    ARIADNE_TEST_THROWS( system.check_mode( q4 ), UnderspecifiedConstraintError );
    ARIADNE_TEST_EXECUTE( system.new_mode( q5,{dot(x)=x} ) );
    ARIADNE_TEST_EXECUTE( system.new_update( q5, e, q5, {prime(x)=y} ) );
    ARIADNE_TEST_THROWS( system.check_mode( q5 ), UnderspecifiedResetError );
}

Void
TestHybridAutomaton::test_underspecified_reset()
{
    HybridAutomaton system;
    StringVariable q("q");
    DiscreteLocation q1(q|"1");
    DiscreteLocation q2(q|"2");
    DiscreteLocation q3(q|"3");
    DiscreteLocation qt(q|"t");
    DiscreteEvent e("e");
    RealVariable x("x");
    RealVariable y("y");
    RealVariable z("z");

    ARIADNE_TEST_EXECUTE( system.new_mode( qt,{let(x)=1},{dot(y)=1} ) );
    ARIADNE_TEST_EXECUTE( system.new_mode( q1,{dot(x)=x} ) );
    ARIADNE_TEST_EXECUTE( system.new_mode( q2,{dot(x)=x} ) );
    ARIADNE_TEST_EXECUTE( system.new_mode( q3,{dot(x)=x} ) );
    ARIADNE_TEST_EXECUTE( system.new_update( q1, e, qt, {next(y)=1} ) );
    ARIADNE_TEST_EXECUTE( system.check_mode( q1 ) );
    ARIADNE_TEST_EXECUTE( system.new_update( q2, e, qt ) );
    ARIADNE_TEST_THROWS( system.check_mode( q2 ), UnderspecifiedResetError );
    ARIADNE_TEST_EXECUTE( system.new_update( q3, e, qt, {next(y)=1,prime(z)=1} ) );
    ARIADNE_TEST_THROWS( system.check_mode( q3 ), UnderspecifiedDynamicError );
}



Void
TestHybridAutomaton::test_build_hybrid_system()
{
    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0_dec);
    RealConstant hmin("hmin",5.5_dec);
    RealConstant hmax("hmax",8.0_dec);
    RealConstant delta("delta",0.05_dec);
    RealConstant lambda("lambda",0.02_dec);
    RealConstant rate("rate",0.3_dec);

    // Declare the system variables
    RealVariable height("height");
    RealVariable alpha("alpha");

    // Create the tank object
    HybridAutomaton tank_automaton;
    DiscreteLocation tank;
    // The water level is always given by the same dynamic
    // The inflow is controlled by the valve alpha, the outflow depends on the
    // pressure, which is proportional to the water height.
    tank_automaton.new_mode(tank,{dot(height)=-lambda*height+rate*alpha});

    ARIADNE_TEST_PRINT(tank_automaton);

    // Describe the valve model

    // Declare the events we use in the valve automaton
    DiscreteEvent start_opening("start_opening");
    DiscreteEvent start_closing("start_closing");
    DiscreteEvent finished_opening("finished_opening");
    DiscreteEvent finished_closing("finished_closing");

    // Declare the locations we use
    String open("open");
    String opening("opening");
    String closed("closed");
    String closing("closing");

    RealVariable beta("beta");

    StringVariable valve("valve");
    HybridAutomaton valve_automaton;

    // Since alpha is a known constant when the valve is open or closed,
    // specify alpha by an algebraic equation.
    valve_automaton.new_mode({valve|open},{let(alpha)=1});
    valve_automaton.new_mode({valve|closed},{let(alpha)=-1});
    // Specify the differential equation for how the valve opens/closes.
    valve_automaton.new_mode({valve|opening},{dot(alpha)=1/T});
    valve_automaton.new_mode({valve|closing},{dot(alpha)=-1/T});

    ARIADNE_TEST_PRINT(valve_automaton);

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    //valve_automaton.new_invariant(open,start_closing,height<=hmax || (height>=hmin && !(height<=hmin+1)));
    valve_automaton.new_invariant({valve|open},height<=hmax,start_closing);
    valve_automaton.new_invariant({valve|opening},height<=hmax,start_closing);
    valve_automaton.new_invariant({valve|opening},alpha<=1,finished_opening);
    valve_automaton.new_invariant({valve|closed},height>=hmin,start_opening);
    valve_automaton.new_invariant({valve|closing},height>=hmin,start_opening);
    valve_automaton.new_invariant({valve|closing},alpha>=0,finished_closing);

    valve_automaton.new_transition({valve|closed},start_opening,{valve|opening},{next(alpha)=alpha},height<=hmin,EventKind::PERMISSIVE);
    valve_automaton.new_transition({valve|closing},start_opening,{valve|opening},{next(alpha)=alpha},height<=hmin,EventKind::PERMISSIVE);
    valve_automaton.new_transition({valve|open},start_closing,{valve|closing},{next(alpha)=alpha},height>=hmax,EventKind::PERMISSIVE);
    valve_automaton.new_transition({valve|opening},start_closing,{valve|closing},{next(alpha)=alpha},height>=hmax,EventKind::PERMISSIVE);

    // Set the transitions for when the valve finished opening.
    // Since alpha is defined by an algebraic equation in the new mode,
    // it may not be specified in the reset.
    valve_automaton.new_transition({valve|opening},finished_opening,{valve|open},alpha>=1,EventKind::PERMISSIVE);
    valve_automaton.new_transition({valve|closing},finished_closing,{valve|closed},alpha<=0,EventKind::PERMISSIVE);

    ARIADNE_TEST_PRINT(valve_automaton);

    CompositeHybridAutomaton watertank_system({tank_automaton,valve_automaton});
    std::cout << "watertank_system:\n" << watertank_system << "\n";

    _system=watertank_system;
    _system.verbosity=0;

}



Void
TestHybridAutomaton::test_static_analysis()
{
    StringVariable valve("valve");
    String opening("opening");
    DiscreteLocation valve_opening={valve|opening};
    DiscreteLocation valve_open={valve|"open"};
    ARIADNE_TEST_PRINT(valve_opening);

    CompositeHybridAutomaton const& system = _system;
    ARIADNE_TEST_PRINT (system.mode(valve_opening));

    HybridAutomaton const& tank_automaton = _system.component(0);
    HybridAutomaton const& valve_automaton = _system.component(1);
    ARIADNE_TEST_PRINT(tank_automaton);
    ARIADNE_TEST_PRINT(valve_automaton);

    ARIADNE_TEST_EXECUTE(_system.discrete_reachability(valve_opening));

    ARIADNE_TEST_EXECUTE(_system.check_mode(valve_opening));
    ARIADNE_TEST_EXECUTE(_system.check_mode(valve_open));
    ARIADNE_TEST_EXECUTE(_system.check_reachable_modes(valve_opening));
    ARIADNE_TEST_EXECUTE(_system.discrete_reachability(valve_opening));
}



Void
TestHybridAutomaton::test_build_intensional_hybrid_automaton()
{
    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0_dec);
    RealConstant hmin("hmin",5.5_dec);
    RealConstant hmax("hmax",8.0_dec);
    RealConstant delta("delta",0.05_dec);
    RealConstant lambda("lambda",0.02_dec);
    RealConstant rate("rate",0.3_dec);

    // Declare the system variables
    RealVariable height("height");
    RealVariable alpha("alpha");

    // Describe the valve model

    // Declare the events we use in the valve automaton
    DiscreteEvent start_opening("start_opening");
    DiscreteEvent start_closing("start_closing");
    DiscreteEvent finished_opening("finished_opening");
    DiscreteEvent finished_closing("finished_closing");

    // Declare the locations we use
    StringVariable valve("valve");
    StringConstant open("open");
    StringConstant opening("opening");
    StringConstant closed("closed");
    StringConstant closing("closing");

    RealVariable beta("beta");

    HybridAutomaton valve_automaton;

    // Since alpha is a known constant when the valve is open or closed,
    // specify alpha by an algebraic equation.
    valve_automaton.new_mode({valve|open},(let(alpha)=1));
    valve_automaton.new_mode({valve|closed},(let(alpha)=-1));
    // Specify the differential equation for how the valve opens/closes.
    valve_automaton.new_mode({valve|opening},{dot(alpha)=1/T});
    valve_automaton.new_mode({valve|closing},{dot(alpha)=-1/T});

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    //valve_automaton.new_invariant(open,start_closing,height<=hmax || (height>=hmin && !(height<=hmin+1)));
    valve_automaton.new_action({valve|open},height<=hmax,start_closing,height>=hmax-delta,EventKind::PERMISSIVE);
    valve_automaton.new_action({valve|opening},height<=hmax,start_closing,height>=hmax-delta,EventKind::PERMISSIVE);
    valve_automaton.new_action({valve|closed},height>=hmin,start_opening,height<=hmin+delta,EventKind::PERMISSIVE);
    valve_automaton.new_action({valve|opening},height>=hmin,start_closing,height<=hmin+delta,EventKind::PERMISSIVE);
    valve_automaton.new_action({valve|closing},finished_closing,alpha>=0,EventKind::URGENT);
    valve_automaton.new_action({valve|opening},finished_opening,alpha<=1,EventKind::URGENT);

    valve_automaton.new_update({valve|closed},start_opening,{valve|opening},{next(alpha)=alpha});
    valve_automaton.new_update({valve|closing},start_opening,{valve|opening},{next(alpha)=alpha});
    valve_automaton.new_update({valve|open},start_closing,{valve|closing},{next(alpha)=alpha});
    valve_automaton.new_update({valve|opening},start_closing,{valve|closing},{next(alpha)=alpha});

    // Set the transitions for when the valve finished opening.
    // Since alpha is defined by an algebraic equation in the new mode,
    // it may not be specified in the reset.
    valve_automaton.new_update({valve|opening},finished_opening,{valve|open});
    valve_automaton.new_update({valve|closing},finished_closing,{valve|closed});

    ARIADNE_TEST_PRINT(valve_automaton);

}



Int main() {
    //ARIADNE_TEST_CALL(TestHybridAutomaton().test_build_atomic_hybrid_automaton());
    ARIADNE_TEST_CALL(TestHybridAutomaton().test());

    //ARIADNE_TEST_CALL(TestHybridSystem().test());
    return ARIADNE_TEST_FAILURES;
}

