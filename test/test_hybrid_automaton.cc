/***************************************************************************
 *            test_hybrid_automaton.cc
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
#include <fstream>

#include "test.h"

#include "expression.h"
#include "valuation.h"
#include "hybrid_automaton.h"
#include <boost/concept_check.hpp>

using namespace std;
using namespace Ariadne;

Set<DiscreteEvent> operator,(const DiscreteEvent e1, const DiscreteEvent e2) { Set<DiscreteEvent> r; r.insert(e1); r.insert(e2); return r; }
Set<DiscreteEvent> operator,(const Set<DiscreteEvent> s, const DiscreteEvent e) { Set<DiscreteEvent> r(s); r.insert(e); return r; }


class TestHybridSystem {
  public:
    void test();
  private:
    void test_build_hybrid_system();
    void test_static_analysis();
  private:
    HybridSystem _system;
};

void
TestHybridSystem::test()
{
    std::clog<<std::boolalpha;
    std::cerr<<std::boolalpha;
    ARIADNE_TEST_CALL(test_build_hybrid_system());
    ARIADNE_TEST_CALL(test_static_analysis());
}

void
TestHybridSystem::test_build_hybrid_system()
{
    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0);
    RealConstant hmin("hmin",5.5);
    RealConstant hmax("hmax",8.0);
    RealConstant delta("delta",0.05);
    RealConstant lambda("lambda",0.02);
    RealConstant rate("rate",0.3);

    // Declare the system variables
    RealVariable height("height");
    RealVariable alpha("alpha");

    // Create the tank object
    HybridSystem tank_component;
    StringVariable q("q"); StringVariable t("t");
    // The water level is always given by the same dynamic
    // The inflow is controlled by the valve alpha, the outflow depends on the
    // pressure, which is proportional to the water height.
    tank_component.new_dynamic(DiscretePredicate(true),(dot(height)=-lambda*height+rate*alpha));
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

    HybridSystem valve_component;

    valve_component.disable_events((valve==open || valve==opening),complement((start_closing,finished_opening)));
    valve_component.disable_events((valve==closed || valve==closing),complement((start_opening,finished_closing)));
    valve_component.disable_events((valve==open),finished_opening);
    valve_component.disable_events((valve==closed),finished_closing);

    // Since alpha is a known constant when the valve is open or closed,
    // specify alpha by an algebraic equation.
    valve_component.new_auxiliary(valve==open,(alpha=1.0));
    valve_component.new_auxiliary(valve==closed,(alpha=-1.0));
    // Specify the differential equation for how the valve opens/closes.
    valve_component.new_dynamic(valve==opening,(dot(alpha)=1.0/T));
    valve_component.new_dynamic(valve==closing,(dot(alpha)=-1.0/T));

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    //valve_automaton.new_invariant(open,start_closing,height<=hmax || (height>=hmin && !(height<=hmin+1)));
    valve_component.new_invariant(valve==open,start_closing,height<=hmax);
    valve_component.new_invariant(valve==opening,start_closing,height<=hmax);
    valve_component.new_invariant(valve==opening,finished_opening,alpha<=1.0);
    valve_component.new_invariant(valve==closed,start_opening,height>=hmin);
    valve_component.new_invariant(valve==closing,start_opening,height>=hmin);
    valve_component.new_invariant(valve==closing,finished_closing,alpha>=0.0);

    valve_component.new_transition(valve==closed,start_opening,(next(valve)=opening),(next(alpha)=alpha),height<=hmin);
    valve_component.new_transition(valve==closing,start_opening,next(valve)=opening,(next(alpha)=alpha),height<=hmin);
    valve_component.new_transition(valve==open,start_closing,(next(valve)=closing),(next(alpha)=alpha),height>=hmax);
    valve_component.new_transition(valve==opening,start_closing,(next(valve)=closing),(next(alpha)=alpha),height>=hmax);

    // Set the transitions for when the valve finished opening.
    // Since alpha is defined by an algebraic equation in the new mode,
    // it may not be specified in the reset.
    valve_component.new_transition(valve==opening,finished_opening,next(valve)=open,alpha>=1.0);
    valve_component.new_transition(valve==closing,finished_closing,next(valve)=closed,alpha<=0.0);

    ARIADNE_TEST_PRINT(valve_component.mode((valve|opening)));
    ARIADNE_TEST_PRINT(valve_component.mode((valve|closing)));
    ARIADNE_TEST_PRINT(valve_component.mode((valve|open)));
    ARIADNE_TEST_PRINT(valve_component.mode((valve|closed)));
    ARIADNE_TEST_PRINT(valve_component);

    HybridSystem watertank_system=compose((tank_component,valve_component));
    std::cout << "watertank_system:\n" << watertank_system << "\n";

    _system=watertank_system;
    _system.verbosity=0;

    Set<DiscreteLocation> reachable = valve_component.reachable_locations((valve|open));
    std::cerr << reachable;
    for(Set<DiscreteLocation>::const_iterator iter=reachable.begin(); iter!=reachable.end(); ++iter) {
        std::cerr<<valve_component.mode(*iter)<<"\n";
        valve_component.check_mode(*iter);
    }

}





void TestHybridSystem::test_static_analysis()
{
    DiscreteLocation valve_opening(StringVariable("valve")|"opening");
    ARIADNE_TEST_EXECUTE(this->_system.check_mode(valve_opening));
    //ARIADNE_TEST_EXECUTE(_system.check_mode((tank|draining,valve|open)));
    //ARIADNE_TEST_EXECUTE(_system.check_reachable_modes((tank|draining,valve|opening)));
    //ARIADNE_TEST_EXECUTE(_system.discrete_reachability((tank|draining,valve|opening)));

    //DiscreteLocation invalid_location=(AtomicDiscreteLocation("nonexistent"),AtomicDiscreteLocation("opening"));
    //ARIADNE_TEST_FAIL(_system.discrete_reachability(invalid_location));
}



class TestHybridAutomaton {
  public:
    void test();

    void test_build_hybrid_system();
    void test_build_atomic_hybrid_automaton();
    void test_build_intensional_hybrid_automaton();
    void test_static_analysis();
  private:
    CompositeHybridAutomaton _system;
};

void
TestHybridAutomaton::test()
{
    std::clog<<std::boolalpha;
    std::cerr<<std::boolalpha;
    ARIADNE_TEST_CALL(test_build_hybrid_system());
    ARIADNE_TEST_CALL(test_build_atomic_hybrid_automaton());
    ARIADNE_TEST_CALL(test_static_analysis());
}


void
TestHybridAutomaton::test_build_hybrid_system()
{
    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0);
    RealConstant hmin("hmin",5.5);
    RealConstant hmax("hmax",8.0);
    RealConstant delta("delta",0.05);
    RealConstant lambda("lambda",0.02);
    RealConstant rate("rate",0.3);

    // Declare the system variables
    RealVariable height("height");
    RealVariable alpha("alpha");

    // Create the tank object
    HybridAutomaton tank_automaton;
    DiscreteLocation tank;
    // The water level is always given by the same dynamic
    // The inflow is controlled by the valve alpha, the outflow depends on the
    // pressure, which is proportional to the water height.
    tank_automaton.new_mode(tank,(dot(height)=-lambda*height+rate*alpha));

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
    valve_automaton.new_mode((valve|open),(alpha=1.0));
    valve_automaton.new_mode((valve|closed),(alpha=-1.0));
    // Specify the differential equation for how the valve opens/closes.
    valve_automaton.new_mode((valve|opening),(dot(alpha)=1.0/T));
    valve_automaton.new_mode((valve|closing),(dot(alpha)=-1.0/T));

    ARIADNE_TEST_PRINT(valve_automaton);

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    //valve_automaton.new_invariant(open,start_closing,height<=hmax || (height>=hmin && !(height<=hmin+1)));
    valve_automaton.new_invariant((valve|open),height<=hmax,start_closing);
    valve_automaton.new_invariant((valve|opening),height<=hmax,start_closing);
    valve_automaton.new_invariant((valve|opening),alpha<=1.0,finished_opening);
    valve_automaton.new_invariant((valve|closed),height>=hmin,start_opening);
    valve_automaton.new_invariant((valve|closing),height>=hmin,start_opening);
    valve_automaton.new_invariant((valve|closing),alpha>=0.0,finished_closing);

    valve_automaton.new_transition((valve|closed),start_opening,(valve|opening),(next(alpha)=alpha),height<=hmin,PERMISSIVE);
    valve_automaton.new_transition((valve|closing),start_opening,(valve|opening),(next(alpha)=alpha),height<=hmin,PERMISSIVE);
    valve_automaton.new_transition((valve|open),start_closing,(valve|closing),(next(alpha)=alpha),height>=hmax,PERMISSIVE);
    valve_automaton.new_transition((valve|opening),start_closing,(valve|closing),(next(alpha)=alpha),height>=hmax,PERMISSIVE);

    // Set the transitions for when the valve finished opening.
    // Since alpha is defined by an algebraic equation in the new mode,
    // it may not be specified in the reset.
    valve_automaton.new_transition((valve|opening),finished_opening,(valve|open),alpha>=1.0,PERMISSIVE);
    valve_automaton.new_transition((valve|closing),finished_closing,(valve|closed),alpha<=0.0,PERMISSIVE);

    ARIADNE_TEST_PRINT(valve_automaton);

    CompositeHybridAutomaton watertank_system((tank_automaton,valve_automaton));
    std::cout << "watertank_system:\n" << watertank_system << "\n";

    _system=watertank_system;
    _system.verbosity=0;

}



void
TestHybridAutomaton::test_static_analysis()
{
    StringVariable valve("valve");
    String opening("opening");
    DiscreteLocation valve_opening(valve|opening);
    DiscreteLocation valve_open(valve|"open");
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

void
TestHybridAutomaton::test_build_atomic_hybrid_automaton()
{
    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0);
    RealConstant hmin("hmin",5.5);
    RealConstant hmax("hmax",8.0);
    RealConstant delta("delta",0.05);
    RealConstant lambda("lambda",0.02);
    RealConstant rate("rate",0.3);

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
    StringConstant open("open");
    StringConstant opening("opening");
    StringConstant closed("closed");
    StringConstant closing("closing");

    RealVariable beta("beta");

    AtomicHybridAutomaton valve_automaton("valve");

    // Since alpha is a known constant when the valve is open or closed,
    // specify alpha by an algebraic equation.
    valve_automaton.new_mode(open,(alpha=1.0));
    valve_automaton.new_mode(closed,(alpha=-1.0));
    // Specify the differential equation for how the valve opens/closes.
    valve_automaton.new_mode(opening,(dot(alpha)=1.0/T));
    valve_automaton.new_mode(closing,(dot(alpha)=-1.0/T));

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    //valve_automaton.new_invariant(open,start_closing,height<=hmax || (height>=hmin && !(height<=hmin+1)));
    valve_automaton.new_invariant(open,height<=hmax,start_closing);
    valve_automaton.new_invariant(opening,height<=hmax,start_closing);
    valve_automaton.new_invariant(opening,alpha<=1.0,finished_opening);
    valve_automaton.new_invariant(closed,height>=hmin,start_opening);
    valve_automaton.new_invariant(closing,height>=hmin,start_opening);
    valve_automaton.new_invariant(closing,alpha>=0.0,finished_closing);

    valve_automaton.new_transition(closed,start_opening,opening,(next(alpha)=alpha),height<=hmin,PERMISSIVE);
    valve_automaton.new_transition((closing),start_opening,(opening),(next(alpha)=alpha),height<=hmin,PERMISSIVE);
    valve_automaton.new_transition(open,start_closing,(closing),(next(alpha)=alpha),height>=hmax,PERMISSIVE);
    valve_automaton.new_transition((opening),start_closing,(closing),(next(alpha)=alpha),height>=hmax,PERMISSIVE);

    // Set the transitions for when the valve finished opening.
    // Since alpha is defined by an algebraic equation in the new mode,
    // it may not be specified in the reset.
    valve_automaton.new_transition((opening),finished_opening,open,alpha>=1.0,PERMISSIVE);
    valve_automaton.new_transition((closing),finished_closing,(closed),alpha<=0.0,PERMISSIVE);

    ARIADNE_TEST_PRINT(valve_automaton);

}


void
TestHybridAutomaton::test_build_intensional_hybrid_automaton()
{
    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0);
    RealConstant hmin("hmin",5.5);
    RealConstant hmax("hmax",8.0);
    RealConstant delta("delta",0.05);
    RealConstant lambda("lambda",0.02);
    RealConstant rate("rate",0.3);

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
    valve_automaton.new_mode((valve|open),(alpha=1.0));
    valve_automaton.new_mode((valve|closed),(alpha=-1.0));
    // Specify the differential equation for how the valve opens/closes.
    valve_automaton.new_mode((valve|opening),(dot(alpha)=1.0/T));
    valve_automaton.new_mode((valve|closing),(dot(alpha)=-1.0/T));

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    //valve_automaton.new_invariant(open,start_closing,height<=hmax || (height>=hmin && !(height<=hmin+1)));
    valve_automaton.new_action((valve|open),height<=hmax,start_closing,height>=hmax-delta,PERMISSIVE);
    valve_automaton.new_action((valve|opening),height<=hmax,start_closing,height>=hmax-delta,PERMISSIVE);
    valve_automaton.new_action((valve|closed),height>=hmin,start_opening,height<=hmin+delta,PERMISSIVE);
    valve_automaton.new_action((valve|opening),height>=hmin,start_closing,height<=hmin+delta,PERMISSIVE);
    valve_automaton.new_action((valve|closing),finished_closing,alpha>=0.0,URGENT);
    valve_automaton.new_action((valve|opening),finished_opening,alpha<=1.0,URGENT);

    valve_automaton.new_update((valve|closed),start_opening,(valve|opening),(next(alpha)=alpha));
    valve_automaton.new_update((valve|closing),start_opening,(valve|opening),(next(alpha)=alpha));
    valve_automaton.new_update((valve|open),start_closing,(valve|closing),(next(alpha)=alpha));
    valve_automaton.new_update((valve|opening),start_closing,(valve|closing),(next(alpha)=alpha));

    // Set the transitions for when the valve finished opening.
    // Since alpha is defined by an algebraic equation in the new mode,
    // it may not be specified in the reset.
    valve_automaton.new_update((valve|opening),finished_opening,(valve|open));
    valve_automaton.new_update((valve|closing),finished_closing,(valve|closed));

    ARIADNE_TEST_PRINT(valve_automaton);

}



int main() {
    ARIADNE_TEST_CALL(TestHybridAutomaton().test_build_atomic_hybrid_automaton());
    ARIADNE_TEST_CALL(TestHybridAutomaton().test());

    //ARIADNE_TEST_CALL(TestHybridSystem().test());
    return ARIADNE_TEST_FAILURES;
}

