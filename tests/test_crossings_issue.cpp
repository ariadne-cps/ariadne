#include "ariadne.hpp"
#include "test.hpp"

using namespace std;
using namespace Ariadne;

class TestCrossingsIssue
{
  public:
	TestCrossingsIssue(unsigned int verbosity);
    void test();
  private:
    unsigned int _verbosity;
    void test_crossings();
};

TestCrossingsIssue::TestCrossingsIssue(unsigned int verbosity) : _verbosity(verbosity) { }

void TestCrossingsIssue::test()
{
    ARIADNE_TEST_CALL(test_crossings());
}

void TestCrossingsIssue::test_crossings()
{
    /// Build the Hybrid System

    RealConstant v("v",0.092_dec);
    RealConstant L("L",0.00025_dec);
    RealConstant M("M",0.0023_dec);

    /// Create a HybridAutomaton object
    AtomicHybridAutomaton automaton("compute_crossings_issue");

    /// Modes

    AtomicDiscreteLocation far("far");
    AtomicDiscreteLocation close("close");

    // Variables

    RealVariable x("x"); // X position
    TimeVariable t;

    automaton.new_mode(far, {dot(x)=v});
    automaton.new_mode(close, {dot(x)=v});

    DiscreteEvent comes("comes");
    DiscreteEvent leaves("leaves");

    automaton.new_transition(far,comes,close,{next(x)=x},sqr(x-M)<=sqr(L));
	automaton.new_transition(close,leaves,far,{next(x)=x},sqr(x-M)>=sqr(L));

    // Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver(automaton);
    evolver.verbosity = _verbosity;

    // Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(0.5);
    evolver.configuration().set_maximum_step_size(0.0001);

    typedef GeneralHybridEvolver::OrbitType OrbitType;

    HybridSet initial_set({automaton|far},{x==0});

    Nat max_n = 3; Real max_t = 0.046_dec;
    HybridTime evolution_time(max_t,max_n);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    std::cout << "done." << std::endl;
/*
    Real max_x = max_t * v;
    HybridBoxSet guard_set(automaton|far,{t.in(0,2*max_t),x.in(M-L,M+L)});
    plot("crossings_issue",Axes2d(0.0,t,1.1*max_t.get_d(), 0.0,x,1.1*max_x.get_d()), Colour(1.0,0.5,0.0), guard_set, Colour(0.0,0.5,1.0), orbit);
*/

    for (auto enclosure : orbit.final()) {
    	ARIADNE_TEST_ASSERT(enclosure.previous_events().size() == 2);
    }
}

int main(Int argc, const char* argv[]) {

    auto verbosity=get_verbosity(argc,argv);

    TestCrossingsIssue(verbosity).test();
    return ARIADNE_TEST_FAILURES;
}

