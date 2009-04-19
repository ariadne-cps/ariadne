#include "../test/test.h"

#include "python_hybrid_evolver.h"
#include "taylor_set.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"

namespace Ariadne {


} // namespace Ariadne

using namespace Ariadne;

int main() {

    DiscreteState q1(1);
    HybridAutomaton system;
    HybridTaylorSet initial(q1,Vector<Interval>(2, -0.125,0.125, 0,0.25));
    HybridTime time(1.0,2);
    PythonHybridEvolver evolver;
    ListSet<HybridTaylorSet> final;


    evolver.evolve(system,initial,time);

}