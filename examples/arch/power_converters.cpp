#include "power_converters.hpp"

namespace Ariadne {

using std::cout;

void test() {

    HybridAutomaton buck = make_buck_system();

    cout << buck << "\n";

    TimeVariable t;

    DiscreteLocation initial_location={swtch|open};
    RealVariablesBox initial_box={V.in(0.0_dec,0.0_dec),I.in(1.0_dec,1.0_dec),tau.in(0.0_dec,0.0_dec)};
    HybridBoxSet initial_set(initial_location,initial_box);

    Integer n=20;
    Real tmax=n*T;
    GeneralHybridEvolver evolver(buck);
    HybridTime evolution_time(tmax,2*n);

    auto orbit=evolver.orbit(initial_set,evolution_time,UPPER_SEMANTICS);

    plot("buck",{0<=V<=1,-1<=I<=1},Colour(.5,.0,.5),orbit.reach());
    plot("buck-V",{0<=t<=tmax,0<=V<=1},Colour(.5,.0,.5),orbit.reach());
    plot("buck-I",{0<=t<=tmax,0<=I<=1},Colour(.5,.0,.5),orbit.reach());

}

} // namespace Ariadne

int main() {
    Ariadne::test();
}
