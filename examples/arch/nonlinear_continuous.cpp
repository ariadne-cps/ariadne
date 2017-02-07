#include "utility/logging.hpp"
#include "nonlinear_continuous.hpp"
#include "dynamics/reachability_analyser.hpp"

namespace Ariadne {

using std::cout;

using Axis = ApproximateDoubleVariableInterval;

int verbosity;

void test() {

//    auto problem = make_fitzhugh_nagumo_problem();
    auto problem = make_dai_problem();

    auto initial_set = problem.initial_set;
    auto system = problem.system;
    auto unsafe_set = problem.unsafe_set;




    cout << "initial_set=" << initial_set << "\n";
    cout << "system=" << system << "\n";

    TimeVariable t;


//    Box<RealInterval> initial_box_set({{-0.625_dec,-0.5_dec},{-2.0_dec,-1.875_dec}});
//    Box<RealInterval> initial_box_set({{1.000_dec,1.001_dec},{-2.000_dec,-1.999_dec}});
//    Box<RealInterval> initial_box_set({{-0.500_dec,-0.499_dec},{-2.000_dec,-1.999_dec}});
    Box<RealInterval> initial_box_set({{-0.1_dec,+0.1_dec},{-2.000_dec,-1.999_dec}});
    BoundedConstraintSet initial_constraint_set(initial_box_set);

    Real tmax=1.00_dec;

    MaximumError max_err=0.01;

    TaylorSeriesIntegrator integrator(max_err);
    integrator.set_maximum_step_size(8);
    auto& factory=integrator.function_factory();
    VectorFieldEvolver evolver(system,integrator);
    evolver.verbosity=1;

    Precision64 pr;

    for(Dyadic rad(1,8u); rad<=4; rad=2*rad) {
        Box<Interval<Float64Value>> init(Box<Interval<Dyadic>>({{-1-rad,-1+rad},{-2-rad,-2+rad}}),pr);
        Float64 step(8.0);
        auto phi=integrator.flow_step(system.function(),init,step);
        FloatError<Precision64>::set_output_places(9);
        std::cerr<<"rad="<<rad.get_d()<<", step="<<step<<"\nphi.errors()="<<phi.errors()<<"\n\n";
        //std::cerr<<"phi="<<phi<<"\n";
    }
    return;

    cout << "initial_box_set=" << initial_box_set << "\n";
    auto orbit=evolver.orbit(evolver.enclosure(initial_box_set),tmax,UPPER_SEMANTICS);

    Figure g; g.set_bounding_box({{-3,+3},{-4,+5}});
    g.set_fill_colour(0.75,0,0.75); g.draw(orbit.reach()); g.set_fill_colour(0.5,0,0.5); g.draw(orbit.final()); g.draw(orbit.initial()); g.write("dai-evolve");

//    Enclosure intermediate_enclosure=orbit.final()[0];
//    intermediate_enclosure.recondition();
//    intermediate_enclosure = Enclosure(cast_exact_box(intermediate_enclosure.bounding_box()),factory);
//    auto orbit2=evolver.orbit(intermediate_enclosure,evolution_time,UPPER_SEMANTICS);

    ReachabilityAnalyser analyser(system,evolver);
    analyser.configuration().set_lock_to_grid_time(0.5_dec);
    analyser.configuration().set_maximum_grid_depth(4);
    analyser.verbosity=verbosity;
    cout << "\nComputing upper reach set...\n";
    auto chain_reach = analyser.upper_reach(initial_constraint_set,2.0_dec);
    cout << "\nplotting...\n";
//    auto chain_reach = analyser.outer_chain_reach(initial_constraint_set);

    Axis taxis(0,t,tmax), xaxis(-2,x,+4), yaxis(-4,y,+3);

    g.clear(); g.set_fill_colour(0.5,0,0.5); g.draw(chain_reach); g.write("dai-chain_reach");
//    plot("fitzhugh_nagumo",{xaxis,yaxis},Colour(0,1,0),HybridBoundedConstraintSet(initial_set),Colour(.5,.0,.5),orbit.reach(),Colour(.25,.0,.25),orbit.final());
//    plot("fitzhugh_nagumo-x",{taxis,xaxis},Colour(.5,.0,.5),orbit.reach());
//    plot("fitzhugh_nagumo-y",{taxis,yaxis},Colour(.5,.0,.5),orbit.reach());

}

} // namespace Ariadne

int main(int argc, const char* argv[]) {
    Ariadne::verbosity=Ariadne::get_verbosity(argc,argv);
    Ariadne::test();
}
