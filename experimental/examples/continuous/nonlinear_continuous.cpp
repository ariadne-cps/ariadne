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

#include "output/logging.hpp"
#include "nonlinear_continuous.hpp"
#include "dynamics/reachability_analyser.hpp"

namespace Ariadne {

inline int get_verbosity(int argc, const char* argv[]) {
    if(argc>1) {
        if(std::strcmp(argv[1],"-v")==0) { if(argc>2) { return std::atoi(argv[2]); } }
        else { std::cerr << "Unrecognised command-line option \"" << argv[1] << "\"\n"; } }
    return 0;
}

Interval<FloatDPValue> over_approximation(Interval<Real> const& ivl) {
    return cast_exact(Interval<FloatDPUpperBound>(ivl,dp));
}
Box<Interval<FloatDPValue>> over_approximation(Box<Interval<Real>> const& bx) {
    Box<Interval<FloatDPValue>> ebx(bx.dimension(),Interval<FloatDPValue>(0,0,dp));
    for(SizeType i=0; i!=ebx.dimension(); ++i) { ebx[i]=over_approximation(bx[i]); }
    return ebx;
}


using std::cout;

using Axis = ApproximateDoubleVariableInterval;

unsigned int verbosity;

void verify_dai() {
    auto problem = make_dai_problem();

    auto initial_set = problem.initial_set;
    auto system = problem.system;
    auto safe_set = problem.safe_set;

    cout << "initial_set=" << initial_set << "\n";
    cout << "system=" << system << "\n";
    cout << "safe_set=" << safe_set << "\n";

    RealExpressionBoundedConstraintSet unsafe_set = {{-0.5_dec<=x<=+0.5_dec,0.5_dec<=y<=1.5_dec},{sqr(x)+sqr(y-1)<=0.09_dec}};
    cout << "system.state_space()=" << system.state_space() << "\n";
    BoundedConstraintSet initial_constraint_set = initial_set.euclidean_set(system.state_space());
    BoundedConstraintSet safe_constraint_set = safe_set.euclidean_set(system.state_space());
    BoundedConstraintSet unsafe_constraint_set = unsafe_set.euclidean_set(system.state_space());

    cout << "initial_constraint_set=" << initial_constraint_set << "\n";
    cout << "safe_constraint_set=" << safe_constraint_set << "\n";

    TimeVariable t;
    Real gtmax=5.00_dec;
    Real tmax=1.00_dec;
    Real rtmax=1.00_dec;


    MaximumError max_err=0.01;
    TaylorSeriesIntegrator integrator(max_err);
    integrator.set_maximum_step_size(8);
    integrator.set_maximum_step_size(1);
    auto& factory=integrator.function_factory();

    VectorFieldEvolver evolver(system,integrator);
    evolver.verbosity=1;

    ContinuousReachabilityAnalyser analyser(system,evolver);
    analyser.configuration().set_transient_time(0.75_dec);
    analyser.configuration().set_lock_to_grid_time(0.75_dec);
    analyser.configuration().set_maximum_grid_depth(5);
    analyser.verbosity=1;


    Grid grid(system.dimension());
    auto paved_safe_set = inner_approximation(safe_constraint_set, grid, 5u);
    std::cerr << "paved_safe_set=" << paved_safe_set << "\n";
    auto paved_unsafe_set = outer_approximation(unsafe_constraint_set, grid, 5u);
    std::cerr << "paved_unsafe_set=" << paved_unsafe_set << "\n";

//    PavingType paved_safe_set=inner_approximation(safe_set, grid, safe_set_bounding_box, this->_configuration->maximum_grid_depth()+3);

//    safe_set.write(std::cerr);

    Figure g; g.set_bounding_box({{-5,+5},{-4,+6}});
    g.set_fill_colour(1.0,0,1.0); g.draw(paved_unsafe_set); g.write("dai-paved_unsafe_set"); g.clear();
    g.set_fill_colour(1.0,0,0.0); g.draw(unsafe_constraint_set); g.set_fill_colour(0.0,0,1.0); g.draw(initial_constraint_set); g.write("dai-unsafe_set"); g.clear();

    Axis taxis(0,t,gtmax), xaxis(-5,x,+5), yaxis(-4,y,+6);
    g.set_fill_colour(0.75,0,0.75);

    cout << "\nComputing safety...\n";
    auto safety = analyser.verify_safety(initial_constraint_set,safe_constraint_set);

    cout << "\nis_safe="<<safety.is_safe<<"\n";

    cout << "\nplotting...\n";
    g.clear(); g.set_fill_colour(1.0,0.75,0.75); g.draw(safety.safe_set); g.set_fill_colour(0.5,0,0.5); g.draw(safety.chain_reach_set); g.set_fill_colour(0.5,0,0.0); g.draw(unsafe_constraint_set); g.write("dai-chain_reach");

    cout << "\nis_safe="<<safety.is_safe<<"\n";
}

void test() {
//    auto problem = make_fitzhugh_nagumo_problem();
    auto problem = make_dai_problem();

    auto initial_set = problem.initial_set;
    auto system = problem.system;
    auto safe_set = problem.safe_set;

    cout << "initial_set=" << initial_set << "\n";
    cout << "system=" << system << "\n";
    cout << "safe_set=" << safe_set << "\n";

    cout << "system.state_space()=" << system.state_space() << "\n";
    BoundedConstraintSet initial_constraint_set = initial_set.euclidean_set(system.state_space());
    BoundedConstraintSet safe_constraint_set = safe_set.euclidean_set(system.state_space());

    cout << "initial_constraint_set=" << initial_constraint_set << "\n";
    cout << "safe_constraint_set=" << safe_constraint_set << "\n";

    TimeVariable t;
    Real gtmax=5.00_dec;
    Real tmax=1.00_dec;
    Real rtmax=1.00_dec;


    MaximumError max_err=0.01;
    TaylorSeriesIntegrator integrator(max_err);
    integrator.set_maximum_step_size(8);
    integrator.set_maximum_step_size(1);
    auto& factory=integrator.function_factory();
    VectorFieldEvolver evolver(system,integrator);
    evolver.verbosity=1;

    Figure g; g.set_bounding_box({{-3,+3},{-4,+5}});
    Axis taxis(0,t,gtmax), xaxis(-2,x,+4), yaxis(-4,y,+3);
    g.set_fill_colour(0.75,0,0.75);

    if (false) {

    //    Box<RealInterval> initial_box_set({{-0.625_dec,-0.5_dec},{-2.0_dec,-1.875_dec}});
    //    Box<RealInterval> initial_box_set({{1.000_dec,1.001_dec},{-2.000_dec,-1.999_dec}});
    //    Box<RealInterval> initial_box_set({{-0.500_dec,-0.499_dec},{-2.000_dec,-1.999_dec}});
        Box<RealInterval> initial_box_set({{-0.10_dec,+0.10_dec},{-2.001_dec,-1.999_dec}});
    //    Box<RealInterval> initial_box_set({{-0.10_dec,+0.10_dec},{-2.0_dec,-1.8_dec}});



        DoublePrecision pr;

        Box<ExactIntervalType> initial_box=over_approximation(initial_box_set);
        auto flow_bounds=integrator.flow_bounds(system.function(),initial_box,1);
        std::cerr<<"initial_box="<<initial_box<<"\n";
        std::cerr<<"flow_bounds="<<flow_bounds<<"\n";

        auto short_orbit=evolver.orbit(evolver.enclosure(initial_box_set),1.0_dec,Semantics::UPPER);
        g.set_fill_colour(1.0,0.5,1.0); g.draw(short_orbit.reach());
        g.set_fill_colour(0.75,0,0.75); g.draw(short_orbit.initial()); g.draw(short_orbit.final());

    /*
        auto short_orbit0=evolver.orbit(evolver.enclosure(initial_box_set),0.5_dec,Semantics::UPPER);
        Pair<Enclosure,Enclosure> intermediate_enclosures=short_orbit0.final()[0].split();
        auto short_orbit1=evolver.orbit(intermediate_enclosures.first,0.5_dec,Semantics::UPPER);
        auto short_orbit2=evolver.orbit(intermediate_enclosures.second,0.5_dec,Semantics::UPPER);
        g.set_fill_colour(1.0,0.5,1.0); g.draw(short_orbit0.reach());
        g.set_fill_colour(0.5,1.0,1.0); g.draw(short_orbit1.reach());
        g.set_fill_colour(1.0,1.0,0.5); g.draw(short_orbit2.reach());
    */
        g.write("dai-short_orbit"); g.clear();
    //    std::cerr<<"\nshort_orbit_final="<<short_orbit.final().sze()<<"\n";
    //    std::cerr<<"short_orbit_reach="<<short_orbit.reach().size()<<"\n";
        return;
        /*
        for(Dyadic rad(1,8u); rad<=4; rad=2*rad) {
            Box<Interval<FloatDPValue>> init(Box<Interval<Dyadic>>({{-1-rad,-1+rad},{-2-rad,-2+rad}}),pr);
            FloatDP step(8.0);
            auto phi=integrator.flow_step(system.function(),init,step);
            FloatError<DoublePrecision>::set_output_places(9);
            std::cerr<<"rad="<<rad.get_d()<<", step="<<step<<"\nphi.errors()="<<phi.errors()<<"\n\n";
            //std::cerr<<"phi="<<phi<<"\n";
        }
    */
        cout << "initial_box_set=" << initial_box_set << "\n";
        auto orbit=evolver.orbit(evolver.enclosure(initial_box_set),tmax,Semantics::UPPER);

        for (auto enclosure : orbit.reach()) { g.draw(enclosure.affine_over_approximation()); } g.write("dai-affine_orbit");

        g.clear();
        for (auto enclosure : orbit.reach()) { g.draw(enclosure.affine_over_approximation()); } g.write("dai-affine_orbit");

        GridTreePaving gts(2);
        g.clear();
        for (auto enclosure : orbit.reach()) { enclosure.adjoin_outer_approximation_to(gts,3); } g.draw(gts); g.write("dai-grid_orbit");

    //    Enclosure intermediate_enclosure=orbit.final()[0];
    //    intermediate_enclosure.recondition();
    //    intermediate_enclosure = Enclosure(cast_exact_box(intermediate_enclosure.bounding_box()),factory);
    //    auto orbit2=evolver.orbit(intermediate_enclosure,evolution_time,Semantics::UPPER);
    }
}

} // namespace Ariadne

int main(int argc, const char* argv[]) {
    Ariadne::verbosity=Ariadne::get_verbosity(argc,argv);
    Ariadne::verify_dai();
}
