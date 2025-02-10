/***************************************************************************
 *            SPRE20.hpp
 *
 *  Copyright  2020  Luca Geretti
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

#include "ariadne.hpp"
#include "utility/stopwatch.hpp"
#include "arch.hpp"

namespace Ariadne {

template<class IS, class SYS, class SS> struct Problem {
    IS initial_set; SYS system; SS safe_set; };
template<class IS, class SYS, class SS> Problem<IS,SYS,SS> make_problem(IS initial_set, SYS system, SS safe_set) {
    return Problem<IS,SYS,SS>{initial_set,system,safe_set}; }

decltype(auto) make_space_rendezvous_problem() {

    // Differential variables
    RealVariable t("t");
    RealVariable x("x");
    RealVariable y("y");
    RealVariable vx("vx");
    RealVariable vy("vy");
    // Algebraic variables
    RealVariable rc("rc");
    RealVariable ux("ux");
    RealVariable uy("uy");
    RealVariable av("av");

    RealConstant mu("mu",3.986e14_dec*3600_dec);
    RealConstant r("r",42164000);
    RealConstant mc("mc",500);
    RealConstant n("n",sqrt(mu/pow(r,3)));

    Matrix<Real> K1({{-28.8287_dec, 0.1005_dec, -1449.9754_dec, 0.0046_dec},{-0.087_dec, -33.2562_dec, 0.00462_dec, -1451.5013_dec}});
    Matrix<Real> K2({{-288.0288_dec, 0.1312_dec, -9614.9898_dec, 0},{-0.1312_dec, -288, 0, -9614.9883_dec}});

    StringVariable spacecraft("spacecraft");
    StringConstant approaching("approaching");
    StringConstant rendezvous("rendezvous");
    StringConstant aborting("aborting");

    DiscreteEvent attempt("attempt");
    DiscreteEvent can_timeout("can_timeout");
    DiscreteEvent must_timeout("must_timeout");

    HybridAutomaton space_rendezvous("space_rendezvous");

//    HybridBoundedConstraintSet initial_set(spacecraft|approaching,{t==0,-925<=x<=-875,-425<=y<=-375,vx==0,vy==0});
    HybridBoundedConstraintSet initial_set(spacecraft|approaching,{t==0,-906.25_dy<=x<=-893.75_dy,-406.25_dy<=y<=-393.75_dy,vx==0,vy==0});

    HybridConstraintSet safe_set;
    safe_set.adjoin(spacecraft|rendezvous,{y>=x*tan(Ariadne::pi/6u),-y>=x*tan(Ariadne::pi/6u),sqrt(sqr(vx)+sqr(vy))<=3.3_dec});
    safe_set.adjoin(spacecraft|aborting,{x<=-0.2_dec,x>=0.2_dec,y<=-0.2_dec,y>=0.2_dec});

    auto rc_dyn = sqrt(sqr(r+x)+sqr(y));
    auto t_dyn = 1;
    auto x_dyn = vx;
    auto y_dyn = vy;
    auto vx_dyn = sqr(n)*x+2*n*vy+mu/sqr(r)-mu/pow(r,3)*(r+x)+ux/mc;
    auto vy_dyn = sqr(n)*y-2*n*vx-mu/pow(rc,3)*y+uy/mc;

    space_rendezvous.new_mode(spacecraft|approaching, {let(rc)=rc_dyn,
                                                       let(ux)=K1[0][0]*x+K1[0][1]*y+K1[0][2]*vx+K1[0][3]*vy,
                                                       let(uy)=K1[1][0]*x+K1[1][1]*y+K1[1][2]*vx+K1[1][3]*vy},
                                                      {dot(t)=t_dyn, dot(x)=vx, dot(y)=vy, dot(vx)=vx_dyn, dot(vy)=vy_dyn});
    space_rendezvous.new_mode(spacecraft|rendezvous, {let(rc)=rc_dyn,
                                                      let(ux)=K2[0][0]*x+K2[0][1]*y+K2[0][2]*vx+K2[0][3]*vy,
                                                      let(uy)=K2[1][0]*x+K2[1][1]*y+K2[1][2]*vx+K2[1][3]*vy},
                                                     {dot(t)=t_dyn, dot(x)=vx, dot(y)=vy, dot(vx)=vx_dyn, dot(vy)=vy_dyn});
    space_rendezvous.new_mode(spacecraft|aborting, {let(rc)=rc_dyn, let(ux)=0, let(uy)=0}, {dot(t)=t_dyn, dot(x)=vx, dot(y)=vy, dot(vx)=vx_dyn, dot(vy)=vy_dyn});

    space_rendezvous.new_transition(spacecraft|approaching,attempt,spacecraft|rendezvous, next({t,x,y,vx,vy})={t,x,y,vx,vy}, x>=-100, EventKind::URGENT);
    space_rendezvous.new_transition(spacecraft|approaching,can_timeout,spacecraft|aborting, next({t,x,y,vx,vy})={t,x,y,vx,vy}, t>=120, EventKind::PERMISSIVE);
    space_rendezvous.new_transition(spacecraft|rendezvous,can_timeout,spacecraft|aborting, next({t,x,y,vx,vy})={t,x,y,vx,vy}, t>=120, EventKind::PERMISSIVE);

    space_rendezvous.new_invariant(spacecraft|approaching, t<=150, must_timeout);
    space_rendezvous.new_invariant(spacecraft|rendezvous, t<=150, must_timeout);



    std::cerr<<space_rendezvous<<"\n";
    std::cerr<<"alg_eqns:"<<space_rendezvous.auxiliary_function({spacecraft|approaching})<<"\n";
    auto dif_eqns=space_rendezvous.dynamic_function({spacecraft|approaching});
    std::cerr<<"dif_eqns: R"<<dif_eqns.argument_size()<<"->R"<<dif_eqns.result_size()<<":"<<dif_eqns<<"\n";
    std::cerr<<"J="<<dif_eqns.jacobian(Vector<FloatDPApproximation>({10,100,1000,100000,1000000},dp))<<"\n";;
    //J=[0.,0.,0.,0.,0.; 0.,0.,0.,1.000,0.; 0.,0.,0.,0.,1.000; 0.,-0.05766,0.0002010,-2.900,0.008760; 0.,-0.0001740,-0.06651,-0.008741,-2.903]
    auto J=Matrix<FloatDPApproximation>({{0.,0.,0.,0.,0.},{ 0.,0.,0.,1.000,0.},{ 0.,0.,0.,0.,1.000},{ 0.,-0.05766,0.0002010,-2.900,0.008760},{ 0.,-0.0001740,-0.06651,-0.008741,-2.903}},dp);

    auto loc=DiscreteLocation({spacecraft|approaching});
    auto spc=space_rendezvous.continuous_state_space(loc);
    //auto init=initial_set.euclidean_set(loc,spc);
    std::cerr<<"init="<<initial_set.continuous_set(loc)<<"\n";
    std::cerr<<"spc="<<spc<<"\n";
    auto init=cast_exact_box(initial_set.euclidean_set(loc,spc).bounding_box());
    init[0]=IntervalDomainType(0,0); init[3]=IntervalDomainType(0,0); init[4]=IntervalDomainType(0,0);
    std::cerr<<"init="<<init<<"\n";
    std::cerr<<"image(init,dif_eqns)="<<image(init,dif_eqns)<<"\n";
    StepSizeType h=1.0_x;

    EulerBounder bnd;
    std::cerr<<"bounds="<<bnd.compute(dif_eqns,init,suggest(h))<<"\n";


    {char c; std::cin>>c;}
    return make_problem(initial_set,space_rendezvous,safe_set);


}

void SPRE20()
{
    ArchBenchmark benchmark("SPRE20");
    auto instance = benchmark.create_instance();
    instance.write();

    auto problem = make_space_rendezvous_problem();

    auto initial_set = problem.initial_set;
    auto system = problem.system;
    HybridConstraintSet safe_set = problem.safe_set;

    CONCLOG_PRINTLN("Space Rendezvous system:");

    DiscreteLocation initial_location = initial_set.location();
    RealSpace initial_space = system.state_space()[initial_location];
    BoundedConstraintSet initial_constraint_set = initial_set.euclidean_set(initial_location,initial_space);

    StepMaximumError max_err=8e-6;
    max_err=1e-4;
    TaylorSeriesIntegrator integrator(max_err,Order(3u));

    GeneralHybridEvolver evolver(system);
    evolver.set_integrator(integrator);
    evolver.configuration().set_maximum_enclosure_radius(12.0);
    evolver.configuration().set_maximum_step_size(1.0);
    evolver.configuration().set_maximum_spacial_error(1e-3);
    evolver.configuration().set_enable_subdivisions(true);

    Stopwatch<Milliseconds> sw;

    CONCLOG_PRINTLN("Computing orbit...");
    HybridTime evolution_time(115,2);
    auto orbit=evolver.orbit(initial_set,evolution_time,Semantics::UPPER);

    StringVariable spacecraft("spacecraft");
    StringConstant approaching("approaching");
    StringConstant rendezvous("rendezvous");
    StringConstant aborting("aborting");

    std::cerr<<"\n\n\n\n\n"<<orbit.final().size()<<orbit.final().front().state_function().errors();
/*
    CONCLOG_PRINTLN("Checking properties...");
    Nat num_ce = 0;
    for (auto reach : orbit.reach()) {
        auto reach_box = HybridExactBox(reach.location(),reach.state_space(),cast_exact_box(reach.bounding_box().euclidean_set()));
        if (reach.location() == DiscreteLocation(spacecraft|rendezvous) and not(definitely(safe_set.covers(reach_box)))) {
            CONCLOG_PRINTLN_AT(1,"Found counterexample in location " << reach.location() << " with bounding box " << reach_box << ", unsafe");
            ++num_ce;
        }
        if (reach.location() == DiscreteLocation(spacecraft|aborting) and not(definitely(safe_set.separated(reach_box)))) {
            CONCLOG_PRINTLN_AT(1,"Found counterexample in location " << reach.location() << " with bounding box " << reach_box << ", unsafe");
            ++num_ce;
        }
    }
    if (num_ce>0) CONCLOG_PRINTLN("Number of counterexamples: " << num_ce);

    sw.click();
    CONCLOG_PRINTLN("Done in " << sw.elapsed_seconds() << " seconds.");
*/
    RealVariable t("t"), x("x"), y("y"), vx("vx"), vy("vy");

    CONCLOG_PRINTLN("Plotting...");
CONCLOG_RUN_MUTED(
    plot("SPRE20",{-1000<=x<=200,-450<=y<=0},Colour(1.0,0.75,0.5),orbit.reach());
)
    CONCLOG_PRINTLN("File SPRE20.png written.");
}

} // namespace Ariadne
