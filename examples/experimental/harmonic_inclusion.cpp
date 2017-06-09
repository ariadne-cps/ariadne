/***************************************************************************
 *            test_differential_inclusions.cpp
 *
 *  Copyright  2008-17  Pieter Collins, Sanja Zivanovic
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

#include "dynamics/differential_inclusion.hpp"

#include "algebra/sweeper.hpp"
#include "solvers/integrator_interface.hpp"
#include "geometry/box.hpp"
#include "ariadne.hpp"

//#include "geometry/zonotope.hpp"
//#include "geometry/polytope.hpp"

namespace Ariadne {

Nat verbosity=0;

template<class F, class S> List<ResultOf<F(S)>> map(F const& f, List<S> const& list) {
    List<ResultOf<F(S)>> result; for(auto item : list) { result.append(f(item)); } return result;
}

ValidatedConstrainedImageSet range(ValidatedVectorFunctionModelType const& fm) {
    return ValidatedConstrainedImageSet(fm.domain(),fm);
}

void damped_harmonic(InclusionIntegratorBase const& integrator, Real evolution_time, Real damping, Vector<Real> noise, Real box_radius) {
    std::cerr<<"evolution_time="<<evolution_time
             <<", damping="<<damping
             <<", noise="<<noise
             <<", box_radius="<<box_radius
             <<"\n";
    SizeType n=2;

    Real T=evolution_time;
    Real d=damping;
    Vector<Real> v=noise;
    Real e=box_radius;

    auto I=IntervalDomainType(-1,+1);
    auto c=EffectiveVectorFunction::constant(2,1);
    auto x=EffectiveVectorFunction::identity(2);
    //auto f=EffectiveVectorFunction({-x[0]*d-x[1],x[0]-x[1]*d});
    auto f=EffectiveVectorFunction({x[0]/2-c,c});
    RealBox Vr({RealInterval(-v[0],+v[0]),RealInterval(-v[1],+v[1])});
    Box<Interval<ValidatedNumber>> Vb(Vr);
    BoxDomainType V=over_approximation(Vr);
    auto starting_set=RealBox({{1-e,1+e},{-e,+e}});
    BoxDomainType X0=over_approximation(starting_set);

    std::cerr<<"Vr="<<Vr<<", Vb="<<Vb<<"\n";
    std::cerr<<"f="<<f<<", V="<<V<<", X0="<<X0<<", T="<<T<<"\n";

    std::cerr<<"Computing evolution... ";
//    auto start_time=time.clock();
    List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,V,X0,T);
//    auto end_time=time.clock() - start_time;
//    ARIADNE_LOG(2,"end_time="<<end_time);
    std::cerr<<"done! Computed "<<flow_functions.size()<<" patches\n";

    List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);

    ValidatedVectorFunctionModelType evolve_function = partial_evaluate(flow_functions.back(),X0.size(),Float64Bounds(T,Precision64()));
    ValidatedConstrainedImageSet evolve_set = range(evolve_function);

    std::cerr<<"Plotting... ";
    auto fig=Figure();
    fig.set_bounding_box(BoxDomainType({{-2,1.5},{-1.3,2}}));
    fig.set_line_colour(0.0,0.0,0.0);
    fig.set_fill_colour(1.0,0.0,0.8);
    for (auto set : reach_sets) { fig.draw(set); }
    fig.set_fill_colour(0.75,0.0,0.6);
    fig.draw(RealBoxSet(starting_set));
    fig.draw(evolve_set);
    fig.write("damped_harmonic-test");
    std::cerr<<"done!\n";
}

void test() {
    // damped_harmonic( evolution_time=2*pi, damping=1.0/4, noise=(0.0,0.1), step_size=8.0/32 );
    // damped_harmonic( evolution_time=2*pi, damping=0.0, noise=(0.0,0.1), delta=0.01, step_size=2*pi/50 );
    ThresholdSweeper64 sweeper(Precision64(),1e-8);
    auto integrator = InclusionIntegrator2ndOrder(sweeper, step_size=1.0/4, number_of_steps_between_simplifications=64, number_of_variables_to_keep=32);

    Real evolution_time=3/4_q;
    Real damping=1/100_q;
    Vector<Real> noise={1/16_q,1/32_q};
    Real box_radius=1/128_q;
    damped_harmonic( integrator, evolution_time, damping, noise, box_radius);
    // damped_harmonic( evolution_time=6.5, damping=1.0/16, noise=1.0/8, step_size=1.0/4 );
}

} // namespace Ariadne;

int main() {
    Ariadne::test();
}
