/***************************************************************************
 *            noisy-utilities.hpp
 *
 *  Copyright  2008-18 Luca Geretti
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

#include "ariadne.hpp"

#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>

namespace Ariadne {

typedef Tuple<String,DottedRealAssignments,RealVariablesBox,RealVariablesBox,Real,double> SystemType;

template<class F, class S> List<ResultOf<F(S)>> map(F const& f, List<S> const& list) {
    List<ResultOf<F(S)>> result; for(auto item : list) { result.append(f(item)); } return result;
}

FloatDP score(ValidatedConstrainedImageSet const& evolve_set) {
    auto bbx = evolve_set.bounding_box();
    return 1.0/std::pow(volume(bbx).get_d(),1.0/bbx.size());
}

ThresholdSweeperDP make_threshold_sweeper(double thr) { return ThresholdSweeperDP(DoublePrecision(),thr); }

template<class C> struct Reverse {
    C const& _c;
    Reverse(C const& c) :  _c(c) {}
    typename C::const_reverse_iterator begin() const{ return _c.rbegin(); }
    typename C::const_reverse_iterator end() const { return _c.rend(); }
};
template<class C> Reverse<C> reverse(C const& c) { return Reverse<C>(c); }


void run_single(String name, DifferentialInclusionIVP const& ivp, Real evolution_time, double step, List<InputApproximation> approximations, SweeperDP sweeper, SizeType freq, int verbosity) {

    typedef typename ValidatedVectorFunctionModelType::NumericType NumericType;
    typedef typename NumericType::PrecisionType PrecisionType;

    PrecisionType prec;

    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step,number_of_steps_between_simplifications=freq,number_of_variables_to_keep=20000);
    integrator.verbosity = verbosity;

    tms start_time, end_time;
    times(&start_time);

    List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(ivp,evolution_time);

    times(&end_time);
    clock_t ticks = end_time.tms_utime - start_time.tms_utime;
    clock_t const hz = sysconf(_SC_CLK_TCK);

    List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return ValidatedConstrainedImageSet(fm.domain(),fm);},flow_functions);
    auto final_set = flow_functions.back();
    ValidatedVectorFunctionModelType evolve_function = partial_evaluate(final_set,final_set.argument_size()-1,NumericType(evolution_time,prec));
    auto evolve_set = ValidatedConstrainedImageSet(evolve_function.domain(),evolve_function);

    std::cout << "score: " << score(evolve_set) << ", time: " << ticks / hz << "." << ticks % hz << " s" << std::endl;
/*
    std::cout << "plotting..." << std::endl;
    Box<FloatDPUpperInterval> graphics_box(f.result_size());
    for (auto set: reach_sets) {
        graphics_box = hull(graphics_box,set.bounding_box());
    }
    for (SizeType i : range(0,f.result_size()-1)) {
        for (SizeType j : range(i+1,f.result_size())) {
            Figure fig=Figure();
            fig.set_bounding_box(graphics_box);
            fig.set_projection(f.result_size(),i,j);
            fig.set_line_colour(0.0,0.0,0.0);
            fig.set_line_style(false);
            fig.set_fill_colour(0.5,0.5,0.5);
            fig.draw(starting_set);
            fig.set_fill_colour(1.0,0.75,0.5);
            for (auto set : reverse(reach_sets)) { fig.draw(set); }
            fig.draw(evolve_set);
            char num_char[7] = "";
            if (f.result_size() > 2)
                sprintf(num_char,"[%lu,%lu]",i,j);
            fig.write((name+num_char).c_str());
        }
    }
*/
}

void run_each_approximation(String name, DifferentialInclusionIVP const& ivp, Real evolution_time, double step, List<InputApproximation> approximations, SweeperDP sweeper, SizeType freq, int verbosity) {

    for (auto appro: approximations) {
        List<InputApproximation> singleapproximation = {appro};
        std::cout << appro << std::endl;
        run_single(name,ivp,evolution_time,step,singleapproximation,sweeper,freq,verbosity);
    }
}

void run_noisy_system(String name, const DottedRealAssignments& dynamics, const RealVariablesBox& inputs,
              const RealVariablesBox& initial, Real evolution_time, double step) {

    DifferentialInclusionIVP ivp(dynamics,inputs,initial);

    SizeType freq=12;
    ThresholdSweeperDP sweeper = make_threshold_sweeper(1e-8);

    int verbosity = 1;

    List<InputApproximation> approximations;
    approximations.append(InputApproximation::ZERO);
    approximations.append(InputApproximation::CONSTANT);
    approximations.append(InputApproximation::AFFINE);
    approximations.append(InputApproximation::SINUSOIDAL);
    approximations.append(InputApproximation::PIECEWISE);

    run_single(name,ivp,evolution_time,step,approximations,sweeper,freq,verbosity);
    //run_each_approximation(name,ivp,evolution_time,step,approximations,sweeper,freq,verbosity);
}

void run_noisy_system(SystemType system) {
    run_noisy_system(std::get<0>(system),std::get<1>(system),std::get<2>(system),std::get<3>(system),std::get<4>(system),std::get<5>(system));
}

}
