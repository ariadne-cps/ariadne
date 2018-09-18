/***************************************************************************
 *            noisy-utilities.hpp
 *
 *  Copyright  2008-18 Luca Geretti
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

namespace Ariadne {

typedef Tuple<String,DottedRealAssignments,RealVariablesBox,RealVariablesBox,Real,double> SystemType;

void run_single(String name, DifferentialInclusionIVP const& ivp, Real evolution_time, double step, List<InputApproximation> approximations, SweeperDP sweeper, SizeType freq, unsigned int verbosity, bool draw);
void run_each_approximation(String name, DifferentialInclusionIVP const& ivp, Real evolution_time, double step, List<InputApproximation> approximations, SweeperDP sweeper, SizeType freq, unsigned int verbosity, bool draw);
void run_noisy_system(String name, const DottedRealAssignments& dynamics, const RealVariablesBox& inputs, const RealVariablesBox& initial, Real evolution_time, double step);
void run_noisy_system(SystemType system);


template<class F, class S> List<ResultOf<F(S)>> map(F const& f, List<S> const& list) {
    List<ResultOf<F(S)>> result; for(auto item : list) { result.append(f(item)); } return result;
}

inline FloatDP score(ValidatedConstrainedImageSet const& evolve_set) {
    auto bbx = evolve_set.bounding_box();
    return 1.0/std::pow(volume(bbx).get_d(),1.0/bbx.size());
}

template<class C> struct Reverse {
    C const& _c;
    Reverse(C const& c) :  _c(c) {}
    typename C::const_reverse_iterator begin() const{ return _c.rbegin(); }
    typename C::const_reverse_iterator end() const { return _c.rend(); }
};
template<class C> Reverse<C> reverse(C const& c) { return Reverse<C>(c); }


void run_single(String name, DifferentialInclusionIVP const& ivp, Real evolution_time, double step, List<InputApproximation> approximations, SweeperDP sweeper, SizeType freq, unsigned int verbosity, bool draw) {

    typedef typename ValidatedVectorMultivariateFunctionModelType::NumericType NumericType;
    typedef typename NumericType::PrecisionType PrecisionType;

    PrecisionType prec;

    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step,number_of_steps_between_simplifications=freq,number_of_variables_to_keep=20000);
    integrator.verbosity = verbosity;

    StopWatch sw;

    List<ValidatedVectorMultivariateFunctionModelType> flow_functions = integrator.flow(ivp,evolution_time);

    sw.click();

    List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorMultivariateFunctionModelType const& fm){return ValidatedConstrainedImageSet(fm.domain(),fm);},flow_functions);
    auto final_set = flow_functions.back();
    ValidatedVectorMultivariateFunctionModelType evolve_function = 
        partial_evaluate(final_set,final_set.argument_size()-1,NumericType(evolution_time,prec));
    auto evolve_set = ValidatedConstrainedImageSet(evolve_function.domain(),evolve_function);

    std::cout << "score: " << score(evolve_set) << ", time: " << sw.elapsed() << " s" << std::endl;

    if (draw) {
        std::cout << "plotting..." << std::endl;
        auto n = ivp.di().num_variables();
        Box<FloatDPUpperInterval> graphics_box(n);
        for (auto set: reach_sets) {
            graphics_box = hull(graphics_box,set.bounding_box());
        }
        for (SizeType i : range(0,n-1)) {
            for (SizeType j : range(i+1,n)) {
                Figure fig=Figure();
                fig.set_bounding_box(graphics_box);
                fig.set_projection(n,i,j);
                fig.set_line_colour(0.0,0.0,0.0);
                fig.set_line_style(false);
                fig.set_fill_colour(0.5,0.5,0.5);
                fig.draw(ivp.X0());
                fig.set_fill_colour(1.0,0.75,0.5);
                for (auto set : reverse(reach_sets)) { fig.draw(set); }
                fig.draw(evolve_set);
                char num_char[7] = "";
                if (n > 2) sprintf(num_char,"[%lu,%lu]",i,j);
                fig.write((name+num_char).c_str());
            }
        }
    }
}

void run_each_approximation(String name, DifferentialInclusionIVP const& ivp, Real evolution_time, double step, List<InputApproximation> approximations, SweeperDP sweeper, SizeType freq, unsigned int verbosity, bool draw) {

    for (auto appro: approximations) {
        List<InputApproximation> singleapproximation = {appro};
        std::cout << appro << std::endl;
        run_single(name,ivp,evolution_time,step,singleapproximation,sweeper,freq,verbosity,draw);
    }
}

void run_noisy_system(String name, const DottedRealAssignments& dynamics, const RealVariablesBox& inputs,
              const RealVariablesBox& initial, Real evolution_time, double step) {

    DifferentialInclusionIVP ivp(dynamics,inputs,initial);

    SizeType freq=12;
    ThresholdSweeperDP sweeper(DoublePrecision(),1e-8);

    unsigned int verbosity = 1;
    bool draw = false;
    DRAWING_METHOD = DrawingMethod::BOX;

    List<InputApproximation> approximations;
    approximations.append(InputApproximation::ZERO);
    approximations.append(InputApproximation::CONSTANT);
    approximations.append(InputApproximation::AFFINE);
    approximations.append(InputApproximation::SINUSOIDAL);
    approximations.append(InputApproximation::PIECEWISE);

    run_single(name,ivp,evolution_time,step,approximations,sweeper,freq,verbosity,draw);
    //run_each_approximation(name,ivp,evolution_time,step,approximations,sweeper,freq,verbosity,draw);
}

void run_noisy_system(SystemType system) {
    run_noisy_system(std::get<0>(system),std::get<1>(system),std::get<2>(system),std::get<3>(system),std::get<4>(system),std::get<5>(system));
}

}
