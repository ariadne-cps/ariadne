/***************************************************************************
 *            test_differential_inclusion.cpp
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
#include "function/function.hpp"
#include "function/formula.hpp"
#include "function/taylor_model.hpp"
#include "algebra/algebra.hpp"
#include "geometry/function_set.hpp"
#include "output/graphics.hpp"

#include "test/test.hpp"
#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>

namespace Ariadne {


template<class F, class S> List<ResultOf<F(S)>> map(F const& f, List<S> const& list) {
    List<ResultOf<F(S)>> result; for(auto item : list) { result.append(f(item)); } return result;
}

ValidatedConstrainedImageSet range(ValidatedVectorFunctionModelType const& fm) {
    return ValidatedConstrainedImageSet(fm.domain(),fm);
}

FloatDP score(ValidatedConstrainedImageSet const& evolve_set) {
    auto bbx = evolve_set.bounding_box();
    return 1.0/pow(volume(bbx).get_d(),1.0/bbx.size());
}

ThresholdSweeperDP make_threshold_sweeper(double thr) { return ThresholdSweeperDP(DoublePrecision(),thr); }
GradedSweeperDP make_graded_sweeper(SizeType deg) { return GradedSweeperDP(DoublePrecision(),deg); }
GradedThresholdSweeperDP make_graded_threshold_sweeper(SizeType deg, double thr) { return GradedThresholdSweeperDP(DoublePrecision(),deg, thr); }

template<class C> struct Reverse {
    C const& _c;
    Reverse(C const& c) :  _c(c) {}
    typename C::const_reverse_iterator begin() const{ return _c.rbegin(); }
    typename C::const_reverse_iterator end() const { return _c.rend(); }
};
template<class C> Reverse<C> reverse(C const& c) { return Reverse<C>(c); }

} // namespace Ariadne

using namespace Ariadne;

class TestInclusionIntegrator {


    Void run_battery_frequency_topdown(String name,
                               ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, RealVector noise_levels,
                               RealBox real_starting_set, Real evolution_time, double step, SweeperDP sweeper, SizeType max_freq) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;

        PrecisionType prec;
        BoxDomainType noise=cast_exact_box(UpperIntervalType(-1,+1)*noise_levels);
        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        SizeType currfreq = max_freq;
        List<SizeType> frequencies;
        while (currfreq >= 1) {
            frequencies.push_back(currfreq);
            currfreq = floor(0.95*currfreq);
        }

        while (!frequencies.empty()) {

            auto freq = frequencies.back();
            frequencies.pop_back();

            SizeType base = 100000;
            List<SharedPointer<DIApproximation>> approximations;
            approximations.append(SharedPointer<DIApproximation>(new PiecewiseDIApproximation(sweeper)));
            approximations.append(SharedPointer<DIApproximation>(new SinusoidalDIApproximation(sweeper)));
            approximations.append(SharedPointer<DIApproximation>(new AffineDIApproximation(sweeper)));
            approximations.append(SharedPointer<DIApproximation>(new ConstantDIApproximation(sweeper)));
            approximations.append(SharedPointer<DIApproximation>(new ZeroDIApproximation(sweeper)));
            auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
            integrator.verbosity = 0;

            tms start_time, end_time;
            times(&start_time);

            List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);

            times(&end_time);
            clock_t ticks = end_time.tms_utime - start_time.tms_utime;
            clock_t const hz = sysconf(_SC_CLK_TCK);

            List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);
            ValidatedVectorFunctionModelType evolve_function = partial_evaluate(flow_functions.back(),flow_functions.back().argument_size()-1,NumericType(evolution_time,prec));
            ValidatedConstrainedImageSet evolve_set = range(evolve_function);

            std::cout << freq << ": " << score(evolve_set) << ", " << ticks / hz << "." << ticks % hz << "s" << std::endl;
        }
    }


    Void run_battery_each_approximation(String name,
                                        ValidatedVectorFunction const &f, Vector<ValidatedVectorFunction> const &g,
                                        RealVector noise_levels,
                                        RealBox real_starting_set, Real evolution_time, double step, SizeType freq) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
        PrecisionType prec;

        BoxDomainType noise=cast_exact_box(UpperIntervalType(-1,+1)*noise_levels);
        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        auto sweeper = make_threshold_sweeper(1e-8);

        for (auto approx: range(0,5)) {

            List<SharedPointer<DIApproximation>> approximations;

            SizeType ppi = 0;
            switch (approx) {
                case 0:
                    approximations.append(SharedPointer<DIApproximation>(new ZeroDIApproximation(sweeper)));
                    ppi = 0;
                    break;
                case 1:
                    approximations.append(SharedPointer<DIApproximation>(new ConstantDIApproximation(sweeper)));
                    ppi = 1;
                    break;
                case 2:
                    approximations.append(SharedPointer<DIApproximation>(new AffineDIApproximation(sweeper)));
                    ppi = 2;
                    break;
                case 3:
                    approximations.append(SharedPointer<DIApproximation>(new SinusoidalDIApproximation(sweeper)));
                    ppi = 2;
                    break;
                case 4:
                    approximations.append(SharedPointer<DIApproximation>(new PiecewiseDIApproximation(sweeper)));
                    ppi = 2;
                    break;
                default:
                    break;
            }

            auto n = f.result_size();
            auto m = noise.size();
            double rho = 6.0;
            SizeType base = 1000;
            auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
            integrator.verbosity = 0;

            std::cout << approximations.at(0)->getKind() << std::endl;

            tms start_time, end_time;
            times(&start_time);

            List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);

            times(&end_time);
            clock_t ticks = end_time.tms_utime - start_time.tms_utime;
            clock_t const hz = sysconf(_SC_CLK_TCK);

            List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);
            ValidatedVectorFunctionModelType evolve_function = partial_evaluate(flow_functions.back(),flow_functions.back().argument_size()-1,NumericType(evolution_time,prec));
            ValidatedConstrainedImageSet evolve_set = range(evolve_function);

            std::cout << "score " << score(evolve_set) << ", " << ticks / hz << "." << ticks % hz << "s" << std::endl;
        }
    }

    Void run_single_test(String name, InclusionIntegratorInterface& integrator,
                  ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, RealVector noise_levels,
                  RealBox real_starting_set, Real evolution_time) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
        PrecisionType prec;

        FloatDPUpperBound ratio(1.0);
        BoxDomainType noise=cast_exact_box(UpperIntervalType(-ratio,+ratio)*noise_levels);
        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        std::cout << "flowing..." << std::endl;

        tms start_time, end_time;
        times(&start_time);

        List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);

        times(&end_time);
        clock_t ticks = end_time.tms_utime - start_time.tms_utime;
        clock_t const hz = sysconf(_SC_CLK_TCK);

        List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);
        auto final_set = flow_functions.back();
        ValidatedVectorFunctionModelType evolve_function = partial_evaluate(final_set,final_set.argument_size()-1,NumericType(evolution_time,prec));
        ValidatedConstrainedImageSet evolve_set = range(evolve_function);

        std::cout << "score: " << score(evolve_set) << ", " << ticks / hz << "." << ticks % hz << "s" << std::endl;
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
                fig.write(("test_differential_inclusion-"+name+num_char).c_str());
            }
        }
*/
    }

    Void run_test(String name, ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, RealVector noise_levels,
                  RealBox real_starting_set, Real evolution_time, double step) const {

        SizeType freq=12;
        ThresholdSweeperDP sweeper = make_threshold_sweeper(1e-8);

        List<SharedPointer<DIApproximation>> approximations;
        approximations.append(SharedPointer<DIApproximation>(new PiecewiseDIApproximation(sweeper)));
        approximations.append(SharedPointer<DIApproximation>(new SinusoidalDIApproximation(sweeper)));
        approximations.append(SharedPointer<DIApproximation>(new AffineDIApproximation(sweeper)));
        approximations.append(SharedPointer<DIApproximation>(new ConstantDIApproximation(sweeper)));
        approximations.append(SharedPointer<DIApproximation>(new ZeroDIApproximation(sweeper)));

        auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=20000);
        integrator.verbosity = 0;

        run_single_test(name,integrator,f,g,noise_levels,real_starting_set,evolution_time);
        //this->run_battery_each_approximation(name,f,g,noise_levels,real_starting_set,evolution_time,step,freq);
    }

  public:
    void test() const;

    void test_wiggins_18_7_3() const;
    void test_order7() const;
    void test_3dsphere() const;
    void test_vinograd() const;
    void test_higgins_selkov() const;
    void test_reactor() const;
    void test_lotka_volterra() const;
    void test_fitzhugh_nagumo() const;
    void test_jet_engine() const;
    void test_pi_controller() const;
    void test_jerk21() const;
    void test_lorenz() const;
    void test_rossler() const;
    void test_jerk16() const;
    void test_DCDC() const;
    void test_harmonic() const;

    void test_van_der_pol() const;
    void test_clock() const;
};


void TestInclusionIntegrator::test_wiggins_18_7_3() const {
    double step=1.0/16;

    RealVector noise_levels={1/100_q,1/100_q};

    auto x = ValidatedVectorFunction::identity(2u);
    auto one = ValidatedScalarFunction::constant(2u,1_z);
    auto zero = ValidatedScalarFunction::constant(2u,0_z);

    auto f = ValidatedVectorFunction({-x[0]+2*x[1]+pow(x[0],2)*x[1]+pow(x[0],4)*pow(x[1],5),
                                      -x[1]-pow(x[0],4)*pow(x[1],6)+pow(x[0],8)*pow(x[1],9)});

    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e=1/10000000_q;
    RealBox starting_set={{Real(1/3_q)-e,Real(1/3_q)+e},{Real(1/3_q)-e,Real(1/3_q)+e}};
    Real evolution_time=90/10_q;

    this->run_test("wiggins_18_7_3",f,g,noise_levels,starting_set,evolution_time,step);
}


void TestInclusionIntegrator::test_order7() const {
    double step=1.0/32;

    RealVector noise_levels={1/100_q};

    auto x = ValidatedVectorFunction::identity(2u);
    auto one = ValidatedScalarFunction::constant(2u,1_z);
    auto zero = ValidatedScalarFunction::constant(2u,0_z);

    auto f = ValidatedVectorFunction({-42*pow(x[0],7) +68*pow(x[0],6)*x[1] -46*pow(x[0],5)*x[1] +256*pow(x[0],4)*x[1] +156*pow(x[0],3)*x[1] +50*pow(x[0],2)*x[1] +20*x[0]*pow(x[1],6) -8*pow(x[1],7),
                                      x[1]*(1110*pow(x[0],6) -220*pow(x[0],5)*x[1] -3182*pow(x[0],4)*x[1] +478*pow(x[0],3)*pow(x[1],3) + 487*pow(x[0],2)*pow(x[1],4) -102*x[0]*pow(x[1],5) -12*pow(x[1],6))});

    Vector<ValidatedVectorFunction> g({{zero,one}});

    Real e=1/10_q;
    RealBox starting_set={{Real(-1.0)-e,Real(-1.0)+e},{Real(1.0)-e,Real(1.0)+e}};
    Real evolution_time=40/8_q;

    this->run_test("order7",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_3dsphere() const {

    double step=1.0/16;

    RealVector noise_levels={4/1000_q,4/1000_q,4/1000_q};

    auto x = ValidatedVectorFunction::identity(3u);
    auto one = ValidatedScalarFunction::constant(3u,1_z);
    auto zero = ValidatedScalarFunction::constant(3u,0_z);

    auto f = ValidatedVectorFunction({pow(x[0],2) - x[0]*(pow(x[0],3)+pow(x[1],3)+pow(x[2],3)),
                                      pow(x[1],2) - x[1]*(pow(x[0],3)+pow(x[1],3)+pow(x[2],3)),
                                      pow(x[2],2) - x[2]*(pow(x[0],3)+pow(x[1],3)+pow(x[2],3))});

    Vector<ValidatedVectorFunction> g({{one,zero,zero},{zero,one,zero},{zero,zero,one}});

    Real e=1/1000000000_q;
    Real x0_i(1/4_q);
    Real x1_i(1/8_q);
    Real x2_i(1/10_q);
    RealBox starting_set={{x0_i-e,x0_i+e},{x1_i-e,x1_i+e},{x2_i-e,x2_i+e}};
    Real evolution_time=90/10_q;

    this->run_test("3dsphere",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_vinograd() const {
    double step=1.0/16;

    RealVector noise_levels={1/1000_q,1/1000_q};

    auto x = ValidatedVectorFunction::identity(2u);
    auto one = ValidatedScalarFunction::constant(2u,1_z);
    auto zero = ValidatedScalarFunction::constant(2u,0_z);
    auto f = ValidatedVectorFunction({pow(x[1],5)+pow(x[0],2)*(-x[0]+x[1]),
                                      pow(x[1],2)*(-2*x[0]+x[1])});

    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e=1/10000000_q;
    RealBox starting_set={{Real(1)-e,Real(1)+e},{-e,+e}};
    Real evolution_time=180/10_q;

    this->run_test("vinograd",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_higgins_selkov() const {

    double step=1.0/50;

    RealVector noise_levels={2/10000_q,2/10000_q,2/10000_q};

    auto x = ValidatedVectorFunction::identity(2u);
    auto one = ValidatedScalarFunction::constant(2u,1_z);
    auto zero = ValidatedScalarFunction::constant(2u,0_z);

    Real k2(1.00001);

    auto f = ValidatedVectorFunction({one-x[0]*x[1]*x[1],x[0]*x[1]*x[1] - k2*x[1]});
    Vector<ValidatedVectorFunction> g({{one,zero},{-x[0]*x[1]*x[1],x[0]*x[1]*x[1]},{zero,-x[1]}});

    Real e=1/100_q;
    Real x0_i(2.0);
    Real x1_i(1.0);
    RealBox starting_set={{x0_i-e,x0_i+e},{x1_i-e,x1_i+e}};
    Real evolution_time=100/10_q;

    this->run_test("higgins-selkov",f,g,noise_levels,starting_set,evolution_time,step);

}

void TestInclusionIntegrator::test_reactor() const {
    double step=1.0/16;

    RealVector noise_levels={1/1000_q,1/1000_q,2/10_q};

    auto x = ValidatedVectorFunction::identity(4u);
    auto one = ValidatedScalarFunction::constant(4u,1_z);
    auto zero = ValidatedScalarFunction::constant(4u,0_z);

    Real U3(30.0);
    Real k2(0.4);
    Real iV(0.05);
    Real ka(0.1);
    Real kb(0.045);

    auto f = ValidatedVectorFunction({-U3*x[0]*x[1]-k2*x[0]*x[2]+iV-ka*x[0],-U3*x[0]*x[1]+kb-ka*x[1],
                                      U3*x[0]*x[1]-k2*x[0]*x[2]-ka*x[2],k2*x[0]*x[2]-ka*x[3]});

    Vector<ValidatedVectorFunction> g({{one*iV,zero,zero,zero},{zero,one*iV,zero,zero},{x[0]*x[1],x[0]*x[1],x[0]*x[1],zero}});

    Real e=1/1000000_q;
    RealBox starting_set={{0,e},{0,e},{0,e},{0,e}};
    Real evolution_time=100/10_q;

    this->run_test("reactor",f,g,noise_levels,starting_set,evolution_time,step);
}


void TestInclusionIntegrator::test_lotka_volterra() const {
    double step=1.0/50;

    RealVector noise_levels={1/100_q,1/100_q};

    auto x = ValidatedVectorFunction::identity(2u);
    auto one = ValidatedScalarFunction::constant(2u,1_z);
    auto zero = ValidatedScalarFunction::constant(2u,0_z);
    auto three = ValidatedScalarFunction::constant(2u,3_z);

    auto f = ValidatedVectorFunction({three*x[0]*(one-x[1]),x[1]*(x[0]-one)});

    Vector<ValidatedVectorFunction> g({{x[0]*(one-x[1]),zero},{zero,x[1]*(x[0]-one)}});

    Real e=1/100000000_q;
    RealBox starting_set={{Real(1.2)-e,Real(1.2)+e},{Real(1.1)-e,Real(1.1)+e}};
    Real evolution_time=100/10_q;

    this->run_test("lotka-volterra",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_fitzhugh_nagumo() const {
    double step=1.0/20;

    RealVector noise_levels={1/10000_q,1/10000_q};

    auto x = ValidatedVectorFunction::identity(2u);
    auto one = ValidatedScalarFunction::constant(2u,1_z);
    auto zero = ValidatedScalarFunction::constant(2u,0_z);

    auto f = ValidatedVectorFunction({x[0] - x[0]*x[0]*x[0] - x[1] + Real(7.0/8),Real(8.0/100)*(x[0] + Real(0.7) - Real(0.8)*x[1])});

    Vector<ValidatedVectorFunction> g({{one,zero},{zero,x[0] + Real(0.7) - Real(0.8)*x[1]}});

    Real e=0/100_q;
    RealBox starting_set={{Real(-1.0)-e,Real(-1.0)+e},{Real(1.0)-e,Real(1.0)+e}};
    Real evolution_time=400/10_q;

    this->run_test("fitzhugh-nagumo",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_jet_engine() const {
    double step=1.0/50;

    RealVector noise_levels={5/1000_q,5/1000_q};

    auto x = ValidatedVectorFunction::identity(2u);
    auto one = ValidatedScalarFunction::constant(2u,1_z);
    auto zero = ValidatedScalarFunction::constant(2u,0_z);

    auto f = ValidatedVectorFunction({-x[1]-Real(1.5)*x[0]*x[0]-Real(0.5)*x[0]*x[0]*x[0]-Real(0.5),3*x[0]-x[1]});

    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e1=5/100_q;
    Real e2=7/100_q;
    RealBox starting_set={{Real(1.0)-e1,Real(1.0)+e1},{Real(1.0)-e2,Real(1.0)+e2}};
    Real evolution_time=40/8_q;

    this->run_test("jet-engine",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_pi_controller() const {
    double step=1.0/32;

    RealVector noise_levels={1/10_q};

    auto x = ValidatedVectorFunction::identity(2u);
    auto one = ValidatedScalarFunction::constant(2u,1_z);
    auto zero = ValidatedScalarFunction::constant(2u,0_z);

    auto f = ValidatedVectorFunction({Real(-0.101)*(x[0]-Real(20.0))+Real(1.3203)*(x[1]-Real(0.1616))-Real(0.01)*x[0]*x[0], Real(-1.0)*(Real(-0.101)*(x[0]-Real(20.0))+Real(1.3203)*(x[1]-Real(0.1616))-Real(0.01)*x[0]*x[0]) + Real(3.0)*(Real(20.0)-x[0])});

    Vector<ValidatedVectorFunction> g({{zero,one}});

    Real e=1/1024_q;
    RealBox starting_set={{Real(5.0),Real(10.0)},{-e,+e}};
    Real evolution_time=40/8_q;

    this->run_test("pi-controller",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_jerk21() const {

    double step=1.0/16;

    RealVector noise_levels={1/1000_q};

    auto x = ValidatedVectorFunction::identity(3u);
    auto one = ValidatedScalarFunction::constant(3u,1_z);
    auto zero = ValidatedScalarFunction::constant(3u,0_z);

    Real A(0.25);

    auto f = ValidatedVectorFunction({x[1],x[2],-x[2]*x[2]*x[2]-x[1]*x[0]*x[0]-A*x[0]});

    Vector<ValidatedVectorFunction> g({{zero,zero,-x[0]}});

    Real e=1/1024_q;
    Real x0_i(0.25);
    Real x1_i(0.0);
    Real x2_i(0.0);
    RealBox starting_set={{x0_i-e,x0_i+e},{x1_i-e,x1_i+e},{x2_i-e,x2_i+e}};
    Real evolution_time=100/10_q;

    this->run_test("jerk21",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_lorenz() const {
    double step=1.0/256;

    RealVector noise_levels={1/100_q};

    auto x = ValidatedVectorFunction::identity(3u);
    auto one = ValidatedScalarFunction::constant(3u,1_z);
    auto zero = ValidatedScalarFunction::constant(3u,0_z);

    Real sigma(10.0);
    Real rho(28.0);
    Real beta(8.0/3);

    auto f = ValidatedVectorFunction({(x[1]-x[0])*sigma,x[0]*(one*rho - x[2]) - x[1],x[0]*x[1] - x[2]*beta});

    Vector<ValidatedVectorFunction> g({{zero,x[0],zero}});

    Real e=1/1024_q;
    Real x0_i(1.0);
    Real x1_i(1.0);
    Real x2_i(1.0);
    RealBox starting_set={{x0_i-e,x0_i+e},{x1_i-e,x1_i+e},{x2_i-e,x2_i+e}};
    Real evolution_time=10/10_q;

    auto f2 = ValidatedVectorFunction({x[0]*2_z,x[1],x[0]*x[2]*x[2]});
    std::cout << "Df: " << f2.differential(cast_singleton(over_approximation(starting_set)),2) << std::endl;

    this->run_test("lorenz",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_rossler() const {
    double step=1.0/128;

    RealVector noise_levels={1/1000_q};

    auto x = ValidatedVectorFunction::identity(3u);
    auto one = ValidatedScalarFunction::constant(3u,1_z);
    auto zero = ValidatedScalarFunction::constant(3u,0_z);

    Real a(0.1);
    Real b(0.1);
    Real c(6.0);

    auto f = ValidatedVectorFunction({-x[1]-x[2],x[0] + x[1]*a,one*b + x[2]*(x[0]-one*c)});

    Vector<ValidatedVectorFunction> g({{zero,zero,one}});

    Real e=1/1024_q;
    Real x0_i(-9.0);
    Real x1_i(0.0);
    Real x2_i(0.01);
    RealBox starting_set={{x0_i-e,x0_i+e},{x1_i-e,x1_i+e},{x2_i-e,x2_i+e}};
    Real evolution_time=120/10_q;

    this->run_test("rossler",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_jerk16() const {
    double step=1.0/16;

    RealVector noise_levels={1/1000_q};

    auto x = ValidatedVectorFunction::identity(3u);
    auto one = ValidatedScalarFunction::constant(3u,1_z);
    auto zero = ValidatedScalarFunction::constant(3u,0_z);

    Real B(0.03);

    auto f = ValidatedVectorFunction({x[1],x[2],-x[1]+x[0]*x[0]-one*B});

    Vector<ValidatedVectorFunction> g({{zero,zero,one}});

    Real e=1/1024_q;
    Real x0_i(0.0);
    Real x1_i(0.0);
    Real x2_i(0.0);
    RealBox starting_set={{x0_i-e,x0_i+e},{x1_i-e,x1_i+e},{x2_i-e,x2_i+e}};
    Real evolution_time=100/10_q;

    this->run_test("jerk16",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_DCDC() const {
    double step=1.0/10;

    RealVector noise_levels={2/1000_q,1/15_q};

    Real k0(0.002987);
    Real fp0(-0.018);//(-11+k0)/600
    Real fp1(-0.066); //(k0-1)/15
    Real fq0(0.071);//(1-k0)/14
    Real fq1(-0.00853);//-k0*20/7
    Real gp0 = 1/600_q;
    Real gp1 = 1/15_q;
    Real gq0 = -1/14_q;
    Real gq1 = -20/7_q;

    auto x = ValidatedVectorFunction::identity(2u);
    auto one = ValidatedScalarFunction::constant(2u,1_z);
    auto zero = ValidatedScalarFunction::constant(2u,0_z);
    auto third = ValidatedScalarFunction::constant(2u,1/3_z);

    auto f = ValidatedVectorFunction({third+x[0]*fp0+x[1]*fp1,x[0]*fq0+x[1]*fq1});

    Vector<ValidatedVectorFunction> g({{gp0*x[0]+gp1*x[1],gq0*x[0]+gq1*x[1]},{one,zero}});

    Real e=1/1000000_q;
    RealBox starting_set={{Real(1)-e,Real(1)+e},{Real(5)-e,Real(5)+e}};
    Real evolution_time=50/10_q;

    this->run_test("DCDC",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_harmonic() const {

    double step=1.0/64;

    RealVector noise_levels={4/100_q};

    auto x = ValidatedVectorFunction::identity(2u);
    auto one = ValidatedScalarFunction::constant(2u,1_z);
    auto zero = ValidatedScalarFunction::constant(2u,0_z);

    auto f=ValidatedVectorFunction({x[1],-x[0]});

    Vector<ValidatedVectorFunction> g({{one,zero}});

    Real e=1/10000000_q;
    RealBox starting_set={{-e,e},{-e,+e}};
    Real evolution_time=314/100_q;

    this->run_test("harmonic",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_van_der_pol() const {

    double step=1.0/8;

    RealVector noise_levels={1/20_q,1/10000_q};

    auto x = ValidatedVectorFunction::identity(2u);
    auto one = ValidatedScalarFunction::constant(2u,1_z);
    auto zero = ValidatedScalarFunction::constant(2u,0_z);

    auto f=ValidatedVectorFunction({x[1],-x[0]+x[1]*(1 - x[0]*x[0])});

    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e=1/1024_q;
    RealBox starting_set={{Real(1.21)-e,Real(1.21)+e},{Real(2.01)-e,Real(2.01)+e}};
    Real evolution_time=8/4_q;

    this->run_test("vanderpol",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test_clock() const {

    double step=1.0/256;

    RealVector noise_levels={1/16_q,1/16_q};

    auto x = ValidatedVectorFunction::identity(2u);
    auto one = ValidatedScalarFunction::constant(2u,1_z);
    auto zero = ValidatedScalarFunction::constant(2u,0_z);

    auto f=ValidatedVectorFunction({one,one});
    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e=1/128_q;
    RealBox starting_set={{-e,e},{-e,+e}};
    Real evolution_time=20/4_q;

    this->run_test("clock",f,g,noise_levels,starting_set,evolution_time,step);
}

void TestInclusionIntegrator::test() const {

    //ARIADNE_TEST_CALL(test_wiggins_18_7_3());
    //ARIADNE_TEST_CALL(test_order7());
    //ARIADNE_TEST_CALL(test_3dsphere());
    //ARIADNE_TEST_CALL(test_vinograd());
    ARIADNE_TEST_CALL(test_higgins_selkov());
    ARIADNE_TEST_CALL(test_reactor());
    ARIADNE_TEST_CALL(test_lotka_volterra());
    ARIADNE_TEST_CALL(test_jet_engine());
    //ARIADNE_TEST_CALL(test_fitzhugh_nagumo());
    ARIADNE_TEST_CALL(test_pi_controller());
    ARIADNE_TEST_CALL(test_jerk21());
    ARIADNE_TEST_CALL(test_lorenz());
    ARIADNE_TEST_CALL(test_rossler());
    ARIADNE_TEST_CALL(test_jerk16());
    ARIADNE_TEST_CALL(test_DCDC());
    ARIADNE_TEST_CALL(test_harmonic());
    //ARIADNE_TEST_CALL(test_van_der_pol());
    //ARIADNE_TEST_CALL(test_clock());
}

int main() {
    TestInclusionIntegrator().test();
    return ARIADNE_TEST_FAILURES;
}
