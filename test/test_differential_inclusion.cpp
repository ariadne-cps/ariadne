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

    Void run_battery_fixed_avg(String name,
                                ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, RealVector noise_levels,
                                RealBox real_starting_set, Real evolution_time, double step, SizeType avg, SizeType min_freq, SizeType max_freq, SizeType ppi) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
        PrecisionType prec;

        BoxDomainType noise=cast_exact_box(UpperIntervalType(-1,+1)*noise_levels);
        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        SizeType n = f.result_size();
        SizeType m = g.size();

        SizeType pps = n+ppi*m;

        for (auto freq : range(min_freq,max_freq+1)) {

            SizeType base = avg-n-pps*(freq-1u)/2u;

            auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);

            tms start_time, end_time;
            times(&start_time);

            List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);

            times(&end_time);
            clock_t ticks = end_time.tms_utime - start_time.tms_utime;
            clock_t const hz = sysconf(_SC_CLK_TCK);

            List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);
            ValidatedVectorFunctionModelType evolve_function = partial_evaluate(flow_functions.back(),starting_set.size(),NumericType(evolution_time,prec));
            ValidatedConstrainedImageSet evolve_set = range(evolve_function);

            FloatDPUpperBound total_diameter(0.0);
            auto ebb = evolve_set.bounding_box();
            for (auto i : range(ebb.size())) {
                total_diameter += ebb[i].width();
            }
            std::cout << freq << " (b: " << base << "): " << total_diameter << ", " << ticks / hz << "." << ticks % hz << "s" << std::endl;
        }
    }

    Void run_battery_fixed_relativebase(String name,
                                ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, RealVector noise_levels,
                                RealBox real_starting_set, Real evolution_time, double step, SizeType ppi, double ratio, SizeType min_freq, SizeType max_freq) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
        PrecisionType prec;

        BoxDomainType noise=cast_exact_box(UpperIntervalType(-1,+1)*noise_levels);
        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        SizeType n = f.result_size();
        SizeType m = g.size();

        SizeType pps = n+ppi*m;

        for (auto freq : range(min_freq,max_freq+1)) {

            SizeType base(round(ratio*pps*freq));

            auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);

            tms start_time, end_time;
            times(&start_time);

            List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);

            times(&end_time);
            clock_t ticks = end_time.tms_utime - start_time.tms_utime;
            clock_t const hz = sysconf(_SC_CLK_TCK);

            List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);
            ValidatedVectorFunctionModelType evolve_function = partial_evaluate(flow_functions.back(),starting_set.size(),NumericType(evolution_time,prec));
            ValidatedConstrainedImageSet evolve_set = range(evolve_function);

            FloatDPUpperBound total_diameter(0.0);
            auto ebb = evolve_set.bounding_box();
            for (auto i : range(ebb.size())) {
                total_diameter += ebb[i].width();
            }
            std::cout << freq << " (b: " << base << "): " << total_diameter << ", " << ticks / hz << "." << ticks % hz << "s" << std::endl;
        }
    }


    Void run_battery_fixed_absolutebase(String name,
                                        ValidatedVectorFunction const &f, Vector<ValidatedVectorFunction> const &g,
                                        RealVector noise_levels,
                                        RealBox real_starting_set, Real evolution_time, double step, SizeType base,
                                        SizeType min_freq, SizeType max_freq) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
        PrecisionType prec;

        BoxDomainType noise=cast_exact_box(UpperIntervalType(-1,+1)*noise_levels);
        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        for (auto freq : range(min_freq,max_freq+1)) {

            auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);

            tms start_time, end_time;
            times(&start_time);

            List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);

            times(&end_time);
            clock_t ticks = end_time.tms_utime - start_time.tms_utime;
            clock_t const hz = sysconf(_SC_CLK_TCK);

            List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);
            ValidatedVectorFunctionModelType evolve_function = partial_evaluate(flow_functions.back(),starting_set.size(),NumericType(evolution_time,prec));
            ValidatedConstrainedImageSet evolve_set = range(evolve_function);

            FloatDPUpperBound total_diameter(0.0);
            auto ebb = evolve_set.bounding_box();
            for (auto i : range(ebb.size())) {
                total_diameter += ebb[i].width();
            }
            std::cout << freq << ": " << total_diameter << ", " << ticks / hz << "." << ticks % hz << "s" << std::endl;
        }
    }

    Void run_test(String name, InclusionIntegratorInterface& integrator,
                  ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, RealVector noise_levels,
                  RealBox real_starting_set, Real evolution_time) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
        PrecisionType prec;

        BoxDomainType noise=cast_exact_box(UpperIntervalType(-1,+1)*noise_levels);
        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        std::cout << "flowing..." << std::endl;

        tms start_time, end_time;
        times(&start_time);

        List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);

        times(&end_time);
        clock_t ticks = end_time.tms_utime - start_time.tms_utime;
        clock_t const hz = sysconf(_SC_CLK_TCK);

        List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);
        ValidatedVectorFunctionModelType evolve_function = partial_evaluate(flow_functions.back(),starting_set.size(),NumericType(evolution_time,prec));
        ValidatedConstrainedImageSet evolve_set = range(evolve_function);

        FloatDPUpperBound total_diameter(0.0);
        auto ebb = evolve_set.bounding_box();
        for (auto i : range(ebb.size())) {
            total_diameter += ebb[i].width();
        }
        std::cout << "total diameter: " << total_diameter << ", " << ticks / hz << "." << ticks % hz << "s" << std::endl;

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
    }

  public:
    void test() const;

    void test_reactor() const;
    void test_lorenz() const;
    void test_rossler() const;
    void test_jerk21() const;
    void test_jerk16() const;
    void test_higgins_selkov() const;
    void test_lotka_volterra() const;
    void test_jet_engine() const;
    void test_pi_controller() const;
    void test_van_der_pol() const;
    void test_circle() const;
    void test_clock() const;
};

void TestInclusionIntegrator::test() const {
    ARIADNE_TEST_CALL(test_reactor());
    //ARIADNE_TEST_CALL(test_lorenz());
    //ARIADNE_TEST_CALL(test_rossler());
    //ARIADNE_TEST_CALL(test_jerk21());
    //ARIADNE_TEST_CALL(test_jerk16());
    //ARIADNE_TEST_CALL(test_higgins_selkov());
    //ARIADNE_TEST_CALL(test_lotka_volterra());
    //ARIADNE_TEST_CALL(test_jet_engine());
    //ARIADNE_TEST_CALL(test_pi_controller());
    //ARIADNE_TEST_CALL(test_van_der_pol());
    //ARIADNE_TEST_CALL(test_circle());
    //ARIADNE_TEST_CALL(test_clock());
}

void TestInclusionIntegrator::test_reactor() const {
    double step=1.0/32;
    SizeType freq=5;
    SizeType base=60;
    SizeType avg=120;

    auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/100_q,1/100_q,1/100_q};

    auto x = EffectiveVectorFunction::identity(4u);
    auto one = EffectiveScalarFunction::constant(4u,1_z);
    auto zero = EffectiveScalarFunction::constant(4u,0_z);

    Real U3(30.0);
    Real k2(0.4);
    Real iV(0.05);
    Real ka(0.1);
    Real kb(0.045);

    auto f = EffectiveVectorFunction({-U3*x[0]*x[1]-k2*x[0]*x[2]+iV-ka*x[0],-U3*x[0]*x[1]+kb-ka*x[1],
                                      U3*x[0]*x[1]-k2*x[0]*x[1]-ka*x[2],k2*x[0]*x[2]-ka*x[3]});

    Vector<ValidatedVectorFunction> g({{one*iV,zero,zero,zero},{zero,one*iV,zero,zero},{x[0]*x[1],x[0]*x[1],x[0]*x[1],zero}});

    Real e=1/1000000_q;
    RealBox starting_set={{0,e},{0,e},{0,e},{0,e}};
    Real evolution_time=80/10_q;

    this->run_test("reactor",integrator,f,g,noise_levels,starting_set,evolution_time);
    //this->run_battery_fixed_absolutebase("reactor",f,g,noise_levels,starting_set,evolution_time,step,base,1,15);
}

void TestInclusionIntegrator::test_lorenz() const {
    double step=1.0/256;
    SizeType freq=5;
    SizeType base=50;

    auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/100_q};

    auto x = EffectiveVectorFunction::identity(3u);
    auto one = EffectiveScalarFunction::constant(3u,1_z);
    auto zero = EffectiveScalarFunction::constant(3u,0_z);

    Real sigma(10.0);
    Real rho(28.0);
    Real beta(8.0/3);

    auto f = EffectiveVectorFunction({(x[1]-x[0])*sigma,x[0]*(one*rho - x[2]) - x[1],x[0]*x[1] - x[2]*beta});

    Vector<ValidatedVectorFunction> g({{zero,x[0],zero}});

    Real e=1/1024_q;
    Real x0_i(1.0);
    Real x1_i(1.0);
    Real x2_i(1.0);
    RealBox starting_set={{x0_i-e,x0_i+e},{x1_i-e,x1_i+e},{x2_i-e,x2_i+e}};
    Real evolution_time=10/10_q;

    this->run_test("lorenz",integrator,f,g,noise_levels,starting_set,evolution_time);
    //this->run_battery_fixed_absolutebase("lorenz",f,g,noise_levels,starting_set,evolution_time,step,base,1,12);
}

void TestInclusionIntegrator::test_rossler() const {
    double step=1.0/128;
    SizeType freq=12;
    SizeType base=20;
    SizeType avg=120;

    auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/1000_q};

    auto x = EffectiveVectorFunction::identity(3u);
    auto one = EffectiveScalarFunction::constant(3u,1_z);
    auto zero = EffectiveScalarFunction::constant(3u,0_z);

    Real a(0.1);
    Real b(0.1);
    Real c(6.0);

    auto f = EffectiveVectorFunction({-x[1]-x[2],x[0] + x[1]*a,one*b + x[2]*(x[0]-one*c)});

    Vector<ValidatedVectorFunction> g({{zero,zero,one}});

    Real e=1/1024_q;
    Real x0_i(-9.0);
    Real x1_i(0.0);
    Real x2_i(0.01);
    RealBox starting_set={{x0_i-e,x0_i+e},{x1_i-e,x1_i+e},{x2_i-e,x2_i+e}};
    Real evolution_time=120/10_q;

    //this->run_test("rossler",integrator,f,g,noise_levels,starting_set,evolution_time);
    this->run_battery_fixed_absolutebase("rossler",f,g,noise_levels,starting_set,evolution_time,step,base,1,25);
}

void TestInclusionIntegrator::test_jerk21() const {

    double step=1.0/16;
    SizeType freq=12;
    SizeType base=30;

    auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/1000_q};

    auto x = EffectiveVectorFunction::identity(3u);
    auto one = EffectiveScalarFunction::constant(3u,1_z);
    auto zero = EffectiveScalarFunction::constant(3u,0_z);

    Real A(0.25);

    auto f = EffectiveVectorFunction({x[1],x[2],-x[2]*x[2]*x[2]-x[1]*x[0]*x[0]-A*x[0]});

    Vector<ValidatedVectorFunction> g({{zero,zero,-x[0]}});

    Real e=1/1024_q;
    Real x0_i(0.25);
    Real x1_i(0.0);
    Real x2_i(0.0);
    RealBox starting_set={{x0_i+e,x0_i+e},{x1_i-e,x1_i+e},{x2_i-e,x2_i+e}};
    Real evolution_time=100/10_q;

    //this->run_test("jerk21",integrator,f,g,noise_levels,starting_set,evolution_time);
    this->run_battery_fixed_absolutebase("jerk21",f,g,noise_levels,starting_set,evolution_time,step,base,1,16);
}

void TestInclusionIntegrator::test_jerk16() const {
    double step=1.0/16;
    SizeType freq=12;
    SizeType base=30;

    auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/1000_q};

    auto x = EffectiveVectorFunction::identity(3u);
    auto one = EffectiveScalarFunction::constant(3u,1_z);
    auto zero = EffectiveScalarFunction::constant(3u,0_z);

    Real B(0.03);

    auto f = EffectiveVectorFunction({x[1],x[2],-x[1]+x[0]*x[0]-one*B});

    Vector<ValidatedVectorFunction> g({{zero,zero,one}});

    Real e=1/1024_q;
    Real x0_i(0.0);
    Real x1_i(0.0);
    Real x2_i(0.0);
    RealBox starting_set={{x0_i+e,x0_i+e},{x1_i-e,x1_i+e},{x2_i-e,x2_i+e}};
    Real evolution_time=100/10_q;

    //this->run_test("jerk16",integrator,f,g,noise_levels,starting_set,evolution_time);
    this->run_battery_fixed_absolutebase("jerk16",f,g,noise_levels,starting_set,evolution_time,step,base,1,12);
}

void TestInclusionIntegrator::test_higgins_selkov() const {

    double step=1.0/200;
    SizeType freq=2;
    SizeType base=60;
    SizeType avg=100;

    auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={2/10000_q,2/10000_q,2/10000_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);

    Real k1(1.00001);

    auto f = EffectiveVectorFunction({one-x[0]*k1*x[1]*x[1],x[0]*k1*x[1]*x[1] - x[1]});

    Vector<ValidatedVectorFunction> g({{one,zero},{-x[0]*x[1]*x[1],x[0]*x[1]*x[1]},{zero,-x[1]}});

    Real e=1/100_q;
    Real x0_i(2.0);
    Real x1_i(1.0);
    RealBox starting_set={{x0_i-e,x0_i+e},{x1_i-e,x1_i+e}};
    Real evolution_time=100/10_q;

    //this->run_test("higgins-selkov",integrator,f,g,noise_levels,starting_set,evolution_time);
    this->run_battery_fixed_absolutebase("higgins-selkov",f,g,noise_levels,starting_set,evolution_time,step,base,1,20);
}

void TestInclusionIntegrator::test_lotka_volterra() const {
    double step=1.0/128;
    SizeType freq=11;
    SizeType base=31;
    SizeType avg=100;

    auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/100_q,1/100_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);
    auto three = EffectiveScalarFunction::constant(2u,3_z);

    auto f = EffectiveVectorFunction({three*x[0]*(one-x[1]),x[1]*(x[0]-one)});


    Vector<ValidatedVectorFunction> g({{x[0]*(one-x[1]),zero},{zero,x[1]*(x[0]-one)}});

    Real e=1/1024_q;
    RealBox starting_set={{Real(1.2)-e,Real(1.2)+e},{Real(1.1)-e,Real(1.1)+e}};
    Real evolution_time=36/10_q;

    //this->run_test("lotka-volterra",integrator,f,g,noise_levels,starting_set,evolution_time);
    this->run_battery_fixed_absolutebase("lotka-volterra",f,g,noise_levels,starting_set,evolution_time,step,base,1,15);
}

void TestInclusionIntegrator::test_jet_engine() const {
    double step=1.0/50;
    SizeType base=30;
    SizeType freq=10;
    SizeType avg=100;

    auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={5/1000_q,5/1000_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);

    auto f = EffectiveVectorFunction({-x[1]-Real(1.5)*x[0]*x[0]-Real(0.5)*x[0]*x[0]*x[0]-Real(0.5),3*x[0]-x[1]});

    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e1=5/100_q;
    Real e2=7/100_q;
    RealBox starting_set={{Real(1.0)-e1,Real(1.0)+e1},{Real(1.0)-e2,Real(1.0)+e2}};
    Real evolution_time=40/8_q;

    //this->run_test("jet-engine",integrator,f,g,noise_levels,starting_set,evolution_time);
    this->run_battery_fixed_absolutebase("jet-engine", f, g, noise_levels, starting_set, evolution_time, step, base, 1, 12);
}

void TestInclusionIntegrator::test_pi_controller() const {
    double step=1.0/32;
    SizeType base=30;
    SizeType freq=10;
    SizeType avg=100;

    auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/10_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);

    auto f = EffectiveVectorFunction({Real(-0.101)*(x[0]-Real(20.0))+Real(1.3203)*(x[1]-Real(0.1616))-Real(0.01)*x[0]*x[0], Real(-1.0)*(Real(-0.101)*(x[0]-Real(20.0))+Real(1.3203)*(x[1]-Real(0.1616))-Real(0.01)*x[0]*x[0]) + Real(3.0)*(Real(20.0)-x[0])});

    Vector<ValidatedVectorFunction> g({{zero,one}});

    Real e=1/1024_q;
    RealBox starting_set={{Real(5.0),Real(10.0)},{-e,+e}};
    Real evolution_time=40/8_q;

    //this->run_test("pi-controller",integrator,f,g,noise_levels,starting_set,evolution_time);
    this->run_battery_fixed_absolutebase("pi-controller", f, g, noise_levels, starting_set, evolution_time, step, base, 26, 35);
}

void TestInclusionIntegrator::test_van_der_pol() const {
    auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=1.0/16, number_of_steps_between_simplifications=16, number_of_variables_to_keep=10);
    integrator.verbosity = 2;

    RealVector noise_levels={1/1024_q,1/1024_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);

    auto f=EffectiveVectorFunction({x[1],-x[0]+x[1]*(1 - x[0]*x[0])});

    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e=1/1024_q;
    RealBox starting_set={{Real(2.01)-e,Real(2.01)+e},{-e,+e}};
    Real evolution_time=27/4_q;

    this->run_test("vanderpol",integrator,f,g,noise_levels,starting_set,evolution_time);
}

void TestInclusionIntegrator::test_circle() const {
    auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=1.0/16, number_of_steps_between_simplifications=11, number_of_variables_to_keep=25);
    integrator.verbosity = 2;

    RealVector noise_levels={1/1024_q,1/1024_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);

    auto f=EffectiveVectorFunction({x[1],-x[0]});

    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e=1/1024_q;
    RealBox starting_set={{1-e,1+e},{-e,+e}};
    Real evolution_time=26/4_q;

    this->run_test("circle",integrator,f,g,noise_levels,starting_set,evolution_time);
}

void TestInclusionIntegrator::test_clock() const {
    auto integrator = InclusionIntegrator(make_threshold_sweeper(1e-8), step_size=1.0/256, number_of_steps_between_simplifications=11, number_of_variables_to_keep=24);
    integrator.verbosity = 2;

    RealVector noise_levels={1/16_q,1/16_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);

    auto f=EffectiveVectorFunction({one,one});
    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e=1/128_q;
    RealBox starting_set={{-e,e},{-e,+e}};
    Real evolution_time=20/4_q;

    this->run_test("clock",integrator,f,g,noise_levels,starting_set,evolution_time);
}

int main() {
    TestInclusionIntegrator().test();
    return ARIADNE_TEST_FAILURES;
}
