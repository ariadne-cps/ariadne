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

namespace Ariadne {


template<class F, class S> List<ResultOf<F(S)>> map(F const& f, List<S> const& list) {
    List<ResultOf<F(S)>> result; for(auto item : list) { result.append(f(item)); } return result;
}

ValidatedConstrainedImageSet range(ValidatedVectorFunctionModelType const& fm) {
    return ValidatedConstrainedImageSet(fm.domain(),fm);
}

ThresholdSweeperDP make_threshold_sweeper(double d) { return ThresholdSweeperDP(DoublePrecision(),d); }

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
    EffectiveVectorFunction x; EffectiveScalarFunction one;

    Void run_test(String name, InclusionIntegratorInterface const& integrator,
                  ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, RealVector noise_levels,
                  RealBox real_starting_set, Real evolution_time) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
        PrecisionType prec;

        BoxDomainType noise=cast_exact_box(UpperIntervalType(-1,+1)*noise_levels);
        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);
        List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);
        ValidatedVectorFunctionModelType evolve_function = partial_evaluate(flow_functions.back(),starting_set.size(),NumericType(evolution_time,prec));
        ValidatedConstrainedImageSet evolve_set = range(evolve_function);

        Figure fig=Figure();

        Box<FloatDPUpperInterval> graphics_box(2);
        for (auto set: reach_sets) {
            graphics_box = hull(graphics_box,set.bounding_box());
        }

        fig.set_bounding_box(graphics_box);
        fig.set_line_colour(0.0,0.0,0.0);
        fig.set_line_style(false);
        fig.set_fill_colour(0.5,0.5,0.5);
        fig.draw(starting_set);
        fig.set_fill_colour(1.0,0.75,0.5);
        for (auto set : reverse(reach_sets)) { fig.draw(set); }
        fig.draw(evolve_set);
        fig.write(("test_differential_inclusion-"+name).c_str());
    }

  public:
    TestInclusionIntegrator();

    void test() const;
    void test_clock() const;
    void test_rotation() const;
    void test_van_der_pol() const;
    void test_jet_engine() const;
    void test_lotka_volterra() const;
};



TestInclusionIntegrator::TestInclusionIntegrator()
    : x(EffectiveVectorFunction::identity(2u)), one(EffectiveScalarFunction::constant(2u,1_z)) { }

void TestInclusionIntegrator::test() const {
    ARIADNE_TEST_CALL(test_clock());
    ARIADNE_TEST_CALL(test_rotation());
    ARIADNE_TEST_CALL(test_van_der_pol());
    ARIADNE_TEST_CALL(test_jet_engine());
    ARIADNE_TEST_CALL(test_lotka_volterra());
}

void TestInclusionIntegrator::test_lotka_volterra() const {
    auto integrator = InclusionIntegratorAffineW(make_threshold_sweeper(1e-8), step_size=1.0/32, number_of_steps_between_simplifications=4, number_of_variables_to_keep=32);
    integrator.verbosity = 0;

    RealVector noise_levels={1/100_q,1/100_q};

    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);
    auto three = EffectiveScalarFunction::constant(2u,3_z);

    auto f = EffectiveVectorFunction({three*x[0]*(one-x[1]),x[1]*(x[0]-one)});


    Vector<ValidatedVectorFunction> g({{x[0]*(one-x[1]),zero},{zero,x[1]*(x[0]-one)}});

    Real e=1/1024_q;
    RealBox starting_set={{Real(1.2)-e,Real(1.2)+e},{Real(1.1)-e,Real(1.1)+e}};
    Real evolution_time=36/10_q;

    this->run_test("lotka-volterra",integrator,f,g,noise_levels,starting_set,evolution_time);
}

void TestInclusionIntegrator::test_jet_engine() const {
    auto integrator = InclusionIntegratorAffineW(make_threshold_sweeper(1e-8), step_size=1.0/32, number_of_steps_between_simplifications=4, number_of_variables_to_keep=8);
    integrator.verbosity = 0;

    RealVector noise_levels={5/1000_q,5/1000_q};

    auto f = EffectiveVectorFunction({-x[1]-Real(1.5)*x[0]*x[0]-Real(0.5)*x[0]*x[0]*x[0]-Real(0.5),3*x[0]-x[1]});

    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);
    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e=1/10_q;
    RealBox starting_set={{Real(1.0)-e,Real(1.0)+e},{Real(1.0)-e,Real(1.0)+e}};
    Real evolution_time=40/8_q;

    this->run_test("jetengine",integrator,f,g,noise_levels,starting_set,evolution_time);
}

void TestInclusionIntegrator::test_van_der_pol() const {
    auto integrator = InclusionIntegratorAffineW(make_threshold_sweeper(1e-8), step_size=1.0/16, number_of_steps_between_simplifications=1, number_of_variables_to_keep=4);
    integrator.verbosity = 0;

    RealVector noise_levels={1/1024_q,1/1024_q};

    auto f=EffectiveVectorFunction({x[1],-x[0]+x[1]*(1 - x[0]*x[0])});

    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);
    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e=1/1024_q;
    RealBox starting_set={{Real(2.01)-e,Real(2.01)+e},{-e,+e}};
    Real evolution_time=27/4_q;

    this->run_test("vanderpol",integrator,f,g,noise_levels,starting_set,evolution_time);
}

void TestInclusionIntegrator::test_rotation() const {
    auto integrator = InclusionIntegratorAffineW(make_threshold_sweeper(1e-8), step_size=1.0/16, number_of_steps_between_simplifications=1, number_of_variables_to_keep=4);
    integrator.verbosity = 0;

    RealVector noise_levels={1/1024_q,1/1024_q};

    auto f=EffectiveVectorFunction({x[1],-x[0]});

    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);
    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e=1/1024_q;
    RealBox starting_set={{1-e,1+e},{-e,+e}};
    Real evolution_time=26/4_q;

    this->run_test("rotation",integrator,f,g,noise_levels,starting_set,evolution_time);
}

void TestInclusionIntegrator::test_clock() const {
    auto integrator = InclusionIntegratorConstantW(make_threshold_sweeper(1e-8), step_size=1.0/4, number_of_steps_between_simplifications=64, number_of_variables_to_keep=32);

    RealVector noise_levels={1/16_q,1/16_q};

    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);

    auto f=EffectiveVectorFunction({one,one});
    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e=1/128_q;
    RealBox starting_set={{-e,e},{-e,+e}};
    Real evolution_time=16/4_q;

    this->run_test("clock",integrator,f,g,noise_levels,starting_set,evolution_time);
}

int main() {
    TestInclusionIntegrator().test();
    return ARIADNE_TEST_FAILURES;
}
