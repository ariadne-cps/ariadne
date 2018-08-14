/***************************************************************************
 *            test_differential_inclusion.cpp
 *
 *  Copyright  2008-18 Luca Geretti, Pieter Collins, Sanja Zivanovic
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
#include "symbolic/expression_set.hpp"
#include "dynamics/vector_field.hpp"

#include "test.hpp"
#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>

namespace Ariadne {


template<class F, class S> List<ResultOf<F(S)>> map(F const& f, List<S> const& list) {
    List<ResultOf<F(S)>> result; for(auto item : list) { result.append(f(item)); } return result;
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


    Void run_each_approximation(String name, DifferentialInclusionIVP const& ivp, Real evolution_time, double step, List<InputApproximation> approximations, SweeperDP sweeper, SizeType freq, int verbosity) const
    {
        for (auto appro: approximations) {
            List<InputApproximation> singleapproximation = {appro};
            std::cout << appro << std::endl;
            run_single_test(name,ivp,evolution_time,step,singleapproximation,sweeper,freq,verbosity);
        }
    }

    Void run_single_test(String name, DifferentialInclusionIVP const& ivp, Real evolution_time, double step, List<InputApproximation> approximations, SweeperDP sweeper, SizeType freq, int verbosity) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
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
                fig.write(("test_differential_inclusion-"+name+num_char).c_str());
            }
        }
*/
    }

    Void run_test(String name, const DottedRealAssignments& dynamics, const RealVariablesBox& inputs,
                  const RealVariablesBox& initial, Real evolution_time, double step) const {

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

        this->run_single_test(name,ivp,evolution_time,step,approximations,sweeper,freq,verbosity);
        //this->run_each_approximation(name,ivp,evolution_time,step,approximations,sweeper,freq,verbosity);
    }

  public:
    void test() const;

    void test_wiggins_18_7_3() const;
    void test_order7() const;
    void test_3dsphere() const;
    void test_vinograd() const;
    void test_laub_loomis() const;
    void test_fitzhugh_nagumo() const;
    void test_van_der_pol() const;
    void test_clock() const;

    void test_higgins_selkov() const;
    void test_reactor() const;
    void test_lotka_volterra() const;
    void test_jet_engine() const;
    void test_pi_controller() const;
    void test_jerk21() const;
    void test_lorenz() const;
    void test_rossler() const;
    void test_jerk16() const;
    void test_DCDC() const;
    void test_harmonic() const;
};


void TestInclusionIntegrator::test_wiggins_18_7_3() const {
    RealVariable x("x"), y("y"), u1("u1"), u2("u2");
    DottedRealAssignments dynamics={dot(x)=-x+2*y+pow(x,2)*y+pow(x,4)*pow(y,5)+u1,
    		                             dot(y)=-y-pow(x,4)*pow(y,6)+pow(x,8)*pow(y,9)+u2};
    RealVariablesBox inputs={-2/100_q<=u1<=2/100_q,-2/100_q<=u2<=2/100_q};

    auto e=1/10000000_q;
    RealVariablesBox initial={{1/3_q-e<=x<=1/3_q+e},{1/3_q-e<=y<=1/3_q+e}};

    Real evolution_time=9;
    double step=1.0/16;

    this->run_test("wiggins_18_7_3",dynamics,inputs,initial,evolution_time,step);
}


void TestInclusionIntegrator::test_order7() const {
    RealVariable x("x"), y("y"), u("u");
    DottedRealAssignments dynamics={dot(x)=-42*pow(x,7)+68*pow(x,6)*y-46*pow(x,5)*y+256*pow(x,4)*y+156*pow(x,3)*y+50*pow(x,2)*y+20*x*pow(y,6)-8*pow(y,7),
                                         dot(y)=y*(1110*pow(x,6)-220*pow(x,5)*y-3182*pow(x,4)*y+478*pow(x,3)*pow(y,3)+487*pow(x,2)*pow(y,4)-102*x*pow(y,5)-12*pow(y,6))+u};
    RealVariablesBox inputs={-1/100_q<=u<=1/100_q};

    auto e=1/10_q;
    RealVariablesBox initial={{-1-e<=x<=-1+e},{1-e<=y<=1+e}};

    Real evolution_time=5;
    double step=1.0/32;

    this->run_test("order7",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_3dsphere() const {
    RealVariable x("x"), y("y"), z("z"), u1("u1"), u2("u2"), u3("u3");
    RealExpression cuberadius(pow(x,3)+pow(y,3)+pow(z,3));
    DottedRealAssignments dynamics={dot(x)=pow(x,2) - u1*x*cuberadius,
                                    dot(y)=pow(y,2) - u2*y*cuberadius,
                                    dot(z)=pow(z,2) - u3*z*cuberadius};
    RealVariablesBox inputs={-1/100_q+1<=u1<=1/100_q+1,-1/100_q+1<=u2<=1/100_q+1,-1/100_q+1<=u3<=1/100_q+1};

    Real e=1/1000000000_q;
    RealVariablesBox initial={{1/4_q-e<=x<=1/4_q+e},{1/8_q-e<=y<=1/8_q+e},{1/10_q-e<=z<=1/10_q+e}};

    Real evolution_time=9;
    double step=1.0/16;

    this->run_test("3dsphere",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_vinograd() const {
    RealVariable x("x"), y("y"), u1("u1"), u2("u2");
    DottedRealAssignments dynamics={dot(x)=pow(y,5)+pow(x,2)*(-x+y)+u1,
                                         dot(y)=pow(y,2)*(-2*x+y)+u2};
    RealVariablesBox inputs={-1/1000_q<=u1<=1/1000_q,-1/1000_q<=u2<=1/1000_q};

    Real e=1/10000000_q;
    RealVariablesBox initial={{1-e<=x<=1+e},{-e<=y<=+e}};

    Real evolution_time=18;
    double step=1.0/16;

    this->run_test("vinograd",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_laub_loomis() const {
    RealVariable x1("x1"), x2("x2"), x3("x3"), x4("x4"), x5("x5"), x6("x6"), x7("x7"), u("u");
    DottedRealAssignments dynamics={dot(x1)=1.4_dec*x3-0.9_dec*x1,
                                    dot(x2)=2.5_dec*x5-1.5_dec*x2,
                                    dot(x3)=0.6_dec*x7-0.8_dec*x2*x3,
                                    dot(x4)=2-1.3_dec*x3*x4,
                                    dot(x5)=0.7_dec*x1-x4*x5,
                                    dot(x6)=0.3_dec*x1-3.1_dec*x6,
                                    dot(x7)=1.8_dec*x6-1.5_dec*x2*x7,
                                    };
    RealVariablesBox inputs={0.0_dec<=u<=0.0_dec};

    Real W0=2/100_q;
    RealVariablesBox initial={{1.2_dec-W0<=x1<=1.2_dec+W0},{1.05_dec-W0<=x2<=1.05_dec+W0},{1.5_dec-W0<=x3<=1.5_dec+W0},
                              {2.4_dec-W0<=x4<=2.4_dec+W0},{1-W0<=x5<=1+W0},{0.1_dec-W0<=x6<=0.1_dec+W0},{0.45_dec-W0<=x7<=0.45_dec+W0}};

    Real evolution_time=10;
    double step=2.0/100;

    this->run_test("laub-loomis",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_fitzhugh_nagumo() const {
    RealVariable x("x"), y("y"), u1("u1"), u2("u2");
    DottedRealAssignments dynamics={dot(x)=x-pow(x,3)-y+u1,dot(y)=x+u2-0.8_dec*y};
    RealVariablesBox inputs={0.874_dec<=u1<=0.876_dec,0.699_dec<=u2<=0.701_dec};

    Real e=1/100_q;
    RealVariablesBox initial={{-1-e<=x<=-1+e},{1-e<=x<=1+e}};

    Real evolution_time=10;
    double step=1.0/20;

    this->run_test("fitzhugh-nagumo",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_van_der_pol() const {
    RealVariable x("x"), y("y"), u1("u1"), u2("u2");
    DottedRealAssignments dynamics={dot(x)=y+u1,dot(y)=-x+y*(1-pow(x,2))+u2};
    RealVariablesBox inputs={-1/20_q<=u1<=1/20_q,-1/10000_q<=u2<=1/10000_q};

    Real e=1/1024_q;
    RealVariablesBox initial={{1.21_dec-e<=x<=1.21_dec+e},{2.01_dec-e<=y<=2.01_dec+e}};

    Real evolution_time=24/4_q;
    double step=1.0/8;

    this->run_test("vanderpol",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_clock() const {
    RealVariable x("x"), y("y"), u1("u1"), u2("u2");
    DottedRealAssignments dynamics={dot(x)=u1,dot(y)=u2};
    RealVariablesBox inputs={0.95_dec<=u1<=1.05_dec,0.95_dec<=u2<=1.05_dec};

    Real e=1/128_q;
    RealVariablesBox initial={{x==0},{y==0}};

    auto evolution_time=5;
    double step=1.0/256;

    this->run_test("clock",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_higgins_selkov() const {
    RealVariable S("S"), P("P"), v0("v0"), k1("k1"), k2("k2");
    DottedRealAssignments dynamics={dot(S)=v0-S*k1*pow(P,2),dot(P)=S*k1*pow(P,2)-k2*P};
    RealVariablesBox inputs={0.9998_dec<=v0<=1.0002_dec,0.9998_dec<=k1<=1.0002_dec,0.99981_dec<=k2<=1.00021_dec};

    Real e=1/100_q;
    RealVariablesBox initial={{2-e<=S<=2+e},{1-e<=P<=1+e}};

    Real evolution_time=10;
    double step=1.0/50;

    this->run_test("higgins-selkov",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_reactor() const {
    RealVariable xA("xA"), xB("xB"), xC("xC"), xD("xD"), u1("u1"), u2("u2"), u3("u3");
    DottedRealAssignments dynamics={dot(xA)=-u3*xA*xB-0.4_dec*xA*xC+0.05_dec*u1-0.1_dec*xA,
                                    dot(xB)=-u3*xA*xB+0.05_dec*u2-0.1_dec*xB,
                                    dot(xC)=u3*xA*xB-0.4_dec*xA*xC-0.1_dec*xC,
                                    dot(xD)=0.4_dec*xA*xC-0.1_dec*xD};
    RealVariablesBox inputs={0.999_dec<=u1<=1.001_dec,0.899_dec<=u2<=0.901_dec,29.8_dec<=u3<=30.2_dec};

    Real e=1/1000000_q;
    RealVariablesBox initial={{0<=xA<=e},{0<=xB<=e},{0<=xC<=e},{0<=xD<=e}};

    Real evolution_time=10;
    double step=1.0/16;

    this->run_test("reactor",dynamics,inputs,initial,evolution_time,step);
}


void TestInclusionIntegrator::test_lotka_volterra() const {
    RealVariable x("x"), y("y"), u1("u1"), u2("u2");
    DottedRealAssignments dynamics={dot(x)=u1*x*(1-y),dot(y)=u2*y*(x-1)};
    RealVariablesBox inputs={2.99_dec<=u1<=3.01_dec,0.99_dec<=u2<=1.01_dec};

    RealVariablesBox initial={{x==1.2_dec},{y==1.1_dec}};

    Real evolution_time=10;
    double step=1.0/50;

    this->run_test("lotka-volterra",dynamics,inputs,initial,evolution_time,step);
}


void TestInclusionIntegrator::test_jet_engine() const {
    RealVariable x("x"), y("y"), u1("u1"), u2("u2");
    DottedRealAssignments dynamics={dot(x)=-y-1.5_dec*pow(x,2)-0.5_dec*pow(x,3)-0.5_dec+u1,dot(y)=3*x-y+u2};
    RealVariablesBox inputs={-5/1000_q<=u1<=5/1000_q,-5/1000_q<=u2<=5/1000_q};

    Real e1=5/100_q; Real e2=7/100_q;
    RealVariablesBox initial={{1-e1<=x<=1+e1},{1-e2<=y<=1+e2}};

    Real evolution_time=5;
    double step=1.0/50;

    this->run_test("jet-engine",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_pi_controller() const {
    RealVariable v("v"), x("x"), u("u");
    RealExpression dynv = -0.101_dec*(v-20)+1.3203_dec*(x-0.1616_dec)-0.01_dec*pow(v,2);
    DottedRealAssignments dynamics={dot(v)=dynv,dot(x)=-dynv + 3*(20-v) + u};
    RealVariablesBox inputs={-1/10_q<=u<=1/10_q};

    Real e=1/1024_q;
    RealVariablesBox initial={{5<=v<=10},{-e<=x<=+e}};

    Real evolution_time=5;
    double step=1.0/32;

    this->run_test("pi-controller",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_jerk21() const {
    RealVariable x("x"), y("y"), z("z"), u("u");
    DottedRealAssignments dynamics={dot(x)=y,dot(y)=z,dot(z)=-pow(z,3)-y*pow(x,2)-u*x};
    RealVariablesBox inputs={0.249_dec<=u<=0.251_dec};

    Real e=1/1024_q;
    RealVariablesBox initial={{0.25_dec-e<=x<=0.25_dec+e},{-e<=y<=e},{-e<=z<=e}};

    Real evolution_time=10;
    double step=1.0/16;

    this->run_test("jerk21",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_lorenz() const {
    RealVariable x("x"), y("y"), z("z"), u("u");
    DottedRealAssignments dynamics={dot(x)=10*(y-x),
    									 dot(y)=x*(28 - z) - y + x*u,
										 dot(z)=x*y - z*8/3_q};
    RealVariablesBox inputs={-1/100_q<=u<=1/100_q};

    Real e=1/1024_q;
    RealVariablesBox initial={{1-e<=x<=1+e},{1-e<=y<=1+e},{1-e<=z<=1+e}};

    Real evolution_time=1;
    double step=1.0/256;

    this->run_test("lorenz",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_rossler() const {
    RealVariable x("x"), y("y"), z("z"), u("u");
    DottedRealAssignments dynamics={dot(x)=-y-z,dot(y)=x+y*0.1_dec,dot(z)=z*(x-6)+u};
    RealVariablesBox inputs={0.099_dec<=u<=0.101_dec};

    Real e=1/1024_q;
    RealVariablesBox initial={{-9-e<=x<=-9+e},{-e<=y<=e},{0.01_dec-e<=z<=0.01_dec+e}};

    Real evolution_time=12;
    double step=1.0/128;

    this->run_test("rossler",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_jerk16() const {
    RealVariable x("x"), y("y"), z("z"), u("u");
    DottedRealAssignments dynamics={dot(x)=y,dot(y)=z,dot(z)=-y+pow(x,2)+u};
    RealVariablesBox inputs={-0.031_dec<=u<=-0.029_dec};

    Real e=1/1024_q;
    RealVariablesBox initial={{-e<=x<=e},{-e<=y<=e},{-e<=z<=e}};

    Real evolution_time=10;
    double step=1.0/16;

    this->run_test("jerk16",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_DCDC() const {
    Real k0(0.002987);
    Real fp0 = (-11+k0)/600;
    Real fp1 = (k0-1)/15_q;
    Real fq0 = (1-k0)/14;
    Real fq1 = -k0*20/7_q;
    Real gp0 = 1/600_q;
    Real gp1 = 1/15_q;
    Real gq0 = -1/14_q;
    Real gq1 = -20/7_q;

    RealVariable x("x"), y("y"), u1("u1"), u2("u2");
    DottedRealAssignments dynamics={dot(x)=x*fp0+y*fp1+u1*(gp0*x+gp1*y)+u2,
    							    dot(y)=x*fq0+y*fq1+u1*(gq0*x+gq1*y)};
    RealVariablesBox inputs={-2/1000_q<=u1<=2/1000_q,4/15_q<=u2<=6/15_q};

    RealVariablesBox initial={{x==1},{y==5}};

    Real evolution_time=5;
    double step=1.0/10;

    this->run_test("DCDC",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_harmonic() const {
    RealVariable x("x"), y("y"), u("u");
    DottedRealAssignments dynamics={dot(x)=y+u,dot(y)=-x};
    RealVariablesBox inputs={-4/100_q<=u<=4/100_q};

    Real e=1/10000000_q;
    RealVariablesBox initial={{-e<=x<=e},{-e<=y<=e}};

    Real evolution_time=3.141592_dec;
    double step=1.0/64;

    this->run_test("harmonic",dynamics,inputs,initial,evolution_time,step);
}


void TestInclusionIntegrator::test() const {

    //ARIADNE_TEST_CALL(test_wiggins_18_7_3());
    //ARIADNE_TEST_CALL(test_order7());
    //ARIADNE_TEST_CALL(test_3dsphere());
    //ARIADNE_TEST_CALL(test_vinograd());
    //ARIADNE_TEST_CALL(test_laub_loomis());
    //ARIADNE_TEST_CALL(test_fitzhugh_nagumo());
    //ARIADNE_TEST_CALL(test_van_der_pol());
    //ARIADNE_TEST_CALL(test_clock());
    /*ARIADNE_TEST_CALL(test_higgins_selkov());
    ARIADNE_TEST_CALL(test_reactor());
    ARIADNE_TEST_CALL(test_lotka_volterra());
    ARIADNE_TEST_CALL(test_jet_engine());
    ARIADNE_TEST_CALL(test_pi_controller());
    ARIADNE_TEST_CALL(test_jerk21());
    ARIADNE_TEST_CALL(test_lorenz());
    ARIADNE_TEST_CALL(test_rossler());
    ARIADNE_TEST_CALL(test_jerk16());*/
    ARIADNE_TEST_CALL(test_DCDC());
    ARIADNE_TEST_CALL(test_harmonic());
}

int main() {
    TestInclusionIntegrator().test();
    return ARIADNE_TEST_FAILURES;
}
