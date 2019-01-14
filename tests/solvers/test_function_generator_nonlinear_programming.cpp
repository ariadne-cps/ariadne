/***************************************************************************
 *            test_function_generator_nonlinear_programming.cpp
 *
 *          Copyright  2018  Nicola Dessì
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

#include <chrono>
#include <fstream>
#include <iostream>

#include "config.hpp"

#include "../test.hpp"

#include "algebra/algebra.hpp"
#include "algebra/vector.hpp"
#include "function/formula.hpp"
#include "function/function.hpp"
#include "function/taylor_model.hpp"
#include "geometry/box.hpp"
#include "numeric/numeric.hpp"
#include "solvers/nonlinear_programming.hpp"

#include "DataStream.hpp"
#include "RandomFunction.hpp"

using namespace std;
using namespace Ariadne;

#define M_PI 3.14159265358979323846 /* pi */
#define RED "\033[1;31m\n"
#define GREEN "\033[1;32m\n"
#define WHITE "\033[0m\n"

#define ci(a) (static_cast<int>(a))
#define cd(a) (static_cast<double>(a))
#define cfdp(a) (static_cast<FloatDP>(a))
class DomainException
{
};

class TestOptimiser
{
  private:
    std::unique_ptr<OptimiserInterface> optimiser;
    // DoublePrecision pr;

  public:
    TestOptimiser(const OptimiserInterface &opt) : optimiser(opt.clone()) {}

    Void test(const std::string &name, unsigned seed)
    {
        test_1(name+"_1",seed);  //!< min F_1
        test_2(name+"_2",seed);  //!< min F_1 s.t. G_1
        test_3(name+"_3",seed);  //!< min F_2 s.t. G_2
        test_4(name + "_4", seed); //!< min F_4 s.t. G_4
        test_5(name + "_5", seed); //!< min F_5 s.t. G_3
        test_X(name+"_X",seed);  //!< min F_1 s.t. G_3
        // performance_test();
        // null_gradient_test();
    }


    Void test_1(const std::string &name, unsigned seed)
    {

        std::cerr << "Testing " << name << "...\n";
        DataStream ds = csv::make_datastream();
        ds.setHeader({"#sequence_" + name, "nVars_" + name, "order_" + name,
                      "correct_" + name, "f_obj_" + name, "time_" + name,
                      "status_" + name});
        ds.open(name + ".csv");

        // FunctionDistribution<EffectiveVectorMultivariateFunction>
        // generator_g(2,2);
        EffectiveScalarMultivariateFunction f;
        FloatBoundsVector x_optimal;

        float elapsed_time;
        FloatDP f_obj;
        int index = 0;
        unsigned maxNVars = 10;
        int maxOrder = 10;
        for (unsigned nVars = 1; nVars < maxNVars; ++nVars) {
            for (int order = 1; order < maxOrder; ++order, ++index) {
                fprintf(stderr, "%3d %% ]\r",
                        ci((static_cast<double>(ci(index) * 100) /
                            (ci(maxNVars) * maxOrder))));
                try {
                    single_test_1(seed, nVars, order, f_obj, elapsed_time, f,
                                  x_optimal);
                    ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order),
                                        cfdp(1), f_obj, cfdp(elapsed_time), 1});
                } catch (SQPStatusLeq1Exception sqpsl1) {
                    ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order),
                                        sqpsl1.st, f_obj, cfdp(elapsed_time),
                                        sqpsl1.st});
                } catch (std::runtime_error e) {
                    ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order),
                                        cfdp(-15), f_obj, cfdp(elapsed_time),
                                        -15});
                }
            }
        }

        // ds.close();
    }

    Void test_2(const std::string &name, unsigned seed)
    {

        EffectiveScalarMultivariateFunction f;
        EffectiveVectorMultivariateFunction g;
        FloatBoundsVector x_optimal;

        float elapsed_time;
        FloatDP f_obj;

        std::cerr << "Testing " << name << "...\n";
        DataStream ds = csv::make_datastream();
        ds.setHeader({"#sequence_" + name, "nVars_" + name, "order_" + name,
                      "correct_" + name, "f_obj_" + name, "time_" + name,
                      "status_" + name});
        ds.open(name + ".csv");

        // FunctionDistribution<EffectiveVectorMultivariateFunction>
        // generator_g(2,2);

        int index = 0;
        unsigned maxNVars = 10;
        int maxOrder = 20;
        for (unsigned nVars = 1; nVars < maxNVars; ++nVars) {
            for (int order = 1; order < maxOrder; ++order, ++index) {
                fprintf(stderr, "%3d %% ]\r",
                        ci((cd(ci(index) * 100) / (ci(maxNVars) * maxOrder))));
                try {
                    single_test_2(seed, nVars, 1, order, f_obj, elapsed_time, f,
                                  g, x_optimal);
                    ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order),
                                        cfdp(1), f_obj, cfdp(elapsed_time), 1});
                } catch (DomainException de) {
                    std::cerr << f << ", " << g << ", " << x_optimal << ", "
                              << f(x_optimal) << ", " << g(x_optimal) << "\n";
                    ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order),
                                        cfdp(-1), f_obj, cfdp(elapsed_time),
                                        0});
                } catch (SQPStatusLeq1Exception sqpsl1) {
                    // std::cerr<<f<<", "<<g<<", "<<x_optimal<<",
                    // "<<f(x_optimal)<<", "<<g(x_optimal)<<"\n";
                    ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order),
                                        sqpsl1.st, f_obj, cfdp(elapsed_time),
                                        sqpsl1.st});
                } catch (std::runtime_error e) {
                    std::cerr << f << ", " << g << ", " << x_optimal << ", "
                              << f(x_optimal) << ", " << g(x_optimal) << "\n";
                    ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order),
                                        cfdp(-15), f_obj, cfdp(elapsed_time),
                                        -15});
                }
            }
        }

        ds.close();
    }

    Void test_3(const std::string &name, unsigned seed)
    {

        // FunctionDistribution<EffectiveVectorMultivariateFunction>
        // generator_g(2,2);
        EffectiveScalarMultivariateFunction f;
        EffectiveVectorMultivariateFunction g;
        FloatBoundsVector x_optimal;

        float elapsed_time;
        FloatDP f_obj;

        std::cerr << "Testing " << name << "...\n";
        DataStream ds = csv::make_datastream();
        ds.setHeader({"#sequence_" + name, "nVars_" + name, "order_" + name,
                      "correct_" + name, "f_obj_" + name, "time_" + name,
                      "status_" + name});
        ds.open(name + ".csv");

        int index = 0;
        unsigned maxNVars = 30;
        for (unsigned nVars = 1; nVars < maxNVars; ++nVars) {
            fprintf(stderr, "%3d %% ]\r",
                    ci((cd(ci(nVars) * 100) / (ci(maxNVars)))));
            try {
                single_test_3(seed, nVars, 1, 0, f_obj, elapsed_time, f, g,
                              x_optimal);
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(1), cfdp(1),
                                    f_obj, cfdp(elapsed_time), 1});
            } catch (SQPStatusLeq1Exception sqpsl1) {
                std::cerr << f << ", " << g << ", " << x_optimal << ", "
                          << f(x_optimal) << ", " << g(x_optimal) << "\n";
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(0),
                                    sqpsl1.st, f_obj, cfdp(elapsed_time),
                                    sqpsl1.st});
            } catch (std::runtime_error e) {
                std::cerr <<"instance: ("<<nVars<<") -> "<<e.what()<<"\n";
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(0),
                                    cfdp(-15), f_obj, cfdp(elapsed_time), -15});
            }
        }

        ds.close();
    }

    Void test_4(const std::string &name, unsigned seed)
    {

        EffectiveScalarMultivariateFunction f;
        EffectiveVectorMultivariateFunction g;
        FloatBoundsVector x_optimal;

        float elapsed_time;
        FloatDP f_obj;

        std::cerr << "Testing " << name << "...\n";
        DataStream ds = csv::make_datastream();
        ds.setHeader({"#sequence_" + name, "nVars_" + name, "order_" + name,
                      "correct_" + name, "f_obj_" + name, "time_" + name,
                      "status_" + name});
        ds.open(name + ".csv");

        int index = 0;
        unsigned maxNVars = 40;
        for (unsigned nVars = 1; nVars < maxNVars; ++nVars, ++index) {

            fprintf(stderr, "%3d %% ]\r",
                    ci((cd(ci(index) * 100) / (ci(maxNVars)))));
            try {
                single_test_4(seed, nVars, f_obj, elapsed_time, f, g,
                              x_optimal);
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
                                    cfdp(1), f_obj, cfdp(elapsed_time), 1});
            } catch (DomainException de) {
                std::cerr << f << ", " << g << ", " << x_optimal << ", "
                          << f(x_optimal) << ", " << g(x_optimal) << "\n";
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
                                    cfdp(-1), f_obj, cfdp(elapsed_time), 0});
            } catch (SQPStatusLeq1Exception sqpsl1) {
                // std::cerr<<f<<", "<<g<<", "<<x_optimal<<",
                // "<<f(x_optimal)<<", "<<g(x_optimal)<<"\n";
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
                                    sqpsl1.st, f_obj, cfdp(elapsed_time),
                                    sqpsl1.st});

            }
            catch(InfeasibleProblemException ipe)
            {
              ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
                cfdp(1), f_obj, cfdp(elapsed_time), 2});
            }
            catch (std::runtime_error e) {
                std::cerr <<"instance: ("<<nVars<<") -> "<<e.what()<<"\n";
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
                                    cfdp(-15), f_obj, cfdp(elapsed_time), -15});
            }
        }

        ds.close();
    }

    Void test_5(const std::string &name, unsigned seed)
    {

      EffectiveScalarMultivariateFunction f;
      EffectiveVectorMultivariateFunction g;
      FloatBoundsVector x_optimal;

      float elapsed_time;
      FloatDP f_obj;

      std::cerr << "Testing " << name << "...\n";
      DataStream ds = csv::make_datastream();
      ds.setHeader({"#sequence_" + name, "nVars_" + name, "order_" + name,
      "correct_" + name, "f_obj_" + name, "time_" + name,
      "status_" + name});
      ds.open(name + ".csv");

      // FunctionDistribution<EffectiveVectorMultivariateFunction>
      // generator_g(2,2);

      int index = 0;
      unsigned maxNVars = 20;
      for (unsigned nVars = 1; nVars < maxNVars; ++nVars, ++index) {

        fprintf(stderr, "%3d %% ]\r",
        ci((cd(ci(index) * 100) / (ci(maxNVars)))));
        try
        {
          single_test_5(seed, nVars, f_obj, elapsed_time, f, g,
            x_optimal);
            ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
              cfdp(1), f_obj, cfdp(elapsed_time), 1});
            }
            catch (DomainException de)
            {
              std::cerr << f << ", " << g << ", " << x_optimal << ", "
              << f(x_optimal) << ", " << g(x_optimal) << "\n";
              ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
                cfdp(-1), f_obj, cfdp(elapsed_time), 0});
              }
              catch (SQPStatusLeq1Exception sqpsl1)
              {
                // std::cerr<<f<<", "<<g<<", "<<x_optimal<<",
                // "<<f(x_optimal)<<", "<<g(x_optimal)<<"\n";
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
                  sqpsl1.st, f_obj, cfdp(elapsed_time),
                  sqpsl1.st});
              }
              catch(InfeasibleProblemException ipe)
              {
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
                  cfdp(1), f_obj, cfdp(elapsed_time), 1});
              }
              catch (std::runtime_error e)
              {
                std::cerr <<"instance: ("<<nVars<<") -> "<<e.what()<<"\n";
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
                  cfdp(-15), f_obj, cfdp(elapsed_time), -15});
                }
              }

                ds.close();
              }

    Void test_X(const std::string &name, unsigned seed)
    {

        EffectiveScalarMultivariateFunction f;
        EffectiveVectorMultivariateFunction g;
        FloatBoundsVector x_optimal;

        float elapsed_time;
        FloatDP f_obj;

        std::cerr << "Testing " << name << "...\n";
        DataStream ds = csv::make_datastream();
        ds.setHeader({"#sequence_" + name, "nVars_" + name, "order_" + name,
                      "correct_" + name, "f_obj_" + name, "time_" + name,
                      "status_" + name});
        ds.open(name + ".csv");

        // FunctionDistribution<EffectiveVectorMultivariateFunction>
        // generator_g(2,2);

        int index = 0;
        unsigned maxNVars = 20;
        for (unsigned nVars = 1; nVars < maxNVars; ++nVars, ++index) {

            fprintf(stderr, "%3d %% ]\r",
                    ci((cd(ci(index) * 100) / (ci(maxNVars)))));
            try {
                single_test_X(seed, nVars, f_obj, elapsed_time, f, g,
                              x_optimal);
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
                                    cfdp(1), f_obj, cfdp(elapsed_time), 1});
            } catch (DomainException de) {
                std::cerr << f << ", " << g << ", " << x_optimal << ", "
                          << f(x_optimal) << ", " << g(x_optimal) << "\n";
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
                                    cfdp(-1), f_obj, cfdp(elapsed_time), 0});
            } catch (SQPStatusLeq1Exception sqpsl1) {
                // std::cerr<<f<<", "<<g<<", "<<x_optimal<<",
                // "<<f(x_optimal)<<", "<<g(x_optimal)<<"\n";
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
                                    sqpsl1.st, f_obj, cfdp(elapsed_time),
                                    sqpsl1.st});
            } catch (std::runtime_error e) {
                std::cerr <<"instance: ("<<nVars<<") -> "<<e.what()<<"\n";
                ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(nVars),
                                    cfdp(-15), f_obj, cfdp(elapsed_time), -15});
            }
        }

        ds.close();
    }


    Void performance_test()
    {
      RandomPolynomialEngine engine(0, 0);
      FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_5>
          generator_f(5);
      FunctionDistribution<EffectiveVectorMultivariateFunction, G_TEST_3>
          generator_g(5);
      ExactBoxType D(5, ExactIntervalType(-10, 11));
      ExactBoxType C(1, ExactIntervalType(-1, 100));

      FloatBoundsVector x_optimal;
      EffectiveScalarMultivariateFunction f = generator_f(engine);
      EffectiveVectorMultivariateFunction g = generator_g(engine);
      unsigned max = 10;

      for(unsigned i=0;i<max;++i)
      {
        fprintf(stderr, "%3d %% ]\r",
                ci((cd(i * 100) / (ci(max)))));
        x_optimal = optimiser->minimise(f, D, g, C);
      }



    }

    Void null_gradient_test()
    {
      EffectiveScalarMultivariateFunction f;
      EffectiveVectorMultivariateFunction g;
      FloatBoundsVector x_optimal;

      float elapsed_time;
      FloatDP f_obj;

      single_test_2(0, 2, 1, 2, f_obj, elapsed_time, f,
                    g, x_optimal);
    }

    Void single_test_1(unsigned seed, unsigned nVars, int max_order,
                       FloatDP &f_obj, float &elapsed_time,
                       EffectiveScalarMultivariateFunction &f,
                       FloatBoundsVector &x_optimal)
    {
        RandomPolynomialEngine engine(seed, max_order);
        FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_1>
            generator(nVars);
        ExactBoxType D(nVars, ExactIntervalType(-10, 10));
        ExactBoxType test(nVars, ExactIntervalType(-10.01, 10.01));
        EffectiveVectorMultivariateFunction g(0u, nVars);
        ExactBoxType C = {};

        f = generator(engine);

        clock_t s_time = clock();
        // run code
        x_optimal = optimiser->minimise(f, D, g, C);
        // End time
        clock_t e_time = clock();
        f_obj = cfdp(f(x_optimal));
        if (!decide(element(x_optimal, test))) {
            throw DomainException();
        }

        elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
        // std::cout << "Elapsed time: " << elapsed_time << " sec\n";
        // std::cout<<"f(x_opt)="<<f(x_optimal)<<",
        // g(x_opt)="<<g(x_optimal)<<"\n";

        // Start time

        // ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
        // FloatDPValue required_accuracy(1e-6);
        // ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
        // end code
    }

    Void single_test_2(unsigned seed, unsigned nVars, unsigned nCons,
                       int max_order, FloatDP &f_obj, float &elapsed_time,
                       EffectiveScalarMultivariateFunction &f,
                       EffectiveVectorMultivariateFunction &g,
                       FloatBoundsVector &x_optimal)
    {
        RandomPolynomialEngine engine(seed, max_order);
        FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_1>
            generator_f(nVars);
        FunctionDistribution<EffectiveVectorMultivariateFunction, G_TEST_1>
            generator_g(nVars, nCons);
        ExactBoxType D(nVars, ExactIntervalType(-10, 10));
        ExactBoxType test(nVars, ExactIntervalType(-10.01, 10.01));
        ExactBoxType C(nCons, ExactIntervalType(-1, 100));
        ExactBoxType testC(nCons, ExactIntervalType(-1.01, 100.01));

        f = generator_f(engine);
        g = generator_g(engine);

        clock_t s_time = clock();
        // run code
        x_optimal = optimiser->minimise(f, D, g, C);
        // End time
        clock_t e_time = clock();
        f_obj = cfdp(f(x_optimal));
        if (!decide(element(x_optimal, test)) ||
            !decide(element(g(x_optimal), testC))) {
            // throw std::runtime_error("Predicato sbagliato\n");
            throw DomainException();
        }

        elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
        // std::cout << "Elapsed time: " << elapsed_time << " sec\n";
        // std::cout<<"f(x_opt)="<<f(x_optimal)<<",
        // g(x_opt)="<<g(x_optimal)<<"\n";

        // Start time

        // ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
        // FloatDPValue required_accuracy(1e-6);
        // ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
        // end code
    }

    Void single_test_3(unsigned seed, unsigned nVars, unsigned nCons,
                       int max_order, FloatDP &f_obj, float &elapsed_time,
                       EffectiveScalarMultivariateFunction &f,
                       EffectiveVectorMultivariateFunction &g,
                       FloatBoundsVector &x_optimal)
    {
        RandomPolynomialEngine engine(seed);
        FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_2>
            generator_f(nVars);
        FunctionDistribution<EffectiveVectorMultivariateFunction, G_TEST_2>
            generator_g(nVars);
        ExactBoxType D(nVars, ExactIntervalType(0, 100));
        ExactBoxType testD(nVars, ExactIntervalType(-0.01, 100.01));
        ExactBoxType C = ExactBoxType{{-1000, 100}, {10, 1000}};
        ExactBoxType testC = ExactBoxType{{-10000.01, 100.01}, {9.99, 1000.01}};

        f = generator_f(engine);
        g = generator_g(engine);

        clock_t s_time = clock();
        // run code
        x_optimal = optimiser->minimise(f, D, g, C);
        // End time
        clock_t e_time = clock();
        f_obj = cfdp(f(x_optimal));
        if (!decide(element(x_optimal, testD)) ||
            !decide(element(g(x_optimal), testC))) {
            throw std::runtime_error("Predicato sbagliato\n");
            // throw DomainException();
        }

        elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
    }

    Void single_test_4(unsigned seed, unsigned nVars, FloatDP &f_obj,
                       float &elapsed_time,
                       EffectiveScalarMultivariateFunction &f,
                       EffectiveVectorMultivariateFunction &g,
                       FloatBoundsVector &x_optimal)
    {
        RandomPolynomialEngine engine(seed, 0);
        FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_4>
            generator_f(nVars);
        FunctionDistribution<EffectiveVectorMultivariateFunction, G_TEST_4>
            generator_g(nVars);
        ExactBoxType D(nVars, ExactIntervalType(0.1, +inf));
        ExactBoxType test(nVars, ExactIntervalType(0.099, +inf));
        ExactBoxType C(nVars, ExactIntervalType(1, 10));
        ExactBoxType testC(nVars, ExactIntervalType(0.99, 10.01));

        f = generator_f(engine);
        g = generator_g(engine);

        clock_t s_time = clock();
        // run code
        x_optimal = optimiser->minimise(f, D, g, C);
        // End time
        clock_t e_time = clock();
        f_obj = cfdp(f(x_optimal));
        if (!decide(element(x_optimal, test)) ||
            !decide(element(g(x_optimal), testC))) {
            // throw std::runtime_error("Predicato sbagliato\n");
            throw DomainException();
        }

        elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
    }

    Void single_test_5(unsigned seed, unsigned nVars, FloatDP &f_obj,
      float &elapsed_time,
      EffectiveScalarMultivariateFunction &f,
      EffectiveVectorMultivariateFunction &g,
      FloatBoundsVector &x_optimal)
      {
        RandomPolynomialEngine engine(seed, 0);
        FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_5>
        generator_f(nVars);
        FunctionDistribution<EffectiveVectorMultivariateFunction, G_TEST_3>
        generator_g(nVars);
        ExactBoxType D(nVars, ExactIntervalType(-10, 11));
        ExactBoxType test(nVars, ExactIntervalType(-10.01, 11.01));
        ExactBoxType C(1, ExactIntervalType(-1, 100));
        ExactBoxType testC(1, ExactIntervalType(-1.01, 100.01));

        f = generator_f(engine);
        g = generator_g(engine);

        clock_t s_time = clock();
        // run code
        x_optimal = optimiser->minimise(f, D, g, C);
        // End time
        clock_t e_time = clock();
        f_obj = cfdp(f(x_optimal));
        if (!decide(element(x_optimal, test)) ||
        !decide(element(g(x_optimal), testC))) {
          // throw std::runtime_error("Predicato sbagliato\n");
          throw DomainException();
        }

        elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
        // std::cout << "Elapsed time: " << elapsed_time << " sec\n";
        // std::cout<<"f(x_opt)="<<f(x_optimal)<<",
        // g(x_opt)="<<g(x_optimal)<<"\n";

        // Start time

        // ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
        // FloatDPValue required_accuracy(1e-6);
        // ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
        // end code
      }

    Void single_test_X(unsigned seed, unsigned nVars, FloatDP &f_obj,
                       float &elapsed_time,
                       EffectiveScalarMultivariateFunction &f,
                       EffectiveVectorMultivariateFunction &g,
                       FloatBoundsVector &x_optimal)
    {
        RandomPolynomialEngine engine(seed, 0);
        FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_1>
            generator_f(nVars);
        FunctionDistribution<EffectiveVectorMultivariateFunction, G_TEST_3>
            generator_g(nVars);
        ExactBoxType D(nVars, ExactIntervalType(-10, 11));
        ExactBoxType test(nVars, ExactIntervalType(-10.01, 11.01));
        ExactBoxType C(1, ExactIntervalType(-1, 100));
        ExactBoxType testC(1, ExactIntervalType(-1.01, 100.01));

        f = generator_f(engine);
        g = generator_g(engine);

        clock_t s_time = clock();
        // run code
        x_optimal = optimiser->minimise(f, D, g, C);
        // End time
        clock_t e_time = clock();
        f_obj = cfdp(f(x_optimal));
        if (!decide(element(x_optimal, test)) ||
            !decide(element(g(x_optimal), testC))) {
            // throw std::runtime_error("Predicato sbagliato\n");
            throw DomainException();
        }

        elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
        // std::cout << "Elapsed time: " << elapsed_time << " sec\n";
        // std::cout<<"f(x_opt)="<<f(x_optimal)<<",
        // g(x_opt)="<<g(x_optimal)<<"\n";

        // Start time

        // ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
        // FloatDPValue required_accuracy(1e-6);
        // ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
        // end code
    }

};

Int main(Int argc, const char *argv[])
{
    Nat optimiser_verbosity = get_verbosity(argc, argv);

    // unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    unsigned seed = 0;

    
    NonlinearSQPOptimiser nlsqp;
    nlsqp.verbosity = optimiser_verbosity;
    TestOptimiser(nlsqp).test("nlsqp", seed);
    // return ARIADNE_TEST_FAILURES;

    NonlinearInteriorPointOptimiser nlipm;
    nlipm.verbosity = optimiser_verbosity;
    TestOptimiser(nlipm).test("nlipm", seed);
    // return ARIADNE_TEST_FAILURES;


    NonlinearInfeasibleInteriorPointOptimiser nliipm;
    nliipm.verbosity = optimiser_verbosity;
    TestOptimiser(nliipm).test("nliipm", seed);
    // return ARIADNE_TEST_FAILURES;


    NonlinearMixedOptimiser nlmop;
    nlmop.verbosity = optimiser_verbosity;
    TestOptimiser(nlmop).test("nlmop",seed);
    return ARIADNE_TEST_FAILURES;
}
