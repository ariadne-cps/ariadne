/***************************************************************************
 *            test_linear_programming.cpp
 *
 *  Copyright  2009  Pieter Collins
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

#include <iostream>
#include <fstream>
#include <chrono>

#include "config.hpp"
#include "../test.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "solvers/quadratic_programming.hpp"

using namespace std;
using namespace Ariadne;

class TestASMQPSolver
{
    ASMQPSolver *optimiser;

  public:
    TestASMQPSolver(ASMQPSolver &o) : optimiser(&o) {}

    Void test()
    {
        // ARIADNE_TEST_CALL(test_optimization1()); // [1.4;1.7], 4.3179 (octave)
        // ARIADNE_TEST_CALL(test_optimization2()); // [0.5;1.0], 9.25 (octave)
        // ARIADNE_TEST_CALL(test_optimization3()); // [0.76250;0.47500], 4.31787 (octave)
        // ARIADNE_TEST_CALL(test_optimization4()); // [-6], -46 (octave)
        // ARIADNE_TEST_CALL(test_optimization5()); // [-50], ... (octave)
        // ARIADNE_TEST_CALL(test_optimization6()); // [-171.25,40.2975], -16567.8 (octave)
        // ARIADNE_TEST_CALL(test_optimization7()); // [0] 0 (octave)
        // ARIADNE_TEST_CALL(test_optimization8()); // [1.42857,4.28571] -8.16327 (octave)
        // ARIADNE_TEST_CALL(test_optimization9()); // [2,2] -14 (octave)
        // ARIADNE_TEST_CALL(test_optimization10());// [2,2] (octave)
        ARIADNE_TEST_CALL(test_optimizationtmp());
    }

    Void test_optimization1()
    {
        FloatMatrix Q = {{2, 0}, {0, 2}};
        RawFloatVector d = {-2, -5};
        FloatMatrix A = {{1, -2}, {-1, -2}, {-1, 2}, {1, 0}, {0, 1}};

        RawFloatVector xl = {0.0, -inf};
        RawFloatVector xu = {4, +inf};

        RawFloatVector wl = {-2, -6, -2, 0, 0};
        RawFloatVector wu = {+inf, +inf, +inf, +inf, +inf};
        RawFloatVector x0 = {0, 0};

        int status = 0;

        std::cerr<<"Testing\n min 0.5 x'"<<Q<<"x + "<<d<<"'x\n\ts.t. "<<wl<<" <= "<<A<<" <= "<<wu<<"\n\t    "<<xl<<" <= x <= "<<xu<<"\n";

        clock_t start_time = clock();
        ARIADNE_TEST_PRINT(optimiser->minimise(Q, d, xl,xu,A, wl,wu,status,x0));
        clock_t end_time = clock();

        float diff = static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC;
        std::cerr << "Elapsed time: " << diff << "\n";
    }
    Void test_optimization2()
    {
        FloatMatrix Q = {{6, 2}, {2, 2}};
        RawFloatVector d = {1, 6};
        FloatMatrix A = {{2, 3}, {1, 0}, {0, 1}};
        RawFloatVector wl = {4, 0, 0};
        RawFloatVector wu = {+inf, +inf, +inf};

        RawFloatVector xl = {-inf, -inf};
        RawFloatVector xu = {+inf, +inf};
        RawFloatVector x0 = {1,1};


        int status = 0;
        std::cerr<<"Testing\n min 0.5 x'"<<Q<<"x + "<<d<<"'x\n\ts.t. "<<wl<<" <= "<<A<<" <= "<<wu<<"\n\t    "<<xl<<" <= x <= "<<xu<<"\n";

        clock_t start_time = clock();
        ARIADNE_TEST_PRINT(optimiser->minimise(Q, d, xl,xu,A, wl,wu,status,x0));
        clock_t end_time = clock();

        float diff = static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC;
        std::cerr << "Elapsed time: " << diff << "\n";
    }
    Void test_optimization3()
    {
        FloatMatrix Q = {{8, 2}, {2, 10}};
        RawFloatVector d = {1.5, -2};
        FloatMatrix A = {{2, 1}, {1, -2}, {1, 0}, {-1, 0}, {0, 1}};
        RawFloatVector wl = {2, -6, 0, -20, 0};
        RawFloatVector wu = {+inf, +inf, +inf, +inf, +inf};

        RawFloatVector xl = {-inf, -inf};
        RawFloatVector xu = {+inf, +inf};
        RawFloatVector x0 = {1,0};

        int status = 0;
        std::cerr<<"Testing\n min 0.5 x'"<<Q<<"x + "<<d<<"'x\n\ts.t. "<<wl<<" <= "<<A<<" <= "<<wu<<"\n\t    "<<xl<<" <= x <= "<<xu<<"\n";


        clock_t start_time = clock();
        ARIADNE_TEST_PRINT(optimiser->minimise(Q, d, xl,xu,A, wl,wu,status,x0));
        clock_t end_time = clock();

        float diff = static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC;
        std::cerr << "Elapsed time: " << diff << "\n";
    }
    Void test_optimization4()
    {
        FloatMatrix Q = {{-2}};
        RawFloatVector d = {0};
        FloatMatrix A = {{1}};

        RawFloatVector xl = {-6.0};
        RawFloatVector xu = {-5.0};

        RawFloatVector wl = {-8};
        RawFloatVector wu = {+inf};
        RawFloatVector x0 = {-5};

        int status = 0;
        std::cerr<<"Testing\n min 0.5 x'"<<Q<<"x + "<<d<<"'x\n\ts.t. "<<wl<<" <= "<<A<<" <= "<<wu<<"\n\t    "<<xl<<" <= x <= "<<xu<<"\n";

        clock_t start_time = clock();
        ARIADNE_TEST_PRINT(optimiser->minimise(Q, d, xl,xu,A, wl,wu,status,x0));
        clock_t end_time = clock();

        float diff = static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC;
        std::cerr << "Elapsed time: " << diff << "\n";
    }
    Void test_optimization5()
    {
        FloatMatrix Q = {{1}};
        RawFloatVector d = {5.18471e+21};
        FloatMatrix A = {{1}};

        RawFloatVector xl = {-50.0};
        RawFloatVector xu = {100.0};

        RawFloatVector wl = {-50};
        RawFloatVector wu = {50};
        RawFloatVector x0 = {0};


        int status = 0;
        std::cerr<<"Testing\n min 0.5 x'"<<Q<<"x + "<<d<<"'x\n\ts.t. "<<wl<<" <= "<<A<<" <= "<<wu<<"\n\t    "<<xl<<" <= x <= "<<xu<<"\n";

        clock_t start_time = clock();
        ARIADNE_TEST_PRINT(optimiser->minimise(Q, d, xl,xu,A, wl,wu,status,x0));
        clock_t end_time = clock();

        float diff = static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC;
        std::cerr << "Elapsed time: " << diff << "\n";
    }
    Void test_optimization6() //interesting case where p_y is almoust 0, and so Z'HZ has negative eigenvalues! algorithm cannot converge (see 16.5 of Nocedal (Numerical optimizaiont 2006)).
    {
        FloatMatrix Q = {{1,0},{0,1}};
        RawFloatVector d = {244.692,244.692};
        FloatMatrix A = {{0.181818, 0.705544},{0, 0.70554}};

        RawFloatVector xl = {-inf,-inf};
        RawFloatVector xu = {inf, inf};

        RawFloatVector wl = {-2.70475,-0.70867};
        RawFloatVector wu = {inf,inf};
        RawFloatVector x0 = {5.5,0};

        int status = 0;
        std::cerr<<"Testing\n min 0.5 x'"<<Q<<"x + "<<d<<"'x\n\ts.t. "<<wl<<" <= "<<A<<" <= "<<wu<<"\n\t    "<<xl<<" <= x <= "<<xu<<"\n";

        clock_t start_time = clock();
        ARIADNE_TEST_PRINT(optimiser->minimise(Q, d, xl,xu,A, wl,wu,status,x0));
        clock_t end_time = clock();
        std::cout<<"Ended with status "<<status<<"\n";

        float diff = static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC;
        std::cerr << "Elapsed time: " << diff << "\n";
    }
    Void test_optimization7()
    {
        FloatMatrix Q = {{0.008}};
        RawFloatVector d = {1};
        FloatMatrix A = {{1}};

        RawFloatVector xl = {-10.0};
        RawFloatVector xu = {10};

        RawFloatVector wl = {0};
        RawFloatVector wu = {10};
        RawFloatVector x0 = {5};

        int status = 0;
        std::cerr<<"Testing\n min 0.5 x'"<<Q<<"x + "<<d<<"'x\n\ts.t. "<<wl<<" <= "<<A<<" <= "<<wu<<"\n\t    "<<xl<<" <= x <= "<<xu<<"\n";

        clock_t start_time = clock();
        ARIADNE_TEST_PRINT(optimiser->minimise(Q, d, xl,xu,A, wl,wu,status,x0));
        clock_t end_time = clock();

        float diff = static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC;
        std::cerr << "Elapsed time: " << diff << "\n";
    }
    Void test_optimization8()
    {
        FloatMatrix Q = {{1,0},{0,-1}};
        RawFloatVector d = {0,0};
        FloatMatrix A = {{1,2},{5,-4}};

        RawFloatVector xl = {-10.0,0.0};
        RawFloatVector xu = {+inf,+inf};

        RawFloatVector wl = {2,-10};
        RawFloatVector wu = {10,10};
        RawFloatVector x0 = {2,0};

        int status = 0;
        std::cerr<<"Testing\n min 0.5 x'"<<Q<<"x + "<<d<<"'x\n\ts.t. "<<wl<<" <= "<<A<<" <= "<<wu<<"\n\t    "<<xl<<" <= x <= "<<xu<<"\n";

        clock_t start_time = clock();
        ARIADNE_TEST_PRINT(optimiser->minimise(Q, d, xl,xu,A, wl,wu,status,x0));
        clock_t end_time = clock();

        float diff = static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC;
        std::cerr << "Elapsed time: " << diff << "\n";
    }
    Void test_optimization9()
    {
      FloatMatrix Q = {{1, 0}, {0, -1}};
      RawFloatVector d = {-2, -5};
      FloatMatrix A = {{1, -1}, {-1, -2}, {-3, 1}};

      RawFloatVector xl = {2.0, -inf};
      RawFloatVector xu = {4, +inf};

      RawFloatVector wl = {0, -6, -6};
      RawFloatVector wu = {10, +inf, +inf};
      RawFloatVector x0 = {0, 0};

      int status = 0;
      std::cerr<<"Testing\n min 0.5 x'"<<Q<<"x + "<<d<<"'x\n\ts.t. "<<wl<<" <= "<<A<<" <= "<<wu<<"\n\t    "<<xl<<" <= x <= "<<xu<<"\n";

      clock_t start_time = clock();
      ARIADNE_TEST_PRINT(optimiser->minimise(Q,d,xl,xu,A,wl,wu,status,x0));
      clock_t end_time = clock();

      float diff = static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC;
      std::cerr << "Elapsed time: " << diff << "\n";
    }
    Void test_optimization10()
    {
        FloatMatrix Q = {{2, 0}, {0, 2}};
        RawFloatVector d = {-2, -5};
        FloatMatrix A = {{1, -1}, {-1, -2}, {-3, 1}};

        RawFloatVector xl = {2.0, -inf};
        RawFloatVector xu = {4, +inf};

        RawFloatVector wl = {0, -6, -6};
        RawFloatVector wu = {10, +inf, +inf};
        RawFloatVector x0 = {1, 0};

        int status = 0;
        std::cerr<<"Testing\n min 0.5 x'"<<Q<<"x + "<<d<<"'x\n\ts.t. "<<wl<<" <= "<<A<<" <= "<<wu<<"\n\t    "<<xl<<" <= x <= "<<xu<<"\n";

        clock_t start_time = clock();
        ARIADNE_TEST_PRINT(optimiser->minimise(Q,d,xl,xu,A,wl,wu,status,x0));
        clock_t end_time = clock();

        float diff = static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC;
        std::cerr << "Elapsed time: " << diff << "\n";
    }
    Void test_optimizationtmp()
    {
        FloatMatrix Q={{601784,6241.68,304271,3347.92,227713,15991.2,-40104.4},{6241.68,2053.05,3595.1,-37.4783,1665.66,231.691,180.486},{304271,3595.1,157071,1128.21,117486,5715.61,-18073.7},{3347.92,-37.4783,1128.21,3803.17,-66.4693,30.1081,-1633.77},{227713,1665.66,117486,-66.4693,89862.6,5235.32,-15490.9},{15991.2,231.691,5715.61,30.1081,5235.32,4078.46,-2639.17},{-40104.4,180.486,-18073.7,-1633.77,-15490.9,-2639.17,10218.7}};
        RawFloatVector d={-5261.66,-9.62588,32265,141.593,1434.62,-388.882,-3607.15};

        FloatMatrix A = {{37.6024,1.41572,-60.6355,-1.49069,4.27012,1.32972,17.5082},{-37.6024,-1.41572,60.6355,1.49069,-4.27012,-1.32972,-17.5082},{1,0,0,0,0,0,0},{0,1,0,0,0,0,0},{0,0,1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,0,1,0,0},{0,0,0,0,0,1,0},{0,0,0,0,0,0,1},{-1,0,0,0,0,0,0},{0,-1,0,0,0,0,0},{0,0,-1,0,0,0,0},{0,0,0,-1,0,0,0},{0,0,0,0,-1,0,0},{0,0,0,0,0,-1,0},{0,0,0,0,0,0,-1}};

        RawFloatVector xl = Vector<FloatDP>();// = {-inf, -inf,-inf,-inf,-inf,-inf,-inf};
        RawFloatVector xu = Vector<FloatDP>();// = {+inf, +inf,+inf,+inf,+inf,+inf,+inf};
        RawFloatVector wl={-194.62,93.6202,-7.27182,-4.68706,-0.254834,-4.47946,-4.89095,-4.68092,-5.83653,-12.7282,-15.3129,-19.7452,-15.5205,-15.1091,-15.3191,-14.1635};
        RawFloatVector wu(wl.size(),+inf);// = {+inf, +inf, +inf,+inf,+inf};
        RawFloatVector x0={-2.72818,-5.31294,-9.74517,-5.52054,-5.10905,-5.31908,-4.16347};

        int status = 0;
        std::cerr<<"Testing\n min 0.5 x'"<<Q<<"x + "<<d<<"'x\n\ts.t. "<<wl<<" <= "<<A<<" <= "<<wu<<"\n\t    "<<xl<<" <= x <= "<<xu<<"\n";

        clock_t start_time = clock();
        ARIADNE_TEST_PRINT(optimiser->minimise(Q,d,xl,xu,A,wl,wu,status,x0));
        clock_t end_time = clock();

        float diff = static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC;
        std::cerr << "Elapsed time: " << diff << "\n";
    }
};

Int main(Int argc, const char *argv[])
{
    auto verbosity = get_verbosity(argc, argv);

    ASMQPSolver asmqp_optimiser;
    asmqp_optimiser.verbosity = verbosity;
    TestASMQPSolver(asmqp_optimiser).test();

    std::cerr << "INCOMPLETE\n";
    return ARIADNE_TEST_FAILURES;
}
