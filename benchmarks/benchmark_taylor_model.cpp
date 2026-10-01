/***************************************************************************
 *            benchmark_taylor_model.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "function/taylor_model.hpp"
#include "utility/stopwatch.hpp"

using namespace Ariadne;

template<class Run>
Void benchmark(Nat repetitions, const char* name, const Run& run)
{
    auto result=run();

    Stopwatch<Microseconds> stopwatch;
    for(Nat i=0u; i!=repetitions; ++i) {
        result=run();
    }
    stopwatch.click();

    const double average_time_in_microseconds=
        static_cast<double>(stopwatch.duration().count())/
        static_cast<double>(repetitions);

    std::cout << std::setw(20) << std::left << name << std::right
              << std::setw(10) << std::fixed << std::setprecision(2)
              << average_time_in_microseconds << " "
              << std::setw(12) << std::scientific << std::setprecision(4)
              << result.error()
              << std::setw(8) << result.number_of_terms()
              << std::endl;
}

Int main(Int argc, const char* argv[])
{
    Nat repetitions=20u;
    if(argc>1) {
        repetitions=static_cast<Nat>(std::strtoul(argv[1],nullptr,10));
        if(repetitions==0u) {
            repetitions=1u;
        }
    }

    DP pr;
    Sweeper<FloatDP> trivial_sweeper{TrivialSweeper<FloatDP>(pr)};
    Sweeper<FloatDP> threshold_sweeper{ThresholdSweeper<FloatDP>(pr,1e-5)};

    using TM=ValidatedTaylorModelDP;

    TM w(3u,threshold_sweeper);
    SizeType i=0u;
    for(MultiIndex a(3u); a.degree()<=9u; ++a) {
        const double di=static_cast<double>(i);
        if(i%7u<3u) {
            w.expansion().append(a,1.0/(1.0+di*di*di*di*di));
        } else if(i%7u<4u) {
            w.expansion().append(a,1.0/(1.0+di));
        }
        ++i;
    }

    const FloatDPBounds interval(0.33_x,0.49_x,pr);
    const FloatDP constant(0.41_x,pr);

    TM x(3u,threshold_sweeper);
    TM y(3u,threshold_sweeper);

    i=0u;
    for(MultiIndex a(3u); a.degree()<=5u; ++a) {
        const double di=static_cast<double>(i);
        if(i%7u<4u) {
            x.expansion().append(a,1.0/(1.0+di));
        }
        if(i%3u<2u) {
            y.expansion().append(a,1.0/(2.0+di));
        }
        ++i;
    }

    TM z(3u,trivial_sweeper);
    z.expansion().append(MultiIndex({0u,0u,0u}),1.0);
    z.expansion().append(MultiIndex({1u,0u,0u}),0.5);
    z.expansion().append(MultiIndex({0u,1u,0u}),-0.25);
    z.expansion().append(MultiIndex({0u,0u,1u}),0.625);

    const TM variable=TM::coordinate(2u,0u,trivial_sweeper);
    const TM one=TM::constant(2u,1,trivial_sweeper);

    std::cout << std::setw(20) << std::left << "name" << std::right
              << std::setw(11) << "time(us)"
              << std::setw(12) << "error"
              << std::setw(8) << "size"
              << std::endl;

    benchmark(repetitions*10000u,"sweep-02",[&]() {
        TM result=w;
        result.sweep();
        return result;
    });
    benchmark(repetitions*10000u,"copy-02",[&]() {
        return TM(w);
    });
    benchmark(repetitions*10000u,"iadd-noinsert-02",[&]() {
        TM result=x;
        result+=interval;
        return result;
    });
    benchmark(repetitions*10000u,"iscal-02",[&]() {
        TM result=x;
        result*=interval;
        return result;
    });
    benchmark(repetitions*10000u,"fscal-02",[&]() {
        TM result=x;
        result*=constant;
        return result;
    });
    benchmark(repetitions*1000u,"isum-02",[&]() {
        TM result=x;
        result+=y;
        return result;
    });
    benchmark(repetitions*1000u,"sum-02",[&]() {
        return x+y;
    });
    benchmark(repetitions*10u,"prod-02",[&]() {
        return x*y;
    });
    benchmark(repetitions,"exp-02",[&]() {
        return exp(z);
    });
    benchmark(repetitions,"exp_cos-01",[&]() {
        return exp(variable)*cos(one);
    });
    benchmark(repetitions,"sigmoid-01",[&]() {
        return exp(-variable/10);
    });

    return 0;
}
