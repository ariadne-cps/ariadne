/***************************************************************************
 *            profile_taylor_variable.cc
 *
 *  Copyright 2008  Pieter Collins
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <string>
using std::string;

#include <boost/timer.hpp>
#include <boost/progress.hpp>

#include "taylor_model.h"
using namespace Ariadne;

TaylorModel sum(const Vector<TaylorModel>& v) {
    return v[0]+v[1];
}

TaylorModel prod(const Vector<TaylorModel>& v) {
    return v[0]*v[1];
}


TaylorModel exp_cos(const Vector<TaylorModel>& v) {
    return exp(v[0])*cos(v[1]);
}

TaylorModel sigmoid(const Vector<TaylorModel>& v) {
    const double a=10;
    return exp(-v[0]/a);
}

typedef TaylorModel(*TaylorFunctionPtr)(const Vector<TaylorModel>&);


void profile(uint ntries, string name, TaylorFunctionPtr fn, const Vector<TaylorModel>& args)
{
    TaylorModel res=fn(args);
    std::cerr<< "\n" << name << "(" << args << ")=\n  " << res << "\n\n";

    boost::timer tm;

    tm.restart();
    for(uint i=0; i!=ntries; ++i) {
        res=fn(args);
    }

    double total_time = tm.elapsed();
    double average_time_in_microseconds = 1000000*(total_time/ntries);
    double error = res.error();
    unsigned int size = res.number_of_nonzeros();

    std::cout << std::setw(20) << std::left << name << std::right
              << std::setw(9) << std::fixed << std::setprecision(1) << average_time_in_microseconds << " "
              << std::setw(12) << std::scientific << std::setprecision(4) << error
              << std::setw(8) << size
              << std::endl;
}

int main(int argc, const char* argv[]) {
    uint ntries=100;
    if(argc>1) { ntries=atoi(argv[1]); }
    std::cerr<<ntries<<std::endl;

    Vector<Float> c(2, 1.0,2.0);

    Vector<TaylorModel> v(2,TaylorModel(2));
    v[0]=TaylorModel::variable(2,0);
    v[1]=TaylorModel::constant(2,1.0);

    Vector<TaylorModel> x(2,TaylorModel(3));
    int i=0;
    std::cerr<<x<<"\n";
    for(MultiIndex a(3); a.degree()<=5; ++a) {
        if(i%7<4) { x[0][a]=1/(1.0+i); }
        if(i%3<2) { x[1][a]=1/(2.0+i); }
        std::cerr<<x<<"\n";
    }
    std::cerr<<"v="<<v<<"\nx="<<x<<"\n";

    std::cout << std::setw(20) << std::left << "name" << std::right
              << std::setw(10) << "time(us)"
              << std::setw(12) << "error"
              << std::setw(8) << "size"
              << std::endl;

    profile(ntries*100,"sum-01",sum,x);
    profile(ntries*10,"prod-01",prod,x);
    profile(ntries,"exp_cos-01",exp_cos,v);
    profile(ntries,"sigmoid-01",sigmoid,v);
}
