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

#include "taylor_variable.h"
using namespace Ariadne;

TaylorVariable sum(const Vector<TaylorVariable>& v) {
    return v[0]+v[1];
}

TaylorVariable prod(const Vector<TaylorVariable>& v) {
    return v[0]*v[1];
}


TaylorVariable exp_cos(const Vector<TaylorVariable>& v) {
    return exp(v[0])*cos(v[1]);
}

TaylorVariable sigmoid(const Vector<TaylorVariable>& v) {
    const double a=10;
    return exp(-v[0]/a);
}

typedef TaylorVariable(*TaylorFunctionPtr)(const Vector<TaylorVariable>&);


void profile(uint ntries, string name, TaylorFunctionPtr fn, const Vector<TaylorVariable>& args)
{
    TaylorVariable res=fn(args);
    std::cerr<< "\n" << name << "(" << args << ")=\n  " << res << "\n\n";

    boost::timer tm; double t=0;

    tm.restart();
    for(uint i=0; i!=ntries; ++i) {
        res=fn(args);
    }
    t=tm.elapsed();
    std::cout << name << ":\n"
              << "  time = "<<std::setprecision(5)<<1000000*(t/ntries)<<"us\n"
              << "  size = "<<res.number_of_nonzeros()<<"\n"
              << "  error = "<<res.error()<<"\n"
              << std::endl;
}

int main(int argc, const char* argv[]) {
    Vector<Interval> d(2,Interval(-1,+1));
    Vector<Float> c(2, 1.0,2.0);

    Vector<TaylorVariable> v(2,2);
    v[0]=TaylorVariable::variable(d,0);
    v[1]=TaylorVariable::constant(d,1.0);

    Vector<TaylorVariable> x(2,2);
    x[0]=TaylorVariable(d,Expansion<Float>(2,3, 1.0,2.0,0.0,4.0,0.0,6.0,0.0,8.0,9.0,10.0, 0.25));
    x[1]=TaylorVariable(d,Expansion<Float>(2,3, 1.0,0.0,3.0,4.0,0.0,6.0,7.0,8.0,0.0,10.0, 0.5));
    std::cerr<<"v="<<v<<"\nx="<<x<<"\n";

    profile(100000,"sum",sum,x);
    profile(10000,"prod",prod,x);
    profile(1000,"exp_cos",exp_cos,v);
    profile(1000,"sigmoid",sigmoid,v);
}
