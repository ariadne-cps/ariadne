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

TaylorVariable exp_cos(const Vector<TaylorVariable>& v) {
    return exp(v[0])*cos(v[1]);
}

TaylorVariable sigmoid(const Vector<TaylorVariable>& v) {
    const double a=10;
    return exp(-v[0]/a);
}

typedef TaylorVariable(*TaylorFunctionPtr)(const Vector<TaylorVariable>&);


void profile(uint ntries, string name, uint nargs, TaylorFunctionPtr fn) {
    
    TaylorVariable x(2,1, 0.0,1.0,0.0, 0.0);
    TaylorVariable y(2,1, 1.0,0.0,0.0, 0.0);
    TaylorVariable z(2);

    Vector<TaylorVariable> a(nargs); a[0]=x; a[1]=y;

    boost::timer tm; double t=0;
 
    tm.restart();
    for(uint i=0; i!=ntries; ++i) {
        z=fn(a);
    }
    t=tm.elapsed();
    std::cout << name << ":\n"
              << "  time = "<<std::setprecision(5)<<1000000*(t/ntries)<<"us\n"
              << "  size = "<<z.nnz()<<"\n"
              << "  error = "<<z.error()<<"\n"
              << std::endl;
}

int main(int argc, const char* argv[]) {
    profile(1000,"exp_cos",2,exp_cos);
    profile(1000,"sigmoid",2,sigmoid);
}
