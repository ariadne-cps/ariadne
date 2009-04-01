/***************************************************************************
 *            profile_taylor_calculus.cc
 *
 *  Copyright 2009  Pieter Collins
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

#include "numeric.h"
#include "taylor_set.h"
#include "taylor_variable.h"
#include "taylor_function.h"
#include "taylor_calculus.h"
#include "function.h"
#include "box.h"
using namespace Ariadne;

typedef Polynomial<Float> FPolynomial;
typedef Vector<Float> FVector;
typedef Vector<Interval> IVector;
FPolynomial p(int as, int j) { return FPolynomial::variable(as,j); }
FVector e(int rs, int i) { return FVector::unit(rs,i); }

Float a=1.5; Float b=0.375;
FPolynomial x=p(2,0); FPolynomial y=p(2,1);
FVector e0=e(2,0); FVector e1=e(2,1);
typedef Vector< Polynomial<Float> > FPolyVec;

PolynomialFunction henon_map =  FPolyVec((a+x*x+b*y)*e0+x*e1);
PolynomialFunction spiral_vf = FPolyVec((1.0+0.5*x-0.75*y)*e0+(0.75*x+0.25*y)*e1);
PolynomialFunction affine_pred =  FPolyVec( (x+y-0.25)*e0);
TaylorFunction spiral_flow_model;
TaylorVariable affine_pred_model;

Box unit_box(2, -1,+1,-1,+1);
TaylorSet unit_box_model=unit_box;

//TaylorFunction henon(Box(2, -1,+1,-1,+1),henon_poly);

struct ProfileReset {
    ProfileReset(TaylorCalculus* c_, FunctionInterface& f_, TaylorSet& s_, int n_=1) : c(c_), f(f_), s(s_), n(n_) { }
    TaylorCalculus* c; FunctionInterface& f; TaylorSet s; int n; typedef TaylorSet Result;
    TaylorSet operator()() const { TaylorSet r=c->reset_step(f,s); for(int i=1; i<n; ++i) { r=c->reset_step(f,s); } return r; }
};

struct ProfileFlow {
    ProfileFlow(TaylorCalculus* c_, FunctionInterface& f_, IVector d_, Float h_, IVector b_) : c(c_), f(f_), d(d_), h(h_), b(b_) { }
    TaylorCalculus* c; FunctionInterface& f; Vector<Interval> d; Float h; Vector<Interval> b; typedef TaylorFunction Result;
    TaylorFunction operator()() const { return c->flow_model(f,d,h,b); }
};

struct ProfileCrossing {
    ProfileCrossing(TaylorCalculus* c_, TaylorVariable& g_, TaylorFunction& f_, TaylorSet s_) : c(c_), g(g_), f(f_), s(s_) { }
    TaylorCalculus* c; TaylorVariable& g; TaylorFunction f; TaylorSet s; typedef TaylorModel Result;
    TaylorModel operator()() const { return c->crossing_time(g,f,s); }
};


template<class T>
double error(const T& f) {
    Float max_error=0.0;
    for(uint i=0; i!=f.size(); ++i) {
        max_error=std::max(max_error,f.models()[i].error());
    }
    return max_error;
}

template<class T>
unsigned int number_of_nonzeros(const T& f) {
    unsigned int nnz=0;
    for(uint i=0; i!=f.size(); ++i) {
        nnz+=f.models()[i].number_of_nonzeros();
    }
    return nnz;
}

template<> double error(const TaylorVariable& t) { return t.model().error(); }
template<> unsigned int number_of_nonzeros(const TaylorVariable& t) { return t.model().number_of_nonzeros(); }
template<> double error(const TaylorModel& t) { return t.error(); }
template<> unsigned int number_of_nonzeros(const TaylorModel& t) { return t.number_of_nonzeros(); }

template<class Test>
void
profile(const char* name, const Test& test, unsigned int tries)
{

    boost::timer tm; double t=0;

    typename Test::Result res=test();
    //TaylorFunction res=static_cast<TaylorFunction>(test(calc,args1,args2));
    unsigned int nnz=number_of_nonzeros(res);
    double err=error(res);

    tm.restart();
    for(uint i=0; i!=tries; ++i) {
        test();
    }
    t=tm.elapsed();
    std::cout << std::setw(17) << std::left << name << std::right
              << std::setw(8) << std::setprecision(5)<<1000000*(t/tries)
              << std::setw(10) << nnz
              << std::setw(10) << std::fixed << err
              << std::setw(8) << tries
              << std::endl;
}

int main(int argc, const char* argv[]) {
    int n=(1<<3);
    if(argc>1) { n=(1<<atoi(argv[1])); }
    TaylorCalculus calculus;

    std::cout<<"name                 time(us)  size     error   tries\n";

    profile("apply",ProfileReset(&calculus,henon_map,unit_box_model),n*100);
    profile("apply*5",ProfileReset(&calculus,henon_map,unit_box_model,5),n*100);
    profile("flow",ProfileFlow(&calculus,spiral_vf,unit_box/2,0.125,unit_box),n*100);
    profile("crossing",ProfileCrossing(&calculus,affine_pred_model,spiral_flow_model,unit_box_model),n*100);
}
