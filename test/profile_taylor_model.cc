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
using StringType;

#include <boost/timer.hpp>
#include <boost/progress.hpp>

#include "config.h"
#include "function/taylor_model.h"
using namespace Ariadne;

template<class A1>
struct InplaceUnaryClosure {
    typedef A1 result_type;
    InplaceUnaryClosure(A1&(*f_)(A1&),const A1& a1_) : f(f_), a1(a1_) { };
    A1& operator()() const { r=a1; return f(r); }
    A1&(*f)(A1&); mutable A1 r; const A1& a1;
};

template<class A1, class A2>
struct InplaceBinaryClosure {
    typedef A1 result_type;
    InplaceBinaryClosure(A1&(*f_)(A1&,const A2&),const A1& a1_,const A2& a2_) : f(f_), a1(a1_), a2(a2_) { };
    A1& operator()() const { r=a1; return f(r,a2); }
    A1&(*f)(A1&,const A2&); mutable A1 r; const A1& a1; const A2& a2;
};

template<class R, class A1>
struct UnaryClosure {
    typedef R result_type;
    UnaryClosure(R(*f_)(const A1&),const A1& a1_) : f(f_), a1(a1_) { };
    R operator()() const { return f(a1); }
    R(*f)(const A1&); const A1& a1;
};

template<class R, class A1, class A2>
struct BinaryClosure {
    typedef R result_type;
    BinaryClosure(R(*f_)(const A1&,const A2&),const A1& a1_, const A2& a2_) : f(f_), a1(a1_), a2(a2_) { };
    R operator()() const { return f(a1,a2); }
    R(*f)(const A1&,const A2&); const A1& a1; const A2& a2;
};

template<class F, class A1>
InplaceUnaryClosure<A1>
inplace_bind(const F& f,const A1& a1) {
    return InplaceUnaryClosure<A1>((A1&(*)(A1&))f,a1);
}

template<class F, class A1, class A2>
InplaceBinaryClosure<A1,A2>
inplace_bind(const F& f,const A1& a1,const A2& a2) {
    return InplaceBinaryClosure<A1,A2>((A1&(*)(A1&,const A2&))f,a1,a2);
}

template<class F, class A1>
UnaryClosure<A1,A1>
bind(const F& f,const A1& a1) {
    return UnaryClosure<A1,A1>((A1(*)(const A1&))f,a1);
}

template<class F, class A1, class A2>
BinaryClosure<A1,A1,A2>
bind(const F& f,const A1& a1,const A2& a2) {
    return BinaryClosure<A1,A1,A2>((A1(*)(const A1&,const A2&))f,a1,a2);
}


namespace Ariadne {
Void _mul_clear(ValidatedTaylorModel64& r, const ValidatedTaylorModel64& x, const ValidatedTaylorModel64& y);
Void _mul_full(ValidatedTaylorModel64& r, const ValidatedTaylorModel64& x, const ValidatedTaylorModel64& y);
}


ValidatedTaylorModel64 copy(const ValidatedTaylorModel64& x) {
    return x;
}

ValidatedTaylorModel64& isweep(ValidatedTaylorModel64& x) {
    x.sweep(); return x;
}

ValidatedTaylorModel64& ivladd(ValidatedTaylorModel64& x, const ExactIntervalType& ivl) {
    return x+=ivl;
}

ValidatedTaylorModel64& ivlscal(ValidatedTaylorModel64& x, const ExactIntervalType& ivl) {
    return x*=ivl;
}

ValidatedTaylorModel64& fscal(ValidatedTaylorModel64& x, const Float64& c) {
    return x*=c;
}

ValidatedTaylorModel64& isum(ValidatedTaylorModel64& x, const ValidatedTaylorModel64& y) {
    return x+=y;
}

ValidatedTaylorModel64 sum(const ValidatedTaylorModel64 x1, const ValidatedTaylorModel64& x2) {
    return x1+x2;
}

ValidatedTaylorModel64 prod(const ValidatedTaylorModel64 x1, const ValidatedTaylorModel64& x2) {
    return x1*x2;
}

ValidatedTaylorModel64 prod_full(const ValidatedTaylorModel64 x1, const ValidatedTaylorModel64& x2) {
    ValidatedTaylorModel64 r(x2.argument_size(),x2.sweeper()); _mul_full(r,x1,x2); return r;
}

ValidatedTaylorModel64 prod_clear(const ValidatedTaylorModel64 x1, const ValidatedTaylorModel64& x2) {
    ValidatedTaylorModel64 r(x2.argument_size(),x2.sweeper()); _mul_clear(r,x1,x2); return r;
}


ValidatedTaylorModel64 exp_cos(const ValidatedTaylorModel64 x, const ValidatedTaylorModel64& y) {
    return exp(x)*cos(y);
}

ValidatedTaylorModel64 sigmoid(const ValidatedTaylorModel64 x, const ValidatedTaylorModel64& y) {
    const double a=10;
    return exp(-x/a);
}

typedef ValidatedTaylorModel64(*TaylorFunctionPtr)(const Vector<ValidatedTaylorModel64>&);


template<class T>
Void profile(Nat ntries, const char* name, const T& run)
{
    typename T::result_type res=run();
    //std::cerr<< "\n" << name << "(" << args << ")=\n  " << res << "\n\n";

    boost::timer tm;

    tm.restart();
    for(Nat i=0; i!=ntries; ++i) {
        res=run();
    }

    double total_time = tm.elapsed();
    double average_time_in_microseconds = 1000000*(total_time/ntries);
    Float64 error = res.error();
    unsigned Int size = res.number_of_nonzeros();

    std::cout << std::setw(20) << std::left << name << std::right
              << std::setw(10) << std::fixed << std::setprecision(2) << average_time_in_microseconds << " "
              << std::setw(12) << std::scientific << std::setprecision(4) << error
              << std::setw(8) << size
              << std::endl;
}

Int main(Int argc, const char* argv[]) {
    Nat ntries=20;
    if(argc>1) { ntries=atoi(argv[1]); }

    Vector<Float64> c={1.0,2.0};
    Int i;

    TrivialSweeper trivial_sweeper;
    Vector<ValidatedTaylorModel64> v(2,ValidatedTaylorModel64(2,trivial_sweeper));
    v[0]=ValidatedTaylorModel64::variable(2,0,trivial_sweeper);
    v[1]=ValidatedTaylorModel64::constant(2,1.0,trivial_sweeper);

    // Use in clean()
    ValidatedTaylorModel64 w(3,trivial_sweeper);
    i=0;
    for(MultiIndex a(3); a.degree()<=9; ++a) {
        if(i%7<3) { w.expansion().append(a,1/(1.0+i*i*i*i*i)); }
        else if(i%7<4) { w.expansion().append(a,1/(1.0+i)); }
        ++i;
    }
    //w.set_maximum_degree(7);
    ThresholdSweeper threshold_sweeper=ThresholdSweeper(1e-5);
    w.set_sweeper(threshold_sweeper);


    ExactIntervalType ivl(0.33,0.49);
    Float64 cnst=0.41;
    ValidatedTaylorModel64 x(3,threshold_sweeper);
    ValidatedTaylorModel64 y(3,threshold_sweeper);

    for(MultiIndex a(3); a.degree()<=7; ++a) {
        if(i%7<4) { x.expansion().append(a,1/(1.0+i)); }
    }
    i=0;
    for(MultiIndex a(3); a.degree()<=5; ++a) {
        if(i%7<4) { x.expansion().append(a,1/(1.0+i)); }
        if(i%3<2) { y.expansion().append(a,1/(2.0+i)); }
        ++i;
    }
    //y.set_maximum_degree(9);
    threshold_sweeper = ThresholdSweeper(1e-3);
    y.set_sweeper(threshold_sweeper);

    ValidatedTaylorModel64 z(Expansion<Float64>({ {{0,0,0},1.0}, {{1,0,0},0.5}, {{0,1,0},-0.25}, {{0,0,1},0.625} }),0.0,threshold_sweeper);

    std::cout << std::setw(20) << std::left << "name" << std::right
              << std::setw(11) << "time(us)"
              << std::setw(12) << "error"
              << std::setw(8) << "size"
              << std::endl;

    typedef ValidatedTaylorModel64 TM;
    //std::cerr<<"\n\nexp("<<z<<")=\n  "<<exp(z)<<"\n\n";
    //std::cerr<<"exp(1)="<<Ariadne::exp(1.0)<<"  exp([1:1])="<<Ariadne::exp(ExactIntervalType(1))<<"\n";

    profile(ntries*10000,"isweep-02",inplace_bind(&isweep,w));
    profile(ntries*10000,"copy-02",bind(&copy,w));
    profile(ntries*10000,"iadd-noinsert-02",inplace_bind(&ivladd,x,ivl));
    profile(ntries*10000,"iadd-insert-02",inplace_bind(&ivladd,x,ivl));
    profile(ntries*10000,"iscal-02",inplace_bind(&ivlscal,x,ivl));
    profile(ntries*10000,"fscal-02",inplace_bind(&fscal,x,cnst));
    profile(ntries*1000,"isum-02",inplace_bind(&isum,x,y));
    profile(ntries*1000,"sum-02",bind(&sum,x,y));
    profile(ntries*10,"prod-02",bind(&prod,x,y));
    profile(ntries*10,"prod_clear-02",bind(&prod_clear,x,y));
    profile(ntries*10,"prod_full-02",bind(&prod_full,x,y));
    profile(ntries*1,"exp-02",bind((TM(*)(const TM&))&Ariadne::exp,z));
    profile(ntries,"exp_cos-01",bind(exp_cos,v[0],v[1]));
    profile(ntries,"sigmoid-01",bind(sigmoid,v[0]));
}
