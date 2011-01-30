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
void _mul_clear(IntervalTaylorModel& r, const IntervalTaylorModel& x, const IntervalTaylorModel& y);
void _mul_full(IntervalTaylorModel& r, const IntervalTaylorModel& x, const IntervalTaylorModel& y);
}


IntervalTaylorModel copy(const IntervalTaylorModel& x) {
    return x;
}

IntervalTaylorModel& isweep(IntervalTaylorModel& x) {
    x.sweep(); return x;
}

IntervalTaylorModel& ivladd(IntervalTaylorModel& x, const Interval& ivl) {
    return x+=ivl;
}

IntervalTaylorModel& ivlscal(IntervalTaylorModel& x, const Interval& ivl) {
    return x*=ivl;
}

IntervalTaylorModel& fscal(IntervalTaylorModel& x, const Float& c) {
    return x*=c;
}

IntervalTaylorModel& isum(IntervalTaylorModel& x, const IntervalTaylorModel& y) {
    return x+=y;
}

IntervalTaylorModel sum(const IntervalTaylorModel x1, const IntervalTaylorModel& x2) {
    return x1+x2;
}

IntervalTaylorModel prod(const IntervalTaylorModel x1, const IntervalTaylorModel& x2) {
    return x1*x2;
}

IntervalTaylorModel prod_full(const IntervalTaylorModel x1, const IntervalTaylorModel& x2) {
    IntervalTaylorModel r(x2.argument_size(),x2.sweeper()); _mul_full(r,x1,x2); return r;
}

IntervalTaylorModel prod_clear(const IntervalTaylorModel x1, const IntervalTaylorModel& x2) {
    IntervalTaylorModel r(x2.argument_size(),x2.sweeper()); _mul_clear(r,x1,x2); return r;
}


IntervalTaylorModel exp_cos(const IntervalTaylorModel x, const IntervalTaylorModel& y) {
    return exp(x)*cos(y);
}

IntervalTaylorModel sigmoid(const IntervalTaylorModel x, const IntervalTaylorModel& y) {
    const double a=10;
    return exp(-x/a);
}

typedef IntervalTaylorModel(*TaylorFunctionPtr)(const Vector<IntervalTaylorModel>&);


template<class T>
void profile(uint ntries, const char* name, const T& run)
{
    typename T::result_type res=run();
    //std::cerr<< "\n" << name << "(" << args << ")=\n  " << res << "\n\n";

    boost::timer tm;

    tm.restart();
    for(uint i=0; i!=ntries; ++i) {
        res=run();
    }

    double total_time = tm.elapsed();
    double average_time_in_microseconds = 1000000*(total_time/ntries);
    Float error = res.error();
    unsigned int size = res.number_of_nonzeros();

    std::cout << std::setw(20) << std::left << name << std::right
              << std::setw(10) << std::fixed << std::setprecision(2) << average_time_in_microseconds << " "
              << std::setw(12) << std::scientific << std::setprecision(4) << error
              << std::setw(8) << size
              << std::endl;
}

int main(int argc, const char* argv[]) {
    uint ntries=20;
    if(argc>1) { ntries=atoi(argv[1]); }

    Vector<Float> c(2, 1.0,2.0);
    int i;

    Vector<IntervalTaylorModel> v(2,IntervalTaylorModel(2));
    v[0]=IntervalTaylorModel::variable(2,0);
    v[1]=IntervalTaylorModel::constant(2,1.0);

    // Use in clean()
    IntervalTaylorModel w(3);
    i=0;
    for(MultiIndex a(3); a.degree()<=9; ++a) {
        if(i%7<3) { w.expansion().append(a,1/(1.0+i*i*i*i*i)); }
        else if(i%7<4) { w.expansion().append(a,1/(1.0+i)); }
        ++i;
    }
    //w.set_maximum_degree(7);
    w.set_sweeper(new ThresholdSweeper(1e-5));


    Interval ivl(0.33,0.49);
    Float cnst=0.41;
    IntervalTaylorModel x(3);
    IntervalTaylorModel y(3);

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
    y.set_sweeper(new ThresholdSweeper(1e-3));

    IntervalTaylorModel z(Expansion<Float>(3,4, 0,0,0,1.0, 1,0,0,0.5, 0,1,0,-0.25, 0,0,1,0.625),0.0);

    std::cout << std::setw(20) << std::left << "name" << std::right
              << std::setw(11) << "time(us)"
              << std::setw(12) << "error"
              << std::setw(8) << "size"
              << std::endl;

    typedef IntervalTaylorModel TM;
    //std::cerr<<"\n\nexp("<<z<<")=\n  "<<exp(z)<<"\n\n";
    //std::cerr<<"exp(1)="<<Ariadne::exp(1.0)<<"  exp([1:1])="<<Ariadne::exp(Interval(1))<<"\n";

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
