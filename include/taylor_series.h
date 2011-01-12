/***************************************************************************
 *            taylor_series.h
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

/*! \file taylor_series.h
 *  \brief Approximate power series in one variable on a bounded domain with a dense
 *  representation as a floating-point polynomial plus an error bound.
 */

#ifndef ARIADNE_TAYLOR_SERIES_H
#define ARIADNE_TAYLOR_SERIES_H

#include "numeric.h"
#include "series.h"
#include "stlio.h"

namespace Ariadne {

typedef Series<Interval>(*series_function_pointer)(uint,const Interval&);

// A deprecated class for computing power series expansions
// The difference between this and the Series class is that a TaylorSeries
// uses floating-point coefficients and an interval remainder, and is only
// valid on a subdomain, while a Series class uses arbitrary coefficients
// and is valid within the entire radius of convergence
class TaylorSeries {
    typedef Series<Interval>(*series_function_pointer)(uint,const Interval&);
  public:
    TaylorSeries(uint d) : expansion(d+1), error(0) { }
    TaylorSeries(uint degree, series_function_pointer function,
                 const Float& centre, const Interval& domain);
    uint degree() const { return expansion.size()-1; }
    Float& operator[](uint i) { return expansion[i]; }
    Array<Float> expansion;
    Interval error;
    void sweep(Float e) {
        for(uint i=0; i<=degree(); ++i) {
            if(abs(expansion[i])<=e) {
                error+=expansion[i]*Interval(-1,1);
                expansion[i]=0; } } }
};

inline
TaylorSeries::TaylorSeries(uint d, series_function_pointer fn,
                           const Float& c, const Interval& r)
    : expansion(d+1), error(0)
{
    Series<Interval> centre_series=fn(d,Interval(c));
    Series<Interval> range_series=fn(d,r);
    Interval p=1;
    Interval e=r-c;
    //std::cerr<<"\nc="<<c<<" r="<<r<<" e="<<e<<"\n";
    //std::cerr<<"centre_series="<<centre_series<<"\nrange_series="<<range_series<<"\n";
    for(uint i=0; i!=d; ++i) {
        this->expansion[i]=midpoint(centre_series[i]);
        this->error+=(centre_series[i]-this->expansion[i])*p;
        p*=e;
    }
    //this->expansion[d]=midpoint(centre_series[d]);
    this->expansion[d]=midpoint(range_series[d]);
    this->error+=(range_series[d]-this->expansion[d])*p;
    //std::cerr<<"expansion="<<this->expansion<<"\nerror="<<this->error<<"\n";
}

inline
std::ostream&
operator<<(std::ostream& os, const TaylorSeries& ts) {
    return os<<"TS("<<ts.expansion<<","<<ts.error<<")";
}


} // namespace Ariadne

#endif // ARIADNE_TAYLOR_SERIES_H
