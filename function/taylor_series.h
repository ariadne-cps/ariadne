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
 *  \brief Approximate power series in one variable on a is_bounded domain with a dense
 *  representation as a floating-point polynomial plus an error bound.
 */

#ifndef ARIADNE_TAYLOR_SERIES_H
#define ARIADNE_TAYLOR_SERIES_H

#include "numeric/numeric.h"
#include "algebra/series.h"
#include "geometry/interval.h"
#include "utility/container.h"

namespace Ariadne {

class AnalyticFunction;

// A univariate Taylor function for computing power series expansions
// The difference between this and the Series class is that a TaylorSeries
// uses floating-point coefficients and an interval remainder, and is only
// valid on a subdomain, while a Series class uses arbitrary coefficients
// and is valid within the entire radius of convergence
template<class X> class TaylorSeries;

template<> class TaylorSeries<ValidatedFloat64> {
    ExactInterval _domain;
    Array<ExactFloat64> _expansion;
    ErrorFloat64 _error;
  public:
    TaylorSeries(const ExactInterval& dom, DegreeType deg) : _domain(dom), _expansion(deg+1), _error(0u) { }

    TaylorSeries(const ExactInterval& domain, const ExactFloat64& centre, DegreeType degree,
                 AnalyticFunction const& function);

    template<class OP> TaylorSeries(OP unary_operator, const ExactInterval& domain, const ExactFloat64& centre, DegreeType degree);

    DegreeType degree() const { return _expansion.size()-1; }
    ExactFloat64 const& operator[](DegreeType i) const { return _expansion[i]; }
    Array<ExactFloat64> expansion() const { return _expansion; }
    ErrorFloat64 error() const { return _error; }
    Void sweep(ExactFloat64 threshold);

    friend OutputStream& operator<<(OutputStream&, TaylorSeries<ValidatedFloat64> const&);
};


template<class OP> inline
TaylorSeries<ValidatedFloat64>::TaylorSeries(OP unary_operator, const ExactInterval& domain, const ExactFloat64& centre, DegreeType degree)
    : _domain(domain), _expansion(degree+1), _error(0u)
{
    Series<ValidatedNumericType> centre_series=Series<ValidatedFloat64>(unary_operator,ValidatedNumericType(centre));
    Series<ValidatedNumericType> range_series=Series<ValidatedFloat64>(unary_operator,ValidatedNumericType(cast_singleton(domain)));
    for(DegreeType i=0; i!=degree; ++i) {
        this->_expansion[i]=centre_series[i].value();
        this->_error+=mag(centre_series[i]-this->_expansion[i]);
    }
    DegreeType d=degree;
    this->_expansion[d]=range_series[d].value();
    this->_error+=mag(range_series[d]-this->_expansion[d]);
}


inline
Void TaylorSeries<ValidatedFloat64>::sweep(ExactFloat64 threshold) {
    for(DegreeType i=0; i<=degree(); ++i) {
        if(definitely(mag(_expansion[i])<=threshold)) {
            _error+=mag(_expansion[i]);
            _expansion[i]=0;
        }
    }
}

inline
OutputStream& operator<<(OutputStream& os, const TaylorSeries<ValidatedFloat64>& ts) {
    return os<<"TS("<<ts._expansion<<"+/-"<<ts._error<<")";
}


} // namespace Ariadne

#endif // ARIADNE_TAYLOR_SERIES_H
