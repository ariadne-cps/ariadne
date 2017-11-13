/***************************************************************************
 *            series.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file series.hpp
 *  \brief Taylor series in a single variable.
 */

#ifndef ARIADNE_SERIES_HPP
#define ARIADNE_SERIES_HPP

#include "utility/array.hpp"
#include "utility/container.hpp"
#include "numeric/operators.hpp"

namespace Ariadne {

template<class X> class Series;

//! \brief Algorithm for computing a power series.
template<class X> class SeriesGeneratorInterface {
  public:
    virtual ~SeriesGeneratorInterface() = default;
    //! \brief Compute the \a d-th term of a power series, given centering at \a x and previous terms \a s.
    virtual X _next(DegreeType d, X const& x, List<X>& s) const = 0;
};

template<class OP, class X> class SeriesGenerator;

template<class X> class SeriesGenerator<Rec,X> : public SeriesGeneratorInterface<X> {
    virtual X _next(DegreeType deg, X const& c, List<X>& y) const final { Int d=deg; return (d==0) ? (1/c) : y[d-1]*(-y[0]); }
};

template<class X> class SeriesGenerator<Sqrt,X> : public SeriesGeneratorInterface<X> {
    virtual X _next(DegreeType deg, X const& c, List<X>& y) const final { Int d=deg; return (d==0) ? sqrt(c) : ((2*d-3)*(-1/c)/2)/d*y[d-1]; }
};

template<class X> class SeriesGenerator<Exp,X> : public SeriesGeneratorInterface<X> {
    virtual X _next(DegreeType deg, X const& c, List<X>& y) const final { Int d=deg; return (d==0) ? exp(c) : y[d-1]/d; }
};

template<class X> class SeriesGenerator<Log,X> : public SeriesGeneratorInterface<X> {
    virtual X _next(DegreeType deg, X const& c, List<X>& y) const final { Int d=deg; return (d==0) ? log(c) : (d==1) ? 1/c : -y[d-1]*y[1]*(d-1)/d; }
};

template<class X> class SeriesGenerator<Sin,X> : public SeriesGeneratorInterface<X> {
    virtual X _next(DegreeType deg, X const& c, List<X>& y) const final { Int d=deg; return (d==0) ? sin(c) : (d==1) ? cos(c) : -y[d-2]/(d*(d-1)); }
};

template<class X> class SeriesGenerator<Cos,X> : public SeriesGeneratorInterface<X> {
    virtual X _next(DegreeType deg, X const& c, List<X>& y) const final { Int d=deg; return (d==0) ? cos(c) : (d==1) ? -sin(c) : -y[d-2]/(d*(d-1)); }
};

template<class X> class SeriesGenerator<Tan,X> : public SeriesGeneratorInterface<X> {
    virtual X _next(DegreeType deg, X const& c, List<X>& y) const final {
        Int d=deg;
        if (d==0) { return tan(c); }
        else if (d==1) { return 1+sqr(y[0]); }
        else if (d==2) { return y[0]*y[1]; }
        else { X r=nul(c); for(DegreeType i=0; i!=d; ++i) { r+=y[i]*y[d-1-i]; } return r/d; }
    }
};

template<class X> class SeriesGenerator<Atan,X> : public SeriesGeneratorInterface<X> {
    virtual X _next(DegreeType deg, X const& c, List<X>& y) const final {
        Int d=deg;
        if (d==0) { return atan(c); }
        else if (d==1) { return rec(1+sqr(c)); }
        else if (d==2) { return -c*sqr(y[1]); }
        else { return -y[1]*((d-2)*y[d-2]+2*c*(d-1)*y[d-1])/d; }
    }
};



//! \brief A power series, centred on a value, with unlimited many coefficients.
template<class X>
class Series
{
    std::shared_ptr<const SeriesGeneratorInterface<X>> _ptr;
    X _centre;
    mutable List<X> _data;
  private:
    friend class AnalyticFunction;
    Series<X>(std::shared_ptr<const SeriesGeneratorInterface<X>> p, X const& c, std::nullptr_t dummy);
    Void _compute(DegreeType n) const;
    OutputStream& _write(OutputStream& os) const;
  public:
    template<class OP> Series(OP op, X const& x);
    const X& operator[](DegreeType n) const;
    const List<X>& coefficients(DegreeType n) const;
    friend OutputStream& operator<<(OutputStream& os, Series<X> const& s) { return s._write(os); }

    static Series<X> rec(const X& x) { return Series<X>(Rec(),x); }
    static Series<X> sqrt(const X& x) { return Series<X>(Sqrt(),x); }
    static Series<X> exp(const X& x) { return Series<X>(Exp(),x); }
    static Series<X> log(const X& x) { return Series<X>(Log(),x); }

    static Series<X> sin(const X& x) { return Series<X>(Sin(),x); }
    static Series<X> cos(const X& x) { return Series<X>(Cos(),x); }
    static Series<X> tan(const X& x) { return Series<X>(Tan(),x); }
    static Series<X> atan(const X& x) { return Series<X>(Atan(),x); }
};


template<class X> template<class OP> Series<X>::Series(OP op, X const& x)
    : Series(std::make_shared<SeriesGenerator<OP,X>>(),x, nullptr) {
}

template<class X> inline List<X> const& Series<X>::coefficients(DegreeType n) const {
    this->_compute(n); return this->_data;
}

template<class X> inline X const& Series<X>::operator[] (DegreeType n) const {
    this->_compute(n); return this->_data[n];
}

template<class X> inline Series<X>::Series(std::shared_ptr<const SeriesGeneratorInterface<X>> p, X const& c, std::nullptr_t)
    : _ptr(p), _centre(c), _data() {
}

template<class X> inline Void Series<X>::_compute(DegreeType n) const {
    while(_data.size()<=n) { _data.append(_ptr->_next(_data.size(),_centre,_data)); }
}

template<class X> inline OutputStream& Series<X>::_write(OutputStream& os) const {
    this->_compute(2); return os << "Series( centre=" << this->_centre << ", terms=" << this->_data << " )"; }

template<class OP, class X> inline Series<X> make_series(OP op, DegreeType d, const X& x) { return Series<X>(op,x); }
template<class OP, class X> inline Series<X> make_series(OP op, const X& x) { return Series<X>(op,x); }


class AnalyticFunction {
    OperatorCode _op_code;
    std::shared_ptr<const SeriesGeneratorInterface<ApproximateNumericType>> _asg_ptr;
    std::shared_ptr<const SeriesGeneratorInterface<ValidatedNumericType>> _vsg_ptr;
  public:
    template<class OP> explicit AnalyticFunction(OP op)
        : _op_code(op.code()), _asg_ptr(new SeriesGenerator<OP,ApproximateNumericType>()), _vsg_ptr(new SeriesGenerator<OP,ValidatedNumericType>()) { }
    Series<ApproximateNumericType> series(ApproximateNumericType c) const {
        return Series<ApproximateNumericType>(_asg_ptr,c,nullptr); }
    Series<ValidatedNumericType> series(ValidatedNumericType c) const {
        return Series<ValidatedNumericType>(_vsg_ptr,c,nullptr); }
    Series<ValidatedNumericType> series(ExactNumericType c) const {
        return Series<ValidatedNumericType>(_vsg_ptr,c,nullptr); }
    friend OutputStream& operator<<(OutputStream& os, AnalyticFunction const& af) {
        return os << "AnalyticFunction( " << af._op_code << " )"; }
};


// Code below is for missing functions
/*

template<class X> Series<X> Series<X>::asin(Nat d, const X& c) {
    if(d==0) { Series<X> y(d,c*0); y[0]=Ariadne::asin(c); return y; }
    Series<X> y(d-1,c*0); y[0]=c; y[1]=X(1);
    y = Ariadne::rec(Ariadne::sqrt(X(1)-Ariadne::sqr(y)));
    return antiderivative(y,Ariadne::asin(c));
}

template<class X> Series<X> Series<X>::acos(Nat d, const X& c) {
    if(d==0) { Series<X> y(d,c*0); y[0]=Ariadne::acos(c); return y; }
    Series<X> y(d-1,c*0); y[0]=c; y[1]=X(1);
    y = Ariadne::neg(Ariadne::rec(Ariadne::sqrt(X(1)-Ariadne::sqr(y))));
    return antiderivative(y,Ariadne::acos(c));
}

template<class X> Series<X> Series<X>::atan(Nat d, const X& c) {
    if(d==0) { Series<X> y(d,c*0); y[0]=Ariadne::atan(c); return y; }
    Series<X> y(d-1,c*0); y[0]=c; y[1]=X(1);
    y = Ariadne::rec(Ariadne::sqrt(X(1)+Ariadne::sqr(y)));
    return antiderivative(y,Ariadne::atan(c));
}

template<class X> Series<X> antiderivative(const Series<X>& x, const X& c) {
    Nat n=x.degree(); Series<X> r(n+1,x[0]*0);
    r[0]=c; for(Nat i=0; i<=n; ++i) { r[i+1]=x[i]/(i+1); }
    return r;
}

*/


} // namespace Ariadne


#endif // ARIADNE_SERIES_HPP

