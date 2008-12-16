/***************************************************************************
 *            taylor_variable.h
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
 
/*! \file taylor_variable.h
 *  \brief Differential algebra variables with a sparse representation.
 */

#ifndef ARIADNE_TAYLOR_VARIABLE_H
#define ARIADNE_TAYLOR_VARIABLE_H

#include <map>

#include "macros.h"
#include "array.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "series.h"
#include "sparse_differential.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Series;
template<class D> class SparseDifferential;

class TaylorSeries;
class TaylorVariable;

TaylorVariable operator+(const TaylorVariable& x);
TaylorVariable operator-(const TaylorVariable& x);
TaylorVariable operator+(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable operator-(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable operator*(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable operator/(const TaylorVariable& x, const TaylorVariable& y);

TaylorVariable operator+(const TaylorVariable& x, const Float& c);
TaylorVariable operator-(const TaylorVariable& x, const Float& c);
TaylorVariable operator*(const TaylorVariable& x, const Float& c);
TaylorVariable operator/(const TaylorVariable& x, const Float& c);
TaylorVariable operator+(const Float& c, const TaylorVariable& x);
TaylorVariable operator-(const Float& c, const TaylorVariable& x);
TaylorVariable operator*(const Float& c, const TaylorVariable& x);
TaylorVariable operator/(const Float& c, const TaylorVariable& x);

TaylorVariable operator+(const TaylorVariable& x, const Interval& c);
TaylorVariable operator-(const TaylorVariable& x, const Interval& c);
TaylorVariable operator*(const TaylorVariable& x, const Interval& c);
TaylorVariable operator/(const TaylorVariable& x, const Interval& c);
TaylorVariable operator+(const Interval& c, const TaylorVariable& x);
TaylorVariable operator-(const Interval& c, const TaylorVariable& x);
TaylorVariable operator*(const Interval& c, const TaylorVariable& x);
TaylorVariable operator/(const Interval& c, const TaylorVariable& x);

TaylorVariable add(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable mul(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable neg(const TaylorVariable& x);
TaylorVariable rec(const TaylorVariable& x);
TaylorVariable pow(const TaylorVariable& x, int n);
TaylorVariable sqr(const TaylorVariable& x);
TaylorVariable sqrt(const TaylorVariable& x);
TaylorVariable exp(const TaylorVariable& x);
TaylorVariable log(const TaylorVariable& x);
TaylorVariable sin(const TaylorVariable& x);
TaylorVariable cos(const TaylorVariable& x);
TaylorVariable tan(const TaylorVariable& x);

pair<TaylorVariable,TaylorVariable> split(const TaylorVariable& x, uint j);
Vector<Interval> evaluate(const TaylorVariable& y, const Vector<Interval>& z);
TaylorVariable compose(const TaylorSeries& x, const TaylorVariable& y);
TaylorVariable derivative(const TaylorVariable& x, uint j);
TaylorVariable antiderivative(const TaylorVariable& x, uint j);






/*! \brief A class representing a quantity depending on other quantities. */
class TaylorVariable
{
    static const Float _zero;
    SparseDifferential<Float> _expansion;
    Interval _error;
  public:
    typedef MultiIndex IndexType;
    typedef Float ValueType;
    typedef Float ScalarType;
    typedef std::map<MultiIndex,Float>::iterator iterator;
    typedef std::map<MultiIndex,Float>::const_iterator const_iterator;

    TaylorVariable() : _expansion(), _error(0) { }
    TaylorVariable(uint as, uint deg) : _expansion(as,deg), _error(0) { }
    TaylorVariable(uint as) : _expansion(as,255), _error(0) { }
    TaylorVariable(const SparseDifferential<Float>& d, const Interval& e) : _expansion(d), _error(e) { }
    template<class XX, class XXX> TaylorVariable(uint as, uint deg, const XX& eps, const XXX* ptr) : _expansion(as,deg,ptr), _error(eps) { }

    TaylorVariable& operator=(const Float& c) { this->_expansion=c; this->_error=0; return *this; }
    TaylorVariable& operator=(const Interval& c) { this->_expansion=c.midpoint(); this->_error=(c-c.midpoint()); return *this; }

    const SparseDifferential<Float>& expansion() const { return this->_expansion; }
    const Interval& error() const { return this->_error; }
    Interval& error() { return this->_error; }
    void set_error(const Interval& e) { this->_error=e; }
    const Float& value() const { return this->_expansion.value(); }
    void set_value(const Float& x) { return this->_expansion.set_value(x); }


    Float& operator[](const uint& j) { return this->_expansion[j]; }
    Float& operator[](const MultiIndex& a) { return this->_expansion[a]; }

    void set_degree(uint d) { this->_expansion.set_degree(d); }

    const_iterator begin() const { return this->_expansion.begin(); }
    const_iterator end() const { return this->_expansion.end(); }
    uint argument_size() const { return this->_expansion.argument_size(); }
    uint degree() const { return this->_expansion.degree(); }
    const Float& operator[](const uint& j) const { return this->_expansion[j]; }
    const Float& operator[](const MultiIndex& a) const { return this->_expansion[a]; }

    static TaylorVariable constant(uint as, const Float& c) {
        TaylorVariable r(as,0u); r._expansion.set_value(c); return r; }
    static TaylorVariable variable(uint as, const Float& x, uint i) {
        TaylorVariable r(as,1u); r._expansion.set_value(x); r._expansion.set_gradient(i,1.0); return r; }
    static Vector<TaylorVariable> variables(const Vector<Float>& x) {
        Vector<TaylorVariable> result(x.size()); for(uint i=0; i!=x.size(); ++i) { 
            result[i]=TaylorVariable::variable(x.size(),x[i],i); } return result; }

    bool operator==(const TaylorVariable& sd) const {
        return this->_expansion==sd._expansion && this->_error == sd._error; }
    bool operator!=(const TaylorVariable& sd) const { 
        return !(*this==sd); }

    TaylorVariable& operator+=(const TaylorVariable& x);
    TaylorVariable& operator-=(const TaylorVariable& x);
    TaylorVariable& operator+=(const Float& c);
    TaylorVariable& operator+=(const Interval& c);
    TaylorVariable& operator-=(const Float& c);
    TaylorVariable& operator-=(const Interval& c);
    TaylorVariable& operator*=(const Float& c);
    TaylorVariable& operator*=(const Interval& c);
    TaylorVariable& operator/=(const Float& c);
    TaylorVariable& operator/=(const Interval& c);

    Vector<Interval> domain() const;
    Interval range() const;
    Interval evaluate(const Vector<Interval>& x) const;

    std::string str() const;

    friend TaylorVariable operator+(const TaylorVariable& x);
    friend TaylorVariable operator-(const TaylorVariable& x);
    friend TaylorVariable operator+(const TaylorVariable& x, const TaylorVariable& y);
    friend TaylorVariable operator-(const TaylorVariable& x, const TaylorVariable& y);
    friend TaylorVariable operator*(const TaylorVariable& x, const TaylorVariable& y);
    friend TaylorVariable operator/(const TaylorVariable& x, const TaylorVariable& y);
    friend TaylorVariable compose(const TaylorSeries& x, const TaylorVariable& y);
    friend TaylorVariable derivative(const TaylorVariable& x, uint i);
    friend TaylorVariable antiderivative(const TaylorVariable& x, uint i);
  public:
    void clean();
    void sweep(const Float& eps);
};




TaylorVariable add(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable mul(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable max(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable min(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable abs(const TaylorVariable& x);
TaylorVariable neg(const TaylorVariable& x);
TaylorVariable rec(const TaylorVariable& x);
TaylorVariable sqr(const TaylorVariable& x);
TaylorVariable pow(const TaylorVariable& x, int n);
TaylorVariable sqrt(const TaylorVariable& x);
TaylorVariable exp(const TaylorVariable& x);
TaylorVariable log(const TaylorVariable& x);
TaylorVariable sin(const TaylorVariable& x);
TaylorVariable cos(const TaylorVariable& x);
TaylorVariable tan(const TaylorVariable& x);

std::ostream& operator<<(std::ostream& os, const TaylorVariable& x);


} // namespace Ariadne

#endif // ARIADNE_TAYLOR_VARIABLE_H
