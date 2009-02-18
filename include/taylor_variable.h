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

class TaylorVariable;
template<> class Vector<TaylorVariable>;

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
pair< Vector<TaylorVariable>, Vector<TaylorVariable> > split(const Vector<TaylorVariable>& x, uint j);

// Scale the variabe by post-composing with an affine map taking the interval \a ivl to the unit interval
TaylorVariable unscale(const TaylorVariable& x, const Interval& ivl);
Vector<TaylorVariable> unscale(const Vector<TaylorVariable>& x, const Vector<Interval>& bx);

// Scale the variabe by post-composing with an affine map taking the unit interval to \a ivl.
TaylorVariable scale(const TaylorVariable& x, const Interval& ivl);
Vector<TaylorVariable> scale(const Vector<TaylorVariable>& x, const Vector<Interval>& bx);

// Evaluate an array of Taylor variables on a vector.
Vector<Interval> evaluate(const Vector<TaylorVariable>& x, const Vector<Interval>& sy);
Interval evaluate(const TaylorVariable& x, const Vector<Interval>& sy);

Vector<Float> value(const Vector<TaylorVariable>& x);
Matrix<Float> jacobian(const Vector<TaylorVariable>& x, const Vector<Float>& sy);


// Compose an array of Taylor variables with another, after scaling by the interval vectors
Vector<TaylorVariable> compose(const Vector<TaylorVariable>& x, const Vector<Interval>& bx, const Vector<TaylorVariable>& y);

// Wrappers for univariate composition
TaylorVariable compose(const TaylorVariable& x, const Vector<Interval>& bx, const Vector<TaylorVariable>& y);
TaylorVariable compose(const TaylorVariable& x, const Interval& b, const TaylorVariable& y);

SparseDifferential<Float> expansion(const TaylorVariable& x, const Vector<Interval>& d);
Vector< SparseDifferential<Float> > expansion(const Vector<TaylorVariable>& x, const Vector<Interval>& d);

TaylorVariable antiderivative(const TaylorVariable& x, const Interval& dk, uint k);
Vector<TaylorVariable> antiderivative(const Vector<TaylorVariable>& x, const Interval& dk, uint k);

pair<Vector<Interval>,Interval> bounds(const Vector<TaylorVariable>& vf, const Vector<Interval>& d, const Interval& maxh);
Vector<TaylorVariable> flow(const Vector<TaylorVariable>& vf, const Vector<Interval>& d, const Interval& h, const Vector<Interval>& b);
Vector<TaylorVariable> implicit(const Vector<TaylorVariable>& x, const Vector<Interval>& d);

TaylorVariable embed(const TaylorVariable& tv, uint as, uint b);
Vector<TaylorVariable> embed(const Vector<TaylorVariable>& tvs, uint as, uint b);

bool refines(const TaylorVariable& tv1, const TaylorVariable& tv2);
bool refines(const Vector<TaylorVariable>& tv1, const Vector<TaylorVariable>& tv2);



/*! \brief A class representing a quantity depending on other quantities. */
class TaylorVariable
{
    static const Float _zero;
    SparseDifferential<Float> _expansion;
    Float _error;
    double _sweep_threshold;
    uint _maximum_degree;
  private:
    static double _default_sweep_threshold;
    static uint _default_maximum_degree;
  public:
    static const double em;
    static const double ec;
  public:
    typedef MultiIndex IndexType;
    typedef Float ValueType;
    typedef Float ScalarType;
    typedef std::map<MultiIndex,Float>::iterator iterator;
    typedef std::map<MultiIndex,Float>::const_iterator const_iterator;

    TaylorVariable() : _expansion(0), _error(0), _sweep_threshold(_default_sweep_threshold), _maximum_degree(_default_maximum_degree) { }
    TaylorVariable(uint as) : _expansion(as), _error(0), _sweep_threshold(_default_sweep_threshold), _maximum_degree(_default_maximum_degree) { }
    TaylorVariable(const SparseDifferential<Float>& d, const Float& e);
    TaylorVariable(uint as, uint deg, const double* ptr, const double& err);
    TaylorVariable(uint as, uint deg, double d0, ...);

    TaylorVariable& operator=(const Float& c) { this->_expansion=c; this->_error=0; return *this; }
    TaylorVariable& operator=(const Interval& c) { this->_expansion=c.midpoint(); this->_error=c.radius(); return *this; }

    const SparseDifferential<Float>& expansion() const { return this->_expansion; }
    SparseDifferential<Float>& expansion() { return this->_expansion; }
    const Float& error() const { return this->_error; }
    Float& error() { return this->_error; }
    const Float& value() const { return this->_expansion.value(); }
    Float& value() { return this->_expansion.data()[MultiIndex::zero(this->argument_size())]; }

    void set_error(const Float& ne) { ARIADNE_ASSERT(ne>=0); this->_error=ne; }

    Float& operator[](uint j) { return this->_expansion[j]; }
    Float& operator[](const MultiIndex& a) { return this->_expansion[a]; }
    const Float& operator[](uint j) const { return this->_expansion[j]; }
    const Float& operator[](const MultiIndex& a) const { return this->_expansion[a]; }

    iterator begin() { return this->_expansion.begin(); }
    iterator end() { return this->_expansion.end(); }
    const_iterator begin() const { return this->_expansion.begin(); }
    const_iterator end() const { return this->_expansion.end(); }
    
    uint argument_size() const { return this->_expansion.argument_size(); }
    uint degree() const { return this->_expansion.degree(); }
    
    static TaylorVariable zero(uint as) {
        TaylorVariable r(as); r._expansion.set_value(0.0); return r; }
    static TaylorVariable constant(uint as, const Float& c) {
        TaylorVariable r(as); r._expansion.set_value(c); return r; }
    static TaylorVariable variable(uint as, const Float& x, uint i) {
        TaylorVariable r(as); r._expansion.set_value(x); r._expansion.set_gradient(i,1.0); return r; }
    static TaylorVariable affine(const Float& x, const Vector<Float>& dx) {
        TaylorVariable r(dx.size()); r._expansion.set_value(x); for(uint j=0; j!=dx.size(); ++j) { r._expansion.set_gradient(j,dx[j]); } return r; }
    static Vector<TaylorVariable> zeroes(uint rs, uint as);
    static Vector<TaylorVariable> constants(uint as, const Vector<Float>& c);
    static Vector<TaylorVariable> variables(const Vector<Float>& x);

    bool operator==(const TaylorVariable& sd) const {
        return this->_expansion==sd._expansion && this->_error == sd._error; }
    bool operator!=(const TaylorVariable& sd) const { 
        return !(*this==sd); }

    TaylorVariable& scal(const Float& c);
    TaylorVariable& scal(const Interval& c);

    TaylorVariable& acc(const Float& c);
    TaylorVariable& acc(const Interval& c);
    TaylorVariable& acc(const TaylorVariable& x);
    TaylorVariable& acc(const TaylorVariable& x, const TaylorVariable& y);
    
    Vector<Interval> domain() const;
    Interval range() const;
    Interval evaluate(const Vector<Interval>& x) const;
    
    template<class XX> Interval evaluate(const Vector<XX>& x) const;

    static void set_default_maximum_degree(uint md) { _default_maximum_degree=md; }
    static void set_default_sweep_threshold(double me) { ARIADNE_ASSERT(me>=0.0); _default_sweep_threshold=me; }
    static uint default_maximum_degree() { return _default_maximum_degree; }
    static double default_sweep_threshold() { return _default_sweep_threshold; }
    void set_maximum_degree(uint md) { this->_maximum_degree=md; }
    void set_sweep_threshold(double me) { ARIADNE_ASSERT(me>=0.0); this->_sweep_threshold=me; }
    uint maximum_degree() const { return this->_maximum_degree; }
    double sweep_threshold() const { return this->_sweep_threshold; }

    std::string str() const;

    friend TaylorVariable& operator+=(TaylorVariable& x, const TaylorVariable& y);
    friend TaylorVariable& operator-=(TaylorVariable& x, const TaylorVariable& y);
    friend TaylorVariable& operator+=(TaylorVariable& x, const Float& c);
    friend TaylorVariable& operator+=(TaylorVariable& x, const Interval& c);
    friend TaylorVariable& operator-=(TaylorVariable& x, const Float& c);
    friend TaylorVariable& operator-=(TaylorVariable& x, const Interval& c);
    friend TaylorVariable& operator*=(TaylorVariable& x, const Float& c);
    friend TaylorVariable& operator*=(TaylorVariable& x, const Interval& c);
    friend TaylorVariable& operator/=(TaylorVariable& x, const Float& c);
    friend TaylorVariable& operator/=(TaylorVariable& x, const Interval& c);

    friend TaylorVariable operator+(const TaylorVariable& x);
    friend TaylorVariable operator-(const TaylorVariable& x);
    friend TaylorVariable operator+(const TaylorVariable& x, const TaylorVariable& y);
    friend TaylorVariable operator-(const TaylorVariable& x, const TaylorVariable& y);
    friend TaylorVariable operator*(const TaylorVariable& x, const TaylorVariable& y);
    friend TaylorVariable operator/(const TaylorVariable& x, const TaylorVariable& y);

    friend TaylorVariable max(const TaylorVariable& x, const TaylorVariable& y);
    friend TaylorVariable min(const TaylorVariable& x, const TaylorVariable& y);
    friend TaylorVariable abs(const TaylorVariable& x);
    friend TaylorVariable neg(const TaylorVariable& x);
    friend TaylorVariable rec(const TaylorVariable& x);
    friend TaylorVariable sqr(const TaylorVariable& x);
    friend TaylorVariable pow(const TaylorVariable& x, int n);
    friend TaylorVariable sqrt(const TaylorVariable& x);
    friend TaylorVariable exp(const TaylorVariable& x);
    friend TaylorVariable log(const TaylorVariable& x);
    friend TaylorVariable sin(const TaylorVariable& x);
    friend TaylorVariable cos(const TaylorVariable& x);
    friend TaylorVariable tan(const TaylorVariable& x);

    friend std::ostream& operator<<(std::ostream& os, const TaylorVariable& x);

  public:
    void clean();
    TaylorVariable& truncate(uint d);
    TaylorVariable& truncate(const MultiIndex& a);
    TaylorVariable& truncate();
    TaylorVariable& sweep(double eps);
    TaylorVariable& sweep();
  public:
    TaylorVariable& clobber();
    TaylorVariable& clobber(uint o);
    TaylorVariable& clobber(uint so, uint to);
};

TaylorVariable max(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable min(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable abs(const TaylorVariable& x);
    

template<> 
class Vector<TaylorVariable>
    : public ublas::vector<TaylorVariable>
{
  public:
    Vector() : ublas::vector<TaylorVariable>() { }
    Vector(uint rs) : ublas::vector<TaylorVariable>(rs) { for(uint i=0; i!=rs; ++i) { (*this)[i]=TaylorVariable(); } }
    Vector(uint rs, uint as) : ublas::vector<TaylorVariable>(rs) { for(uint i=0; i!=rs; ++i) { (*this)[i]=TaylorVariable(as); } }
    Vector(uint rs, const TaylorVariable& x) : ublas::vector<TaylorVariable>(rs) { for(uint i=0; i!=rs; ++i) { (*this)[i]=x; } }
    template<class E> Vector(const ublas::vector_expression<E> &ve) : ublas::vector<TaylorVariable>(ve) { }
    uint result_size() const { return this->size(); }
    uint argument_size() const { 
        ARIADNE_ASSERT(this->size()>0); 
        for(uint i=1; i!=this->size(); ++i) { ARIADNE_ASSERT((*this)[0].argument_size()==(*this)[i].argument_size()); } 
        return (*this)[0].argument_size(); }
    
    void check() const { for(uint i=0; i!=this->size(); ++i) { ARIADNE_ASSERT((*this)[0].argument_size()==(*this)[i].argument_size()); } }
};


inline Vector<TaylorVariable> TaylorVariable::zeroes(uint rs, uint as) {
    Vector<TaylorVariable> result(rs); for(uint i=0; i!=rs; ++i) { 
        result[i]=TaylorVariable::zero(as); } return result; }
inline Vector<TaylorVariable> TaylorVariable::constants(uint as, const Vector<Float>& c) {
    Vector<TaylorVariable> result(c.size()); for(uint i=0; i!=c.size(); ++i) { 
        result[i]=TaylorVariable::constant(as,c[i]); } return result; }
inline Vector<TaylorVariable> TaylorVariable::variables(const Vector<Float>& x) {
    Vector<TaylorVariable> result(x.size()); for(uint i=0; i!=x.size(); ++i) { 
        result[i]=TaylorVariable::variable(x.size(),x[i],i); } return result; }


template<class XX> Interval TaylorVariable::evaluate(const Vector<XX>& x) const {
    return this->evaluate(Vector<Interval>(x));
}

} // namespace Ariadne

#endif // ARIADNE_TAYLOR_VARIABLE_H
