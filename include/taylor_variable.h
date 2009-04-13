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
 *  \brief Approximate functions on a bounded domain with a sparse representation.
 */

#ifndef ARIADNE_TAYLOR_VARIABLE_H
#define ARIADNE_TAYLOR_VARIABLE_H

#include <map>

#include "macros.h"
#include "array.h"
#include "vector.h"
#include "multi_index.h"
#include "taylor_model.h"

namespace Ariadne {

template<class T1, class T2> class Product;
template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Expansion;
template<class X> class Polynomial;
class TaylorModel;
class TaylorVariable;
class ExpressionInterface;

// Restrict to a smaller domain. REQUIRED
TaylorVariable restrict(const TaylorVariable& x, const Vector<Interval>& d);

// Test if a variable refines another
bool refines(const TaylorVariable& tv1, const TaylorVariable& tv2);

// Evaluate an array of Taylor variables on a vector.
Interval evaluate(const TaylorVariable& x, const Vector<Interval>& sy);

// Split the variable over two domains, subdividing along the independent variable j.
pair<TaylorVariable,TaylorVariable> split(const TaylorVariable& x, uint j);


// Embed the variable in a space of higher dimension
TaylorVariable embed(const TaylorVariable& tv1, const Interval& d2);
TaylorVariable embed(const Vector<Interval>& d1, const TaylorVariable& tv2);

// Antidifferentiation operator
TaylorVariable antiderivative(const TaylorVariable& x, uint k);

// Implicit function solver
TaylorVariable implicit(const TaylorVariable& f);


// Combine two functions over different domains
Vector<TaylorVariable> combine(const TaylorVariable& x1, const TaylorVariable& x2);


/*! \brief A class representing a quantity depending on other quantities.
 *  Based on a power series expansion, scaled to the unit box.
 *
 * See also Expansion, TaylorModel, TaylorFunction, TaylorSet.
 */
class TaylorVariable
{
    typedef Vector<Interval> DomainType;
    typedef TaylorModel ModelType;
    typedef Expansion<Float> ExpansionType;
    typedef Float ErrorType;
    static const Float _zero;
    DomainType _domain;
    ModelType _model;
  public:
    //! \brief An iterator through the (index,coefficient) pairs of the expansion expansion.
    typedef ExpansionType::iterator iterator;
    //! \brief A constant iterator through the (index,coefficient) pairs of the expansion expansion.
    typedef ExpansionType::const_iterator const_iterator;

    //@{
    /*! \name Constructors and destructors. */
    //! \brief Default constructor.
    explicit TaylorVariable();
    //! \brief Construct a TaylorVariable over the domain \a d.
    explicit TaylorVariable(const DomainType& d);
    //! \brief Construct a TaylorVariable over the domain \a d, based on the scaled model \a m.
    explicit TaylorVariable(const DomainType& d, const TaylorModel& m);
    //! \brief Construct a TaylorVariable over the domain \a d, with scaled power series expansion \a f and error \a e.
    explicit TaylorVariable(const DomainType& d, const ExpansionType& f, const ErrorType& e=0);
    //! \brief Construct a TaylorVariable over the domain \a d from the polynomial \a p.
    template<class X> explicit TaylorVariable(const DomainType& d, const Polynomial<X>& p);
    //! \brief Construct a TaylorVariable over the domain \a d from the polynomial \a p.
    explicit TaylorVariable(const DomainType& d, const ExpressionInterface& f);
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorVariable& operator=(const Float& c) { this->_model=c; return *this; }
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorVariable& operator=(const Interval& c) { this->_model=c; return *this; }
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorVariable constant(const DomainType& d, const Float& c);
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorVariable constant(const DomainType& d, const Interval& c);
    //! \brief Construct the quantity \f$x_j\f$ over the domain \a d.
    static TaylorVariable variable(const DomainType& d, unsigned int j);
    //! \brief Construct the quantity \f$c+\sum g_jx_j\f$ over the domain \a d.
    static TaylorVariable affine(const DomainType& d, const Float& c, const Vector<Float>& g);
    //! \brief Construct the quantity \f$c+\sum g_jx_j \pm e\f$ over domain \a d.
    static TaylorVariable affine(const DomainType& d, const Float& x, const Vector<Float>& g, const Float& e) ;

    //! \brief Return the vector of constants with values \a c over domain \a d.
    static Vector<TaylorVariable> constants(const DomainType& d, const Vector<Float>& c);
    //! \brief Return the vector of constants with interval values \a c over domain \a d.
    static Vector<TaylorVariable> constants(const DomainType& d, const Vector<Interval>& c);
    //! \brief Return the vector of variables with values \a x over domain \a d.
    static Vector<TaylorVariable> variables(const DomainType& d);
    //@}

    //@{
    /*! \name Data access */
    //! \brief The domain of the quantity.
    const DomainType& domain() const { return this->_domain; }
    //! \brief The scaled expansion over a unit box with error bound.
    const ModelType& model() const { return this->_model; }
    //! \brief The scaled expansion over a unit box.
    const ExpansionType& expansion() const { return this->_model._expansion; }
    //! \brief The error of the expansion over the domain.
    const ErrorType& error() const { return this->_model._error; }
    //! \brief A reference to the scaled expansion over a unit box with error bound.
    ModelType& model() { return this->_model; }
    //! \brief A reference to the expansion.
    ExpansionType& expansion() { return this->_model._expansion; }
    //! \brief A reference to the error of the expansion over the domain.
    ErrorType& error() { return this->_model._error; }
    //! \brief A reference to the constant term in the expansion.
    Float& value() { return this->_model.value(); }
    //! \brief The constant term in the expansion.
    const Float& value() const { return this->_model.value(); }

    //! \brief Set the error of the expansion.
    void set_error(const Float& ne) { this->_model.set_error(ne); }
    //! \brief Set the constant term in the expansion.
    void set_value(const Float& c) { this->_model.set_value(c); }

    //! \brief The coefficient of the term in $x^a$.
    const Float& operator[](const MultiIndex& a) const { return this->_model[a]; }
    //! \brief A read/write reference to the coefficient of the term in $x^a$.
    Float& operator[](const MultiIndex& a) { return this->_model[a]; }

    //! \brief An iterator to the first term in the expansion expansion.
    iterator begin() { return this->_model.begin(); }
    //! \brief A constant iterator to the first term in the expansion expansion.
    const_iterator begin() const { return this->_model.begin(); }
    //! \brief An iterator to the end of the expansion expansion.
    iterator end() { return this->_model.end(); }
    //! \brief A constant iterator to the end of the expansion expansion.
    const_iterator end() const { return this->_model.end(); }
    //! \brief An iterator to the term with index \a.
    iterator find(const MultiIndex& a) { return this->_model.find(a); }
    //! \brief A constant iterator to the term with index \a.
    const_iterator find(const MultiIndex& a) const { return this->_model.find(a); }

    //! \brief The number of variables in the argument of the quantity.
    uint argument_size() const { return this->_model.argument_size(); }
    //! \brief The maximum degree of terms in the expansion expansion.
    uint degree() const { return this->_model.degree(); }
    //! \brief The number of nonzero terms in the expansion expansion.
    uint number_of_nonzeros() const { return this->_model.number_of_nonzeros(); }
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    bool operator==(const TaylorVariable& tv) const;
    //! \brief Inequality operator.
    bool operator!=(const TaylorVariable& tv) const { return !(*this==tv); }
    //@}

    //@{
    /*! \name Function operations. */
    //! \brief An over-approximation to the range of the quantity.
    Interval range() const { return this->_model.range(); }
    //! \brief Evaluate the quantity at the point \a x.
    Interval evaluate(const Vector<Float>& x) const;
    //! \brief Evaluate the quantity over the interval of points \a x.
    Interval evaluate(const Vector<Interval>& x) const;
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Truncate to the default maximum degree of the quantity.
    TaylorVariable& truncate() { this->_model.truncate(); return *this; }
    //! \brief Truncate to degree \a deg.
    TaylorVariable& truncate(uint deg) { this->_model.truncate(deg); return *this; }
    //! \brief Truncate all terms with any coefficient higher than \a a.
    TaylorVariable& truncate(const MultiIndex& a) { this->_model.truncate(a); return *this; }
    //! \brief Truncate all terms with any coefficient higher than those given by \a b.
    TaylorVariable& truncate(const MultiIndexBound& b) { this->_model.truncate(b); return *this; }
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorVariable& sweep() { this->_model.sweep(); return *this; }
    //! \brief Remove all terms whose coefficient has magnitude less than \a eps.
    TaylorVariable& sweep(double eps) { this->_model.sweep(eps); return *this; }
    //! \brief Remove all terms whose degree is higher than \a deg or
    //! whose coefficient has magnitude less than \a eps.
    TaylorVariable& clean(const TaylorModel::Accuracy& acc) { this->_model.clean(acc); return *this; }
    //! \brief Remove all terms which have high degree or small magnitude.
    TaylorVariable& clean() { this->_model.clean(); return *this; }
    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \brief .
    void set_maximum_index(MultiIndexBound md) { this->_model.set_maximum_index(md); }
    //! \brief .
    void set_maximum_degree(uint md) { this->_model.set_maximum_degree(md); }
    //! \brief .
    void set_sweep_threshold(double me) { ARIADNE_ASSERT(me>=0.0); this->_model.set_sweep_threshold(me); }
    //! \brief .
    MultiIndexBound maximum_index() const { return this->_model.maximum_index(); }
    //! \brief .
    uint maximum_degree() const { return this->_model.maximum_degree(); }
    //! \brief .
    double sweep_threshold() const { return this->_model.sweep_threshold(); }
    //@}

    //@{
    /*! \name Non-arithmetic operations. */
    //! \brief Test if the quantity is a better approximation than \a t throughout the domain.
    friend bool refines(const TaylorVariable& x1, const TaylorVariable& x2);
    //! \brief Restrict to a subdomain.
    friend TaylorVariable restrict(const TaylorVariable& x, const DomainType& d);
    //@}

    /*! \name Arithmetic operations. */
    //! \brief Inplace addition of another variable.
    friend TaylorVariable& operator+=(TaylorVariable& x, const TaylorVariable& y);
    //! \brief Inplace subtraction of another variable.
    friend TaylorVariable& operator-=(TaylorVariable& x, const TaylorVariable& y);
    //! \brief Inplace addition of a product of two variables.
    friend TaylorVariable& operator+=(TaylorVariable& x, const Product<TaylorVariable,TaylorVariable>& y);
    //! \brief Inplace addition of an exact floating-point constant.
    friend TaylorVariable& operator+=(TaylorVariable& x, const Float& c);
    //! \brief Inplace addition of an interval constant.
    friend TaylorVariable& operator+=(TaylorVariable& x, const Interval& c);
    //! \brief Inplace subtraction of an exact floating-point constant.
    friend TaylorVariable& operator-=(TaylorVariable& x, const Float& c);
    //! \brief Inplace subtraction of an interval constant.
    friend TaylorVariable& operator-=(TaylorVariable& x, const Interval& c);
    //! \brief Inplace multiplication by an exact scalar.
    friend TaylorVariable& operator*=(TaylorVariable& x, const Float& c);
    //! \brief Inplace multiplication by an approximate scalar.
    friend TaylorVariable& operator*=(TaylorVariable& x, const Interval& c);
    //! \brief Inplace division by an exact scalar.
    friend TaylorVariable& operator/=(TaylorVariable& x, const Float& c);
    //! \brief Inplace division by an approximate scalar.
    friend TaylorVariable& operator/=(TaylorVariable& x, const Interval& c);

    //! \brief Unary plus.
    friend TaylorVariable operator+(const TaylorVariable& x);
    //! \brief Unary minus.
    friend TaylorVariable operator-(const TaylorVariable& x);
    //! \brief Addition.
    friend TaylorVariable operator+(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Subtraction.
    friend TaylorVariable operator-(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Multiplication.
    friend TaylorVariable operator*(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Division.
    friend TaylorVariable operator/(const TaylorVariable& x, const TaylorVariable& y);

    //! \brief Addition of a scakar.
    friend TaylorVariable operator+(const TaylorVariable& x, const Float& c);
    //! \brief Subtraction of a scakar.
    friend TaylorVariable operator-(const TaylorVariable& x, const Float& c);
    //! \brief Multiplication by a scakar.
    friend TaylorVariable operator*(const TaylorVariable& x, const Float& c);
    //! \brief Division by a scakar.
    friend TaylorVariable operator/(const TaylorVariable& x, const Float& c);
    //! \brief Addition of a scakar.
    friend TaylorVariable operator+(const TaylorVariable& x, const Interval& c);
    //! \brief Subtraction of a scakar.
    friend TaylorVariable operator-(const TaylorVariable& x, const Interval& c);
    //! \brief Multiplication by a scakar.
    friend TaylorVariable operator*(const TaylorVariable& x, const Interval& c);
    //! \brief Division by a scakar.
    friend TaylorVariable operator/(const TaylorVariable& x, const Interval& c);
    //! \brief Addition of a scakar.
    friend TaylorVariable operator+(const Float& c, const TaylorVariable& x);
    //! \brief Subtraction from a scakar.
    friend TaylorVariable operator-(const Float& c, const TaylorVariable& x);
    //! \brief Multiplication by a scakar.
    friend TaylorVariable operator*(const Float& c, const TaylorVariable& x);
    //! \brief Division through a scalar.
    friend TaylorVariable operator/(const Float& c, const TaylorVariable& x);
    //! \brief Addition of a scakar.
    friend TaylorVariable operator+(const Interval& c, const TaylorVariable& x);
    //! \brief Subtraction from a scakar.
    friend TaylorVariable operator-(const Interval& c, const TaylorVariable& x);
    //! \brief Multiplication by a scakar.
    friend TaylorVariable operator*(const Interval& c, const TaylorVariable& x);
    //! \brief Division through a scalar.
    friend TaylorVariable operator/(const Interval& c, const TaylorVariable& x);

    //! \brief Maximum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend TaylorVariable max(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Minimum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend TaylorVariable min(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Addition.
    friend TaylorVariable add(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Multiplication.
    friend TaylorVariable mul(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Absolute value. Throws an error if one variable is neither greater than
    //! zero over the entire domain, nor less than zero over the entire domain.
    friend TaylorVariable abs(const TaylorVariable& x);
    //! \brief Negation.
    friend TaylorVariable neg(const TaylorVariable& x);
    //! \brief Reciprocal.
    friend TaylorVariable rec(const TaylorVariable& x);
    //! \brief Square.
    friend TaylorVariable sqr(const TaylorVariable& x);
    //! \brief Power.
    friend TaylorVariable pow(const TaylorVariable& x, int n);
    //! \brief Square root.
    friend TaylorVariable sqrt(const TaylorVariable& x);
    //! \brief Natural exponent.
    friend TaylorVariable exp(const TaylorVariable& x);
    //! \brief Natural logarithm.
    friend TaylorVariable log(const TaylorVariable& x);
    //! \brief Sine.
    friend TaylorVariable sin(const TaylorVariable& x);
    //! \brief Cosine.
    friend TaylorVariable cos(const TaylorVariable& x);
    //! \brief Tangent.
    friend TaylorVariable tan(const TaylorVariable& x);
    //! \brief Inverse sine.
    friend TaylorVariable asin(const TaylorVariable& x);
    //! \brief Inverse cosine.
    friend TaylorVariable acos(const TaylorVariable& x);
    //! \brief Inverse tangent.
    friend TaylorVariable atan(const TaylorVariable& x);
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const TaylorVariable& x);
    //@}

  public:
    TaylorVariable& clobber() { this->_model.clobber(); return *this; }
    TaylorVariable& clobber(uint o) { this->_model.clobber(o); return *this; }
    TaylorVariable& clobber(uint so, uint to) { this->_model.clobber(so,to); return *this; }
};


inline tribool operator>(const TaylorVariable& x, const Float& c) {
    Interval r=x.range(); if(r.l>c) { return true; } else if(r.u<=c) { return false; } else { return indeterminate; } }
inline tribool operator<(const TaylorVariable& x, const Float& c) {
    Interval r=x.range(); if(r.l<c) { return true; } else if(r.u>=c) { return false; } else { return indeterminate; } }

inline tribool operator>(const TaylorVariable& x, const TaylorVariable& y) { return (x-y)>0; }
inline tribool operator<(const TaylorVariable& x, const TaylorVariable& y) { return (x-y)<0; }

inline TaylorVariable& operator+=(TaylorVariable& x, const Float& c) {
    x._model+=c; return x; }
inline TaylorVariable& operator-=(TaylorVariable& x, const Float& c) {
    x._model-=c; return x; }
inline TaylorVariable& operator*=(TaylorVariable& x, const Float& c) {
    x._model*=c; return x; }
inline TaylorVariable& operator/=(TaylorVariable& x, const Float& c) {
    x._model*=c; return x; }

inline TaylorVariable& operator+=(TaylorVariable& x, const Interval& c) {
    x._model+=c; return x; }
inline TaylorVariable& operator-=(TaylorVariable& x, const Interval& c) {
    x._model-=c; return x; }
inline TaylorVariable& operator*=(TaylorVariable& x, const Interval& c) {
    x._model*=c; return x; }
inline TaylorVariable& operator/=(TaylorVariable& x, const Interval& c) {
    x._model*=c; return x; }


inline TaylorVariable operator+(const TaylorVariable& x) {
    return TaylorVariable(x._domain,x._model); }
inline TaylorVariable operator-(const TaylorVariable& x) {
    return TaylorVariable(x._domain,-x._model); }

inline TaylorVariable operator+(const TaylorVariable& x, const Float& c) {
    return TaylorVariable(x._domain,x._model+c); }
inline TaylorVariable operator-(const TaylorVariable& x, const Float& c) {
    return TaylorVariable(x._domain,x._model-c); }
inline TaylorVariable operator*(const TaylorVariable& x, const Float& c) {
    return TaylorVariable(x._domain,x._model*c); }
inline TaylorVariable operator/(const TaylorVariable& x, const Float& c) {
    return TaylorVariable(x._domain,x._model/c); }
inline TaylorVariable operator+(const TaylorVariable& x, const Interval& c) {
    return TaylorVariable(x._domain,x._model+c); }
inline TaylorVariable operator-(const TaylorVariable& x, const Interval& c) {
    return TaylorVariable(x._domain,x._model-c); }
inline TaylorVariable operator*(const TaylorVariable& x, const Interval& c) {
    return TaylorVariable(x._domain,x._model*c); }
inline TaylorVariable operator/(const TaylorVariable& x, const Interval& c) {
    return TaylorVariable(x._domain,x._model/c); }
inline TaylorVariable operator+(const Float& c, const TaylorVariable& x) {
    return TaylorVariable(x._domain,c+x._model); }
inline TaylorVariable operator-(const Float& c, const TaylorVariable& x) {
    return TaylorVariable(x._domain,c-x._model); }
inline TaylorVariable operator*(const Float& c, const TaylorVariable& x) {
    return TaylorVariable(x._domain,c*x._model); }
inline TaylorVariable operator/(const Float& c, const TaylorVariable& x) {
    return TaylorVariable(x._domain,c/x._model); }
inline TaylorVariable operator+(const Interval& c, const TaylorVariable& x) {
    return TaylorVariable(x._domain,c+x._model); }
inline TaylorVariable operator-(const Interval& c, const TaylorVariable& x) {
    return TaylorVariable(x._domain,c-x._model); }
inline TaylorVariable operator*(const Interval& c, const TaylorVariable& x) {
    return TaylorVariable(x._domain,c*x._model); }
inline TaylorVariable operator/(const Interval& c, const TaylorVariable& x) {
    return TaylorVariable(x._domain,c/x._model); }


inline TaylorVariable abs(const TaylorVariable& x) {
    return TaylorVariable(x._domain,abs(x._model)); }
inline TaylorVariable neg(const TaylorVariable& x) {
    return TaylorVariable(x._domain,abs(x._model)); }
inline TaylorVariable rec(const TaylorVariable& x) {
    return TaylorVariable(x._domain,rec(x._model)); }
inline TaylorVariable sqr(const TaylorVariable& x) {
    return TaylorVariable(x._domain,sqr(x._model)); }
inline TaylorVariable pow(const TaylorVariable& x, int n) {
    return TaylorVariable(x._domain,pow(x._model,n)); }
inline TaylorVariable sqrt(const TaylorVariable& x) {
    return TaylorVariable(x._domain,sqrt(x._model)); }
inline TaylorVariable exp(const TaylorVariable& x) {
    return TaylorVariable(x._domain,exp(x._model)); }
inline TaylorVariable log(const TaylorVariable& x) {
    return TaylorVariable(x._domain,log(x._model)); }
inline TaylorVariable sin(const TaylorVariable& x) {
    return TaylorVariable(x._domain,sin(x._model)); }
inline TaylorVariable cos(const TaylorVariable& x) {
    return TaylorVariable(x._domain,cos(x._model)); }
inline TaylorVariable tan(const TaylorVariable& x) {
    return TaylorVariable(x._domain,tan(x._model)); }
inline TaylorVariable asin(const TaylorVariable& x) {
    return TaylorVariable(x._domain,asin(x._model)); }
inline TaylorVariable acos(const TaylorVariable& x) {
    return TaylorVariable(x._domain,acos(x._model)); }
inline TaylorVariable atan(const TaylorVariable& x) {
    return TaylorVariable(x._domain,atan(x._model)); }


inline Interval evaluate(const TaylorVariable& tv, const Vector<Interval>& x) {
    return tv.evaluate(x); }

inline TaylorVariable antiderivative(const TaylorVariable& x, uint k) {
    Interval sf=rad_ivl(x.domain()[k]);
    return TaylorVariable(x.domain(),antiderivative(x.model(),k)*sf); }

inline TaylorVariable embed(const TaylorVariable& tv1, const Interval& dom2) {
    return TaylorVariable(join(tv1.domain(),dom2),embed(tv1.model(),1u)); }
inline TaylorVariable embed(const TaylorVariable& tv1, const Vector<Interval>& dom2) {
    return TaylorVariable(join(tv1.domain(),dom2),embed(tv1.model(),dom2.size())); }
inline TaylorVariable embed(const Vector<Interval>& dom1, const TaylorVariable& tv2) {
    return TaylorVariable(join(dom1,tv2.domain()),embed(dom1.size(),tv2.model())); }



} // namespace Ariadne

#endif // ARIADNE_TAYLOR_VARIABLE_H
