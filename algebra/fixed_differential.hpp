/***************************************************************************
 *            fixed_differential.hpp
 *
 *  Copyright 2011-17  Pieter Collins
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

/*! \file fixed_differential.hpp
 *  \brief First- and second-order multivariate differentials.
 */

#ifndef ARIADNE_FIXED_DIFFERENTIAL_HPP
#define ARIADNE_FIXED_DIFFERENTIAL_HPP

#include <map>

#include "utility/macros.hpp"
#include "utility/array.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/series.hpp"
#include "algebra/expansion.hpp"
#include <boost/concept_check.hpp>

namespace Ariadne {

class Float64;
class ExactIntervalType;

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Series;

template<class X, class T> struct EnableIfScalar { typedef T Type; };
template<class X, class T> struct EnableIfScalar<Vector<X>,T> { };
template<class X, class T> struct EnableIfScalar<Matrix<X>,T> { };
template<class X, class T> struct EnableIfVector { };
template<class X, class T> struct EnableIfVector<Vector<X>,T> { typedef T Type; };
template<class X, class T> struct EnableIfMatrix { };
template<class X, class T> struct EnableIfMatrix<Matrix<X>,T> { typedef T Type; };



template<class X> class SecondDifferential;

template<class V1, class V2>
struct SymmetricOuterProduct
    : public MatrixExpression< SymmetricOuterProduct<V1,V2> >
{
    typedef typename V1::value_type value_type;
    const V1& _v1; const V2& _v2;
    SymmetricOuterProduct<V1,V2>(const V1& v1, const V2& v2) : _v1(v1), _v2(v2) { }
    value_type operator() (Nat i, Nat j) { return _v1[i]*_v2[j]+_v1[j]*_v2[i]; }
};
template<class V1, class V2> inline
SymmetricOuterProduct<V1,V2>
symmetric_outer_product(const V1& v1, const V2& v2) {
    return SymmetricOuterProduct<V1,V2>(v1,v2);
}

template<class V1, class V2>
struct OuterProduct
    : public MatrixExpression< OuterProduct<V1,V2> >
{
    typedef typename V1::value_type value_type;
    const V1& _v1; const V2& _v2;
    OuterProduct<V1,V2>(const V1& v1, const V2& v2) : _v1(v1), _v2(v2) { }
    value_type operator() (Nat i, Nat j) { return _v1[i]*_v2[j]; }
};
template<class V1, class V2> inline
OuterProduct<V1,V2>
outer_product(const V1& v1, const V2& v2) {
    return OuterProduct<V1,V2>(v1,v2);
}

template<class X>
Matrix<X> outer(const Vector<X>& v1, const Vector<X>& v2) {
    Matrix<X> r(v1.size(),v2.size());
    for(Nat i1=0; i1!=v1.size(); ++i1) {
        for(Nat i2=0; i2!=v2.size(); ++i2) {
            r[i1][i2]=v1[i1]*v2[i2];
        }
    }
}

template<class X>
Matrix<X> outer(const Vector<X>& v) {
    Matrix<X> r(v.size(),v.size());
    for(Nat i1=0; i1!=v.size(); ++i1) {
        for(Nat i2=0; i2!=v.size(); ++i2) {
            r[i1][i2]=v[i1]*v[i2];
        }
    }
}

//! \ingroup DifferentiationModule
//! \brief A class representing the partial derivatives of a scalar quantity
//! depending on multiple arguments.
//!
//! Based on a power series Expansion, centred on the point at which the partial derivatives are
//! being evaluated.
//!
//! \invariant The expansion is sorted using graded_sort(). In particular, the
//! total degree of the terms is increasing, and the linear terms appear in coordinate order.
template<class X>
class FirstDifferential
{
  public:
    X _value;
    Covector<X> _gradient;
    static const X _zero;
  public:
    //! \brief The type of used to represent numbers in the coefficient.
    typedef typename X::NumericType NumericType;
    //! \brief The type of used to represent the coefficients.
    typedef X ValueType;

    //! \brief Constructs a differential in \a as variables.
    explicit FirstDifferential(Nat as) : _value(_zero), _gradient(as,_zero) { }

    //! \brief Constructs a differential with degree zero in \a as variables, initialising to the zero value z.
    explicit FirstDifferential(Nat as, const X& z) : _value(z), _gradient(as,z) { }

    //! \brief Constructs a differential with degree zero in \a as variables. (Deprecated)
    explicit FirstDifferential(const X& v, const Covector<X>& g) : _value(v), _gradient(g) { }

    //! \brief Conversion constructor from a different numerical type.
    template<class XX> FirstDifferential(const FirstDifferential<XX>& x)
        : _value(x._value), _gradient(x._gradient) { }

    //! \brief Set the differential equal to a constant, without changing the degree or number of arguments.
    FirstDifferential<X>& operator=(const X& c) { _value=c; for(Nat i=0; i!=_gradient.size(); ++i) { _gradient[i]=_zero; } return *this; }

    //! \brief A constant differential of degree \a deg in \a as arguments with value \a c.
    static FirstDifferential<X> constant(Nat as, const X& c) {
        FirstDifferential<X> r(as); r._value=c; return r; }
    //! \brief A differential of degree \a deg in \a as arguments representing the quantity \f$v+x_j\f$.
    static FirstDifferential<X> variable(Nat as, const X& v, Nat j) {
        FirstDifferential<X> r(as); r._value=v; r._gradient[j]=1; return r; }
    //! \brief \brief A vector of differentials of degree \a deg in \a as arguments with values \f$c_i+x_i\f$.
    static Vector< FirstDifferential<X> > variables(const Vector<X>& x) {
        Vector< FirstDifferential<X> > result(x.size(),FirstDifferential<X>(x.size()));
        for(Nat j=0; j!=x.size(); ++j) { result[j]._value=x[j]; result[j]._gradient[j]=1; }
        return result; }

    //! \brief Equality operator.
    Bool operator==(const FirstDifferential<X>& other) const {
        return this->_value==other._value && this->_gradient==other._gradient; }

    //! \brief Inequality operator.
    Bool operator!=(const FirstDifferential<X>& other) const {
        return !(*this==other); }

    //! \brief The number of independent variables.
    Nat argument_size() const { return this->_gradient.size(); }
    //! \brief The maximum degree of the stored terms.
    Nat degree() const { return 1u; }
    //! \brief The value of the differential i.e. the coefficient of \f$1\f$.
    const X& value() const { return this->_value; }
    //! \brief The vector of coefficients of \f$x_j\f$.
    const Covector<X>& gradient() const { return this->_gradient; }
    //! \brief The vector of coefficients of \f$x_j\f$.
    const X& gradient(Nat j) const { return this->_gradient[j]; }


    //! \brief A reference to the coefficient of \f$x_j\f$.
    X& operator[](const Nat& j) { return this->_gradient[j]; }
    //! \brief The coefficient of \f$x_j\f$.
    const X& operator[](const Nat& j) const { return this->_gradient[j]; }
    //! \brief The coefficient of \f$x^a\f$ i.e. \f$\prod x_j^{a_j}\f$.
    const X& operator[](const MultiIndex& a) const {
        ARIADNE_ASSERT_MSG(a.number_of_variables()==this->argument_size()," d="<<*this<<", a="<<a);
        if(a.degree()==0) { return this->_value; }
        else if(a.degree()==1) { for(Nat j=0; j!=this->argument_size(); ++j) { if(a[j]==1) { return this->_gradient[j]; } } }
        else { return _zero; } }

    //! \brief Set the coefficient of the constant term to \a c.
    Void set_value(const X& c) { this->_value=c; }

    //! \brief Set all coefficients to zero.
    Void clear() { *this=_zero; }

};

template<class X>
const X FirstDifferential<X>::_zero=X(0);


template<class X, class R>
EnableIfNumericType<R,FirstDifferential<X>&>
operator+=(FirstDifferential<X>& x, const R& c)
{
    x._value+=static_cast<X>(c);
    return x;
}

template<class X, class R>
EnableIfNumericType<R,FirstDifferential<X>&>
operator-=(FirstDifferential<X>& x, const R& c)
{
    x._value-=static_cast<X>(c);
    return x;
}

template<class X, class R>
EnableIfNumericType<R,FirstDifferential<X>&>
operator*=(FirstDifferential<X>& x, const R& c)
{
    x._value*=static_cast<X>(c);
    x._gradient*=static_cast<X>(c);
    return x;
}


template<class X, class R>
EnableIfNumericType<R,FirstDifferential<X>&>
operator/=(FirstDifferential<X>& x, const R& c)
{
    x._value/=static_cast<X>(c);
    x._gradient/=static_cast<X>(c);
    return x;
}

template<class X, class R>
EnableIfNumericType<R,FirstDifferential<X> >
operator+(const FirstDifferential<X>& x, const R& c)
{
    FirstDifferential<X> r(x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,FirstDifferential<X> >
operator+(const R& c, const FirstDifferential<X>& x)
{
    FirstDifferential<X> r(x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,FirstDifferential<X> >
operator-(const FirstDifferential<X>& x, const R& c)
{
    FirstDifferential<X> r(x); r-=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,FirstDifferential<X> >
operator-(const R& c, const FirstDifferential<X>& x)
{
    FirstDifferential<X> r(-x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,FirstDifferential<X> >
operator*(const FirstDifferential<X>& x, const R& c)
{
    FirstDifferential<X> r(x); r*=X(c); return r;
}

template<class X>
FirstDifferential<X>
operator*(const FirstDifferential<X>& x, const Int& c)
{
    FirstDifferential<X> r(x); r*=X(0); return r;
}

template<class X, class R>
EnableIfNumericType<R,FirstDifferential<X> >
operator*(const R& c, const FirstDifferential<X>& x)
{
    FirstDifferential<X> r(x); r*=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,FirstDifferential<X> >
operator/(const FirstDifferential<X>& x, const R& c)
{
    FirstDifferential<X> r(x); r/=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,FirstDifferential<X> >
operator/(const R& c, const FirstDifferential<X>& x)
{
    FirstDifferential<X> r(rec(x)); r*=X(c); return r;
}


template<class X>
FirstDifferential<X>& operator+=(FirstDifferential<X>& x, const FirstDifferential<X>& y)
{
    x._value += y._value;
    x._gradient += y._gradient;
    return x;
}

template<class X>
FirstDifferential<X>& operator-=(FirstDifferential<X>& x, const FirstDifferential<X>& y)
{
    x._value -= y._value;
    x._gradient -= y._gradient;
    return x;
}

template<class X>
FirstDifferential<X>& operator*=(FirstDifferential<X>& x, const FirstDifferential<X>& y)
{
    x._gradient *= y._value;
    x._gradient += x._value * y._gradient;
    x._value *= y._value;
    return x;
}

template<class X>
FirstDifferential<X>& operator/=(FirstDifferential<X>& x, const FirstDifferential<X>& y)
{
    x._value /= y._value;
    x._gradient -= x._value * y._gradient;
    x._gradient /= y._value;
    return x;
}


template<class X>
FirstDifferential<X> operator+(const FirstDifferential<X>& x)
{
    return FirstDifferential<X>(+x._value,+x._gradient);
}

template<class X>
FirstDifferential<X> operator-(const FirstDifferential<X>& x)
{
    return FirstDifferential<X>(-x._value,-x._gradient);
}


template<class X>
FirstDifferential<X> operator+(const FirstDifferential<X>& x, const FirstDifferential<X>& y)
{
    return FirstDifferential<X>(x._value+y._value,x._gradient+y._gradient);
}

template<class X>
FirstDifferential<X> operator-(const FirstDifferential<X>& x, const FirstDifferential<X>& y)
{
    return FirstDifferential<X>(x._value-y._value,x._gradient-y._gradient);
}

template<class X>
FirstDifferential<X> operator*(const FirstDifferential<X>& x, const FirstDifferential<X>& y)
{
    return FirstDifferential<X>(x._value*y._value,x._value*y._gradient+y._value*x._gradient);
}

template<class X>
FirstDifferential<X> operator/(const FirstDifferential<X>& x, const FirstDifferential<X>& y)
{
    return FirstDifferential<X>(x._value/y._value,((x._value/y._value)*y._gradient+x._gradient)/y._value);
}

template<class X> FirstDifferential<X> create_zero(const FirstDifferential<X>& c) {
    return FirstDifferential<X>(c.argument_size(),create_zero(c.value())); }








template<class X>
FirstDifferential<X>
min(const FirstDifferential<X>& x1, const FirstDifferential<X>& x2)
{
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"min(FirstDifferential<X> x1, FirstDifferential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()<x2.value() ? x1 : x2;
}


template<class X>
FirstDifferential<X>
max(const FirstDifferential<X>& x1,const FirstDifferential<X>& x2)
{
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"max(FirstDifferential<X> x1, FirstDifferential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()>x2.value() ? x1 : x2;
}

template<class X>
FirstDifferential<X>
abs(const FirstDifferential<X>& x)
{
    if(x.value()==0) {
        ARIADNE_THROW(std::runtime_error,"abs(FirstDifferential<X> x)","x[0]==0");
    }
    return x.value()>0 ? pos(x) : neg(x);
}


template<class X>
FirstDifferential<X>
pos(const FirstDifferential<X>& x)
{
    return x;
}

template<class X>
FirstDifferential<X>
neg(const FirstDifferential<X>& x)
{
    return -x;
}

template<class X>
FirstDifferential<X> rec(const FirstDifferential<X>& x)
{
    return FirstDifferential<X>( rec(x._value), x._gradient * (neg(sqr(rec(x._value)))) );
}

template<class X>
FirstDifferential<X> sqr(const FirstDifferential<X>& x)
{
    return FirstDifferential<X>( sqr(x._value), (2*x._value)*x._gradient );
}

template<class X>
FirstDifferential<X> pow(const FirstDifferential<X>& x, Int n)
{
    return FirstDifferential<X>( pow(x._value,n), (n*pow(x._value,n-1))*x._gradient );
}

template<class X>
FirstDifferential<X> sqrt(const FirstDifferential<X>& x)
{
    X sqrt_val = sqrt(x._value);
    return FirstDifferential<X>( sqrt_val, rec(2*sqrt_val)*x._gradient );
}

template<class X>
FirstDifferential<X> exp(const FirstDifferential<X>& x)
{
    X exp_val = exp(x._value);
    return FirstDifferential<X>( exp_val, exp_val*x._gradient );
}

template<class X>
FirstDifferential<X> log(const FirstDifferential<X>& x)
{
    X log_val = log(x._value);
    X rec_val = rec(x._value);
    return FirstDifferential<X>( log_val, rec_val*x._gradient );
}

template<class X>
FirstDifferential<X> sin(const FirstDifferential<X>& x)
{
    X sin_val = sin(x._value);
    X cos_val = cos(x._value);
    return FirstDifferential<X>( sin_val, cos_val*x._gradient );
}

template<class X>
FirstDifferential<X> cos(const FirstDifferential<X>& x)
{
    X cos_val = cos(x._value);
    X neg_sin_val = neg(sin(x._value));
    return FirstDifferential<X>( cos_val, neg_sin_val*x._gradient );
}

template<class X>
FirstDifferential<X> tan(const FirstDifferential<X>& x)
{
    X tan_val = tan(x._value);
    X sqr_sec_val = sqr(rec(cos(x._value)));
    return FirstDifferential<X>( tan_val, sqr_sec_val*x._gradient );
}

template<class X>
FirstDifferential<X> asin(const FirstDifferential<X>& x)
{
    X asin_val = asin(x._value);
    X d_asin_val = rec(sqrt(1.0-sqr(x._value)));
    return FirstDifferential<X>( asin_val, d_asin_val*x._gradient );
}

template<class X>
FirstDifferential<X> acos(const FirstDifferential<X>& x)
{
    X acos_val = acos(x._value);
    X d_acos_val = neg(rec(sqrt(1.0-sqr(x._value))));
    return FirstDifferential<X>( acos_val, d_acos_val*x._gradient );
}

template<class X>
FirstDifferential<X> atan(const FirstDifferential<X>& x)
{
    X atan_val = atan(x._value);
    X d_atan_val = rec(1+sqr(x._value));
    return FirstDifferential<X>( atan_val, d_atan_val*x._gradient );
}


template<class X>
Bool
operator>=(const FirstDifferential<X>& x, const FirstDifferential<X>& y)
{
    return x._value>=y._value;
}


template<class X>
Bool
operator<=(const FirstDifferential<X>& x, const FirstDifferential<X>& y)
{
    return x._value<=y._value;
}


template<class X, class R>
EnableIfNumericType<R,Bool>
operator>=(const FirstDifferential<X>& x, const R& c)
{
    return x._value>=static_cast<X>(c);
}


template<class X, class R>
EnableIfNumericType<R,Bool>
operator<=(const FirstDifferential<X>& x, const R& c)
{
    return x._value<=static_cast<X>(c);
}



template<class X>
OutputStream& operator<<(OutputStream& os, const FirstDifferential<X>& x)
{
    os << "D<R"<<x.argument_size()<<","<<x.degree()<<">{ ";
    os << x._value;
    for(Nat j=0; j!=x.argument_size(); ++j) {
        os << (j==0?"; ":",") << j << ":" << x[j];
    }
    return os << " }";
}






/*! \brief A class representing the derivatives of a vector quantity depending on multiple arguments. */
template<class X>
class Vector< FirstDifferential<X> >
    : public VectorExpression< Vector<FirstDifferential<X> > >
{
    //BOOST_CONCEPT_ASSERT((DifferentialVectorConcept<DifferentialVector<X> >));
    Array< FirstDifferential<X> > _ary;
  public:
    // The type of the class
    typedef Vector< FirstDifferential<X> > SelfType;
    // The type used for accessing elements
    typedef Nat IndexType;
    // The value stored in the vector.
    typedef FirstDifferential<X> ValueType;
    // The type used for scalars.
    typedef X ScalarType;

    Vector(Nat rs, Nat as) : _ary(rs, FirstDifferential<X>(as)) { }
    Vector(Nat rs, const FirstDifferential<X>& d) : _ary(rs,d) {
        for(Nat i=0; i!=rs; ++i) { (*this)[i]=d; } }
    Vector(Nat rs, const FirstDifferential<X>* p) : _ary(rs,FirstDifferential<X>(p[0].argument_size())) {
        for(Nat i=0; i!=rs; ++i) { (*this)[i]=p[i]; } }
    template<class E> Vector(const VectorExpression<E>& ve);
    template<class E> Vector< FirstDifferential<X> >& operator=(const VectorExpression<E>& ve);

    Nat result_size() const { return this->size(); }
    Nat argument_size() const { return (this->size()==0) ? 0 : (*this)[0].argument_size(); }
    Nat degree() const { return 1u; }

    const FirstDifferential<X>& get(Nat i) const { return _ary[i]; }
    FirstDifferential<X>& at(Nat i) { return _ary[i]; }
    Void set(Nat i, const FirstDifferential<X>& d) const { _ary[i]=d; }
    const FirstDifferential<X>& operator[](Nat i) const { return _ary[i]; }
    FirstDifferential<X>& operator[](Nat i) { return _ary[i]; }
    const FirstDifferential<X> zero_element() const { return FixedDifferential(this->argument_size()); }

    Vector<X> value() const {
        Vector<X> r(this->result_size());
        for(Nat i=0; i!=r.size(); ++i) { r[i]=(*this)[i].value(); } return r; }
    Matrix<X> jacobian() const { Matrix<X> r(this->result_size(),this->argument_size());
        for(Nat i=0; i!=r.row_size(); ++i) { for(Nat j=0; j!=r.column_size(); ++j) { r[i][j]=(*this)[i][j]; } } return r; }

};





//! \ingroup DifferentiationModule
//! \brief A class representing the partial derivatives of a scalar quantity
//! depending on multiple arguments.
//!
//! Based on a power series Expansion, centred on the point at which the partial derivatives are
//! being evaluated.
//!
//! \invariant The expansion is sorted using graded_sort(). In particular, the
//! total degree of the terms is increasing, and the linear terms appear in coordinate order.
template<class X>
class SecondDifferential
{
  public:
    X _value;
    Covector<X> _gradient;
    Matrix<X> _hessian;
    static const X _zero;
  public:
    //! \brief The type of used to represent numbers in the coefficient.
    typedef typename X::NumericType NumericType;
    //! \brief The type of used to represent the coefficients.
    typedef X ValueType;

    //! \brief Constructs a differential with degree zero in \a as variables. (Deprecated)
    explicit SecondDifferential(Nat as) : _value(_zero), _gradient(as), _hessian(as,as) { }

    //! \brief Constructs a differential with degree zero in \a as variables. (Deprecated)
    explicit SecondDifferential(const X& v, const Covector<X>& g, const Matrix<X>& h) : _value(v), _gradient(g), _hessian(h) {
        ARIADNE_ASSERT_MSG(h.row_size()==g.size() && h.column_size()==g.size(), "SecondDifferential(v,g,h): v="<<v<<", g="<<g<<", h="<<h<<"\n"); }

    //! \brief Conversion constructor from a different numerical type.
    template<class XX> SecondDifferential(const SecondDifferential<XX>& x)
        : _value(x._value), _gradient(x._gradient), _hessian(x._hessian) { }

    //! \brief Set the differential equal to a constant, without changing the degree or number of arguments.
    SecondDifferential<X>& operator=(const X& c) {
        _value=c; for(Nat j=0; j!=this->argument_size(); ++j) { _gradient[j]=_zero;
            for(Nat k=0; k!=this->argument_size(); ++k) { _hessian[j][k]=_zero; } } return *this; }

    //! \brief A constant differential of degree \a deg in \a as arguments with value \a c.
    static SecondDifferential<X> constant(Nat as, const X& c) {
        SecondDifferential<X> r(as); r._value=c; return r; }
    //! \brief A differential of degree \a deg in \a as arguments representing the quantity \f$v+x_j\f$.
    static SecondDifferential<X> variable(Nat as, const X& v, Nat j) {
        SecondDifferential<X> r(as); r._value=v; r._gradient[j]=1; return r; }
    //! \brief \brief A vector of differentials of degree \a deg in \a as arguments with values \f$c_i+x_i\f$.
    static Vector< SecondDifferential<X> > variables(const Vector<X>& x) {
        Vector< SecondDifferential<X> > result(x.size(),SecondDifferential<X>(x.size()));
        for(Nat j=0; j!=x.size(); ++j) { result[j]._value=x[j]; result[j]._gradient[j]=1; }
        return result; }

    //! \brief Equality operator.
    Bool operator==(const SecondDifferential<X>& other) const {
        return this->_value==other._value && this->_gradient==other._gradient  && this->_hessian==other._hessian; }

    //! \brief Inequality operator.
    Bool operator!=(const SecondDifferential<X>& other) const {
        return !(*this==other); }

    //! \brief The number of independent variables.
    Nat argument_size() const { return this->_gradient.size(); }
    //! \brief The maximum degree of the stored terms.
    Nat degree() const { return 2u; }
    //! \brief The value of the differential i.e. the coefficient of \f$1\f$.
    const X& value() const { return this->_value; }
    //! \brief The vector of coefficients of \f$x_j\f$.
    const Covector<X>& gradient() const { return this->_gradient; }
    //! \brief The Hessian matrix.
    //! \note Note the the components of the Hessian matrix are \em half those of the values indexed by the differential.
    //! This is because the differential stores the coefficients of the Taylor expansion, rather than the derivatives themselves.
    Matrix<X> hessian() const { return this->_hessian; }


    //! \brief A reference to the coefficient of \f$x_j\f$.
    X& operator[](const Nat& j) { return this->_gradient[j]; }
    //! \brief The coefficient of \f$x_j\f$.
    const X& operator[](const Nat& j) const { return this->_gradient[j]; }
    //! \brief The coefficient of \f$x^a\f$ i.e. \f$\prod x_j^{a_j}\f$.
    const X& operator[](const MultiIndex& a) const {
        ARIADNE_NOT_IMPLEMENTED; }
    //! \brief Set the coefficient of the constant term to \a c.
    Void set_value(const X& c) { this->_value=c; }

    //! \brief Set all coefficients to zero.
    Void clear() { *this=_zero; }

};

template<class X>
const X SecondDifferential<X>::_zero=X(0);




template<class X, class R>
EnableIfNumericType<R,SecondDifferential<X>&>
operator+=(SecondDifferential<X>& x, const R& c)
{
    x._value+=static_cast<X>(c);
}

template<class X, class R>
EnableIfNumericType<R,SecondDifferential<X>&>
operator-=(SecondDifferential<X>& x, const R& c)
{
    x._value-=static_cast<X>(c);
}

template<class X, class R>
EnableIfNumericType<R,SecondDifferential<X>&>
operator*=(SecondDifferential<X>& x, const R& c)
{
    x._value*=static_cast<X>(c);
    x._gradient*=static_cast<X>(c);
    x._hessian*=static_cast<X>(c);
}


template<class X, class R>
EnableIfNumericType<R,SecondDifferential<X>&>
operator/=(SecondDifferential<X>& x, const R& c)
{
    x._value/=static_cast<X>(c);
    x._gradient/=static_cast<X>(c);
    x._hessian/=static_cast<X>(c);
}

template<class X, class R>
EnableIfNumericType<R,SecondDifferential<X> >
operator+(const SecondDifferential<X>& x, const R& c)
{
    SecondDifferential<X> r(x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,SecondDifferential<X> >
operator+(const R& c, const SecondDifferential<X>& x)
{
    SecondDifferential<X> r(x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,SecondDifferential<X> >
operator-(const SecondDifferential<X>& x, const R& c)
{
    SecondDifferential<X> r(x); r-=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,SecondDifferential<X> >
operator-(const R& c, const SecondDifferential<X>& x)
{
    SecondDifferential<X> r(-x); r+=X(c); return r;
}

template<class X, class R>

EnableIfNumericType<R,SecondDifferential<X> >
operator*(const SecondDifferential<X>& x, const R& c)
{
    SecondDifferential<X> r(x); r*=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,SecondDifferential<X> >
operator*(const R& c, const SecondDifferential<X>& x)
{
    SecondDifferential<X> r(x); r*=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,SecondDifferential<X> >
operator/(const SecondDifferential<X>& x, const R& c)
{
    SecondDifferential<X> r(x); r/=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,SecondDifferential<X> >
operator/(const R& c, const SecondDifferential<X>& x)
{
    SecondDifferential<X> r(rec(x)); r*=X(c); return r;
}



template<class X>
SecondDifferential<X>& operator+=(SecondDifferential<X>& x, const SecondDifferential<X>& y)
{
    x._value += y._value;
    x._gradient += y._gradient;
    x._hessian += y._hessian;
    return x;
}

template<class X>
SecondDifferential<X>& operator-=(SecondDifferential<X>& x, const SecondDifferential<X>& y)
{
    x._value -= y._value;
    x._gradient -= y._gradient;
    x._hessian -= y._hessian;
    return x;
}

template<class X>
SecondDifferential<X>& operator*=(SecondDifferential<X>& x, const SecondDifferential<X>& y)
{
    x._hessian *= y._value;
    x._hessian += x._value * y._hessian;
    for(Nat j=0; j!=x.argument_size(); ++j) { for(Nat k=0; k!=x.argument_size(); ++k) {
            x._hessian[j][k]+=(x._gradient[j]*y._gradient[k]+x._gradient[k]*y._gradient[j])/2;
    } }
    x._gradient *= y._value;
    x._gradient += x._value * y._gradient;
    x._value *= y._value;
    return x;
}


template<class X>
SecondDifferential<X> operator+(const SecondDifferential<X>& x)
{
    return SecondDifferential<X>(+x._value,+x._gradient,+x._hessian);
}

template<class X>
SecondDifferential<X> operator-(const SecondDifferential<X>& x)
{
    return SecondDifferential<X>(-x._value,-x._gradient,-x._hessian);
}


template<class X>
SecondDifferential<X> operator+(const SecondDifferential<X>& x, const SecondDifferential<X>& y)
{
    return SecondDifferential<X>(x._value+y._value,x._gradient+y._gradient,x._hessian+y._hessian);
}

template<class X>
SecondDifferential<X> operator-(const SecondDifferential<X>& x, const SecondDifferential<X>& y)
{
    return SecondDifferential<X>(x._value-y._value,x._gradient-y._gradient,x._hessian-y._hessian);
}

template<class X>
SecondDifferential<X> operator*(const SecondDifferential<X>& x, const SecondDifferential<X>& y)
{
    Matrix<X> hessian(x.argument_size(),y.argument_size());
    for(Nat j=0; j!=x.argument_size(); ++j) { for(Nat k=0; k!=x.argument_size(); ++k) {
        hessian[j][k]=(x._gradient[j]*y._gradient[k]+x._gradient[k]*y._gradient[j])/2;
    } }
    hessian+=x._value*y._hessian;
    hessian+=y._value*x._hessian;
    return SecondDifferential<X>(x._value*y._value,x._value*y._gradient+y._value*x._gradient,hessian);
}

template<class X>
SecondDifferential<X> operator/(const SecondDifferential<X>& x, const SecondDifferential<X>& y)
{
    return x*rec(y);
}







template<class X>
SecondDifferential<X>
min(const SecondDifferential<X>& x1, const SecondDifferential<X>& x2)
{
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"min(SecondDifferential<X> x1, SecondDifferential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()<x2.value() ? x1 : x2;
}


template<class X>
SecondDifferential<X>
max(const SecondDifferential<X>& x1,const SecondDifferential<X>& x2)
{
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"max(SecondDifferential<X> x1, SecondDifferential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()>x2.value() ? x1 : x2;
}

template<class X>
SecondDifferential<X>
abs(const SecondDifferential<X>& x)
{
    if(x.value()==0) {
        ARIADNE_THROW(std::runtime_error,"abs(SecondDifferential<X> x)","x[0]==0");
    }
    return x.value()>0 ? pos(x) : neg(x);
}

template<class X>
SecondDifferential<X>
pos(const SecondDifferential<X>& x)
{
    return x;
}

template<class X>
SecondDifferential<X>
neg(const SecondDifferential<X>& x)
{
    return -x;
}

template<class X>
SecondDifferential<X> rec(const SecondDifferential<X>& x)
{
    //return SecondDifferential<X>( rec(x._value), x._gradient * (neg(sqr(rec(x._value)))) );
}

template<class X>
SecondDifferential<X> sqr(const SecondDifferential<X>& x)
{
    return SecondDifferential<X>( sqr(x._value), (2*x._value)*x._gradient, (2*x._value)*x._hessian + outer(x._gradient) );
}

template<class X>
SecondDifferential<X> pow(const SecondDifferential<X>& x, Int n)
{
    return SecondDifferential<X>( pow(x._value,n), (n*pow(x._value,n-1))*x._gradient, n*(n-1)*pow(x._value,n-2)*x._hessian );
}

//ddf(y)/dxdx = d/dx ( f'(y) dy/dx) = f''(y) dy/dx dy/dx + f'(y) ddy/dxdx

template<class X>
SecondDifferential<X> sqrt(const SecondDifferential<X>& x)
{
    X sqrt_val = sqrt(x._value);
    X rec_dbl_sqrt_val = rec(2*sqrt_val);
    X neg_rec_quad_pow_val = neg(rec(4*sqrt_val*x._value));
    return SecondDifferential<X>( sqrt_val, rec_dbl_sqrt_val*x._gradient, rec_dbl_sqrt_val*x._hessian + neg_rec_quad_pow_val * outer(x._gradient) );
}

template<class X>
SecondDifferential<X> exp(const SecondDifferential<X>& x)
{
    X exp_val = exp(x._value);
    return SecondDifferential<X>( exp_val, exp_val*x._gradient, exp_val*x._hessian+exp_val*outer(x._gradient) );
}

template<class X>
SecondDifferential<X> log(const SecondDifferential<X>& x)
{
    X log_val = log(x._value);
    X rec_val = rec(x._value);
    X neg_sqr_rec_val = neg(sqr(rec_val));
    return SecondDifferential<X>( log_val, rec_val*x._gradient, rec_val*x._hessian+neg_sqr_rec_val*outer(x._gradient) );
}

template<class X>
SecondDifferential<X> sin(const SecondDifferential<X>& x)
{
    X sin_val = sin(x._value);
    X cos_val = cos(x._value);
    X neg_sin_val = neg(sin_val);
    return SecondDifferential<X>( sin_val, cos_val*x._gradient, cos_val*x._hessian+neg_sin_val*outer(x._gradient) );
}

template<class X>
SecondDifferential<X> cos(const SecondDifferential<X>& x)
{
    X cos_val = cos(x._value);
    X neg_sin_val = neg(sin(x._value));
    X neg_cos_val = neg(cos_val);
    return SecondDifferential<X>( cos_val, neg_sin_val*x._gradient, neg_sin_val*x._hessian+neg_cos_val*outer(x._gradient) );
}

template<class X>
SecondDifferential<X> tan(const SecondDifferential<X>& x)
{
    X tan_val = tan(x._value);
    X sqr_sec_val = sqr(rec(cos(x._value)));
    X dbl_tan_sqr_sec_val = 2*tan_val*sqr_sec_val;
    return SecondDifferential<X>( tan_val, sqr_sec_val*x._gradient, sqr_sec_val*x._hessian+dbl_tan_sqr_sec_val*outer(x._gradient) );
}




template<class X>
OutputStream& operator<<(OutputStream& os, const SecondDifferential<X>& x)
{
    Expansion<X> e=x.expansion();
    //e.graded_sort();
    os << "D<R"<<x.argument_size()<<","<<x.degree()<<">{";
    os << x._value;
    for(Nat j=0; j!=x.argument_size(); ++j) {
        os << "," << j << ":" << x[j];
    }
    for(Nat j=0; j!=x.argument_size(); ++j) {
        for(Nat k=0; k!=x.argument_size(); ++k) {
            os << "," << j << "," << k << ":" << x._hessian[j][k];
        }
    }
    return os << " }";
}






/*! \brief A class representing the derivatives of a vector quantity depending on multiple arguments. */
template<class X>
class Vector< SecondDifferential<X> >
{
    //BOOST_CONCEPT_ASSERT((DifferentialVectorConcept<DifferentialVector<X> >));
    Array< SecondDifferential<X> > _ary;
  public:
    // The type of the class
    typedef Vector< SecondDifferential<X> > SelfType;
    // The type used for accessing elements
    typedef Nat IndexType;
    // The value stored in the vector.
    typedef SecondDifferential<X> ValueType;
    // The type used for scalars.
    typedef X ScalarType;

    Vector(Nat rs, Nat as) : _ary(rs, SecondDifferential<X>(as)) { }
    Vector(Nat rs, const SecondDifferential<X>& d) : _ary(rs) {
        for(Nat i=0; i!=rs; ++i) { (*this)[i]=d; } }
    Vector(Nat rs, const SecondDifferential<X>* p) : _ary(rs) {
        for(Nat i=0; i!=rs; ++i) { (*this)[i]=p[i]; } }
    template<class E> Vector(const VectorExpression<E>& ve);
    template<class E> Vector< SecondDifferential<X> >& operator=(const VectorExpression<E>& ve);


    Nat result_size() const { return this->size(); }
    Nat argument_size() const { return (this->size()==0) ? 0 : (*this)[0].argument_size(); }
    Nat degree() const { return 1u; }

    const FirstDifferential<X>& get(Nat i) const { return _ary[i]; }
    FirstDifferential<X>& at(Nat i) { return _ary[i]; }
    Void set(Nat i, const FirstDifferential<X>& d) const { _ary[i]=d; }
    const FirstDifferential<X>& operator[](Nat i) const { return _ary[i]; }
    FirstDifferential<X>& operator[](Nat i) { return _ary[i]; }

    Vector<X> value() const {
        Vector<X> r(this->result_size());
        for(Nat i=0; i!=r.size(); ++i) { r[i]=(*this)[i].value(); } return r; }
    Matrix<X> jacobian() const { Matrix<X> r(this->result_size(),this->argument_size());
        for(Nat i=0; i!=r.row_size(); ++i) { for(Nat j=0; j!=r.column_size(); ++j) { r[i][j]=(*this)[i].gradient(j); } } return r; }

};





} //namespace Ariadne

#endif // ARIADNE_FIXED_DIFFERENTIAL_HPP
