/***************************************************************************
 *            vector.h
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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

/*! \file vector.h
 *  \brief Vectors in Euclidean space.
 */

#ifndef ARIADNE_VECTOR_H
#define ARIADNE_VECTOR_H

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <string>
#include <sstream>
#include <cstdarg>

#include "macros.h"
#include "numeric.h"
#include "stlio.h"

using namespace boost::numeric;

namespace Ariadne {

//! \brief A vector over a field. See also \link Ariadne::Matrix \c Matrix<X> \endlink.
template<class X>
class Vector
    : public ublas::vector<X>
{
  public:
    //@{
    //! \name Constructors

    //! \brief Default constructor constructs a vector with no elements.
    Vector()
        : ublas::vector<X>() { }
    //! \brief Construct a vector of size \a n, with elements initialised to zero.
    Vector(size_t n)
        : ublas::vector<X>(n) { for(size_t i=0; i!=this->size(); ++i) { (*this)[i]=0; } }
    //! \brief Construct a vector of size \a n, with elements initialised to \a t.
    Vector(size_t n, const X& t)
        : ublas::vector<X>(n) { for(size_t i=0; i!=this->size(); ++i) { (*this)[i]=t; } }
    //! \brief Construct a vector of size \a n, with values initialised from the C-style array beginning at \a ptr.
    template<class XX> Vector(size_t n, const XX* ptr)
        : ublas::vector<X>(n) { for(size_t i=0; i!=this->size(); ++i) { (*this)[i]=ptr[i]; } }
    //! \brief Construct a vector of size \a n, with values initialised from a variadic argument list. WARNING: The values in the list must all be double-precision type; in particular, constants must be floating-point values \c 2.0 rather integer values \c 2 .
    Vector(size_t n, const double& t0, const double& t1, ...);
    //! \brief Construct a matrix from a string literal, with entries enclosed in square braces and separated by commass. e.g. <tt>"[1, 2.3, 4.2]"</tt>.
    Vector(const std::string& str)
        : ublas::vector<X>() { std::stringstream ss(str); ss >> *this; }
    //! \brief Copy constructor allows conversion from a vector using another numerical type.
    template<class XX> Vector(const Vector<XX>& v)
        : ublas::vector<X>(v) { }
#ifdef DOXYGEN
    //! \brief Copy assignement allows conversion from a vector using another numerical type.
    template<class XX> Vector<X>& operator=(const Vector<XX> &v);
#endif
    template<class E> Vector(const ublas::vector_expression<E> &ve)
        : ublas::vector<X>(ve) { }
    template<class E> Vector<X>& operator=(const ublas::vector_expression<E> &ve) {
        this->ublas::vector<X>::operator=(ve); return *this; }
    //@}

    //@{
    //! \name Static constructors

    //! \brief The zero vector of size \a n.
    static Vector<X> zero(size_t n) { return Vector<Float>(n,0.0); }
    //! \brief The vector of size \a n with all entries equal to one.
    static Vector<X> one(size_t n) { return Vector<Float>(n,1.0); }
    //! \brief The unit vector \f$e_i\f$ with value one in the \a i<sup>th</sup> entry, and zero otherwise.
    static Vector<X> unit(size_t n,size_t i) {
        ARIADNE_ASSERT(i<n); Vector<X> result(n,0.0); result[i]=1.0; return result; }
    static Vector<X> unit_box(size_t n) {
        Vector<X> result(n,Interval(-1,1)); return result; }
    //! \brief The unit vector \f$e_i\f$ with value one in the \a i<sup>th</sup> entry, and zero otherwise.
    static array< Vector<X> > basis(size_t n) {
        array< Vector<X> > result(n); for(uint i=0; i!=n; ++i) { result[i]=unit(n,i); } return result; }
    //@}

    //@{
    //! \name Data access

#ifdef DOXYGEN
    //! \brief The number of elements of the vector.
    size_t size() const;
#endif
    //! \brief Get the value stored in the \a i<sup>th</sup> element.
    const X& get(size_t i) const { return (*this)[i]; }
    //! \brief Set the value stored in the \a i<sup>th</sup> element to \a x.
    template<class T> void set(size_t i, const T& x) { (*this)[i] = x; }
#ifdef DOXYGEN
    //! \brief C-style subscripting operator.
    X& operator[](size_t i);
    //! \brief C-style constant subscripting operator.
    const X& operator[](size_t i) const;
#endif
    //@}

#ifdef DOXYGEN
    //! \brief Equality operator.
    friend template<class X1, class X2> bool operator==(const Vector<X1>& v1, const Vector<X2>& v2);
    //! \brief Inequality operator.
    friend template<class X1, class X2> bool operator!=(const Vector<X1>& v1, const Vector<X2>& v2);

     //! \brief %Vector negation.
    friend template<class X> Vector<X> operator-(const Vector<X>& v);
    //! \brief %Vector addition.
    friend template<class X> Vector<X> operator+(const Vector<X>& v1, const Vector<X>& v2);
    //! \brief %Vector subtraction.
    friend template<class X> Vector<X> operator-(const Vector<X>& v1, const Vector<X>& v2);
    //! \brief %Scalar multiplication.
    friend template<class X> Vector<X> operator*(const X& s, const Vector<X>& v);
    //! \brief %Scalar multiplication.
    friend template<class X> Vector<X> operator*(const Vector<X>& v, const X& s);
    //! \brief %Scalar division.
    friend template<class X> Vector<X> operator/(const Vector<X>& v, const X& s);

    //! \brief The supremum norm.
    friend template<class X> X norm(const Vector<X>& v);
    //! \brief The inner product.
    friend template<class X> X dot(const Vector<X>& v1, const Vector<X>& v2);

    //! \brief Join (catenate, make the direct sum of) two vectors.
    friend template<class X> Vector<X> join(const Vector<X>& v1, const Vector<X>& v2);
    //! \brief Join a vector and a scalar.
    friend template<class X> Vector<X> join(const Vector<X>& v1, const X& s2);
    //! \brief Join a scalar and a vector.
    friend template<class X> Vector<X> join(const X& s1, const Vector<X>& v2);
    //! \brief Join two scalars.
    friend template<class X> Vector<X> join(const X& s1, const X& s2);

    //! \brief Write to an output stream.
    friend template<class X> std::ostream& operator<<(std::ostream& os, const Vector<X>& v);
    //! \brief Read from an output stream.
    friend template<class X> std::istream& operator>>(std::istream& is, Vector<X>& v);
#endif // DOXYGEN
};

typedef ublas::slice Slice;
typedef ublas::range Range;

inline Range range(uint start, uint stop) { return Range(start,stop); }
inline Slice slice(uint size, uint start, uint stride) { return Slice(size,start,stride); }


template<class X>
X sup_norm(const Vector<X>& v)
{
    X r=0;
    for(size_t i=0; i!=v.size(); ++i) {
        r=max(r,static_cast<X>(abs(v[i])));  // Need static cast for expression templates
    }
    return r;
}

template<class X>
X norm(const Vector<X>& v)
{
    return Ariadne::sup_norm(v);
}

template<class X>
X dot(const Vector<X>& v1, const Vector<X>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    X r=0;
    for(size_t i=0; i!=v1.size(); ++i) {
        r+=v1[i]*v2[i];
    }
    return r;
}

template<class X>
Vector<X> join(const Vector<X>& v1, const Vector<X>& v2)
{
    size_t n1=v1.size();
    size_t n2=v2.size();
    Vector<X> r(n1+n2);
    ublas::project(r,range(0,n1))=v1;
    ublas::project(r,range(n1,n1+n2))=v2;
    return r;
}

template<class X>
Vector<X> join(const Vector<X>& v1, const X& s2)
{
    size_t n1=v1.size();
    Vector<X> r(n1+1);
    ublas::project(r,range(0,n1))=v1;
    r[n1]=s2;
    return r;
}

template<class X>
Vector<X> join(const X& s1, const Vector<X>& v2)
{
    size_t n2=v2.size();
    Vector<X> r(1+n2);
    r[0]=s1;
    ublas::project(r,range(1,n2+1))=v2;
    return r;
}

template<class X>
Vector<X> join(const X& s1, const X& s2)
{
    Vector<X> r(2);
    r[0]=s1;
    r[1]=s2;
    return r;
}



template<class X1, class X2>
bool operator==(const Vector<X1>& v1, const Vector<X2>& v2)
{
    if(v1.size()!=v2.size()) { return false; }
    for(size_t i=0; i!=v1.size(); ++i) {
        if(v1[i]!=v2[i]) { return false; }
    }
    return true;
}


template<class X1, class X2> inline
bool operator!=(const Vector<X1>& v1, const Vector<X2>& v2)
{
    return !(v1==v2);
}


template<class X1, class X2>
bool operator<(const Vector<X1>& v1, const Vector<X2>& v2)
{
    if(v1.size()!=v2.size()) { return false; }
    for(size_t i=0; i!=v1.size(); ++i) {
        if(v1[i]<v2[i]) { return true; }
        else if(v1[i]>v2[i]) { return false; }
    }
    return true;
}


template<class X>
bool operator<=(const Vector<X>& v, const X& c)
{
    for(size_t i=0; i!=v.size(); ++i) {
        if(v[i]>c) { return false; }
    }
    return true;
}


template<class X> std::ostream& operator<<(std::ostream& os, const Vector<X>& v) {
    if(v.size()==0) { os << '['; }
    for(size_t i=0; i!=v.size(); ++i) {
        os << (i==0 ? '[' : ',') << v[i]; }
    return os << ']';
}

template<class X> std::istream& operator>>(std::istream& is, Vector<X>& v) {
    std::vector<X> vec;
    is >> vec;
    X* ptr=&vec[0];
    v=Vector<X>(vec.size(),ptr);
    return is;
}

bool contains(const Vector<Interval>& v1, const Vector<Float>& v2);
bool intersect(const Vector<Interval>& v1, const Vector<Interval>& v2);
bool subset(const Vector<Interval>& v1, const Vector<Interval>& v2);

bool disjoint(const Vector<Interval>& v1, const Vector<Interval>& v2);
bool overlap(const Vector<Interval>& v1, const Vector<Interval>& v2);
bool inside(const Vector<Interval>& v1, const Vector<Interval>& v2);
bool covers(const Vector<Interval>& v1, const Vector<Interval>& v2);
bool empty(const Vector<Interval>& v);

Vector<Interval> hull(const Vector<Interval>& v1, const Vector<Interval>& v2);
Vector<Interval> intersection(const Vector<Interval>& v1, const Vector<Interval>& v2);
Vector<Float> midpoint(const Vector<Interval>& v);
Vector<Float> lower(const Vector<Interval>& v);
Vector<Float> upper(const Vector<Interval>& v);
Float radius(const Vector<Interval>& z);
Float volume(const Vector<Interval>& z);

inline Vector<Float> add_approx(const Vector<Float>& v1, const Vector<Float>& v2) { return v1+v2; }
inline Vector<Float> sub_approx(const Vector<Float>& v1, const Vector<Float>& v2) { return v1-v2; }



template<class X, class XX> inline Vector<X> vector(size_t d, const XX* ptr) {
    return Vector<Float>(d,ptr); }
inline Vector<Float> point(size_t d, Float* ptr) {
    return Vector<Float>(d,ptr); }
inline Vector<Interval> box(size_t d, Float* ptr) {
    Vector<Interval> bx(d);
    for(size_t i=0; i!=d; ++i) {
        bx[i]=Interval(ptr[2*i],ptr[2*i+1]); }
    return bx;
}



} // namespace Ariadne

#endif // ARIADNE_VECTOR_H

