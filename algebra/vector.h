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

#include <string>
#include <sstream>
#include <initializer_list>
#include <vector>

#include "utility/macros.h"
#include "utility/stlio.h"
#include "utility/metaprogramming.h"
#include "numeric/numeric.h"
#include "utility/array.h"

#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>

namespace Ariadne {

template<class V> struct VectorExpression {
    const V& operator()() const { return static_cast<const V&>(*this); }
};

template<class V> struct VectorContainer : public VectorExpression<V> { };

template<class T> struct IsScalar : IsNumeric<T> { };


//! \ingroup LinearAlgebraModule
//! \brief A vector over a field. See also \link Ariadne::Matrix \c Matrix<X> \endlink.
//!
//! \par Python interface
//!
//! In the Python interface, classes \c RawFloatVector and \c ExactIntervalVector are defined.
//! Further, Ariadne vectors can be constructed from literals in the form of a Python list: <br><br>
//! <code> RawFloatVector([1.1,2.3,4.2,5]) # Create a RawFloatVector from a list of Python \c int and \c float types. <br>
//!        ExactIntervalVector([{1:2.1},[-3,4],2.3,5,{-1.1:2.2}]) # Create an ExactIntervalVector from a list of Python types convertible to Interal</code>
template<class X>
class Vector
    : public VectorContainer< Vector<X> >
{
    Array<X> _ary;
  public:
    //@{
    //! \name Type definitions

    //! \brief The type of the scalar element.
    typedef X ValueType;

    //@}

    //@{
    //! \name Constructors

    //! \brief Default constructor constructs a vector with no elements.
    Vector()
        : _ary() { }
    //! \brief Construct a vector of size \a n, with elements initialised to zero.
    explicit Vector(size_t n) : _ary(n) { for(size_t i=0; i!=this->size(); ++i) { (*this)[i]=0; } }
    //! \brief Construct a vector of size \a n, with elements initialised to \a t.
    Vector(size_t n, const X& t) : _ary(n,t) {  }
    //! \brief Construct a vector of size \a n, with values initialised from the C-style Array beginning at \a ptr.
    template<class XX> Vector(size_t n, const XX* ptr) : _ary(ptr,ptr+n) { }
    //! \brief Convert from an initializer list of the same type.
    Vector(std::initializer_list<X> lst) : _ary(lst.begin(),lst.end()) { }
    template<class XX, typename std::enable_if<std::is_constructible<X,XX>::value and not std::is_convertible<XX,X>::value,int>::type=0>
        explicit Vector(std::initializer_list<XX> lst) : _ary(lst.begin(),lst.end()) { }
    //! \brief Convert from a list of the same type.
    Vector(const std::vector<X>& lst) : _ary(lst.begin(),lst.end()) { }
    //! \brief Construct from a list of a possibly different type.
    template<class XX> explicit Vector(const std::vector<XX>& lst) : _ary(lst.begin(),lst.end()) { }
    //! \brief Construct a matrix from a string literal, with entries enclosed in square braces and separated by commass. e.g. <tt>"[1, 2.3, 4.2]"</tt>.
    explicit Vector(const std::string& str)
        : _ary() { std::stringstream ss(str); ss >> *this; }
     //! \brief Copy constructor.
    Vector(const Vector<X>& v)
        : _ary(v.size()) { for(size_t i=0; i!=this->size(); ++i) { this->_ary[i]=v[i]; } }
    //! \brief Copy assignment.
    Vector<X>& operator=(const Vector<X>& v) {
        if(this!=&v) { this->_ary = v._ary; } return *this; }
#ifdef DOXYGEN
     //! \brief Copy constructor allows conversion from a vector using another numerical type.
    template<class XX> Vector(const Vector<XX>& v)
        : _ary(v.size()) { for(size_t i=0; i!=this->size(); ++i) { this->_ary[i]=v[i]; } }
   //! \brief Copy assignement allows conversion from a vector using another numerical type.
    template<class XX> Vector<X>& operator=(const Vector<XX> &v);
#endif
    template<class XX, typename std::enable_if<std::is_convertible<XX,X>::value,int>::type=0>
        Vector(const Vector<XX>& v) : _ary(v.size()) { for(size_t i=0; i!=this->size(); ++i) { this->_ary[i]=v[i]; } }
    template<class XX, typename std::enable_if<std::is_constructible<X,XX>::value and not std::is_convertible<XX,X>::value,int>::type=0>
        explicit Vector(const Vector<XX>& v) : _ary(v.size()) { for(size_t i=0; i!=this->size(); ++i) { this->_ary[i]=X(v[i]); } }
    template<class V> Vector(const VectorExpression<V>& ve, typename EnableIf<Or<IsSame<typename V::ValueType,X>, IsNumericCastable<typename V::ValueType,X> >,Void>::Type* = 0)
        : _ary(ve().size()) { for(size_t i=0; i!=this->size(); ++i) { this->_ary[i]=static_cast<X>(ve()[i]); } }
    template<class V> typename EnableIf<IsSafelyConvertible<typename V::ValueType,X>,Vector<X>&>::Type operator=(const VectorExpression<V>& ve) {
        this->resize(ve().size()); for(size_t i=0; i!=this->size(); ++i) { this->_ary[i]=ve()[i]; } return *this; }
    //@}

    //@{
    //! \name Static constructors

    //! \brief The zero vector of size \a n.
    static Vector<X> zero(size_t n) { return Vector<Float>(n,0.0); }
    //! \brief The vector of size \a n with all entries equal to one.
    static Vector<X> one(size_t n) { return Vector<Float>(n,1.0); }
    //! \brief The unit vector \f$e_i\f$ with value one in the \a i<sup>th</sup> entry, and zero otherwise.
    static Vector<X> unit(size_t n,size_t i) {
        ARIADNE_ASSERT(i<n); Vector<X> result(n,static_cast<X>(0.0)); result[i]=static_cast<X>(1.0); return result; }
    static Vector<X> unit_box(size_t n) {
        Vector<X> result(n,ExactInterval(-1,1)); return result; }
    //! \brief The unit vector \f$e_i\f$ with value one in the \a i<sup>th</sup> entry, and zero otherwise.
    static Array< Vector<X> > basis(size_t n) {
        Array< Vector<X> > result(n); for(uint i=0; i!=n; ++i) { result[i]=unit(n,i); } return result; }
    //@}

    //@{
    //! \name Data access

    //! \brief Resize to hold \a n elements.
    //! The previous values need not be preserved, and the new values need not be initialised.
    void resize(size_t n) { this->_ary.resize(n); }
    //! \brief The number of elements of the vector.
    size_t size() const { return this->_ary.size(); }
    //! \brief A reference to the value stored in the \a i<sup>th</sup> element.
    X& at(size_t i) { ARIADNE_PRECONDITION_MSG(i<this->size(),*this<<"["<<i<<"]"); return (*this)[i]; }
    //! \brief Get the value stored in the \a i<sup>th</sup> element.
    const X& get(size_t i) const { ARIADNE_PRECONDITION_MSG(i<this->size(),*this<<"["<<i<<"]"); return (*this)[i]; }
    //! \brief Set the value stored in the \a i<sup>th</sup> element to \a x.
    template<class T> void set(size_t i, const T& x) { ARIADNE_PRECONDITION_MSG(i<this->size(),*this<<"["<<i<<"]"); (*this)[i] = x; }
    //! \brief C-style subscripting operator.
    X& operator[](size_t i) { ARIADNE_PRECONDITION_MSG(i<this->size(),*this<<"["<<i<<"]"); return this->_ary[i]; }
    //! \brief C-style constant subscripting operator.
    const X& operator[](size_t i) const { ARIADNE_PRECONDITION_MSG(i<this->size(),*this<<"["<<i<<"]"); return this->_ary[i]; }
    //! \brief The zero of the ring containing the Vector's elements. This may be dependent on class parameters.
    const X zero_element() const { if(this->size()!=0) { return (*this)[0]*0; } else { return X(); } }
    //@}

#ifdef DOXYGEN
    //! \brief Equality operator.
    friend template<class X1, class X2> bool operator==(const Vector<X1>& v1, const Vector<X2>& v2);
    //! \brief Inequality operator.
    friend template<class X1, class X2> bool operator!=(const Vector<X1>& v1, const Vector<X2>& v2);

     //! \brief %Vector unary plus.
    friend template<class X> Vector<X> operator+(const Vector<X>& v);
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
    //! \brief Join (catenate, make the direct sum of) three vectors.
    friend template<class X> Vector<X> join(const Vector<X>& v1, const Vector<X>& v2, const Vector<X>& v3);
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

template<class X> std::ostream& operator<<(std::ostream& os, const Vector<X>& v);
template<class X> std::istream& operator>>(std::istream& is, Vector<X>& v);

template<class V> std::ostream& operator<<(std::ostream& os, const VectorExpression<V>& ve);


class Range {
    size_t _start; size_t _stop;
  public:
    Range(size_t start, size_t stop) : _start(start), _stop(stop) { }
    size_t size() const { return this->_stop-this->_start; }
    size_t start() const { return this->_start; }
    size_t stride() const { return 1u; }
    size_t stop() const { return this->_stop; }
};

class Slice {
    size_t _size; size_t _start; size_t _stride;
  public:
    Slice(size_t size, size_t start, size_t stride) : _size(size), _start(start), _stride(stride) { }
    size_t size() const { return this->_size; }
    size_t start() const { return this->_start; }
    size_t stride() const { return this->_stride; }
    size_t stop() const { return this->_start+this->_size*this->_stride; }
};

inline Range range(size_t start, size_t stop) { return Range(start,stop); }
inline Slice slice(size_t size, size_t start, size_t stride) { return Slice(size,start,stride); }

template<class V> struct VectorRange
    : public VectorContainer< VectorRange<V> >
{
    typedef typename V::ValueType ValueType;
    const V& _v; Range _rng;
    VectorRange(const V& v, Range rng) : _v(v), _rng(rng) { }
    size_t size() const { return _rng.size(); }
    ValueType operator[](size_t i) const { return _v[i+_rng.start()]; }
};

template<class V> struct VectorContainerRange
    : public VectorExpression< VectorContainerRange<V> >
{
    typedef typename V::ValueType ValueType;
    V& _v; Range _rng;
    VectorContainerRange(V& v, Range rng) : _v(v), _rng(rng) { }
    size_t size() const { return _rng.size(); }
    ValueType operator[](size_t i) const { return _v[i+_rng.start()]; }
    ValueType& operator[](size_t i) { return _v[i+_rng.start()]; }
    const ValueType zero_element() const { return _v.zero_element(); }
    void set(size_t i, const ValueType& x) { _v[i+_rng.start()]=x; }
    template<class VE> VectorContainerRange<V>& operator=(const VectorExpression<VE>& ve) {
        ARIADNE_PRECONDITION(this->size()==ve().size());
        for(size_t i=0; i!=this->size(); ++i) { (*this)[i]=ve()[i]; } return *this; }
};

template<class V> struct VectorNegation
    : public VectorExpression< VectorNegation<V> >
{
    const V& _v;
    VectorNegation(const V& v) : _v(v) { }
    typedef typename V::ValueType ValueType;
    size_t size() const { return _v.size(); }
    ValueType operator[](size_t i) const { return -_v[i]; }
    const ValueType zero_element() const { return _v.zero_element(); }
};

template<class V1, class V2> struct VectorSum
    : public VectorExpression< VectorSum<V1,V2> >
{
    const V1& _v1; const V2& _v2;
    VectorSum(const V1& v1, const V2& v2) : _v1(v1), _v2(v2) { }
    typedef typename Arithmetic<typename V1::ValueType, typename V2::ValueType>::ResultType ValueType;
    size_t size() const { return _v1.size(); }
    ValueType operator[](size_t i) const { return _v1[i]+_v2[i]; }
    const ValueType zero_element() const { return _v1.zero_element(); }
};

template<class V1, class V2> struct VectorDifference
    : public VectorExpression< VectorDifference<V1,V2> >
{
    const V1& _v1; const V2& _v2;
    VectorDifference(const V1& v1, const V2& v2) : _v1(v1), _v2(v2) { }
    typedef typename Arithmetic<typename V1::ValueType, typename V2::ValueType>::ResultType ValueType;
    size_t size() const { return _v1.size(); }
    ValueType operator[](size_t i) const { return _v1[i]-_v2[i]; }
    const ValueType zero_element() const { return _v1.zero_element(); }
};

template<class V1, class X2> struct VectorScalarProduct
    : public VectorExpression< VectorScalarProduct<V1,X2> >
{
    const V1& _v1; const X2& _x2;
    VectorScalarProduct(const V1& v1, const X2& x2) : _v1(v1), _x2(x2) { }
    typedef typename Arithmetic<typename V1::ValueType, X2>::ResultType ValueType;
    size_t size() const { return _v1.size(); }
    ValueType operator[](size_t i) const { return _v1[i]*_x2; }
    const ValueType zero_element() const { return ValueType(_v1.zero_element()*_x2); }
};

template<class V1, class X2> struct VectorScalarQuotient
    : public VectorExpression< VectorScalarQuotient<V1,X2> >
{
    const V1& _v1; const X2& _x2;
    VectorScalarQuotient(const V1& v1, const X2& x2) : _v1(v1), _x2(x2) { }
    typedef typename Arithmetic<typename V1::ValueType, X2>::ResultType ValueType;
    size_t size() const { return _v1.size(); }
    ValueType operator[](size_t i) const { return _v1[i]/_x2; }
    const ValueType zero_element() const { return ValueType(_v1.zero_element()/_x2); }
};

template<class V> inline VectorRange<V> project(const VectorExpression<V>& v, Range rng) { return VectorRange<V>(v(),rng); }
template<class X> inline VectorContainerRange< Vector<X> > project(Vector<X>& v, Range rng) { return VectorContainerRange< Vector<X> >(v,rng); }

template<class V> inline
const V&
operator+(const VectorExpression<V>& ve) { return ve(); }

template<class V> inline
VectorNegation<V>
operator-(const VectorExpression<V>& ve) { return VectorNegation<V>(ve()); }

// The code below is simpler, but illegal arithmetical operations are not caught until later
// template<class V1, class V2>
// VectorSum< V1, V2 >
// operator+(const VectorExpression<V1>& v1, const VectorExpression<V2>& v2) { return VectorSum<V1,V2>(v1(),v2()); }


// The code below does not require the use of VectorExpression.
//template<class V1, class V2> inline
//typename EnableIf< And< IsVector<V1>, IsVector<V2>, IsDefined<typename Arithmetic<typename V1::ValueType,typename V2::ValueType>::ResultType> >,
//                   VectorSum< V1, V2 > >::Type
//operator+(const V1& v1, const V2& v2) { return VectorSum<V1,V2>(v1,v2); }


template<class V1, class V2> inline
typename EnableIf< IsDefined<typename Arithmetic<typename V1::ValueType,typename V2::ValueType>::ResultType>, VectorSum< V1, V2 > >::Type
operator+(const VectorExpression<V1>& v1, const VectorExpression<V2>& v2) { return VectorSum<V1,V2>(v1(),v2()); }

template<class V1, class V2> inline
typename EnableIf< IsDefined<typename Arithmetic<typename V1::ValueType,typename V2::ValueType>::ResultType>, VectorDifference< V1, V2 > >::Type
operator-(const VectorExpression<V1>& v1, const VectorExpression<V2>& v2) { return VectorDifference<V1,V2>(v1(),v2()); }

template<class X1, class V2> inline
typename EnableIf< IsDefined<typename Arithmetic<X1,typename V2::ValueType>::ResultType>, VectorScalarProduct< V2, X1 > >::Type
operator*(const X1& x1, const VectorExpression<V2>& v2) { return VectorScalarProduct<V2,X1>(v2(),x1); }

template<class V1, class X2> inline
typename EnableIf< IsDefined<typename Arithmetic<typename V1::ValueType,X2>::ResultType>, VectorScalarProduct< V1, X2 > >::Type
operator*(const VectorExpression<V1>& v1, const X2& x2) { return VectorScalarProduct<V1,X2>(v1(),x2); }

template<class V1, class X2> inline
typename EnableIf< IsDefined<typename Arithmetic<typename V1::ValueType,X2>::ResultType>, VectorScalarQuotient< V1, X2 > >::Type
operator/(const VectorExpression<V1>& v1, const X2& x2) { return VectorScalarQuotient<V1,X2>(v1(),x2); }

template<class X, class V> inline Vector<X>& operator+=(Vector<X>& r, const VectorExpression<V>& ve) {
    const V& v=ve();
    ARIADNE_PRECONDITION(r.size()==v.size());
    for(size_t i=0; i!=r.size(); ++i) { r[i]+=v[i]; }
    return r;
}

template<class X, class V> inline Vector<X>& operator-=(Vector<X>& r, const VectorExpression<V>& ve) {
    const V& v=ve();
    ARIADNE_PRECONDITION(r.size()==v.size());
    for(size_t i=0; i!=r.size(); ++i) { r[i]-=v[i]; }
    return r;
}

template<class X1, class X2> inline typename EnableIfNumeric<X2,Vector<X1>&>::Type operator*=(Vector<X1>& v, const X2& x) {
    for(size_t i=0; i!=v.size(); ++i) { v[i]*=x; }
    return v;
}

template<class X1, class X2> inline typename EnableIfNumeric<X2,Vector<X1>&>::Type operator/=(Vector<X1>& v, const X2& x) {
    for(size_t i=0; i!=v.size(); ++i) { v[i]/=x; }
    return v;
}



template<class V>
typename V::ValueType sup_norm(const VectorExpression<V>& ve)
{
    const V& v=ve();
    typename V::ValueType r=0;
    for(size_t i=0; i!=v.size(); ++i) {
        typename V::ValueType absvi=abs(v[i]);
        // NOTE: The arguments must be this way round to propagate a nan row_sum
        r=max(absvi,r);
    }
    return r;
}

template<class V>
typename V::ValueType norm(const VectorExpression<V>& ve)
{
    const V& v=ve();
    return Ariadne::sup_norm(v);
}

template<class V1, class V2>
typename Arithmetic<typename V1::ValueType,typename V2::ValueType>::ResultType
dot(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2)
{
    const V1& v1=ve1(); const V2& v2=ve2();
    ARIADNE_ASSERT(v1.size()==v2.size());
    typename V1::ValueType r=0;
    for(size_t i=0; i!=v1.size(); ++i) {
        r+=v1[i]*v2[i];
    }
    return r;
}



template<class V1, class V2>
Vector<typename V1::ValueType>
join(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2,
     typename EnableIf< IsSame<typename V1::ValueType,typename V2::ValueType> >::Type* = 0)
{
    const V1& v1=ve1(); const V2& v2=ve2();
    if(v1.size()==0) { return v2; }
    if(v2.size()==0) { return v1; }
    size_t n1=v1.size();
    size_t n2=v2.size();
    Vector<typename V1::ValueType> r(n1+n2,v1[0]);
    project(r,range(0,n1))=v1;
    project(r,range(n1,n1+n2))=v2;
    return r;
}

template<class V1>
Vector<typename V1::ValueType> join(const VectorExpression<V1>& ve1, const typename V1::ValueType& s2)
{
    const V1& v1=ve1();
    size_t n1=v1.size();
    Vector<typename V1::ValueType> r(n1+1,s2);
    project(r,range(0,n1))=v1;
    return r;
}

template<class V2>
Vector<typename V2::ValueType> join(const typename V2::ValueType& s1, const VectorExpression<V2>& ve2)
{
    const V2& v2=ve2();
    size_t n2=v2.size();
    Vector<typename V2::ValueType> r(1+n2,s1);
    project(r,range(1,n2+1))=v2;
    return r;
}

template<class X>
Vector<X> join(const X& s1, const X& s2, typename EnableIf<IsScalar<X>,Void>::Type* = 0)
{
    Vector<X> r(2,s1);
    r[1]=s2;
    return r;
}

template<class V1, class V2, class V3>
Vector<typename V1::ValueType>
join(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2, const VectorExpression<V3>& ve3,
     typename EnableIf< And< IsSame<typename V1::ValueType,typename V2::ValueType>,IsSame<typename V1::ValueType,typename V3::ValueType> > >::Type* = 0)
{
    const V1& v1=ve1(); const V2& v2=ve2(); const V3& v3=ve3();
    typedef typename V1::ValueType X;
    if(v1.size()==0) { return join(v2,v3); }
    size_t n1=v1.size();
    size_t n2=v2.size();
    size_t n3=v3.size();
    Vector<X> r(n1+n2+n3,v1[0]);
    project(r,range(0,n1))=v1;
    project(r,range(n1,n1+n2))=v2;
    project(r,range(n1+n2,n1+n2+n3))=v3;
    return r;
}

template<class V1, class V2> Vector<typename V1::ValueType>
join(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2, const typename V1::ValueType& s3,
     typename EnableIf< IsSame<typename V1::ValueType,typename V2::ValueType> >::Type* = 0)
{
    const V1& v1=ve1(); const V2& v2=ve2();
    typedef typename V1::ValueType X;
    Vector<X> r(v1.size()+v2.size()+1u,s3);
    for(uint i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(uint i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    return r;
}

template<class V1, class V2, class V3> Vector<typename V1::ValueType>
join(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2, const VectorExpression<V3>& ve3, const typename V1::ValueType& s4,
     typename EnableIf< And< IsSame<typename V1::ValueType,typename V2::ValueType>,IsSame<typename V1::ValueType,typename V3::ValueType> > >::Type* = 0)
{
    const V1& v1=ve1(); const V2& v2=ve2(); const V3& v3=ve3();
    typedef typename V1::ValueType X;
    Vector<X> r(v1.size()+v2.size()+v3.size()+1u);
    for(uint i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(uint i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    for(uint i=0; i!=v3.size(); ++i) { r[v1.size()+v2.size()+i]=v3[i]; }
    r[r.size()-1]=s4;
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


template<class V1, class V2>
bool operator==(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2)
{
    return Vector<typename V1::ValueType>(ve1) == Vector<typename V2::ValueType>(ve2);
}


template<class X1, class X2> inline
bool operator!=(const Vector<X1>& v1, const Vector<X2>& v2)
{
    return !(v1==v2);
}


template<class X1, class X2>
bool operator<(const Vector<X1>& v1, const Vector<X2>& v2)
{
    if(v1.size()!=v2.size()) { return v1.size()<v2.size(); }
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


template<class V> std::ostream& operator<<(std::ostream& os, const VectorExpression<V>& ve) {
    return os << Vector<typename V::ValueType>(ve);
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



template<class X> bool refines(const Vector<X>& v1, const Vector<X>& v2) {
    assert(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!refines(v1[i],v2[i])) { return false; }
    }
    return true;
}

Vector<ExactFloat>const& make_exact(const Vector<ApproximateFloat>& av);
Vector<ValidatedFloat> make_bounds(const Vector<ErrorFloat>& ev);
Vector<ExactFloat> midpoint(const Vector<ValidatedFloat>& x);
PositiveUpperFloat sup_error(const Vector<ValidatedFloat>& x);

bool models(const Vector<ValidatedFloat>& x1, const Vector<ExactFloat>& x2);
bool consistent(const Vector<ValidatedFloat>& x1, const Vector<ValidatedFloat>& x2);
bool inconsistent(const Vector<ValidatedFloat>& x1, const Vector<ValidatedFloat>& x2);
bool refines(const Vector<ValidatedFloat>& x1, const Vector<ValidatedFloat>& x2);
Vector<ValidatedFloat> refinement(const Vector<ValidatedFloat>& x1, const Vector<ValidatedFloat>& x2);


} // namespace Ariadne

#endif // ARIADNE_VECTOR_H

