/***************************************************************************
 *            algebra/vector.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file algebra/vector.h
 *  \brief Vectors over a scalar (number or algebra).
 */

#ifndef ARIADNE_VECTOR_H
#define ARIADNE_VECTOR_H

#define SIMPLE_VECTOR_OPERATORS

#include "utility/macros.h"
#include "utility/metaprogramming.h"
#include "utility/container.h"
#include "utility/declarations.h"

namespace Ariadne {

//! \defgroup LinearAlgebraSubModule Linear Algebra Sub-Module
//! \ingroup AlgebraModule
//! \brief %Vector and matrix classes for linear algebra.

/************ Vector *********************************************************/

template<class V> struct VectorExpression { const V& operator()() const { return static_cast<const V&>(*this); } };
template<class V> struct VectorContainer : public VectorExpression<V> { };

template<class X> class Vector;
template<class X> class Covector;
template<class X> class Matrix;

template<class V> struct IsVector : False { };
template<class V> struct IsVectorExpression : IsVector<V> { };

template<class V> struct IsCovector : False { };
template<class V> struct IsCovectorExpression : IsVector<V> { };

template<class M> struct IsMatrix : False { };
template<class M> struct IsMatrixExpression : IsMatrix<M> { };

template<class X> struct IsNumber;
template<class A> struct IsAlgebra;

template<class X> struct IsScalar { static const Bool value = !IsVector<X>::value && !IsCovector<X>::value && !IsMatrix<X>::value; };

template<class V> using ScalarType=typename V::ScalarType;


template<class X, class = Fallback> struct HasCreateZero : False { };
template<class X> struct HasCreateZero<X, EnableIf<IsSame<decltype(declval<X>().create_zero()),X>,Fallback>> : True { };

template<class X, EnableIf<HasCreateZero<X>> = dummy> X create_zero(const X& x) { return x.create_zero(); }
// FIXME: Below should use a non-integral numeric type to prevent constructor of zero-sized object.
template<class X, DisableIf<HasCreateZero<X>> = dummy> X create_zero(const X& x) { return static_cast<X>(0u); }

template<class T> inline T zero_element(Matrix<T> const& m) { return m.zero_element(); }
template<class T> inline T zero_element(Covector<T> const& u) { return u.zero_element(); }
template<class T> inline T zero_element(Vector<T> const& v) { return v.zero_element(); }
template<class T> inline T zero_element(Scalar<T> const& s) { return create_zero(s); }


//! \ingroup LinearAlgebraSubModule
//! \brief Vectors over some type \a X.
//! Corresponds to elements of a \em module over a mathematical \em ring, or a <em>vector space</em> over a field.
//! May also be used if \a X is an \em algebra \a A over another field.
//! It must be possible to add and multiply any two elements of the vector.
template<class X>
class Vector
    : public VectorContainer<Vector<X>>
{
    Array<X> _ary;
  public:
    //@{
    //! \name Type definitions

    //! \brief The type used to index the elements.
    typedef SizeType IndexType;

    //! \brief The type of the scalar element.
    typedef X ScalarType;
    typedef X ValueType;

    //@}

    //@{
    //! \name Constructors

    //! \brief Default constructor constructs a vector with no elements.
    Vector() : _ary() { }
    //! \brief Construct a vector of size \a n, with elements initialised to the default value.
    explicit Vector(SizeType n) : _ary(n,X()) { static_assert(IsDefaultConstructible<X>::value,""); }
    //! \brief Construct a vector of size \a n, with elements initialised to \a t.
    explicit Vector(SizeType n, const X& t) : _ary(n,t) {  }
    //! \brief Construct from an array of the same type.
    explicit Vector(const Array<X>& ary) : _ary(ary) { }
    explicit Vector(Array<X>&& ary) : _ary(ary) { }
    //! \brief Construct from a list of the same type.
    explicit Vector(const List<X>& lst) : _ary(lst.begin(),lst.end()) { }
    //! \brief Convert from an initializer list of the same type.
    Vector(InitializerList<X> lst) : _ary(lst.begin(),lst.end()) { }

    //! \brief Convert from an %VectorExpression of a different type.
    template<class VE, EnableIf<IsConvertible<typename VE::ScalarType,X>> =dummy>
    Vector(VectorExpression<VE> const& ve) : _ary(ve().size(),ve().zero_element()) {
            for(SizeType i=0; i!=this->size(); ++i) { this->_ary[i]=ve()[i]; } }

    //! \brief Construct from an %VectorExpression of a different type.
    template<class VE, EnableIf<IsConstructible<X,typename VE::ScalarType>> =dummy, DisableIf<IsConvertible<typename VE::ScalarType,X>> =dummy>
    explicit Vector(VectorExpression<VE> const& ve) : _ary(ve().size(),X(ve().zero_element())) {
            for(SizeType i=0; i!=this->size(); ++i) { this->_ary[i]=X(ve()[i]); } }


    //! \brief Copy constructor.
    Vector(const Vector<X>& v) = default;
    //! \brief Move constructor.
    Vector(Vector<X>&& v) = default;
    //! \brief Copy assignment.
    Vector<X>& operator=(const Vector<X>& v) = default;
    //@}

    //@{
    //! \name Static constructors

    //! \brief The zero vector of size \a n.
    static Vector<X> zero(SizeType n) { return Vector<X>(n,static_cast<X>(0)); }
    //! \brief The vector of size \a n with all entries equal to one.
    static Vector<X> one(SizeType n) { return Vector<X>(n,static_cast<X>(1)); }
    //! \brief The unit vector \f$e_i\f$ with value one in the \a i<sup>th</sup> entry, and zero otherwise.
    static Vector<X> unit(SizeType n,SizeType i) {
        ARIADNE_ASSERT(i<n); Vector<X> result(n,static_cast<X>(0)); result[i]=static_cast<X>(1); return result; }
    //! \brief The unit vector \f$e_i\f$ with value one in the \a i<sup>th</sup> entry, and zero otherwise.
    static Array< Vector<X> > basis(SizeType n) {
        Array< Vector<X> > result(n); for(Nat i=0; i!=n; ++i) { result[i]=unit(n,i); } return result; }
    //@}

    //@{
    //! \name Data access

    //! \brief Resize to hold \a n elements.
    //! The previous values need not be preserved, and the new values need not be initialised.
    Void resize(SizeType n) { this->_ary.resize(n); }
    //! \brief The number of elements of the vector.
    SizeType size() const { return this->_ary.size(); }
    //! \brief A reference to the value stored in the \a i<sup>th</sup> element.
    X& at(SizeType i) { ARIADNE_PRECONDITION_MSG(i<this->size(),*this<<"["<<i<<"]"); return (*this)[i]; }
    //! \brief A constant reference to the value stored in the \a i<sup>th</sup> element.
    const X& at(SizeType i) const { ARIADNE_PRECONDITION_MSG(i<this->size(),*this<<"["<<i<<"]"); return (*this)[i]; }
    //! \brief Get the value stored in the \a i<sup>th</sup> element.
    const X& get(SizeType i) const { ARIADNE_PRECONDITION_MSG(i<this->size(),*this<<"["<<i<<"]"); return (*this)[i]; }
    //! \brief Set the value stored in the \a i<sup>th</sup> element to \a x.
    Void set(SizeType i, const X& x) { ARIADNE_PRECONDITION_MSG(i<this->size(),*this<<"["<<i<<"]"); (*this)[i] = x; }
    //! \brief C-style subscripting operator.
    X& operator[](SizeType i) { ARIADNE_PRECONDITION_MSG(i<this->size(),*this<<"["<<i<<"]"); return this->_ary[i]; }
    //! \brief C-style constant subscripting operator.
    const X& operator[](SizeType i) const { ARIADNE_PRECONDITION_MSG(i<this->size(),*this<<"["<<i<<"]"); return this->_ary[i]; }
    //! \brief The zero of the ring containing the Vector's elements. This may be dependent on class parameters.
    const X zero_element() const { if(this->size()!=0) { return create_zero((*this)[0]); } else { return X(); } }
    //! \brief The raw data array.
    Array<X> const& array() const { return _ary; }
    //@}

#ifdef DOXYGEN
    //! \brief Equality operator.
    friend template<class X1, class X2> decltype(declval<X1>()==declval<X2>()) operator==(const Vector<X1>& v1, const Vector<X2>& v2);
    //! \brief Inequality operator.
    friend template<class X1, class X2> decltype(declval<X1>()!=declval<X2>()) operator!=(const Vector<X1>& v1, const Vector<X2>& v2);

     //! \brief %Vector unary plus.
    friend template<class X> Vector<decltype(+declval<X>())> operator+(const Vector<X>& v);
     //! \brief %Vector negation.
    friend template<class X> Vector<decltype(-declval<X>())> operator-(const Vector<X>& v);
    //! \brief %Vector addition.
    friend template<class X1, class X2> Vector<decltype(declval<X1>()+declval<X2>())> operator+(const Vector<X1>& v1, const Vector<X2>& v2);
    //! \brief %Vector subtraction.
    friend template<class X1, class X2> Vector<decltype(declval<X1>()+declval<X2>())>(const Vector<X1>& v1, const Vector<X>& v2);
    //! \brief %Scalar multiplication.
    friend template<class X1, class X2> Vector<decltype(declval<X1>()+declval<X2>())>(const X1& s, const Vector<X2>& v);
    //! \brief %Scalar multiplication.
    friend template<class X1, class X2> Vector<decltype(declval<X1>()+declval<X2>())>(const Vector<X1>& v, const X2& s);
    //! \brief %Scalar division.
    friend template<class X1, class X2> Vector<decltype(declval<X1>()+declval<X2>())>(const Vector<X1>& v, const X2& s);

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
    //! \brief Join two scalars. // FIXME: Removed due to poor detection of scalar types
    // friend template<class X, EnableIf<IsScalar<X>> =dummy> Vector<X> join(const X& s1, const X& s2);

    //! \brief Write to an output stream.
    friend template<class X> OutputStream& operator<<(OutputStream& os, const Vector<X>& v);
    //! \brief Read from an output stream.
    friend template<class X> InputStream& operator>>(InputStream& is, Vector<X>& v);
#endif // DOXYGEN
};

template<class X> struct IsVector<Vector<X>> : True { };

class Range {
    SizeType _start; SizeType _stop;
  public:
    Range(SizeType start, SizeType stop) : _start(start), _stop(stop) { }
    SizeType size() const { return this->_stop-this->_start; }
    SizeType start() const { return this->_start; }
    SizeType stride() const { return 1u; }
    SizeType stop() const { return this->_stop; }
};

class Slice {
    SizeType _size; SizeType _start; SizeType _stride;
  public:
    Slice(SizeType size, SizeType start, SizeType stride) : _size(size), _start(start), _stride(stride) { }
    SizeType size() const { return this->_size; }
    SizeType start() const { return this->_start; }
    SizeType stride() const { return this->_stride; }
    SizeType stop() const { return this->_start+this->_size*this->_stride; }
};

inline Range range(SizeType start, SizeType stop) { return Range(start,stop); }
inline Slice slice(SizeType size, SizeType start, SizeType stride) { return Slice(size,start,stride); }


template<class V> struct VectorRange
    : public VectorContainer< VectorRange<V> >
{
    typedef typename V::ScalarType ScalarType;
    typedef typename V::ValueType ValueType;
    const V& _v; Range _rng;
    VectorRange(const V& v, Range rng) : _v(v), _rng(rng) { }
    SizeType size() const { return _rng.size(); }
    ScalarType zero_element() const { return _v.zero_element(); }
    ScalarType operator[](SizeType i) const { return _v[i+_rng.start()]; }
};

template<class V> struct VectorContainerRange
    : public VectorExpression< VectorContainerRange<V> >
{
    typedef typename V::ScalarType ScalarType;
    typedef typename V::ValueType ValueType;
    V& _v; Range _rng;
    VectorContainerRange(V& v, Range rng) : _v(v), _rng(rng) { }
    SizeType size() const { return _rng.size(); }
    ValueType operator[](SizeType i) const { return _v[i+_rng.start()]; }
    ValueType& operator[](SizeType i) { return _v[i+_rng.start()]; }
    const ScalarType zero_element() const { return _v.zero_element(); }
    Void set(SizeType i, const ValueType& x) { _v[i+_rng.start()]=x; }
    template<class VE> VectorContainerRange<V>& operator=(const VectorExpression<VE>& ve) {
        ARIADNE_PRECONDITION(this->size()==ve().size());
        for(SizeType i=0; i!=this->size(); ++i) { (*this)[i]=ve()[i]; } return *this; }
};

template<class V> inline VectorRange<V> project(const VectorExpression<V>& v, Range rng) {
    return VectorRange<V>(v(),rng); }
template<class X> inline VectorContainerRange< Vector<X> > project(Vector<X>& v, Range rng) {
    return VectorContainerRange< Vector<X> >(v,rng); }




template<class X> OutputStream& operator<<(OutputStream& os, Vector<X> const& v) {
    if(v.size()==0) { os << "{"; }
    for(SizeType i=0; i!=v.size(); ++i) { os << (i==0u?"{":",") << v[i]; }
    return os << "}";
}

template<class V, EnableIf<IsVectorExpression<V>> =dummy> OutputStream& operator<<(OutputStream& os, const V& v) {
    typedef decltype(v[0]) X;
    return os << Vector<X>(v);
}


template<class X1, class X2>
auto operator==(const Vector<X1>& v1, const Vector<X2>& v2) -> decltype(v1[0]==v2[0]) {
    decltype(v1[0]==v2[0]) r=true;
    for(SizeType i=0; i!=v1.size(); ++i) {
        r = r && (v1[i]==v2[i]);
    }
    return r;
}

template<class X1, class X2>
auto operator!=(const Vector<X1>& v1, const Vector<X2>& v2) -> decltype(v1[0]!=v2[0]) {
    decltype(v1[0]!=v2[0]) r=false;
    for(SizeType i=0; i!=v1.size(); ++i) {
        r = r || (v1[i]==v2[i]);
    }
    return r;
}


template<class X, class XX, EnableIf<IsConvertible<decltype(declval<X>()+declval<XX>()),X>> =dummy> inline
Vector<X>& operator+=(Vector<X>& v1, const Vector<XX>& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    for(SizeType i=0; i!=v1.size(); ++i) { v1[i]+=v2[i]; } return v1;
}

template<class X, class XX, EnableIf<IsConvertible<decltype(declval<X>()-declval<XX>()),X>> =dummy> inline
Vector<X>& operator-=(Vector<X>& v1, const Vector<XX>& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    for(SizeType i=0; i!=v1.size(); ++i) { v1[i]-=v2[i]; } return v1;
}

template<class X, class XX, EnableIf<IsConvertible<decltype(declval<X>()*declval<XX>()),X>> =dummy> inline
Vector<X>& operator*=(Vector<X>& v, const XX& s) {
    for(SizeType i=0; i!=v.size(); ++i) { v[i]*=s; } return v;
}

template<class X, class XX, EnableIf<IsConvertible<decltype(declval<X>()/declval<XX>()),X>> =dummy> inline
Vector<X>& operator/=(Vector<X>& v, const XX& s) {
    for(SizeType i=0; i!=v.size(); ++i) { v[i]/=s; } return v;
}



#ifdef SIMPLE_VECTOR_OPERATORS

template<class X> Vector<decltype(+declval<X>())> operator+(Vector<X> const& v) {
    Vector<decltype(+declval<X>())> r(v.size(),+v.zero_element());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=+v[i]; }
    return std::move(r);
}

template<class X> Vector<X> operator-(Vector<X> const& v) {
    Vector<decltype(-declval<X>())> r(v.size(),-v.zero_element());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=-v[i]; }
    return std::move(r);
}

template<class X1, class X2> Vector<decltype(declval<X1>()+declval<X2>())> operator+(Vector<X1> const& v1, Vector<X2> const& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    Vector<decltype(declval<X1>()+declval<X2>())> r(v1.size(),v1.zero_element()+v2.zero_element());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=v1[i]+v2[i]; }
    return std::move(r);
}

template<class X1, class X2> Vector<decltype(declval<X1>()-declval<X2>())> operator-(Vector<X1> const& v1, Vector<X2> const& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    Vector<decltype(declval<X1>()-declval<X2>())> r(v1.size(),v1.zero_element()-v2.zero_element());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=v1[i]-v2[i]; }
    return std::move(r);
}

template<class X1, class X2, EnableIf<IsScalar<X1>> = dummy> Vector<decltype(declval<X1>()*declval<X2>())> operator*(X1 const& x1, Vector<X2> const& v2) {
    Vector<decltype(declval<X1>()*declval<X2>())> r(v2.size(),x1*v2.zero_element());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=x1*v2[i]; }
    return std::move(r);
}

template<class X1, class X2, EnableIf<IsScalar<X2>> = dummy> Vector<decltype(declval<X1>()*declval<X2>())> operator*(Vector<X1> const& v1, X2 const& x2) {
    Vector<decltype(declval<X1>()*declval<X2>())> r(v1.size(),v1.zero_element()*x2);
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=v1[i]*x2; }
    return std::move(r);
}

template<class X1, class X2, EnableIf<IsScalar<X2>> = dummy> Vector<decltype(declval<X1>()/declval<X2>())> operator/(Vector<X1> const& v1, X2 const& x2) {
    Vector<decltype(declval<X1>()/declval<X2>())> r(v1.size(),v1.zero_element()/x2);
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=v1[i]/x2; }
    return std::move(r);
}

#else

template<class V, EnableIf<IsVectorExpression<V>> =dummy> inline
const V& operator+(const V& v) {
    return v; }


template<class V> struct VectorNegation {
    typedef NegationType<typename V::ScalarType> ScalarType;
    const V& _v;
    VectorNegation(const V& v) : _v(v) { }
    SizeType size() const { return _v.size(); }
    ScalarType zero_element() const { return -_v.zero_element(); }
    ScalarType operator[](SizeType i) const { return -_v[i]; }
};
template<class V> struct IsVectorExpression<VectorNegation<V>> : True { };

template<class V, EnableIf<IsVectorExpression<V>> =dummy> inline
VectorNegation<V> operator-(const V& v) {
    return VectorNegation<V>(v); }


template<class V1, class V2> struct VectorSum  {
    typedef SumType<typename V1::ScalarType,typename V2::ScalarType> ScalarType;
    const V1& _v1; const V2& _v2;
    VectorSum(const V1& v1, const V2& v2) : _v1(v1), _v2(v2) { }
    SizeType size() const { return _v1.size(); }
    ScalarType zero_element() const { return _v1.zero_element()+_v2.zero_element(); }
    ScalarType operator[](SizeType i) const { return _v1[i]+_v2[i]; }
};
template<class V1, class V2> struct IsVectorExpression<VectorSum<V1,V2>> : True { };

template<class V1, class V2, EnableIf<And<IsVectorExpression<V1>,IsVectorExpression<V2>>> =dummy> inline
VectorSum<V1,V2> operator+(const V1& v1, const V2& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    return VectorSum<V1,V2>(v1,v2); }


template<class V1, class V2> struct VectorDifference {
    typedef DifferenceType<typename V1::ScalarType,typename V2::ScalarType> ScalarType;
    const V1& _v1; const V2& _v2;
    VectorDifference(const V1& v1, const V2& v2) : _v1(v1), _v2(v2) { }
    SizeType size() const { return _v1.size(); }
    ScalarType zero_element() const { return _v1.zero_element()-_v2.zero_element(); }
    ScalarType operator[](SizeType i) const { return _v1[i]-_v2[i]; }
};
template<class V1, class V2> struct IsVectorExpression<VectorDifference<V1,V2>> : True { };

template<class V1, class V2, EnableIf<And<IsVectorExpression<V1>,IsVectorExpression<V2>>> =dummy> inline
VectorDifference<V1,V2> operator-(const V1& v1, const V2& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    return VectorDifference<V1,V2>(v1,v2); }


template<class V1, class X2> struct VectorScalarProduct {
    typedef ProductType<typename V1::ScalarType,X2> ScalarType;
    const V1& _v1; const X2& _x2;
    VectorScalarProduct(const V1& v1, const X2& x2) : _v1(v1), _x2(x2) { }
    SizeType size() const { return _v1.size(); }
    ScalarType zero_element() const { return _v1.zero_element()*_x2; }
    ScalarType operator[](SizeType i) const { return _v1[i]*_x2; }
};
template<class V1, class X2> struct IsVectorExpression<VectorScalarProduct<V1,X2>> : True { };

template<class X1, class V2, EnableIf<And<IsScalar<X1>,IsVectorExpression<V2>>> =dummy> inline
VectorScalarProduct<V2,X1> operator*(const X1& x1, const V2& v2) {
    return VectorScalarProduct<V2,X1>(v2,x1); }

template<class V1, class X2, EnableIf<And<IsVectorExpression<V1>,IsScalar<X2>>> =dummy> inline
VectorScalarProduct<V1,X2> operator*(const V1& v1, const X2& x2) {
    return VectorScalarProduct<V1,X2>(v1,x2); }

template<class V1, class X2> struct VectorScalarQuotient {
    typedef QuotientType<typename V1::ScalarType,X2> ScalarType;
    const V1& _v1; const X2& _x2;
    VectorScalarQuotient(const V1& v1, const X2& x2) : _v1(v1), _x2(x2) { }
    SizeType size() const { return _v1.size(); }
    ScalarType zero_element() const { return _v1.zero_element()/_x2; }
    ScalarType operator[](SizeType i) const { return _v1[i]/_x2; }
};
template<class V1, class X2> struct IsVectorExpression<VectorScalarQuotient<V1,X2>> : True { };

template<class V1, class X2, EnableIf<And<IsVectorExpression<V1>,IsScalar<X2>>> =dummy> inline
VectorScalarQuotient<V1,X2> operator/(const V1& v1, const X2& x2) {
    return VectorScalarQuotient<V1,X2>(v1,x2); }

template<class V> using ScalarType = typename V::ScalarType;

#endif

template<class X> inline auto norm(const Vector<X>& v) -> decltype(abs(declval<X>())) {
    decltype(abs(declval<X>())) r=0u;
    for(SizeType i=0; i!=v.size(); ++i) {
        r=max(r,abs(v[i]));
    }
    return r;
}

template<class X> inline auto sup_norm(const Vector<X>& v) -> decltype(mag(declval<X>())) {
    decltype(mag(declval<X>())) r=0u;
    for(SizeType i=0; i!=v.size(); ++i) {
        r=max(r,mag(v[i]));
    }
    return r;
}

template<class X1, class X2> inline auto dot(const Vector<X1>& v1, const Vector<X2>& v2) -> decltype(v1[0]*v2[0]+v1[0]*v2[0]) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    decltype(declval<X1>()*declval<X2>()+declval<X1>()*declval<X2>()) r(0u);
    for(SizeType i=0; i!=v1.size(); ++i) {
        r+=v1[i]*v2[i];
    }
    return r;
}



template<class X>
Vector<X> join(const Vector<X>& v1, const Vector<X>& v2)
{
    if(v1.size()==0) { return v2; }
    if(v2.size()==0) { return v1; }
    SizeType n1=v1.size();
    SizeType n2=v2.size();
    Vector<X> r(n1+n2,v1[0]);
    for(SizeType i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(SizeType i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    return r;
}


template<class X>
Vector<X> join(const Vector<X>& v1, const Vector<X>& v2, const Vector<X>& v3)
{
    Vector<X> r(v1.size()+v2.size()+v3.size());
    for(SizeType i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(SizeType i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    for(SizeType i=0; i!=v3.size(); ++i) { r[v1.size()+v2.size()+i]=v3[i]; }
    return std::move(r);
}

template<class X>
Vector<X> join(const typename Vector<X>::ScalarType& x1, const Vector<X>& v2)
{
    Vector<X> r(1u+v2.size());
    r[0u]=x1;
    for(SizeType i=0; i!=v2.size(); ++i) { r[1u+i]=v2[i]; }
    return std::move(r);
}

template<class X>
Vector<X> join(const Vector<X>& v1, const typename Vector<X>::ScalarType& x2)
{
    Vector<X> r(v1.size()+1u);
    for(SizeType i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    r[v1.size()]=x2;
    return std::move(r);
}

//template<class X, EnableIf<IsScalar<X>> =dummy>
//Vector<X> join(const X& x1, const X& x2)
//{
//    Vector<X> r(2u);
//    r[0u]=x1;
//    r[1u]=x2;
//    return std::move(r);
//}

template<class X>
Vector<X> join(const Vector<X>& v1, const Vector<X>& v2, const typename Vector<X>::ScalarType& x3)
{
    Vector<X> r(v1.size()+v2.size()+1u);
    for(SizeType i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(SizeType i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    r[v1.size()+v2.size()]=x3;
    return std::move(r);
}


template<class X> inline decltype(error(declval<X>())) error(const Vector<X>& v) {
    decltype(error(declval<X>())) r=0u;
    for(SizeType i=0; i!=v.size(); ++i) {
        r=max(r,error(v[i]));
    }
    return r;
}

template<class X> inline decltype(error(declval<X>())) sup_error(const Vector<X>& v) {
    decltype(error(declval<X>())) r=0u;
    for(SizeType i=0; i!=v.size(); ++i) {
        r=max(r,error(v[i]));
    }
    return r;
}




template<class X> inline Vector<decltype(refinement(declval<X>(),declval<X>()))> refinement(const Vector<X>& v1, const Vector<X>& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    Vector<X> r(v1.size(),v1.zero_element());
    for(SizeType i=0; i!=v1.size(); ++i) {
        r[i]=refinement(v1[i],v2[i]);
    }
    return std::move(r);
}

template<class X> inline decltype(refines(declval<X>(),declval<X>())) refines(const Vector<X>& v1, const Vector<X>& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    for(SizeType i=0; i!=v1.size(); ++i) {
        if(!refines(v1[i],v2[i])) { return false; }
    }
    return true;
}

template<class VX, class EX> inline decltype(models(declval<VX>(),declval<EX>())) models(const Vector<VX>& v1, const Vector<EX>& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    for(SizeType i=0; i!=v1.size(); ++i) {
        if(!models(v1[i],v2[i])) { return false; }
    }
    return true;
}

template<class VX, class EX> inline decltype(represents(declval<VX>(),declval<EX>())) represents(const Vector<VX>& v1, const Vector<EX>& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    for(SizeType i=0; i!=v1.size(); ++i) {
        if(!represents(v1[i],v2[i])) { return false; }
    }
    return true;
}

template<class X1, class X2> inline decltype(inconsistent(declval<X1>(),declval<X2>())) inconsistent(const Vector<X1>& v1, const Vector<X2>& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    for(SizeType i=0; i!=v1.size(); ++i) {
        if(inconsistent(v1[i],v2[i])) { return true; }
    }
    return false;
}

template<class X1, class X2> inline decltype(consistent(declval<X1>(),declval<X2>())) consistent(const Vector<X1>& v1, const Vector<X2>& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    for(SizeType i=0; i!=v1.size(); ++i) {
        if(!consistent(v1[i],v2[i])) { return false; }
    }
    return true;
}


template<class X> using MidpointType = RemoveConst<decltype(midpoint(declval<X>()))>;

template<class X> inline Vector<MidpointType<X>> midpoint(const Vector<X>& v) {
    Vector<MidpointType<X>> r(v.size(),midpoint(v.zero_element()));
    for(SizeType i=0; i!=v.size(); ++i) {
        r[i]=midpoint(v[i]);
    }
    return r;
}

template<class X> using SingletonType = decltype(cast_singleton(declval<X>()));

template<class X> inline Vector<SingletonType<X>> cast_singleton(const Vector<X>& v) {
    Vector<SingletonType<X>> r(v.size(),cast_singleton(v.zero_element()));
    for(SizeType i=0; i!=v.size(); ++i) {
        r[i]=cast_singleton(v[i]);
    }
    return r;
}

template<class X> using BoundsType = decltype(make_bounds(declval<X>()));

template<class X> inline Vector<BoundsType<X>> make_bounds(const Vector<X>& v) {
    Vector<BoundsType<X>> r(v.size(),make_bounds(v.zero_element()));
    for(SizeType i=0; i!=v.size(); ++i) {
        r[i]=make_bounds(v[i]);
    }
    return std::move(r);
}


template<class X> using ExactType = RemoveConst<RemoveReference<decltype(cast_exact(declval<X>()))>>;

template<class X> inline Vector<ExactType<X>> cast_exact(const Vector<X>& v) {
    Vector<ExactType<X>> r(v.size(),cast_exact(v.zero_element()));
    for(SizeType i=0; i!=v.size(); ++i) {
        r[i]=cast_exact(v[i]);
    }
    return std::move(r);
}

template<class P, class PR> class Float;
class Precision64;
class Exact; class Approximate;
typedef Float<Exact,Precision64> ExactFloat64;
typedef Float<Approximate,Precision64> ApproximateFloat64;
inline Vector<ExactFloat64>const& cast_exact(Vector<ApproximateFloat64>const& v) {
    return reinterpret_cast<Vector<ExactFloat64>const&>(v);
}


} // namespace Ariadne

#endif
