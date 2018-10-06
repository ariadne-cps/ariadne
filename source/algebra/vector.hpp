/***************************************************************************
 *            algebra/vector.hpp
 *
 *  Copyright 2013-17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file algebra/vector.hpp
 *  \brief Vectors over a scalar (number or algebra).
 */

#ifndef ARIADNE_VECTOR_HPP
#define ARIADNE_VECTOR_HPP

#define SIMPLE_VECTOR_OPERATORS

#include "../utility/macros.hpp"
#include "../utility/metaprogramming.hpp"
#include "../utility/container.hpp"
#include "../utility/declarations.hpp"
#include "../numeric/builtin.hpp"

namespace Ariadne {

/************ Vector *********************************************************/

template<class X> class Vector;
template<class X> class Covector;
template<class X> class Matrix;

template<class V> struct IsVector : False { };
template<class V> struct IsVectorExpression : IsVector<V> { };

template<class V> struct IsCovector : False { };
template<class V> struct IsCovectorExpression : IsVector<V> { };

template<class M> struct IsMatrix : False { };
template<class M> struct IsMatrixExpression : IsMatrix<M> { };

template<class X> struct IsNumericType;
template<class A> struct IsAlgebra;

template<class X> struct IsScalar { static const Bool value = !IsVector<X>::value && !IsCovector<X>::value && !IsMatrix<X>::value; };

template<class V> using ScalarType=typename V::ScalarType;

template<class X> struct HasCreateZero {
    template<class XX, class=decltype(std::declval<XX>().create_zero())> static std::true_type test(int);
    template<class XX> static std::false_type test(...);
    static const bool value = decltype(test<X>(1))::value;
};
template<class X> struct HasNul {
    template<class XX, class=decltype(nul(std::declval<XX>()))> static std::true_type test(int);
    template<class XX> static std::false_type test(...);
    static const bool value = decltype(test<X>(1))::value;
};

template<class X, EnableIf<HasCreateZero<X>> = dummy> X create_zero(const X& x) { return x.create_zero(); }
template<class X, DisableIf<HasCreateZero<X>> = dummy, EnableIf<HasNul<X>> = dummy> X create_zero(const X& x) { return nul(x); }
// FIXME: Below should use a non-integral numeric type to prevent constructor of zero-sized object.
template<class X, DisableIf<HasCreateZero<X>> = dummy, DisableIf<HasNul<X>> = dummy> X create_zero(const X& x) { return static_cast<X>(0u); }

template<class T> inline T zero_element(Matrix<T> const& m) { return m.zero_element(); }
template<class T> inline T zero_element(Covector<T> const& u) { return u.zero_element(); }
template<class T> inline T zero_element(Vector<T> const& v) { return v.zero_element(); }
template<class T> inline T zero_element(Scalar<T> const& s) { return create_zero(s); }

class Range;
template<class X> class VectorRange;

#ifdef SIMPLE_VECTOR_OPERATORS

struct DeclareVectorOperations {
    template<class X> friend Vector<X> operator+(Vector<X> const& v);
    template<class X> friend Vector<NegationType<X>> operator-(Vector<X> const& v);
    template<class X1, class X2> friend Vector<SumType<X1,X2>> operator+(Vector<X1> const& v1, Vector<X2> const& v2);
    template<class X1, class X2> friend Vector<DifferenceType<X1,X2>> operator-(Vector<X1> const& v1, Vector<X2> const& v2);
    template<class X1, class X2> friend Vector<ProductType<Scalar<X1>,X2>> operator*(X1 const& x1, Vector<X2> const& v2);
    template<class X1, class X2> friend Vector<ProductType<X1,Scalar<X2>>> operator*(Vector<X1> const& v1, X2 const& x2);
    template<class X1, class X2> friend Vector<QuotientType<X1,Scalar<X2>>> operator/(Vector<X1> const& v1, X2 const& x2);
    template<class X1, class X2> friend Vector<InplaceSumType<X1,X2>>& operator+=(Vector<X1>& v1, const Vector<X2>& v2);
    template<class X1, class X2> friend Vector<InplaceDifferenceType<X1,X2>>& operator-=(Vector<X1>& v1, const Vector<X2>& v2);
    template<class X1, class X2> friend Vector<InplaceProductType<X1,X2>>& operator*=(Vector<X1>& v1, X2 const& x2);
    template<class X1, class X2> friend Vector<InplaceQuotientType<X1,X2>>& operator/=(Vector<X1>& v1, X2 const& x2);
    template<class X> friend decltype(abs(declval<X>())) norm(Vector<X> const& v);
    template<class X> friend decltype(mag(declval<X>())) sup_norm(Vector<X> const& v);
    template<class X> friend decltype(sqrt(sqr(declval<X>()))) two_norm(Vector<X> const& v);
    template<class X1, class X2> friend ArithmeticType<X1,X2> dot(Vector<X1> const& v1, Vector<X2> const& v2);
    template<class X1, class X2> friend EqualsType<X1,X2> operator==(Vector<X1> const& v1, Vector<X2> const& v2);
    template<class X1, class X2> friend decltype(declval<X1>()!=declval<X2>()) operator!=(Vector<X1> const& v1, Vector<X2> const& v2);
    template<class X> friend OutputStream& operator<<(OutputStream& os, Vector<X> const& v);
};

#else // SIMPLE_VECTOR_OPERATORS

struct DeclareVectorOperations { };

#endif // SIMPLE_VECTOR_OPERATORS

template<class V> struct VectorExpression : public DeclareVectorOperations { const V& operator()() const { return static_cast<const V&>(*this); } };
template<class V> struct VectorContainer : public VectorExpression<V> { };

//! \ingroup LinearAlgebraModule
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
    //! Construct a vector from parameters of \a X.
    template<class... PRS, EnableIf<IsConstructible<X,PRS...>> =dummy> explicit Vector(SizeType n, PRS... prs) : Vector(n,X(prs...)) { }
    //! \brief Construct from an array of the same type.
    explicit Vector(const Array<X>& ary) : _ary(ary) { }
    explicit Vector(Array<X>&& ary) : _ary(ary) { }
    //! \brief Construct from a list of the same type.
    explicit Vector(const List<X>& lst) : _ary(lst.begin(),lst.end()) { }
    //! \brief Convert from an initializer list of the same type.
    Vector(InitializerList<X> lst) : _ary(lst.begin(),lst.end()) { }

    //! \brief Convert from an initializer list of generic type and a precision parameter.
    template<class PR, EnableIf<IsConstructible<X,ExactDouble,PR>> =dummy> Vector(InitializerList<double> const& lst, PR pr);
    template<class PR, EnableIf<IsConstructible<X,Real,PR>> =dummy> Vector(InitializerList<Real> const& lst, PR pr);
    template<class PR, EnableIf<IsConstructible<X,ExactDouble,ExactDouble,PR>> =dummy> Vector(InitializerList<Pair<double,double>> const& lst, PR pr);
    //! \brief Convert from an array of generic type and a precision parameter.
    template<class Y, class PR, EnableIf<IsConstructible<X,Y,PR>> =dummy> Vector(Array<Y> const& ary, PR pr) : _ary(ary,pr) { }
    //! \brief Convert from an vector of generic type and a precision parameter.
    template<class Y, class PR, EnableIf<IsConstructible<X,Y,PR>> =dummy> Vector(Vector<Y> const& v, PR pr) : _ary(v.array(),pr) { }
    //! \brief Convert from an %VectorExpression of a different type.
    template<class VE, EnableIf<IsConvertible<typename VE::ScalarType,X>> =dummy>
    Vector(VectorExpression<VE> const& ve) : _ary(ve().size(),ve().zero_element()) {
            for(SizeType i=0; i!=this->size(); ++i) { this->_ary[i]=ve()[i]; } }

    /*! \brief Generate from a function (object) \a g of type \a G mapping an index to a value. */
    template<class G, EnableIf<IsInvocableReturning<X,G,SizeType>> =dummy>
    Vector(SizeType n, G const& g) : _ary(n,g) { }

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
    const X& get(SizeType i) const { return this->_ary[i]; }
    //! \brief Set the value stored in the \a i<sup>th</sup> element to \a x.
    Void set(SizeType i, const X& x) { this->_ary[i] = x; }
    //! \brief Subscripting operator. Unchecked access.
    X& operator[](SizeType i) { return this->_ary[i]; }
    //! \brief Constant subscripting operator.
    const X& operator[](SizeType i) const { return this->_ary[i]; }
    //! \brief Range subscripting operator.
    VectorRange<Vector<X>> operator[](Range rng); // { ARIADNE_PRECONDITION_MSG(rng.stop()<this->size(),*this<<"["<<r<<"]"); return VectorRange<X>(*this,rng); }
    //! \brief Constant range subscripting operator.
    VectorRange<const Vector<X>> operator[](Range rng) const;
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
    friend template<class X> Vector<X> operator+(const Vector<X>& v);
     //! \brief %Vector negation.
    friend template<class X> Vector<NegationType<X>> operator-(const Vector<X>& v);
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
    SizeType operator[](SizeType i) const { return _start+i; }
    Range(SizeType start, SizeType stop) : _start(start), _stop(stop) { }
    SizeType size() const { return this->_stop-this->_start; }
    SizeType start() const { return this->_start; }
    SizeType stride() const { return 1u; }
    SizeType stop() const { return this->_stop; }
};
inline Range range(SizeType stop) { return Range(0u,stop); }
inline Range range(SizeType start, SizeType stop) { return Range(start,stop); }

struct RangeIterator {
    explicit inline RangeIterator(SizeType i) : _i(i) { }
    inline RangeIterator& operator++() { ++this->_i; return *this; }
    inline SizeType operator*() const { return this->_i; }
    friend inline bool operator!=(RangeIterator iter1, RangeIterator iter2) { return iter1._i != iter2._i; }
  private:
    SizeType _i;
};
inline RangeIterator begin(Range rng) { return RangeIterator(rng.start()); }
inline RangeIterator end(Range rng) { return RangeIterator(rng.stop()); }

class Slice {
    SizeType _size; SizeType _start; SizeType _stride;
  public:
    Slice(SizeType size, SizeType start, SizeType stride) : _size(size), _start(start), _stride(stride) { }
    SizeType operator[](SizeType i) { return _start+i*_stride; }
    SizeType size() const { return this->_size; }
    SizeType start() const { return this->_start; }
    SizeType stride() const { return this->_stride; }
    SizeType stop() const { return this->_start+this->_size*this->_stride; }
};
inline Slice slice(SizeType size_, SizeType start, SizeType stride) { return Slice(size_,start,stride); }


template<class V> class VectorRange
    : public VectorContainer< VectorRange<V> >
{
    const V& _v; Range _rng;
  public:
    typedef typename V::ScalarType ScalarType;
    VectorRange(const V& v, Range rng) : _v(v), _rng(rng) { }
    SizeType size() const { return _rng.size(); }
    ScalarType zero_element() const { return _v.zero_element(); }
    ScalarType operator[](SizeType i) const { return _v[i+_rng.start()]; }

    Void set(SizeType i, const ScalarType& x) { _v[i+_rng.start()]=x; }
    ScalarType& operator[](SizeType i) { return _v[i+_rng.start()]; }
    template<class VE> VectorRange<V>& operator=(const VectorExpression<VE>& ve) {
        ARIADNE_PRECONDITION(this->size()==ve().size());
        for(SizeType i=0; i!=this->size(); ++i) { (*this)[i]=ve()[i]; } return *this; }
};

template<class V> inline VectorRange<const V> project(const VectorExpression<V>& v, Range rng) {
    return VectorRange<const V>(v(),rng); }
template<class V> inline VectorRange<V> project(VectorContainer<V>& v, Range rng) {
    return VectorRange<V>(v(),rng); }
template<class X> inline VectorRange<Vector<X>> Vector<X>::operator[](Range rng) {
    return project(*this,rng); }
template<class X> inline VectorRange<const Vector<X>> Vector<X>::operator[](Range rng) const {
    return project(*this,rng); }

template<class V, EnableIf<IsVectorExpression<V>> =dummy> OutputStream& operator<<(OutputStream& os, const V& v) {
    typedef decltype(v[0]) X;
    return os << Vector<X>(v);
}



#ifdef SIMPLE_VECTOR_OPERATORS

struct ProvideVectorOperations {
    template<class X> friend Vector<X> operator+(Vector<X> const& v) {
        return Vector<X>( v.size(), [&v](SizeType i){return +v[i];} ); }

    template<class X> friend Vector<NegationType<X>> operator-(Vector<X> const& v) {
        return Vector<NegationType<X>>( v.size(), [&v](SizeType i){return -v[i];} ); }

    template<class X1, class X2> friend Vector<SumType<X1,X2>> operator+(Vector<X1> const& v1, Vector<X2> const& v2) {
        ARIADNE_PRECONDITION(v1.size()==v2.size());
        return Vector<SumType<X1,X2>>( v1.size(), [&v1,&v2](SizeType i){return v1[i]+v2[i];} ); }

    template<class X1, class X2> friend Vector<DifferenceType<X1,X2>> operator-(Vector<X1> const& v1, Vector<X2> const& v2) {
        ARIADNE_PRECONDITION(v1.size()==v2.size());
        return Vector<DifferenceType<X1,X2>>( v1.size(), [&v1,&v2](SizeType i){return v1[i]-v2[i];} ); }

    template<class X1, class X2> friend Vector<ProductType<Scalar<X1>,X2>> operator*(X1 const& x1, Vector<X2> const& v2) {
        return Vector<ProductType<Scalar<X1>,X2>>( v2.size(), [&x1,&v2](SizeType i){return x1*v2[i];} ); }

    template<class X1, class X2> friend Vector<ProductType<X1,Scalar<X2>>> operator*(Vector<X1> const& v1, X2 const& x2) {
        return Vector<ProductType<X1,Scalar<X2>>>( v1.size(), [&v1,&x2](SizeType i){return v1[i]*x2;} ); }

    template<class X1, class X2> friend Vector<QuotientType<X1,Scalar<X2>>> operator/(Vector<X1> const& v1, X2 const& x2) {
        return Vector<QuotientType<X1,Scalar<X2>>>( v1.size(), [&v1,&x2](SizeType i){return v1[i]/x2;} ); }

    template<class X1,class X2> friend Vector<InplaceSumType<X1,X2>>& operator+=(Vector<X1>& v1, Vector<X2> const& v2) {
        ARIADNE_PRECONDITION(v1.size()==v2.size());
        for(SizeType i=0; i!=v1.size(); ++i) { v1[i]+=v2[i]; } return v1;
    }

    template<class X1,class X2> friend Vector<InplaceDifferenceType<X1,X2>>& operator-=(Vector<X1>& v1, Vector<X2> const& v2) {
        ARIADNE_PRECONDITION(v1.size()==v2.size());
        for(SizeType i=0; i!=v1.size(); ++i) { v1[i]-=v2[i]; } return v1;
    }

    template<class X1,class X2> friend Vector<InplaceProductType<X1,X2>>& operator*=(Vector<X1>& v1, X2 const& s2) {
        for(SizeType i=0; i!=v1.size(); ++i) { v1[i]*=s2; } return v1;
    }

    template<class X1,class X2> friend Vector<InplaceQuotientType<X1,X2>>& operator/=(Vector<X1>& v1, X2 const& s2) {
        for(SizeType i=0; i!=v1.size(); ++i) { v1[i]/=s2; } return v1;
    }

    template<class X1, class X2> friend ArithmeticType<X1,X2> dot(const Vector<X1>& v1, const Vector<X2>& v2) {
        ARIADNE_PRECONDITION(v1.size()==v2.size());
        ArithmeticType<X1,X2> r=v1.zero_element()*v2.zero_element();
        for(SizeType i=0; i!=v1.size(); ++i) {
            r+=v1[i]*v2[i];
        }
        return r;
    }

    template<class X> friend decltype(abs(declval<X>())) norm(const Vector<X>& v) {
        decltype(abs(declval<X>())) r=abs(v.zero_element());
        for(SizeType i=0; i!=v.size(); ++i) {
            r=max(r,abs(v[i]));
        }
        return r;
    }

    template<class X> friend decltype(mag(declval<X>())) sup_norm(const Vector<X>& v) {
        decltype(mag(declval<X>())) r=abs(v.zero_element());
        for(SizeType i=0; i!=v.size(); ++i) {
            r=max(r,mag(v[i]));
        }
        return r;
    }

    template<class X> friend decltype(sqrt(sqr(declval<X>()))) two_norm(const Vector<X>& v) {
        decltype(sqr(declval<X>())) s=sqr(v.zero_element());
        for(SizeType i=0; i!=v.size(); ++i) {
            s=add(s,sqr(v[i]));
        }
        return sqrt(s);
    }

    template<class X1, class X2> friend EqualsType<X1,X2> operator==(const Vector<X1>& v1, const Vector<X2>& v2) {
        decltype(v1[0]==v2[0]) r(true);
        for(SizeType i=0; i!=v1.size(); ++i) {
            r = r && (v1[i]==v2[i]);
        }
        return r;
    }

    template<class X1, class X2> friend decltype(declval<X1>()!=declval<X2>()) operator!=(const Vector<X1>& v1, const Vector<X2>& v2) {
        decltype(v1[0]!=v2[0]) r(false);
        for(SizeType i=0; i!=v1.size(); ++i) {
            r = r || (v1[i]!=v2[i]);
        }
        return r;
    }

    template<class X> friend OutputStream& operator<<(OutputStream& os, Vector<X> const& v) {
        if(v.size()==0) { os << "["; }
        for(SizeType i=0; i!=v.size(); ++i) { os << (i==0u?"[":",") << v[i]; }
        return os << "]";
    }

};

template<class X1, class X2> EqualsType<X1,X2> operator==(const Vector<X1>& v1, const Vector<X2>& v2);
template<class X1, class X2> decltype(declval<X1>()!=declval<X2>()) operator!=(const Vector<X1>& v1, const Vector<X2>& v2);



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
-Woverloaded-virtual
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



template<class X>
Vector<X> join(const Vector<X>& v1, const Vector<X>& v2)
{

    Array<X> ra(v1.size()+v2.size(),Uninitialised());
    X* rp=ra.begin();
    for(X const* vp=v1.array().begin(); vp!=v1.array().end(); ++rp, ++vp) { new (rp) X(*vp); }
    for(X const* vp=v2.array().begin(); vp!=v2.array().end(); ++rp, ++vp) { new (rp) X(*vp); }
    assert(rp==ra.end());
    return Vector<X>(std::move(ra));

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


template<class X> inline decltype(error(declval<X>())) sup_error(const Vector<X>& v) {
    decltype(error(declval<X>())) r=error(v.zero_element());
    for(SizeType i=0; i!=v.size(); ++i) {
        r=max(r,error(v[i]));
    }
    return r;
}

template<class X> inline decltype(error(declval<X>())) error(const Vector<X>& v) {
    return sup_error(v);
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

template<class G> decltype(auto) generate_vector(SizeType n, G const& g) {
    typedef ResultOf<G(SizeType)> R;
    return Vector<R>(n,g);
}

template<class F, class X> Vector<ResultOf<F(X)>> elementwise(F const& f, Vector<X> const& v) {
    typedef ResultOf<F(X)> R;
    return Vector<R>(v.size(),[&](SizeType i){return f(v[i]);});
}

template<class F, class X1, class X2> Vector<ResultOf<F(X1,X2)>> elementwise(F const& f, Vector<X1> const& v1, Vector<X2> const& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    typedef ResultOf<F(X1,X2)> R;
    return Vector<R>(v1.size(),[&](SizeType i){return f(v1[i],v2[i]);});
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

template<class X, class PR> using ConcreteSingletonType = decltype(cast_singleton(declval<X>(),declval<PR>()));

template<class X> inline decltype(auto) cast_singleton(const Vector<X>& v) {
    return elementwise([&](X const& x){return cast_singleton(x);},v);
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

template<class X> template<class PR, EnableIf<IsConstructible<X,ExactDouble,PR>>>
Vector<X>::Vector(InitializerList<double> const& lst, PR pr)
    : _ary(Array<ExactDouble>(Array<double>(lst)),pr)
{
}

template<class X> template<class PR, EnableIf<IsConstructible<X,Real,PR>>>
Vector<X>::Vector(InitializerList<Real> const& lst, PR pr)
    : _ary(Array<Real>(lst),pr)
{
}

template<class X> template<class PR, EnableIf<IsConstructible<X,ExactDouble,ExactDouble,PR>>> Vector<X>::Vector(InitializerList<Pair<double,double>> const& lst, PR pr)
    : _ary(lst.size(),X(pr))
{
    SizeType i=0;
    for(auto iter=lst.begin(); iter!=lst.end(); ++iter) {
        _ary[i]=X(ExactDouble(iter->first),ExactDouble(iter->second),pr); ++i;
    }
}

} // namespace Ariadne

#include "../numeric/float.decl.hpp"
namespace Ariadne {
inline Vector<FloatDPValue>const& cast_exact(Vector<FloatDPApproximation>const& v) {
    return reinterpret_cast<Vector<FloatDPValue>const&>(v); }
} // namespace Ariadne

#endif
