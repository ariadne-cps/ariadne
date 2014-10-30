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
 *  \brief
 */



#ifndef ARIADNE_VECTOR_H
#define ARIADNE_VECTOR_H

#define SIMPLE_VECTOR_OPERATORS

#include "utility/metaprogramming.h"
#include "utility/container.h"

namespace Ariadne {

//! \defgroup LinearAlgebraSubModule Linear Algebra Sub-Module
//! \ingroup AlgebraModule
//! \brief %Vector and matrix classes for linear algebra.

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

template<class X> struct IsScalar { static const bool value = not (IsVectorExpression<X>::value or IsCovectorExpression<X>::value or IsMatrixExpression<X>::value); };

template<class V> using ScalarType=typename V::ScalarType;


template<class X, class = Fallback> struct HasCreateZero : False { };
template<class X> struct HasCreateZero<X, EnableIf<IsSame<decltype(declval<X>().create_zero()),X>,Fallback>> : True { };

template<class X, EnableIf<HasCreateZero<X>> = dummy> X create_zero(const X& x) { return x.create_zero(); }
template<class X, DisableIf<HasCreateZero<X>> = dummy> X create_zero(const X& x) { return static_cast<X>(0u); }

//! \ingroup LinearAlgebraSubModule
//! \brief Vectors over some type \a X.
//! Corresponds to elements of a \em module over a mathematical \em ring, or a <em>vector space</em> over a field.
//! May also be used if \a A is an \em algebra over another field.
//! It must be possible to add and multiply any two elements of the vector.
template<class X> class Vector {
    X _zero;
    Array<X> _ary;
  public:
    typedef SizeType IndexType;
    typedef X ScalarType;
    Vector() : _zero(), _ary() { }
    Vector(SizeType n) : _zero(), _ary(n) { }
    Vector(SizeType n, const X& x) : _zero(create_zero(x)), _ary(n,x) { }
    Vector(Array<X> ary) : _zero(create_zero(ary[0])), _ary(std::move(ary)) { }
    Vector(InitializerList<X> lst) : _zero(create_zero(*lst.begin())), _ary(lst.begin(),lst.end()) { }
    template<class... ARGS> Vector(SizeType n, ARGS const& ...args) : _zero(args...), _ary(n,_zero) { }
    template<class V, EnableIf<And<IsVectorExpression<V>,IsConvertible<Ariadne::ScalarType<V>,X>>> =dummy> Vector(const V& v) : _zero(v.zero_element()), _ary(v.size(),_zero) {
        for(SizeType i=0; i!=v.size(); ++i) { this->_ary[i]=v[i]; } }
    template<class V, EnableIf<  And< IsVectorExpression<V>, And<Not<IsConvertible<Ariadne::ScalarType<V>,X>>,IsConstructible<X,Ariadne::ScalarType<V>>> >  > =dummy>
        explicit Vector(const V& v) : _zero(X(v.zero_element())), _ary(v.size(),_zero) {
            for(SizeType i=0; i!=v.size(); ++i) { this->_ary[i]=X(v[i]); } }
    template<class V, EnableIf<IsVectorExpression<V>> =dummy> Vector<X>& operator=(const V& v) {
        this->_ary.resize(v.size()); for(SizeType i=0; i!=v.size(); ++i) { this->_ary[i]=v[i]; } return *this; }
    static Vector<X> unit(SizeType n, SizeType i) { Vector r(n,X(0)); r._ary[i]=1; return std::move(r); }
    SizeType size() const { return _ary.size(); }
    const X& operator[](SizeType i) const { return _ary.at(i); }
    X& operator[](SizeType i) { return _ary.at(i); }
    const X& at(SizeType i) const { return _ary.at(i); }
    X& at(SizeType i) { return _ary.at(i); }
    const X& get(SizeType i) const { return _ary.at(i); }
    void set(SizeType i, X const& x) { x+_zero; _ary.at(i)=x; }
    X zero_element() const { return _zero; }
    Array<X> const& array() const { return _ary; }
};
template<class X> struct IsVector<Vector<X>> : True { };
template<class X> OutputStream& operator<<(OutputStream& os, Vector<X> const& v) {
    if(v.size()==0) { os << "{"; }
    for(SizeType i=0; i!=v.size(); ++i) { os << (i==0u?"{":",") << v[i]; }
    return os << "}";
}

template<class V, EnableIf<IsVectorExpression<V>> =dummy> OutputStream& operator<<(OutputStream& os, const V& v) {
    typedef decltype(v[0]) X;
    return os << Vector<X>(v);
}


template<class V1, class V2, EnableIf<And<IsVectorExpression<V1>,IsVectorExpression<V2>>> =dummy>
auto operator==(const V1& v1, const V2& v2) -> decltype(v1[0]==v2[0]){
    decltype(v1[0]==v2[0]) r=true;
    for(SizeType i=0; i!=v1.size(); ++i) {
        r = r && (v1[i]==v2[i]);
    }
    return r;
}

template<class V1, class V2, EnableIf<And<IsVectorExpression<V1>,IsVectorExpression<V2>>> =dummy> auto operator!=(const V1& v1, const V2& v2) -> decltype(!(v1==v2)) {
    return !(v1==v2);
}

template<class X, class XX> inline
Vector<X>& operator+=(Vector<X>& v1, const Vector<XX>& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    for(SizeType i=0; i!=v1.size(); ++i) { v1[i]+=v2[i]; } return v1;
}

template<class X, class XX> inline
Vector<X>& operator-=(Vector<X>& v1, const Vector<XX>& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    for(SizeType i=0; i!=v1.size(); ++i) { v1[i]-=v2[i]; } return v1;
}

template<class X> inline
Vector<X>& operator*=(Vector<X>& v, const typename Vector<X>::ScalarType& s) {
    for(SizeType i=0; i!=v.size(); ++i) { v[i]*=s; } return v;
}

template<class X> inline
Vector<X>& operator/=(Vector<X>& v, const typename Vector<X>::ScalarType& s) {
    for(SizeType i=0; i!=v.size(); ++i) { v[i]/=s; } return v;
}




#ifdef SIMPLE_VECTOR_OPERATORS

template<class X> Vector<X> operator+(Vector<X> const& v) {
    return v;
}

template<class X> Vector<X> operator-(Vector<X> v) {
    for(SizeType i=0; i!=v.size(); ++i) { v[i]=-v[i]; }
    return std::move(v);
}

template<class X1, class X2> Vector<decltype(declval<X1>()+declval<X2>())> operator+(Vector<X1> const& v1, Vector<X2> const& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    Vector<decltype(declval<X1>()+declval<X2>())> r(v1.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=v1[i]+v2[i]; }
    return std::move(r);
}

template<class X1, class X2> Vector<decltype(declval<X1>()-declval<X2>())> operator-(Vector<X1> const& v1, Vector<X2> const& v2) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    Vector<decltype(declval<X1>()-declval<X2>())> r(v1.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=v1[i]-v2[i]; }
    return std::move(r);
}

template<class X1, class X2, EnableIf<IsScalar<X1>> = dummy> Vector<decltype(declval<X1>()*declval<X2>())> operator*(X1 const& x1, Vector<X2> const& v2) {
    Vector<decltype(declval<X1>()*declval<X2>())> r(v2.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=x1*v2[i]; }
    return std::move(r);
}

template<class X1, class X2, EnableIf<IsScalar<X2>> = dummy> Vector<decltype(declval<X1>()*declval<X2>())> operator*(Vector<X1> const& v1, X2 const& x2) {
    Vector<decltype(declval<X1>()*declval<X2>())> r(v1.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=v1[i]*x2; }
    return std::move(r);
}

template<class X1, class X2, EnableIf<IsScalar<X2>> = dummy> Vector<decltype(declval<X1>()/declval<X2>())> operator/(Vector<X1> const& v1, X2 const& x2) {
    Vector<decltype(declval<X1>()/declval<X2>())> r(v1.size());
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

template<class X1, class X2> inline auto dot(const Vector<X1>& v1, const Vector<X2>& v2) -> decltype(v1[0]*v2[0]) {
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    decltype(declval<X1>()*declval<X2>()+declval<X1>()*declval<X2>()) r=0u;
    for(SizeType i=0; i!=v1.size(); ++i) {
        r+=v1[i]*v2[i];
    }
    return r;
}

template<class V1, class V2> inline auto dot(const V1& v1, const V2& v2) -> decltype(v1[0]*v2[0])
{
    ARIADNE_PRECONDITION(v1.size()==v2.size());
    decltype(v1[0]*v2[0]) r=0u;
    for(SizeType i=0; i!=v1.size(); ++i) {
        r+=v1[i]*v2[i];
    }
    return r;
}



/*
template<class V1, class V2> inline auto join(const V1& v1, const V2& v2)
    -> Vector<decltype(max(v1[0],v2[0]))>
{
    const V2& v2=v2;
    typedef ScalarType<V1> X1; typedef ScalarType<V2> X2;
    typedef decltype(max(declval<X1>(),declval<X2>())) X0;
    Vector<X0> r(v1.size()+v2.size());
    for(SizeType i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(SizeType i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    return std::move(r);
}
*/

template<class V1, class V2, EnableIf<IsSame<ScalarType<V1>,ScalarType<V2>>> = dummy>
Vector<ScalarType<V1>> join(const V1& v1, const V2& v2)
{
    if(v1.size()==0) { return v2; }
    if(v2.size()==0) { return v1; }
    SizeType n1=v1.size();
    SizeType n2=v2.size();
    Vector<ScalarType<V1>> r(n1+n2,v1[0]);
    for(SizeType i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(SizeType i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    return r;
}


template<class V1, class V2, class V3, EnableIf<AreSame<ScalarType<V1>,ScalarType<V2>,ScalarType<V3>>> =dummy>
Vector<ScalarType<V1>> join(const V1& v1, const V2& v2, const V3& v3)
{

//    typedef ScalarType<V1> X1; typedef ScalarType<V2> X2; typedef ScalarType<V3> X3;
    typedef ScalarType<V1> X0;
    Vector<X0> r(v1.size()+v2.size()+v3.size());
    for(SizeType i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(SizeType i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    for(SizeType i=0; i!=v3.size(); ++i) { r[v1.size()+v2.size()+i]=v3[i]; }
    return std::move(r);
}

template<class X1, class V2, EnableIf<IsSame<X1,ScalarType<V2>>> =dummy>
Vector<ScalarType<V2>> join(const X1& x1, const V2& v2)
{
    typedef X1 X0;
    Vector<X0> r(1u+v2.size());
    r[0u]=x1;
    for(SizeType i=0; i!=v2.size(); ++i) { r[1u+i]=v2[i]; }
    return std::move(r);
}

template<class V1, class X2, EnableIf<IsSame<ScalarType<V1>,X2>> =dummy>
Vector<ScalarType<V1>> join(const V1& v1, const X2& x2)
{
    //typedef ScalarType<V1> X1;
    typedef X2 X0;
    Vector<X0> r(v1.size()+1u);
    for(SizeType i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    r[v1.size()]=x2;
    return std::move(r);
}

template<class X, EnableIf<IsScalar<X>> =dummy>
Vector<X> join(const X& x1, const X& x2)
{
    Vector<X> r(2u);
    r[0u]=x1;
    r[1u]=x2;
    return std::move(r);
}

template<class V1, class V2, class X3, EnableIf<AreSame<ScalarType<V1>,ScalarType<V2>,X3>> =dummy>
Vector<ScalarType<V1>> join(const V1& v1, const V2& v2, const X3& x3)
{
    //typedef ScalarType<V1> X1; typedef ScalarType<V2> X2;
    typedef X3 X0;
    Vector<X0> r(v1.size()+v2.size()+1u);
    for(SizeType i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(SizeType i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    r[v1.size()+v2.size()]=x3;
    return std::move(r);
}


template<class X> inline Vector<decltype(make_exact(declval<X>()))> make_exact(const Vector<X>& v) {
    Vector<decltype(make_exact(declval<X>()))> r(v.size(),make_exact(v.zero_element()));
    for(SizeType i=0; i!=v.size(); ++i) {
        r[i]=make_exact(v[i]);
    }
    return std::move(r);
}

template<class X> inline decltype(error(declval<X>())) error(const Vector<X>& v) {
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




} // namespace Ariadne

#endif
