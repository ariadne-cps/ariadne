/***************************************************************************
 *            algebra/covector.h
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

/*! \file algebra/covector.h
 *  \brief
 */



#ifndef ARIADNE_COVECTOR_H
#define ARIADNE_COVECTOR_H

#include "utility/metaprogramming.h"
#include "utility/container.h"
#include "algebra/vector.h"

namespace Ariadne {

template<class X> struct Covector {
    Array<X> _ary;
  public:
    template<class XX, EnableIf<IsConvertible<XX,X>> =dummy>
        Covector(Covector<XX> const& u) : _ary(u._ary) { }
    template<class XX, EnableIf<IsConstructible<X,XX>> =dummy, DisableIf<IsConvertible<XX,X>> =dummy>
        explicit Covector(Covector<XX> const& u) : _ary(u._ary) { }
    explicit Covector() : _ary() { }
    explicit Covector(SizeType n) : _ary(n) { }
    explicit Covector(SizeType n, const X& x) : _ary(n,x) { }
    explicit Covector(Array<X> ary) : _ary(std::move(ary)) { }
    static Covector<X> unit(SizeType n, SizeType j) { Covector<X> r(n); r[j]=1; return r; }
    SizeType size() const { return _ary.size(); }
    Void resize(SizeType n) { _ary.resize(n); }
    const X& operator[](SizeType j) const { return _ary.at(j); }
    X& operator[](SizeType j) { return _ary.at(j); }
    const X& at(SizeType j) const { return _ary.at(j); }
    X& at(SizeType j) { return _ary.at(j); }
    X zero_element() const { ARIADNE_DEBUG_ASSERT(not _ary.empty()); return create_zero(_ary[0]); }
    OutputStream& write(OutputStream& os) const;
};

template<class X> inline Covector<X> const& transpose(Vector<X> const& v) {
    return reinterpret_cast<Covector<X> const&>(v); }

template<class X> inline Vector<X> const& transpose(Covector<X> const& u) {
    return reinterpret_cast<Vector<X> const&>(u); }


template<class X1, class X2> auto operator==(const Covector<X1>& u1, const Covector<X2>& u2) -> decltype(u1[0]==u2[0]) {
    if(u1.size()!=u2.size()) { return false; }
    decltype(u1[0]==u2[0]) r=true; for(SizeType i=0; i!=u1.size(); ++i) { r = r && (u1[i]==u2[i]); } return r; }

template<class X1,class X2> Covector<SumType<X1,X2>> operator+(Covector<X1> const& u1, Covector<X2> const& u2) {
    ARIADNE_PRECONDITION(u1.size()==u2.size()); Covector<SumType<X1,X2>> r(u1.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=u1[i]+u2[i]; } return std::move(r); }

template<class X1,class X2> Covector<ProductType<X1,X2>> operator*(Covector<X1> const& u1, X2 const& s2) {
    Covector<ProductType<X1,X2>> r(u1.size()); for(SizeType i=0; i!=r.size(); ++i) { r[i]=u1[i]*s2; } return std::move(r); }

template<class X1, class X2> inline ArithmeticType<X1,X2> operator*(const Covector<X1>& u1, const Vector<X2>& v2) {
    ARIADNE_PRECONDITION(u1.size()==v2.size());
    ArithmeticType<X1,X2> r=0u;
    for(SizeType i=0; i!=v2.size(); ++i) {
        r+=u1[i]*v2[i];
    }
    return r;
}

template<class X> OutputStream& Covector<X>::write(OutputStream& os) const {
    if(size()==0) { os << "{"; }
    for(SizeType i=0; i!=size(); ++i) { os << (i==0u?"{":",") << this->_ary[i]; }
    return os << "}";
}

template<class X> OutputStream& operator<<(OutputStream& os, const Covector<X>& u) {
    return u.write(os);
}


} // namespace Ariadne

#endif
