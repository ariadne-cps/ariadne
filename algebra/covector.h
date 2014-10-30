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

#define SIMPLE_VECTOR_OPERATORS

#include "utility/metaprogramming.h"
#include "utility/container.h"

namespace Ariadne {

template<class X> struct Covector {
    X _zero;
    Array<X> _ary;
  public:
    Covector(SizeType n) : _zero(), _ary(n) { }
    Covector(SizeType n, const X& x) : _zero(x*0), _ary(n,x) { }
    SizeType size() const { return _ary.size(); }
    const X& operator[](SizeType j) const { return _ary.at(j); }
    X& operator[](SizeType j) { return _ary.at(j); }
    const X& at(SizeType j) const { return _ary.at(j); }
    X& at(SizeType j) { return _ary.at(j); }
    X zero_element() const { return _zero; }
    OutputStream& write(OutputStream& os) const;
};

template<class X> inline Covector<X> const& transpose(Vector<X> const& v) {
    return reinterpret_cast<Covector<X> const&>(v); }

template<class X> inline Vector<X> const& transpose(Covector<X> const& u) {
    return reinterpret_cast<Vector<X> const&>(u); }


template<class X1,class X2> Covector<SumType<X1,X2>> operator+(Covector<X1> const& u1, Covector<X2> const& u2) {
    Vector<SumType<X1,X2>> v0=transpose(u1)+transpose(u2); return transpose(v0); }

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
template<class X> OutputStream& operator<<(OutputStream& os, const Covector<X>& v) {
    return v.write(os);
}


} // namespace Ariadne

#endif
