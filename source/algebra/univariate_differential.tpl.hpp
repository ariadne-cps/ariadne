/***************************************************************************
 *            univariate_differential.tpl.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file univariate_differential.tpl.hpp
 *  \brief Differentials with respect to a single variable.
 */

#include "algebra/univariate_differential.hpp"

namespace Ariadne {



template<class X> UnivariateDifferential<X>::UnivariateDifferential()
    : _ary(1u,X(0)) { }

template<class X> UnivariateDifferential<X>::UnivariateDifferential(DegreeType d)
    : _ary(d+1,X(0)) { }

template<class X> UnivariateDifferential<X>::UnivariateDifferential(DegreeType d, X const& c)
    : _ary(d+1,nul(c)) { _ary[0]=c; }

template<class X> UnivariateDifferential<X>::UnivariateDifferential(DegreeType d, InitializerList<X> lst)
    : _ary(d+1,X(0))
{
    ARIADNE_PRECONDITION(lst.size()==d+1u);
    std::copy(lst.begin(),lst.end(),_ary.begin());
}

template<class X> UnivariateDifferential<X>::UnivariateDifferential(DegreeType d, Series<X> const& s)
    : _ary(d+1u) { for(SizeType i=0; i<=d; ++i) { this->_ary[i]=s[i]; }
}

template<class X> UnivariateDifferential<X> UnivariateDifferential<X>::constant(DegreeType d, X const& c) {
    UnivariateDifferential r(d);
    r[0]=c;
    return std::move(r);
}

template<class X> UnivariateDifferential<X> UnivariateDifferential<X>::variable(DegreeType d, X const& c) {
    UnivariateDifferential r(d);
    r[0]=c;
    if(d>=1) { r[1]=1; }
    return std::move(r);
}

template<class X> DegreeType UnivariateDifferential<X>::degree() const {
    return this->_ary.size()-1u;
}

template<class X> X const& UnivariateDifferential<X>::operator[](SizeType k) const {
    return this->_ary[k];
}

template<class X> X& UnivariateDifferential<X>::operator[](SizeType k) {
    return this->_ary[k];
}

template<class X> UnivariateDifferential<X>& UnivariateDifferential<X>::operator+=(X const& c) {
    this->_ary[0]+=c;
    return *this;
}

template<class X> UnivariateDifferential<X>& UnivariateDifferential<X>::operator*=(X const& c) {
    for(DegreeType i=0; i<=this->degree(); ++i) {
        this->_ary[i]*=c;
    }
    return *this;
}


template<class X>
UnivariateDifferential<X> UnivariateDifferential<X>::_compose(const Series<X>& x, const UnivariateDifferential<X>& y)
{
    DegreeType d=y.degree();

    X y0 = y[0];
    X z = Ariadne::create_zero(y0);
    const_cast<X&>(y[0]) = z;
    UnivariateDifferential<X> r(d, z);
    r[0]=x[d];
    for(DegreeType n=1; n<=d; ++n) {
        r=r*y;
        r[0]+=x[d-n];
    }
    const_cast<X&>(y[0])=y0;
    return r;
}

template<class X>
UnivariateDifferential<X>
UnivariateDifferential<X>::_derivative(const UnivariateDifferential<X>& x)
{
    DegreeType n=x.degree(); DegreeType one=1u;
    UnivariateDifferential<X> r(std::min(n,one)-one);
    if(n==0) { r[0]=x[0]*0; }
    for(DegreeType i=0; i<n; ++i) {
        r[i]=(i+1)*x[i+1];
    }
    return r;
}

template<class X>
UnivariateDifferential<X>
UnivariateDifferential<X>::_antiderivative(const UnivariateDifferential<X>& x)
{
    DegreeType n=x.degree();
    UnivariateDifferential<X> r(n+1,Ariadne::create_zero(x[0]));
    for(DegreeType i=0; i<=n; ++i) {
        r[i+1]=x[i]/(i+1);
    }
    return r;
}

template<class X>
UnivariateDifferential<X>
UnivariateDifferential<X>::_antiderivative(const UnivariateDifferential<X>& x, const X& c)
{
    DegreeType n=x.degree();
    UnivariateDifferential<X> r(n+1,Ariadne::create_zero(x[0]));
    r[0]=c;
    for(DegreeType i=0; i<=n; ++i) {
        r[i+1]=x[i]/(i+1);
    }
    return r;
}



template<class X>
OutputStream& UnivariateDifferential<X>::write(OutputStream& os) const
{
    UnivariateDifferential<X> const& x=*this;
    os << "D<"<<x.degree()<<">";
    for(DegreeType i=0; i<=x.degree(); ++i) {
        os << (i==0 ? '[' : ',') << x[i];
    }
    return os << "]";
}

/*
template<class X> OutputStream& UnivariateDifferential<X>::write(OutputStream& os) const {
    os << this->_ary[0];
    for(DegreeType i=1; i<=this->degree(); ++i) {
        os <<" ";
        if(decide(this->_ary[i]>=X(0))) { os << "+"; }
        os << this->_ary[i] << "*dx^" << (uint)i;
    }
    return os;
}
*/

} // namespace Ariadne
