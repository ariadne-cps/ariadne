/***************************************************************************
 *            expansion.tcc
 *
 *  Copyright 2008-15  Pieter Collins
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


#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/iterator.hpp>
#include <boost/iterator_adaptors.hpp>

#include "algebra/vector.h"
#include "algebra/multi_index.h"
#include "algebra/expansion.h"



namespace Ariadne {

template<class X, class Y> Y horner_evaluate(const Expansion<X>& e, const Vector<Y>& x)
{
    typedef typename Expansion<X>::const_iterator const_iterator;
    const uint n=e.argument_size();
    const Y z=x.zero_element(); // The zero element of the ring Y
    if(e.number_of_nonzeros()==0) { return z; }

    Array< Y > r(e.argument_size(),z); // An Array of "registers" containing working p(x[0],...,x[k])
    const_iterator iter=e.begin();
    const_iterator end=e.end();
    uint k=n;   // The current working register
    const uchar* na=iter->key().begin(); // The values of the next multi-index
    uint j=k;   // The lowest register containing a non-zero value
    X c=iter->data();
    Y t=z;
    const uchar* a=na;
    ++iter;
    while(iter!=end) {
        na=iter->key().begin();
        k=n-1;
        while(a[k]==na[k]) { --k; }
        // Since terms are ordered reverse-lexicographically,
        // previous index must have higher kth value
        assert(a[k]>na[k]);
        // Set r[k]=(((c+r[0])*x[0]^a[0]+r[1])*x[1]^a[1]+...+r[k])*x[k]^(a[k]-na[k])
        // Omit zero terms where possible
        t=numeric_cast<typename Y::NumericType>(c);
        for(uint i=0; i!=min(j,k); ++i) {
            for(uint ii=0; ii!=a[i]; ++ii) {
                t=t*x[i];
            }
        }
        for(uint i=min(j,k); i!=k; ++i) {
            t=t+r[i];
            for(uint ii=0; ii!=a[i]; ++ii) {
                t=t*x[i];
            }
            r[i]=z;
        }
        if(j<=k) {
            t=t+r[k];
        }
        for(uint ii=na[k]; ii!=a[k]; ++ii) {
            t=t*x[k];
        }
        r[k]=t;
        //std::cerr<<"a="<<MultiIndex(n,a)<<" c="<<c<<" k="<<k<<" r="<<r<<"\n";
        j=k;
        c=iter->data();
        a=na;
        ++iter;
    }
    // Set r=(((c+r[0])*x[0]^a[0]+r[1])*x[1]^a[1]+...+r[n-1])*x[n-1]^(a[n-1])
    t=numeric_cast<typename Y::NumericType>(c);
    for(uint i=0; i!=j; ++i) {
        for(uint ii=0; ii!=a[i]; ++ii) {
            t=t*x[i];
        }
    }
    for(uint i=j; i!=n; ++i) {
        t=t+r[i];
        for(uint ii=0; ii!=a[i]; ++ii) {
            t=t*x[i];
        }
    }
    //std::cerr<<"a="<<MultiIndex(n,a)<<" c="<<c<<" k="<<n<<"\n";
    //std::cerr<<"  r="<<t<<"\n";
    return t;
}

template<class X, class Y>
Y power_evaluate(const Expansion<X>& e, const Vector<Y>& y)
{
    Y zero = y.zero_element();
    Y one = zero; one+=1;

    Y r=zero;
    Y t=zero;
    for(typename Expansion<X>::const_iterator iter=e.begin();
        iter!=e.end(); ++iter)
    {
        const MultiIndex& j=iter->key();
        const X& c=iter->data();
        t=one;
        for(uint k=0; k!=e.argument_size(); ++k) {
            for(uint l=0; l!=j[k]; ++l) {
                t=t*y[k];
            }
        }
        t*=c;
        r+=t;
    }

    return r;
}


template<class X, class Y>
Y evaluate(const Expansion<X>& e, const Vector<Y>& y)
{
    return power_evaluate(e,y);
}

template<class X, class Y>
Y simple_evaluate(const Expansion<X>& e, const Vector<Y>& y)
{
    return power_evaluate(e,y);
}

template<class X, class Y>
Vector<Y> evaluate(const Vector< Expansion<X> >& x, const Vector<Y>& y)
{
    Vector<Y> r(x.size(),y.zero_element());
    for(unsigned int i=0; i!=x.size(); ++i) {
        r[i]=evaluate(x[i],y);
    }
    return r;
}

}

