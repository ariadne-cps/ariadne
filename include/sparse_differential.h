/***************************************************************************
 *            sparse_differential.h
 *
 *  Copyright 2008  Pieter Collins
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
 
/*! \file sparse_differential.h
 *  \brief Differential algebra variables with a sparse representation.
 */

#ifndef ARIADNE_SPARSE_DIFFERENTIAL_H
#define ARIADNE_SPARSE_DIFFERENTIAL_H

#include <map>

#include "macros.h"
#include "array.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "series.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Series;

template<class X> class SparseDifferential;
template<class D> class DifferentialVector;


template<class X> SparseDifferential<X> operator+(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> operator-(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> operator+(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
template<class X> SparseDifferential<X> operator-(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
template<class X> SparseDifferential<X> operator*(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
template<class X> SparseDifferential<X> operator/(const SparseDifferential<X>& x, const SparseDifferential<X>& y);

template<class X, class R> SparseDifferential<X> operator+(const SparseDifferential<X>& x, const R& c);
template<class X, class R> SparseDifferential<X> operator-(const SparseDifferential<X>& x, const R& c);
template<class X, class R> SparseDifferential<X> operator*(const SparseDifferential<X>& x, const R& c);
template<class X, class R> SparseDifferential<X> operator/(const SparseDifferential<X>& x, const R& c);
template<class X, class R> SparseDifferential<X> operator+(const R& c, const SparseDifferential<X>& x);
template<class X, class R> SparseDifferential<X> operator-(const R& c, const SparseDifferential<X>& x);
template<class X, class R> SparseDifferential<X> operator*(const R& c, const SparseDifferential<X>& x);
template<class X, class R> SparseDifferential<X> operator/(const R& c, const SparseDifferential<X>& x);

template<class X> SparseDifferential<X> operator+(const SparseDifferential<X>& x, const X& c);
template<class X> SparseDifferential<X> operator-(const SparseDifferential<X>& x, const X& c);
template<class X> SparseDifferential<X> operator*(const SparseDifferential<X>& x, const X& c);
template<class X> SparseDifferential<X> operator/(const SparseDifferential<X>& x, const X& c);
template<class X> SparseDifferential<X> operator+(const X& c, const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> operator-(const X& c, const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> operator*(const X& c, const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> operator/(const X& c, const SparseDifferential<X>& x);

template<class X> SparseDifferential<X> pow(const SparseDifferential<X>& x, int n);
template<class X> SparseDifferential<X> sqr(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> sqrt(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> exp(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> log(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> sin(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> cos(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> tan(const SparseDifferential<X>& x);

template<class X, class Y> Y evaluate(const SparseDifferential<X>& y, const Vector<Y>& z);
template<class X> SparseDifferential<X> compose(const Series<X>& x, const SparseDifferential<X>& y);
template<class X> SparseDifferential<X> translate(const SparseDifferential<X>& x, const Vector<X>& c);
template<class X> SparseDifferential<X> derivative(const SparseDifferential<X>& x, uint i);
template<class X> SparseDifferential<X> antiderivative(const SparseDifferential<X>& x, uint i);






/*! \brief A class representing the derivatives of a scalar variable depending on multiple arguments. */
template<class X>
class SparseDifferential
{
    static const uint MAX_DEGREE=65535;
    static const X _zero;
  public:
    typedef MultiIndex IndexType;
    typedef X ValueType;
    typedef X ScalarType;
    typedef DifferentialVector< SparseDifferential<X> > VectorType;
    typedef typename std::map<MultiIndex,X>::iterator iterator;
    typedef typename std::map<MultiIndex,X>::const_iterator const_iterator;

    explicit SparseDifferential() : _as(0), _deg(0), _data() { _data[MultiIndex(0)]=0; }
    explicit SparseDifferential(uint as) : _as(as), _deg(MAX_DEGREE), _data() { _data[MultiIndex(as)]=0; }
    explicit SparseDifferential(uint as, uint deg) : _as(as), _deg(deg), _data() { _data[MultiIndex(as)]=0; }
    template<class XX> SparseDifferential(uint as, uint deg, const XX* ptr) : _as(as), _deg(deg), _data() { 
        _data[MultiIndex(as)]=0; for(MultiIndex j(as); j.degree()<=deg; ++j) { if(*ptr!=0) { _data[j]=*ptr; } ++ptr; } }
    template<class XX> SparseDifferential(const SparseDifferential<XX>& x) : _as(x.argument_size()), _deg(x.degree()), _data() { 
        _data[MultiIndex(x.argument_size())]=0; for(typename SparseDifferential<XX>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) {
            if(iter->second!=0) { this->_data[iter->first]=X(iter->second); } } }

    SparseDifferential<X>& operator=(const X& c) { this->_data.clear(); this->_data[MultiIndex(this->argument_size())]=c; return *this; }

    SparseDifferential<X>& operator+=(const SparseDifferential<X>& x);
    SparseDifferential<X>& operator-=(const SparseDifferential<X>& x);
    template<class R> SparseDifferential<X>& operator+=(const R& c);
    template<class R> SparseDifferential<X>& operator-=(const R& c);
    template<class R> SparseDifferential<X>& operator*=(const R& c);
    template<class R> SparseDifferential<X>& operator/=(const R& c);

    void set_degree(uint d) { this->_deg = d; }

    X& operator[](const uint& j) { return this->_data[MultiIndex::unit(this->_as,j)]; }
    X& operator[](const MultiIndex& a) { ARIADNE_ASSERT(a.number_of_variables()==this->argument_size()); return this->_data[a]; }

    void set_value(const X& c) { this->operator[](MultiIndex(this->_as))=c; }
    void set_gradient(uint j, const X& d) { this->operator[](MultiIndex::unit(this->_as,j))=d; }

    const_iterator begin() const { return this->_data.begin(); }
    const_iterator end() const { return this->_data.end(); }
    uint argument_size() const { return this->_as; }
    uint degree() const { return this->_deg; }
    const X& operator[](const uint& j) const { 
        return this->operator[](MultiIndex::unit(this->_as,j)); }
    const X& operator[](const MultiIndex& a) const { 
        ARIADNE_ASSERT(a.number_of_variables()==this->argument_size()); 
        const_iterator iter=this->_data.find(a); 
        if(iter==this->_data.end()) { return _zero; } else { return iter->second; } }
    const X& value() const { return this->operator[](MultiIndex(this->_as)); }
    const X& gradient(uint j) const { return this->operator[](MultiIndex::unit(this->_as,j)); }

    static SparseDifferential<X> constant(uint as, uint d, const X& c) {
        SparseDifferential<X> r(as,d); r.set_value(c); return r; }
    static SparseDifferential<X> variable(uint as, uint d, const X& x, uint i) {
        SparseDifferential<X> r(as,d); r._data[MultiIndex::zero(as)]=x; r._data[MultiIndex::unit(as,i)]=1.0; return r; }

    static Vector< SparseDifferential<X> > constants(uint rs, uint as, uint d, const Vector<X>& c) {
        ARIADNE_ASSERT(c.size()==rs);
        Vector< SparseDifferential<X> > result(rs,SparseDifferential(as,d));
        for(uint i=0; i!=rs; ++i) { result[i]=c[i]; }
        return result;
    }
    static Vector< SparseDifferential<X> > variables(uint rs, uint as, uint d, const Vector<X>& x) {
        ARIADNE_ASSERT(x.size()==rs);  ARIADNE_ASSERT(as==x.size());
        Vector< SparseDifferential<X> > result(rs,SparseDifferential<X>(as,d));
        for(uint i=0; i!=rs; ++i) { result[i]=x[i]; result[i][i]=X(1.0); }
        return result; 
    }

    bool operator==(const SparseDifferential<X>& sd) const {
        if(this->argument_size()!=sd.argument_size()) { return false; }
        for(MultiIndex j(this->argument_size()); j.degree()<=std::max(this->degree(),sd.degree()); ++j) {
            if((*this)[j]!=sd[j]) { return false; } }
        return true;
    }
    bool operator!=(const SparseDifferential<X>& sd) const { return !(*this==sd); }

    friend SparseDifferential<X> operator+<>(const SparseDifferential<X>& x);
    friend SparseDifferential<X> operator-<>(const SparseDifferential<X>& x);
    friend SparseDifferential<X> operator+<>(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
    friend SparseDifferential<X> operator-<>(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
    friend SparseDifferential<X> operator*<>(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
    friend SparseDifferential<X> operator/<>(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
    friend SparseDifferential<X> compose<>(const Series<X>& x, const SparseDifferential<X>& y);
    friend SparseDifferential<X> derivative<>(const SparseDifferential<X>& x, uint i);
    friend SparseDifferential<X> antiderivative<>(const SparseDifferential<X>& x, uint i);
  public:
    iterator begin() { return this->_data.begin(); }
    iterator end() { return this->_data.end(); }
    std::map<MultiIndex,X>& data() { return this->_data; }
    const std::map<MultiIndex,X>& data() const { return this->_data; }
  public:
    void cleanup() { 
        typedef typename std::map<MultiIndex,X>::iterator iterator;
        iterator iter=this->_data.begin(); 
        while(iter!=this->_data.end()) {
            if(iter->second==0) { this->_data.erase(iter++); } else { ++iter; } } 
        this->_data[MultiIndex(this->_as)]; }
  private:
    uint _as;
    uint _deg;
    std::map<MultiIndex,X> _data;
  private:
    //BOOST_CONCEPT_ASSERT((DifferentialConcept< SparseDifferential<X> >));
};

template<class X>
const X SparseDifferential<X>::_zero=X(0);

template<class X>
SparseDifferential<X>& SparseDifferential<X>::operator+=(const SparseDifferential<X>& x)
{
    for(const_iterator iter=x._data.begin(); iter!=x._data.end(); ++iter) {
        this->_data[iter->first]+=iter->second;
    }
    return *this;
}

template<class X>
SparseDifferential<X>& SparseDifferential<X>::operator-=(const SparseDifferential<X>& x)
{
    for(const_iterator iter=x._data.begin(); iter!=x._data.end(); ++iter) {
        this->_data[iter->first]-=iter->second;
    }
    return *this;
}

template<class X> template<class R>
SparseDifferential<X>& SparseDifferential<X>::operator+=(const R& c)
{
    this->_data[MultiIndex(this->_as)]+=c; return *this;
}

template<class X> template<class R>
SparseDifferential<X>& SparseDifferential<X>::operator-=(const R& c)
{
    this->_data[MultiIndex(this->_as)]-=c; return *this;
}

template<class X> template<class R>
SparseDifferential<X>& SparseDifferential<X>::operator*=(const R& c)
{
    if(c==0) {
        X zero=this->_data.begin()->second; zero*=0;
        this->_data.clear();
        this->_data[MultiIndex(this->_as)]=zero;
    } else {
        for(iterator iter=this->_data.begin(); iter!=this->_data.end(); ++iter) {
            iter->second*=c;
        }
    }
    return *this;
}


template<class X> template<class R>
SparseDifferential<X>& SparseDifferential<X>::operator/=(const R& c)
{
    for(iterator iter=this->_data.begin(); iter!=this->_data.end(); ++iter) {
        iter->second/=c;
    }
    return *this;
}

template<class X>
SparseDifferential<X> operator-(const SparseDifferential<X>& x)
{
    SparseDifferential<X> r(x.argument_size(),x.degree()); r-=x; return r; 
}


template<class X, class R>
SparseDifferential<X> operator+(const SparseDifferential<X>& x, const R& c)
{
    SparseDifferential<X> r(x); r+=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator+(const R& c, const SparseDifferential<X>& x)
{
    SparseDifferential<X> r(x); r+=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator-(const SparseDifferential<X>& x, const R& c)
{
    SparseDifferential<X> r(x); r-=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator-(const R& c, const SparseDifferential<X>& x)
{
    SparseDifferential<X> r(-x); r+=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator*(const SparseDifferential<X>& x, const R& c)
{
    SparseDifferential<X> r(x); r*=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator*(const R& c, const SparseDifferential<X>& x)
{
    SparseDifferential<X> r(x); r*=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator/(const SparseDifferential<X>& x, const R& c)
{
    SparseDifferential<X> r(x); r/=X(c); return r; 
}

template<class X, class R>
SparseDifferential<X> operator/(const R& c, const SparseDifferential<X>& x)
{
    SparseDifferential<X> r=reX(c)(x); r*=c; return r; 
}


template<class X>
SparseDifferential<X> operator+(const SparseDifferential<X>& x, const X& c)
{
    SparseDifferential<X> r(x); r+=c; return r; 
}

template<class X>
SparseDifferential<X> operator+(const X& c, const SparseDifferential<X>& x)
{
    SparseDifferential<X> r(x); r+=c; return r; 
}

template<class X>
SparseDifferential<X> operator-(const SparseDifferential<X>& x, const X& c)
{
    SparseDifferential<X> r(x); r-=c; return r; 
}

template<class X>
SparseDifferential<X> operator-(const X& c, const SparseDifferential<X>& x)
{
    SparseDifferential<X> r(-x); r+=c; return r; 
}

template<class X>
SparseDifferential<X> operator*(const SparseDifferential<X>& x, const X& c)
{
    SparseDifferential<X> r(x); r*=c; return r; 
}

template<class X>
SparseDifferential<X> operator*(const X& c, const SparseDifferential<X>& x)
{
    SparseDifferential<X> r(x); r*=c; return r; 
}

template<class X>
SparseDifferential<X> operator/(const SparseDifferential<X>& x, const X& c)
{
    SparseDifferential<X> r(x); r/=c; return r; 
}

template<class X>
SparseDifferential<X> operator/(const X& c, const SparseDifferential<X>& x)
{
    SparseDifferential<X> r=rec(x); r*=c; return r; 
}






template<class X>
SparseDifferential<X> operator+(const SparseDifferential<X>& x, const SparseDifferential<X>& y)
{
    SparseDifferential<X> r(x); r+=y; r.cleanup(); return r; 
}

template<class X>
SparseDifferential<X> operator-(const SparseDifferential<X>& x, const SparseDifferential<X>& y)
{
    SparseDifferential<X> r(x); r-=y; r.cleanup(); return r; 
}

template<class X>
SparseDifferential<X> operator*(const SparseDifferential<X>& x, const SparseDifferential<X>& y)
{
    typedef typename SparseDifferential<X>::const_iterator const_iterator;
    assert(x._as==y._as);
    SparseDifferential<X> r(x._as,std::min(x._deg,y._deg));
    for(const_iterator xiter=x._data.begin(); xiter!=x._data.end(); ++xiter) {
        if(xiter->first.degree()>r.degree()) { break; }
        for(const_iterator yiter=y._data.begin(); yiter!=y._data.end(); ++yiter) {
            if(xiter->first.degree()+yiter->first.degree()>r.degree()) { break; }
            r._data[xiter->first+yiter->first]+=(xiter->second*yiter->second);
        }
    }
    r.cleanup();
    return r;
}

template<class X>
SparseDifferential<X> operator/(const SparseDifferential<X>& x, const SparseDifferential<X>& y)
{
    return x*rec(y);
}

template<class X>
SparseDifferential<X> 
min(const SparseDifferential<X>& x1, const SparseDifferential<X>& x2) 
{
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"min(SparseDifferential<X> x1, SparseDifferential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()<x2.value() ? x1 : x2;
}

  
template<class X>
SparseDifferential<X> 
max(const SparseDifferential<X>& x1,const SparseDifferential<X>& x2) 
{
    if(x1.value()==x2.value()) { 
        ARIADNE_THROW(std::runtime_error,"max(SparseDifferential<X> x1, SparseDifferential<X> x2)","x1[0]==x2[0]"); 
    }
    return x1.value()>x2.value() ? x1 : x2;
}

template<class X>
SparseDifferential<X> 
abs(const SparseDifferential<X>& x) 
{
    if(x.value()==0) { 
        ARIADNE_THROW(std::runtime_error,"abs(SparseDifferential<X> x)","x[0]==0"); 
    }
    return x.value()>0 ? pos(x) : neg(x); 
}

 
template<class X>
SparseDifferential<X> 
pos(const SparseDifferential<X>& x)
{
    return x;
}

template<class X>
SparseDifferential<X> 
neg(const SparseDifferential<X>& x)
{
    SparseDifferential<X> y(x.argument_size(),x.degree());
    for(typename SparseDifferential<X>::const_iterator iter=x.begin();
        iter!=x.end(); ++iter)
        {
            y[iter->first] = -iter->second;
        }
    return y;
}

template<class X>
SparseDifferential<X> rec(const SparseDifferential<X>& x)
{
    return compose(Series<X>::rec(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> sqr(const SparseDifferential<X>& x)
{
    return pow(x,2);
}

template<class X>
SparseDifferential<X> pow(const SparseDifferential<X>& x, int n)
{
    return compose(Series<X>::pow(x.degree(),x.value(),n),x);
}

template<class X>
SparseDifferential<X> sqrt(const SparseDifferential<X>& x)
{
    return compose(Series<X>::sqrt(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> exp(const SparseDifferential<X>& x)
{
    return compose(Series<X>::exp(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> log(const SparseDifferential<X>& x)
{
    return compose(Series<X>::log(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> sin(const SparseDifferential<X>& x)
{
    return compose(Series<X>::sin(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> cos(const SparseDifferential<X>& x)
{
    return compose(Series<X>::cos(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> tan(const SparseDifferential<X>& x)
{
    return compose(Series<X>::tan(x.degree(),x.value()),x);
}



template<class X, class Y>
Y
evaluate(const SparseDifferential<X>& y, const Vector<Y>& x) 
{
    //std::cerr<<ARIADNE_PRETTY_FUNCTION<<std::endl;
    ARIADNE_ASSERT(y.argument_size()==x.size());
    //std::cerr << "y=" << y << std::endl;
    //std::cerr << "x=" << x << std::endl;
    uint d=y.degree();
    uint ms=x.size();
    ARIADNE_ASSERT(d>=1);

    Y zero = x[0]; zero*=0;
    Y one = zero; one+=1;

    // Use inefficient brute-force approach with lots of storage...
    array< array< Y > > val(ms, array< Y >(d+1));
    for(uint j=0; j!=ms; ++j) {
        val[j][0]=one;
        val[j][1]=x[j];
        for(uint k=2; k<=d; ++k) {
            val[j][k]=val[j][k-1]*x[j];
        }
    }

    Y r(zero);
    for(typename SparseDifferential<X>::const_iterator iter=y.begin();
        iter!=y.end(); ++iter)
        {
            const MultiIndex& j=iter->first;
            const X& c=iter->second;
            Y t=one;
            for(uint k=0; k!=ms; ++k) {
                t=t*val[k][j[k]];
            }
            t*=c;
            r+=t;
            //std::cerr<<" j="<<j<<" c="<<c<<" r="<<r<<std::endl;
        }
    return r;
}


template<class X>
SparseDifferential<X> compose(const Series<X>& x, const SparseDifferential<X>& y)
{
    uint as=y.argument_size();
    uint d=std::min(x.degree(),y.degree());

    SparseDifferential<X> w=y;
    w[MultiIndex(as)]=0;
    SparseDifferential<X> r(as,d);
    r[MultiIndex(as)]=x[d];
    for(uint n=1; n<=d; ++n) {
        r=r*w; r+=x[d-n];
    }
    return r;
}


template<class X>
SparseDifferential<X> derivative(const SparseDifferential<X>& x, uint i)
{
    if(x.degree()==0) { return SparseDifferential<X>(x.argument_size(),0u); }
    SparseDifferential<X> r(x.argument_size(), x.degree()-1); 
    MultiIndex da=MultiIndex::unit(x.argument_size(),i); 
    MultiIndex ai=MultiIndex::unit(x.argument_size(),i);
    for(typename SparseDifferential<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) 
        {
            const MultiIndex& a=iter->first;
            if(a[i]!=0) { 
                da=a-ai;
                r[da]=x[a]*a[i];
            }
        }
    return r;
}

template<class X>
SparseDifferential<X> antiderivative(const SparseDifferential<X>& x, uint i)
{
    SparseDifferential<X> r(x.argument_size(), x.degree()+1); 
    MultiIndex da=MultiIndex::zero(x.argument_size()); 
    MultiIndex ai=MultiIndex::unit(x.argument_size(),i);
    for(typename SparseDifferential<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) 
        {
            const MultiIndex& a=iter->first;
            const X& xj=x[a];
            da=a+ai;
            uint dai=da[i]; r[da]=xj/dai;
            //r[da]=x[a]/da[i];
        }
    return r;
}


template<class X>
std::ostream& operator<<(std::ostream& os, const SparseDifferential<X>& x)
{
    for(typename SparseDifferential<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) {
        if(iter==x.begin()) { os << "SD("<<x.argument_size()<<","<<x.degree()<<"){"; } else { os << ","; }
        os << iter->first << ":" << iter->second ;
    }
    return os << "}";
}

} // namespace Ariadne




#include "differential_vector.h"

namespace Ariadne {

//! Translate the polynomial given by \a x to one with centre \a v.
template<class X> 
SparseDifferential<X>  
translate(const SparseDifferential<X>& x, const Vector<X>& v)
{
    ARIADNE_ASSERT(x.argument_size()==v.size());
    uint as=x.argument_size();
    uint d=x.degree();
    DifferentialVector< SparseDifferential<X> > t
        =DifferentialVector< SparseDifferential<X> >::variable(as,as,d,v);
    return evaluate(x,t);
}

} //namespace Ariadne

#endif // ARIADNE_SPARSE_DIFFERENTIAL_H
