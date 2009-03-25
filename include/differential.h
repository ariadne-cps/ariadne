/***************************************************************************
 *            differential.h
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
 
/*! \file differential.h
 *  \brief Differential algebra variables with a sparse representation.
 */

#ifndef ARIADNE_DIFFERENTIAL_H
#define ARIADNE_DIFFERENTIAL_H

#include <map>

#include "macros.h"
#include "array.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "series.h"
#include "expansion.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Series;

template<class X> class Expansion;
template<class X> class Differential;


template<class X> Differential<X> operator+(const Differential<X>& x);
template<class X> Differential<X> operator-(const Differential<X>& x);
template<class X> Differential<X> operator+(const Differential<X>& x, const Differential<X>& y);
template<class X> Differential<X> operator-(const Differential<X>& x, const Differential<X>& y);
template<class X> Differential<X> operator*(const Differential<X>& x, const Differential<X>& y);
template<class X> Differential<X> operator/(const Differential<X>& x, const Differential<X>& y);

template<class X, class R> Differential<X> operator+(const Differential<X>& x, const R& c);
template<class X, class R> Differential<X> operator-(const Differential<X>& x, const R& c);
template<class X, class R> Differential<X> operator*(const Differential<X>& x, const R& c);
template<class X, class R> Differential<X> operator/(const Differential<X>& x, const R& c);
template<class X, class R> Differential<X> operator+(const R& c, const Differential<X>& x);
template<class X, class R> Differential<X> operator-(const R& c, const Differential<X>& x);
template<class X, class R> Differential<X> operator*(const R& c, const Differential<X>& x);
template<class X, class R> Differential<X> operator/(const R& c, const Differential<X>& x);

template<class X> Differential<X> operator+(const Differential<X>& x, const X& c);
template<class X> Differential<X> operator-(const Differential<X>& x, const X& c);
template<class X> Differential<X> operator*(const Differential<X>& x, const X& c);
template<class X> Differential<X> operator/(const Differential<X>& x, const X& c);
template<class X> Differential<X> operator+(const X& c, const Differential<X>& x);
template<class X> Differential<X> operator-(const X& c, const Differential<X>& x);
template<class X> Differential<X> operator*(const X& c, const Differential<X>& x);
template<class X> Differential<X> operator/(const X& c, const Differential<X>& x);

template<class X> Differential<X> pow(const Differential<X>& x, int n);
template<class X> Differential<X> sqr(const Differential<X>& x);
template<class X> Differential<X> sqrt(const Differential<X>& x);
template<class X> Differential<X> exp(const Differential<X>& x);
template<class X> Differential<X> log(const Differential<X>& x);
template<class X> Differential<X> sin(const Differential<X>& x);
template<class X> Differential<X> cos(const Differential<X>& x);
template<class X> Differential<X> tan(const Differential<X>& x);

template<class X, class Y> Y evaluate(const Differential<X>& y, const Vector<Y>& z);
template<class X> Differential<X> compose(const Series<X>& x, const Differential<X>& y);
template<class X> Differential<X> derivative(const Differential<X>& x, uint i);
template<class X> Differential<X> antiderivative(const Differential<X>& x, uint i);

template<class X> Differential<X> compose(const Differential<X>&, const Vector< Differential<X> >&);
template<class X> Vector< Differential<X> > compose(const Vector< Differential<X> >&, const Vector< Differential<X> >&);



/*! \brief A class representing the derivatives of a scalar variable depending on multiple arguments. */
template<class X>
class Differential
{
    static const uint MAX_DEGREE=65535;
    static const X _zero;
  public:
    typedef MultiIndex IndexType;
    typedef X ValueType;
    typedef X ScalarType;
    typedef Vector< Differential<X> > VectorType;
    typedef typename Expansion<X>::iterator iterator;
    typedef typename Expansion<X>::const_iterator const_iterator;

    explicit Differential() : _expansion(0), _degree(0) { _expansion[MultiIndex(0)]=0; }
    explicit Differential(uint as) : _expansion(as), _degree(MAX_DEGREE) { _expansion[MultiIndex(as)]=0; }
    explicit Differential(uint as, uint deg) : _expansion(as),_degree(deg) { _expansion[MultiIndex(as)]=0; }
    explicit Differential(const std::map<MultiIndex,X>& map) : _expansion(map),_degree() { }
    template<class XX> Differential(uint as, uint deg, const XX* ptr) : _expansion(as), _degree(deg) {
        for(MultiIndex j(as); j.degree()<=deg; ++j) { if(*ptr!=0) { _expansion[j]=*ptr; } ++ptr; } }
    template<class XX> Differential(const Differential<XX>& x)
        : _expansion(x.expansion()), _degree(x.degree()) { }

    Differential<X>& operator=(const X& c) {
        this->_expansion.clear(); this->_expansion[MultiIndex(this->argument_size())]=c; return *this; }

    Differential<X>& operator+=(const Differential<X>& x);
    Differential<X>& operator-=(const Differential<X>& x);
    template<class R> Differential<X>& operator+=(const R& c);
    template<class R> Differential<X>& operator-=(const R& c);
    template<class R> Differential<X>& operator*=(const R& c);
    template<class R> Differential<X>& operator/=(const R& c);

    void set_degree(uint d) { this->_degree = d; }

    X& operator[](const uint& j) {
        return this->_expansion[MultiIndex::unit(this->argument_size(),j)]; }
    X& operator[](const MultiIndex& a) {
        ARIADNE_ASSERT(a.number_of_variables()==this->argument_size()); return this->_expansion[a]; }

    void set_value(const X& c) { this->operator[](MultiIndex(this->argument_size()))=c; }
    void set_gradient(uint j, const X& d) { this->operator[](MultiIndex::unit(this->argument_size(),j))=d; }

    const_iterator begin() const { return this->_expansion.begin(); }
    const_iterator end() const { return this->_expansion.end(); }
    uint argument_size() const { return this->_expansion.argument_size(); }
    uint degree() const { return this->_degree; }
    const Expansion<X>& expansion() const { return this->_expansion; }
    Expansion<X>& expansion() { return this->_expansion; }
    const X& operator[](const uint& j) const {
        return this->operator[](MultiIndex::unit(this->argument_size(),j)); }
    const X& operator[](const MultiIndex& a) const { 
        ARIADNE_ASSERT(a.number_of_variables()==this->argument_size()); 
        const_iterator iter=this->_expansion.find(a);
        if(iter==this->_expansion.end()) { return _zero; } else { return iter->second; } }
    const X& value() const { return this->operator[](MultiIndex(this->argument_size())); }
    const X& gradient(uint j) const { return this->operator[](MultiIndex::unit(this->argument_size(),j)); }

    static Differential<X> constant(uint as, uint d, const X& c) {
        Differential<X> r(as,d); r.set_value(c); return r; }
    static Differential<X> variable(uint as, uint d, const X& x, uint i) {
        Differential<X> r(as,d); r._expansion[MultiIndex::zero(as)]=x; r._expansion[MultiIndex::unit(as,i)]=1.0; return r; }

    static Vector< Differential<X> > constants(uint rs, uint as, uint d, const Vector<X>& c) {
        ARIADNE_ASSERT(c.size()==rs);
        Vector< Differential<X> > result(rs,Differential(as,d));
        for(uint i=0; i!=rs; ++i) { result[i]=c[i]; }
        return result;
    }
    static Vector< Differential<X> > variables(uint rs, uint as, uint d, const Vector<X>& x) {
        ARIADNE_ASSERT(x.size()==rs);  ARIADNE_ASSERT(as==x.size());
        Vector< Differential<X> > result(rs,Differential<X>(as,d));
        for(uint i=0; i!=rs; ++i) { result[i]=x[i]; result[i][i]=X(1.0); }
        return result; 
    }

    static Vector< Differential<X> > variables(uint d, const Vector<X>& x) {
        Vector< Differential<X> > result(x.size(),Differential<X>(x.size(),d));
        for(uint i=0; i!=x.size(); ++i) { result[i]=x[i]; result[i][i]=X(1.0); }
        return result; 
    }

    bool operator==(const Differential<X>& sd) const {
        if(this->argument_size()!=sd.argument_size()) { return false; }
        for(MultiIndex j(this->argument_size()); j.degree()<=std::max(this->degree(),sd.degree()); ++j) {
            if((*this)[j]!=sd[j]) { return false; } }
        return true;
    }
    bool operator!=(const Differential<X>& sd) const { return !(*this==sd); }

    friend Differential<X> operator+<>(const Differential<X>& x);
    friend Differential<X> operator-<>(const Differential<X>& x);
    friend Differential<X> operator+<>(const Differential<X>& x, const Differential<X>& y);
    friend Differential<X> operator-<>(const Differential<X>& x, const Differential<X>& y);
    friend Differential<X> operator*<>(const Differential<X>& x, const Differential<X>& y);
    friend Differential<X> operator/<>(const Differential<X>& x, const Differential<X>& y);
    friend Differential<X> compose<>(const Series<X>& x, const Differential<X>& y);
    friend Differential<X> derivative<>(const Differential<X>& x, uint i);
    friend Differential<X> antiderivative<>(const Differential<X>& x, uint i);
  public:
    iterator begin() { return this->_expansion.begin(); }
    iterator end() { return this->_expansion.end(); }
  public:
    void cleanup() { this->_expansion.cleanup(); }
  private:
    Expansion<X> _expansion;
    uint _degree;
  private:
    //BOOST_CONCEPT_ASSERT((DifferentialConcept< Differential<X> >));
};

template<class X>
const X Differential<X>::_zero=X(0);

template<class X>
Differential<X>& Differential<X>::operator+=(const Differential<X>& x)
{
    for(const_iterator iter=x._expansion.begin(); iter!=x._expansion.end(); ++iter) {
        this->_expansion[iter->first]+=static_cast<const X&>(iter->second);
    }
    return *this;
}

template<class X>
Differential<X>& Differential<X>::operator-=(const Differential<X>& x)
{
    for(const_iterator iter=x._expansion.begin(); iter!=x._expansion.end(); ++iter) {
        this->_expansion[iter->first]-=static_cast<const X&>(iter->second);
    }
    return *this;
}

template<class X> template<class R>
Differential<X>& Differential<X>::operator+=(const R& c)
{
    this->_expansion[MultiIndex(this->argument_size())]+=c; return *this;
}

template<class X> template<class R>
Differential<X>& Differential<X>::operator-=(const R& c)
{
    this->_expansion[MultiIndex(this->argument_size())]-=c; return *this;
}

template<class X> template<class R>
Differential<X>& Differential<X>::operator*=(const R& c)
{
    if(c==0) {
        X zero=this->_expansion.begin()->second; zero*=0;
        this->_expansion.clear();
        this->_expansion[MultiIndex(this->argument_size())]=zero;
    } else {
        for(iterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
            static_cast<X&>(iter->second)*=c;
        }
    }
    return *this;
}


template<class X> template<class R>
Differential<X>& Differential<X>::operator/=(const R& c)
{
    for(iterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        static_cast<X&>(iter->second)/=c;
    }
    return *this;
}

template<class X>
Differential<X> operator-(const Differential<X>& x)
{
    Differential<X> r(x.argument_size(),x.degree()); r-=x; return r;
}


template<class X, class R>
Differential<X> operator+(const Differential<X>& x, const R& c)
{
    Differential<X> r(x); r+=X(c); return r;
}

template<class X, class R>
Differential<X> operator+(const R& c, const Differential<X>& x)
{
    Differential<X> r(x); r+=X(c); return r;
}

template<class X, class R>
Differential<X> operator-(const Differential<X>& x, const R& c)
{
    Differential<X> r(x); r-=X(c); return r;
}

template<class X, class R>
Differential<X> operator-(const R& c, const Differential<X>& x)
{
    Differential<X> r(-x); r+=X(c); return r;
}

template<class X, class R>
Differential<X> operator*(const Differential<X>& x, const R& c)
{
    Differential<X> r(x); r*=X(c); return r;
}

template<class X, class R>
Differential<X> operator*(const R& c, const Differential<X>& x)
{
    Differential<X> r(x); r*=X(c); return r;
}

template<class X, class R>
Differential<X> operator/(const Differential<X>& x, const R& c)
{
    Differential<X> r(x); r/=X(c); return r;
}

template<class X, class R>
Differential<X> operator/(const R& c, const Differential<X>& x)
{
    Differential<X> r=reX(c)(x); r*=c; return r;
}


template<class X>
Differential<X> operator+(const Differential<X>& x, const X& c)
{
    Differential<X> r(x); r+=c; return r;
}

template<class X>
Differential<X> operator+(const X& c, const Differential<X>& x)
{
    Differential<X> r(x); r+=c; return r;
}

template<class X>
Differential<X> operator-(const Differential<X>& x, const X& c)
{
    Differential<X> r(x); r-=c; return r;
}

template<class X>
Differential<X> operator-(const X& c, const Differential<X>& x)
{
    Differential<X> r(-x); r+=c; return r;
}

template<class X>
Differential<X> operator*(const Differential<X>& x, const X& c)
{
    Differential<X> r(x); r*=c; return r;
}

template<class X>
Differential<X> operator*(const X& c, const Differential<X>& x)
{
    Differential<X> r(x); r*=c; return r;
}

template<class X>
Differential<X> operator/(const Differential<X>& x, const X& c)
{
    Differential<X> r(x); r/=c; return r;
}

template<class X>
Differential<X> operator/(const X& c, const Differential<X>& x)
{
    Differential<X> r=rec(x); r*=c; return r;
}






template<class X>
Differential<X> operator+(const Differential<X>& x, const Differential<X>& y)
{
    Differential<X> r(x); r+=y; r.cleanup(); return r;
}

template<class X>
Differential<X> operator-(const Differential<X>& x, const Differential<X>& y)
{
    Differential<X> r(x); r-=y; r.cleanup(); return r;
}

template<class X>
Differential<X> operator*(const Differential<X>& x, const Differential<X>& y)
{
    typedef typename Differential<X>::const_iterator const_iterator;
    assert(x.argument_size()==y.argument_size());
    Differential<X> r(x.argument_size(),std::min(x._degree,y._degree));
    MultiIndex a(x.argument_size());
    X c(0);
    for(const_iterator xiter=x._expansion.begin(); xiter!=x._expansion.end(); ++xiter) {
        if(xiter->first.degree()>r.degree()) { break; }
        for(const_iterator yiter=y._expansion.begin(); yiter!=y._expansion.end(); ++yiter) {
            if(xiter->first.degree()+yiter->first.degree()>r.degree()) { break; }
            a=xiter->first+yiter->first;
            c=static_cast<const X&>(xiter->second)*static_cast<const X&>(yiter->second);
            r._expansion[a]+=c;
        }
    }
    r.cleanup();
    return r;
}

template<class X>
Differential<X> operator/(const Differential<X>& x, const Differential<X>& y)
{
    return x*rec(y);
}

template<class X>
Differential<X>
min(const Differential<X>& x1, const Differential<X>& x2)
{
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"min(Differential<X> x1, Differential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()<x2.value() ? x1 : x2;
}

  
template<class X>
Differential<X>
max(const Differential<X>& x1,const Differential<X>& x2)
{
    if(x1.value()==x2.value()) { 
        ARIADNE_THROW(std::runtime_error,"max(Differential<X> x1, Differential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()>x2.value() ? x1 : x2;
}

template<class X>
Differential<X>
abs(const Differential<X>& x)
{
    if(x.value()==0) { 
        ARIADNE_THROW(std::runtime_error,"abs(Differential<X> x)","x[0]==0");
    }
    return x.value()>0 ? pos(x) : neg(x); 
}

 
template<class X>
Differential<X>
pos(const Differential<X>& x)
{
    return x;
}

template<class X>
Differential<X>
neg(const Differential<X>& x)
{
    Differential<X> y(x.argument_size(),x.degree());
    for(typename Differential<X>::const_iterator iter=x.begin();
        iter!=x.end(); ++iter)
        {
            y[iter->first] = -iter->second;
        }
    return y;
}

template<class X>
Differential<X> rec(const Differential<X>& x)
{
    return compose(Series<X>::rec(x.degree(),x.value()),x);
}

template<class X>
Differential<X> sqr(const Differential<X>& x)
{
    return pow(x,2);
}

template<class X>
Differential<X> pow(const Differential<X>& x, int n)
{
    return compose(Series<X>::pow(x.degree(),x.value(),n),x);
}

template<class X>
Differential<X> sqrt(const Differential<X>& x)
{
    return compose(Series<X>::sqrt(x.degree(),x.value()),x);
}

template<class X>
Differential<X> exp(const Differential<X>& x)
{
    return compose(Series<X>::exp(x.degree(),x.value()),x);
}

template<class X>
Differential<X> log(const Differential<X>& x)
{
    return compose(Series<X>::log(x.degree(),x.value()),x);
}

template<class X>
Differential<X> sin(const Differential<X>& x)
{
    return compose(Series<X>::sin(x.degree(),x.value()),x);
}

template<class X>
Differential<X> cos(const Differential<X>& x)
{
    return compose(Series<X>::cos(x.degree(),x.value()),x);
}

template<class X>
Differential<X> tan(const Differential<X>& x)
{
    return compose(Series<X>::tan(x.degree(),x.value()),x);
}



template<class X, class Y>
Y
evaluate(const Differential<X>& y, const Vector<Y>& x)
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
    for(typename Differential<X>::const_iterator iter=y.begin();
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
Differential<X> compose(const Series<X>& x, const Differential<X>& y)
{
    uint as=y.argument_size();
    uint d=std::min(x.degree(),y.degree());

    Differential<X> w=y;
    w[MultiIndex(as)]=0;
    Differential<X> r(as,d);
    r[MultiIndex(as)]=x[d];
    for(uint n=1; n<=d; ++n) {
        r=r*w; r+=x[d-n];
    }
    return r;
}


template<class X>
Differential<X> derivative(const Differential<X>& x, uint i)
{
    if(x.degree()==0) { return Differential<X>(x.argument_size(),0u); }
    Differential<X> r(x.argument_size(), x.degree()-1);
    MultiIndex a(x.argument_size());
    for(typename Differential<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) {
        a=iter->first;
        unsigned int n=a[i];
        if(n!=0) {
            const X& xc=x[a];
            --a[i];
            X& rc=r[a];
            rc=xc*n;
        }
    }
    return r;
}

template<class X>
Differential<X> antiderivative(const Differential<X>& x, uint i)
{
    Differential<X> r(x.argument_size(), x.degree()+1);
    MultiIndex a(x.argument_size());
    MultiIndex ra=MultiIndex(x.argument_size());
    MultiIndex ai=MultiIndex::unit(x.argument_size(),i);
    for(typename Differential<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) {
        a=iter->first;
        const X& xc=x[a];
        ++a[i];
        unsigned int n=a[i];
        X& rc=r[a];
        rc=xc/n;
    }
    return r;
}


template<class X>
std::ostream& operator<<(std::ostream& os, const Differential<X>& x)
{
    for(typename Differential<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) {
        if(iter==x.begin()) { os << "SD("<<x.argument_size()<<","<<x.degree()<<"){"; } else { os << ","; }
        os << MultiIndex(iter->first) << ":" << X(iter->second) ;
    }
    return os << "}";
}


//! Embed the arguments in a space of dimension \a size, starting at position \a start.
template<class X>
Differential<X>
embed(const Differential<X>& x, 
      uint size, uint start)
{  
    assert(start+x.argument_size()<=size);
    Differential<X> r(size,x.degree());
    MultiIndex jr(size);
    for(typename Differential<X>::const_iterator iter=x.begin();
        iter!=x.end(); ++iter)
        {
            const MultiIndex& jx=iter->first;
            for(uint k=0; k!=x.argument_size(); ++k) {
                jr[start+k]=jx[k];
            }
            r[jr]=iter->second;
        }
    return r;
}






/*! \brief A class representing the derivatives of a vector quantity depending on multiple arguments. */
template<class X>
class Vector< Differential<X> >
    : public ublas::vector< Differential<X> >
{
    //BOOST_CONCEPT_ASSERT((DifferentialVectorConcept<DifferentialVector<X> >));
  public:
    // The type used for accessing elements
    typedef uint IndexType;
    // The value stored in the vector.
    typedef Differential<X> ValueType;
    // The type used for scalars.
    typedef typename Differential<X>::ScalarType ScalarType;

    Vector() : ublas::vector< Differential<X> >(0) { }
    Vector(uint rs) : ublas::vector< Differential<X> >(rs) { }
    Vector(uint rs, uint as, uint d) : ublas::vector< Differential<X> >(rs) {
        for(uint i=0; i!=rs; ++i) { (*this)[i]=Differential<X>(as,d); } }
    Vector(uint rs, const Differential<X>& sd) : ublas::vector< Differential<X> >(rs) {
        for(uint i=0; i!=rs; ++i) { (*this)[i]=sd; } }
    template<class XX> Vector(const Vector< Differential<XX> > dv);
    template<class XX> Vector(uint rs, uint as, uint d, const XX* ptr);
    Vector(uint rs, uint as, uint d,const Vector<X>& v, const Matrix<X>& A);
    template<class E> Vector(const ublas::vector_expression<E>& ve)
        : ublas::vector< Differential<X> >(ve) { }
    template<class E> Vector< Differential<X> >& operator=(const ublas::vector_expression<E>& ve) {
        ublas::vector< Differential<X> >::operator=(ve); return *this; }


    uint result_size() const { return this->size(); }
    uint argument_size() const { return (*this)[0].argument_size(); }
    uint degree() const { return (*this)[0].degree(); }

    Vector<X> value() const {
        Vector<X> r(this->result_size());
        for(uint i=0; i!=r.size(); ++i) { r[i]=(*this)[i].value(); } return r; }
    Matrix<X> jacobian() const { Matrix<X> r(this->result_size(),this->argument_size());
        for(uint i=0; i!=r.row_size(); ++i) { for(uint j=0; j!=r.column_size(); ++j) { r[i][j]=(*this)[i].gradient(j); } } return r; }

    void set_value(const Vector<X>& c) {
        ARIADNE_ASSERT(this->result_size()==c.size());
        for(uint i=0; i!=c.size(); ++i) { (*this)[i].set_value(c[i]); } }

    static Vector< Differential<X> > constant(uint rs, uint as, uint d, const Vector<X>& c) {
        ARIADNE_ASSERT(c.size()==rs);
        Vector< Differential<X> > result(rs,as,d);
        for(uint i=0; i!=rs; ++i) { result[i]=c[i]; }
        return result;
    }

    static Vector< Differential<X> > variable(uint rs, uint as, uint d, const Vector<X>& x) {
        ARIADNE_ASSERT(x.size()==rs);
        Vector< Differential<X> > result(rs,as,d);
        for(uint i=0; i!=rs; ++i) { result[i]=x[i]; result[i][i]=X(1.0); }
        return result;
    }

    static Vector< Differential<X> > affine(uint rs, uint as, uint d, const Vector<X>& b, const Matrix<X>& A) {
        ARIADNE_ASSERT(b.size()==rs);
        ARIADNE_ASSERT(A.row_size()==rs);
        ARIADNE_ASSERT(A.column_size()==as);
        Vector< Differential<X> > result(rs,as,d);
        for(uint i=0; i!=rs; ++i) { 
            result[i]=b[i]; 
            for(uint j=0; j!=as; ++j) {
                result[i][j]=A[i][j]; 
            }
        }
        return result;
    }

};


template<class X> template<class XX>
Vector< Differential<X> >::Vector(uint rs, uint as, uint d, const XX* ptr)
    : ublas::vector< Differential<X> >(rs)
{ 
    for(uint i=0; i!=rs; ++i) {
    (*this)[i]=Differential<X>(as,d);
        for(MultiIndex j(as); j.degree()<=d; ++j) {
            if(*ptr!=0) { (*this)[i][j]=*ptr; } ++ptr; } }
}

template<class X>
Vector< Differential<X> >::Vector(uint rs, uint as, uint d,
                                  const Vector<X>& v, const Matrix<X>& A)
    :  ublas::vector< Differential<X> >(rs,Differential<X>())
{
    ARIADNE_ASSERT(rs==v.size());
    ARIADNE_ASSERT(rs==A.row_size());
    ARIADNE_ASSERT(as==A.column_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        (*this)[i]=v[i];
        for(uint j=0; j!=this->argument_size(); ++j) {
            const X& x=A[i][j];
            if(x!=0) { (*this)[i][j]=x; }
        }
    }
}


/*! \brief Compose the series of \a x and the series of \a y, assuming that \a x is centred at the value of \a y. The value of \a y is therefore unused by this method. */
template<class X, class Y>
Vector<Y>
evaluate(const Vector< Differential<X> >& x,
         const Vector<Y>& y)
{  
    Vector<Y> r(x.result_size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=evaluate(x[i],y);
    }
    return r;
}

/*! \brief Compose the series of \a x and the series of \a y, assuming that \a x is centred at the value of \a y. The value of \a y is therefore unused by this method. */
template<class X>
Differential<X> 
compose(const Differential<X>& x, 
        const Vector< Differential<X> >& y)
{  
    Vector<X> yv=y.value();
    Vector< Differential<X> >& ync=const_cast<Vector< Differential<X> >&>(y);
    for(uint i=0; i!=ync.result_size(); ++i) { ync[i].set_value(0); }
    Differential<X> r=evaluate(x,ync);
    ync+=yv;
    return r;
}


/*! \brief Compose the series of \a x and the series of \a y, assuming that \a x is centred at the value of \a y. The value of \a y is therefore unused by this method. */
template<class X>
Vector< Differential<X> > 
compose(const Vector< Differential<X> >& x,
        const Vector< Differential<X> >& y)
{  
    //std::cerr<<"compose(DV x, DV y)\n x="<<x<<"\n y="<<y<<std::endl;
    Vector<X> yv=y.value();
    Vector< Differential<X> >& ync=const_cast<Vector< Differential<X> >&>(y);
    for(uint i=0; i!=ync.result_size(); ++i) { ync[i].set_value(0); }
    Vector< Differential<X> > r=evaluate(x,ync);
    //std::cerr<<"r="<<r<<"\n"<<std::endl;
    ync+=yv;
    return r;
}

template<class X>
Vector<X> 
value(const Vector< Differential<X> >& x)
{
    Vector<X> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].value();
    }
    return r;
}

template<class X>
Matrix<X> 
jacobian(const Vector< Differential<X> >& x)
{
    if(x.size()==0) { return Matrix<X>(); }
    for(uint i=1; i!=x.size(); ++i) {
        ARIADNE_ASSERT(x[i].argument_size()==x[0].argument_size());
    }

    Matrix<X> r(x.size(),x[0].argument_size());
    for(uint i=0; i!=x.size(); ++i) {
        for(uint j=0; j!=x[0].argument_size(); ++j) {
            r[i][j]=x[i].gradient(j);

        }
    }
    return r;
}






} //namespace Ariadne

#endif // ARIADNE_DIFFERENTIAL_H
