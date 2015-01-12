/***************************************************************************
 *            differential.tcc
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

#include <map>

#include "utility/macros.h"
#include "utility/array.h"
#include "numeric/float.decl.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/multi_index.h"
#include "algebra/series.h"
#include "algebra/expansion.h"
#include <boost/concept_check.hpp>

namespace Ariadne {

class ExactInterval;

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Series;

template<class X> class Expansion;
template<class X> class Differential;
template<class X> class Vector<Differential<X>>;

typedef Differential<Float> FloatDifferential;
typedef Differential<ExactInterval> ExactIntervalDifferential;
typedef Differential<UpperInterval> UpperIntervalDifferential;
typedef Vector< Differential<Float> > FloatDifferentialVector;
typedef Vector< Differential<ExactInterval> > ExactIntervalDifferentialVector;
typedef Vector< Differential<UpperInterval> > UpperIntervalDifferentialVector;

template<class X> Differential<X>& operator+=(Differential<X>& x, const Differential<X>& y);
template<class X> Differential<X>& operator-=(Differential<X>& x, const Differential<X>& y);

template<class X, class R> Differential<X>& operator+=(Differential<X>& x, const R& c);
template<class X, class R> Differential<X>& operator-=(Differential<X>& x, const R& c);
template<class X, class R> Differential<X>& operator*=(Differential<X>& x, const R& c);
template<class X, class R> Differential<X>& operator/=(Differential<X>& x, const R& c);

template<class X> Differential<X> operator+(const Differential<X>& x);
template<class X> Differential<X> operator-(const Differential<X>& x);
template<class X> Differential<X> operator+(const Differential<X>& x, const Differential<X>& y);
template<class X> Differential<X> operator-(const Differential<X>& x, const Differential<X>& y);
template<class X> Differential<X> operator*(const Differential<X>& x, const Differential<X>& y);
template<class X> Differential<X> operator/(const Differential<X>& x, const Differential<X>& y);

template<class X> Differential<X> neg(const Differential<X>& x);
template<class X> Differential<X> rec(const Differential<X>& x);
template<class X> Differential<X> pow(const Differential<X>& x, Int n);
template<class X> Differential<X> sqr(const Differential<X>& x);
template<class X> Differential<X> sqrt(const Differential<X>& x);
template<class X> Differential<X> exp(const Differential<X>& x);
template<class X> Differential<X> log(const Differential<X>& x);
template<class X> Differential<X> sin(const Differential<X>& x);
template<class X> Differential<X> cos(const Differential<X>& x);
template<class X> Differential<X> tan(const Differential<X>& x);

template<class X, class Y> Y evaluate(const Differential<X>& y, const Vector<Y>& z);
template<class X> Differential<X> compose(const Series<X>& x, const Differential<X>& y);
template<class X> Differential<X> derivative(const Differential<X>& x, SizeType i);
template<class X> Differential<X> antiderivative(const Differential<X>& x, SizeType i);

template<class X> Differential<X> compose(const Differential<X>&, const Vector<Differential<X>>&);
template<class X> Vector<Differential<X>> compose(const Vector<Differential<X>>&, const Vector<Differential<X>>&);

template<class X> Differential<X>::Differential() : _expansion(0), _degree(0) { }

template<class X> Differential<X>::Differential(SizeType as, DegreeType deg) : _expansion(as),_degree(deg) { }

template<class X> Differential<X>::Differential(const Map<MultiIndex,X>& map, DegreeType deg) : _expansion(0u) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class X> Differential<X>::Differential(const Expansion<X>& e, DegreeType deg) : _expansion(e.argument_size()),_degree(deg) {
    for(typename Expansion<X>::ConstIterator iter=e.begin(); iter!=e.end(); ++iter) {
        if(iter->key().degree()<=deg) { this->_expansion.append(iter->key(),iter->data()); }
    }
    this->cleanup();
}


template<class X> Differential<X>& Differential<X>::operator=(const X& c) {
    this->_expansion.clear();
    this->_expansion.append(MultiIndex(this->argument_size()),c);
    return *this;
}


template<class X> Differential<X> Differential<X>::constant(SizeType as, DegreeType deg, const X& c) {
    Differential<X> r(as,deg);
    r._expansion.append(MultiIndex(as),c);
    return r;
}

template<class X> Differential<X> Differential<X>::variable(SizeType as, DegreeType deg, const X& v, SizeType j) {
    Differential<X> r(as,deg);
    MultiIndex a(as);
    r._expansion.append(a,v);
    a[j]=1;
    r._expansion.append(a,_one);
    return r;
}



template<class X> Vector<Differential<X>> Differential<X>::constants(SizeType rs, SizeType as, DegreeType deg, const Vector<X>& c) {
    ARIADNE_ASSERT(c.size()==rs);
    Vector<Differential<X>> result(rs,Differential(as,deg));
    for(SizeType i=0; i!=rs; ++i) { result[i]=c[i]; }
    return result;
}


template<class X> Vector<Differential<X>> Differential<X>::variables(SizeType rs, SizeType as, DegreeType deg, const Vector<X>& x) {
    ARIADNE_ASSERT(x.size()==rs);  ARIADNE_ASSERT(as==x.size());
    return variables(deg,x);
}


template<class X> Vector<Differential<X>> Differential<X>::variables(DegreeType deg, const Vector<X>& x) {
    Vector<Differential<X>> result(x.size(),Differential<X>(x.size(),deg));
    MultiIndex a(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        result[i]._expansion.append(a,x[i]);
        a[i]=1; result[i]._expansion.append(a,_one); a[i]=0;
    }
    return result;
}


template<class X> InequalityType<X> Differential<X>::operator!=(const Differential<X>& other) const {
    return !(*this==other);
}

template<class X> SizeType Differential<X>::argument_size() const {
    return this->_expansion.argument_size();
}

template<class X> DegreeType Differential<X>::degree() const {
    return this->_degree;
}

template<class X> const Expansion<X>& Differential<X>::expansion() const {
    return this->_expansion;
}

template<class X> Expansion<X>& Differential<X>::expansion() {
    return this->_expansion;
}

template<class X> const X& Differential<X>::value() const {
    return this->operator[](MultiIndex(this->argument_size()));
}

template<class X> const X& Differential<X>::gradient(SizeType j) const {
    return this->operator[](MultiIndex::unit(this->argument_size(),j));
}

template<class X> Vector<X> Differential<X>::gradient() const {
    Vector<X> g(this->argument_size());
    for(SizeType j=0; j!=g.size(); ++j) { g[j]=this->gradient(j); }
    return g;
}



template<class X> X& Differential<X>::operator[](const SizeType& j) {
    return this->operator[](MultiIndex::unit(this->argument_size(),j));
}

template<class X> const X& Differential<X>::operator[](const SizeType& j) const {
    return this->operator[](MultiIndex::unit(this->argument_size(),j));
}

template<class X> X& Differential<X>::operator[](const MultiIndex& a) {
    ARIADNE_ASSERT_MSG(a.number_of_variables()==this->argument_size()," d="<<*this<<", a="<<a);
    return this->_expansion.at(a);
}

template<class X> const X& Differential<X>::operator[](const MultiIndex& a) const {
    ARIADNE_ASSERT_MSG(a.number_of_variables()==this->argument_size()," d="<<*this<<", a="<<a);
    ConstIterator iter=this->_expansion.find(a);
    if(iter==this->_expansion.end()) { return _zero; }
    else { return iter->data(); }
}


template<class X> Void Differential<X>::set_degree(DegreeType deg) {
    this->_degree = deg;
}

template<class X> Void Differential<X>::set_value(const X& c) {
    this->operator[](MultiIndex(this->argument_size()))=c;
}

template<class X> Void Differential<X>::set_gradient(SizeType j, const X& d) {
    this->operator[](MultiIndex::unit(this->argument_size(),j))=d;
}


template<class X> typename Differential<X>::Iterator Differential<X>::begin() {
    return this->_expansion.begin();
}

template<class X> typename Differential<X>::Iterator Differential<X>::end() {
    return this->_expansion.end();
}

template<class X> typename Differential<X>::ConstIterator Differential<X>::begin() const {
    return this->_expansion.begin();
}

template<class X> typename Differential<X>::ConstIterator Differential<X>::end() const {
    return this->_expansion.end();
}


template<class X> Differential<X> Differential<X>::create() const {
    return Differential<X>(this->argument_size(),this->degree());
}

template<class X> Differential<X> Differential<X>::create_zero() const {
    return Differential<X>(this->argument_size(),this->degree());
}

template<class X> Void Differential<X>::clear() {
    this->_expansion.clear();
}

template<class X> Void Differential<X>::cleanup() {
    this->_expansion.graded_sort();
    this->_expansion.combine_terms();
    this->_expansion.remove_zeros();
}

//template<class X> Void Differential<X>::check() const;


template<class X>
const X Differential<X>::_zero=X(0);

template<class X>
const X Differential<X>::_one=X(1);

template<class X>
Differential<X>::Differential(SizeType as, DegreeType deg,
                              InitializerList< PairType<InitializerList<Int>,X> > lst)
    : _expansion(Expansion<X>(as,lst)), _degree(deg)
{
    this->cleanup();
}

template<class X>
Void
Differential<X>::check() const
{
    for(typename Differential<X>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        ARIADNE_ASSERT_MSG(iter->key().degree()<=this->degree(), *this);
        typename Differential<X>::ConstIterator next = iter; ++next;
        ARIADNE_ASSERT_MSG(graded_less(iter->key(),next->key()),"Error in ordering Differential "<<this->expansion());
    }
}

/*
template<class X>
Differential<X>& Differential<X>::operator+=(const Differential<X>& x)
{
    for(ConstIterator iter=x._expansion.begin(); iter!=x._expansion.end(); ++iter) {
        this->_expansion[iter->key()]+=static_cast<const X&>(iter->data());
    }
    return *this;
}

template<class X>
Differential<X>& Differential<X>::operator-=(const Differential<X>& x)
{
    for(ConstIterator iter=x._expansion.begin(); iter!=x._expansion.end(); ++iter) {
        this->_expansion[iter->key()]-=static_cast<const X&>(iter->data());
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
        X zero=this->_expansion.begin()->data(); zero*=0;
        this->_expansion.clear();
        this->_expansion[MultiIndex(this->argument_size())]=zero;
    } else {
        for(Iterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
            static_cast<X&>(iter->data())*=c;
        }
    }
    return *this;
}


template<class X> template<class R>
Differential<X>& Differential<X>::operator/=(const R& c)
{
    for(Iterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        static_cast<X&>(iter->data())/=c;
    }
    return *this;
}
*/


template<class X> Matrix<X> Differential<X>::hessian() const {
    ARIADNE_PRECONDITION(this->degree()>=2);
    Matrix<X> H(this->argument_size(),this->argument_size());
    ConstIterator iter=this->begin();
    while(iter!=this->end() && iter->key().degree()<=1) { ++iter; }
    SizeType i=0; SizeType j=1;
    while(iter!=this->end() && iter->key().degree()<=2) {
        const MultiIndex& a=iter->key(); const X& c=iter->data();
        while(a[i]==0) { ++i; j=i+1; }
        if(a[i]==2) { H[i][i]=c*2; }
        else { while(a[j]==0) { ++j; } H[i][j]=c*2; H[j][i]=c*2; }
        ++iter;
    }
    return H;
}

template<class X> EqualityType<X> Differential<X>::operator==(const Differential<X>& other) const {
    Differential<X> const& self=*this;
    if(self.argument_size()!=other.argument_size()) { return false; }
    EqualityType<X> result;
    ComparisonType less;
    ConstIterator self_iter=self.begin(); ConstIterator other_iter=other.begin();
    while(self_iter!=self.end() && other_iter!=other.end()) {
        if(self_iter->key()==other_iter->key()) {
            result = result && (self_iter->data()==other_iter->data());
            ++self_iter; ++other_iter;
        } else if (less(self_iter->key(), other_iter->key())) {
            result = result && (self_iter->data()==_zero);
            ++self_iter;
        } else { // (self_iter->key() > other_iter->key())
            result = result && (other_iter->data()==_zero);
            ++other_iter;
        }
    }
    while(self_iter!=self.end()) {
        result = result && (self_iter->data()==_zero);
        ++self_iter;
    }
    while(other_iter!=other.end()) {
        result = result && (other_iter->data()==_zero);
        ++other_iter;
    }
    return result;
}



template<class X> Differential<X> operator+(const Differential<X>& x, const Differential<X>& y);
template<class X> Differential<X> operator-(const Differential<X>& x, const Differential<X>& y);



template<class X>
Differential<X>& operator+=(Differential<X>& x, const Differential<X>& y)
{
    x=x+y;
    return x;
}

template<class X>
Differential<X>& operator-=(Differential<X>& x, const Differential<X>& y)
{
    x=x-y;
    return x;
}

template<class X, class R>
Differential<X>& operator+=(Differential<X>& x, const R& c)
{
    MultiIndex a(x.argument_size());
    if(x.expansion().empty()) {
        x.expansion().append(a,c);
    } else if(x.begin()->key()!=a) {
        x.expansion().prepend(a,c);
    } else {
        x.begin()->data()+=c;
    }
    return x;
}

template<class X, class R>
Differential<X>& operator-=(Differential<X>& x, const R& c)
{
    MultiIndex a(x.argument_size());
    if(x.expansion().empty()) {
        x.expansion().append(a,-c);
    } else if(x.begin()->key()!=a) {
        x.expansion().prepend(a,-c);
    } else {
        x.begin()->data()-=c;
    }
    return x;
}

template<class X, class R>
Differential<X>& operator*=(Differential<X>& x, const R& c)
{
    typedef typename Differential<X>::Iterator Iterator;
    if(c==static_cast<X>(0)) {
        x.clear();
    } else {
        for(Iterator iter=x.begin(); iter!=x.end(); ++iter) {
            static_cast<X&>(iter->data())*=c;
        }
    }
    return x;
}


template<class X, class R>
Differential<X>& operator/=(Differential<X>& x, const R& c)
{
    typedef typename Differential<X>::Iterator Iterator;
    for(Iterator iter=x.begin(); iter!=x.end(); ++iter) {
        static_cast<X&>(iter->data())/=static_cast<X>(c);
    }
    return x;
}


template<class X>
Differential<X> operator+(const Differential<X>& x)
{
    Differential<X> r(x.argument_size(),x.degree());
    r.expansion().reserve(x.expansion().number_of_nonzeros());
    for(typename Differential<X>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter) {
        r.expansion().append(iter->key(), +iter->data());
    }
    return r;
}

template<class X>
Differential<X> operator-(const Differential<X>& x)
{
    Differential<X> r(x.argument_size(),x.degree());
    r.expansion().reserve(x.expansion().number_of_nonzeros());
    for(typename Differential<X>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter) {
        r.expansion().append(iter->key(), -iter->data());
    }
    return r;
}


template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator+(const Differential<X>& x, const R& c)
{
    Differential<X> r(x); r+=X(c); return r;
}


template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator+(const R& c, const Differential<X>& x)
{
    Differential<X> r(x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator-(const Differential<X>& x, const R& c)
{
    Differential<X> r(x); r-=X(c); return r;
}

template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator-(const R& c, const Differential<X>& x)
{
    Differential<X> r(-x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator*(const Differential<X>& x, const R& c)
{
    Differential<X> r(x); r*=X(c); return r;
}

template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator*(const R& c, const Differential<X>& x)
{
    Differential<X> r(x); r*=X(c); return r;
}

template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator/(const Differential<X>& x, const R& c)
{
    Differential<X> r(x); r/=X(c); return r;
}

template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator/(const R& c, const Differential<X>& x)
{
    Differential<X> r(rec(x)); r*=X(c); return r;
}








template<class X>
Differential<X> operator+(const Differential<X>& x, const Differential<X>& y)
{
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<" y="<<y);
    Differential<X> r(x.argument_size(),std::min(x.degree(),y.degree()));
    typename Differential<X>::ConstIterator xiter=x.begin();
    typename Differential<X>::ConstIterator yiter=y.begin();
    // No need to check if maximum degree has been reached below,
    // since if one Iterator is above the maximum degree, the other is at the end.
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()==yiter->key()) {
            r.expansion().append(xiter->key(),xiter->data()+yiter->data());
            ++xiter; ++yiter;
        } else if(graded_less(xiter->key(),yiter->key())) {
            r.expansion().append(xiter->key(),xiter->data());
            ++xiter;
        } else {
            r.expansion().append(yiter->key(),yiter->data());
            ++yiter;
        }
    }
    while(xiter!=x.end() && xiter->key().degree()<=r.degree()) {
        r.expansion().append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y.end() && yiter->key().degree()<=r.degree()) {
        r.expansion().append(yiter->key(),yiter->data());
        ++yiter;
    }
    //std::cerr<<"x="<<x<<" y="<<y<<" x+y="<<r<<"\n";
    return r;
}

template<class X>
Differential<X> operator-(const Differential<X>& x, const Differential<X>& y)
{
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<" y="<<y);
    Differential<X> r(x.argument_size(),std::min(x.degree(),y.degree()));
    typename Differential<X>::ConstIterator xiter=x.begin();
    typename Differential<X>::ConstIterator yiter=y.begin();
    // No need to check if maximum degree has been reached below,
    // since if one Iterator is above the maximum degree, the other is at the end.
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()==yiter->key()) {
            r.expansion().append(xiter->key(),xiter->data()-yiter->data());
            ++xiter; ++yiter;
        } else if(graded_less(xiter->key(),yiter->key())) {
            r.expansion().append(xiter->key(),xiter->data());
            ++xiter;
        } else {
            r.expansion().append(yiter->key(),-yiter->data());
            ++yiter;
        }
    }
    while(xiter!=x.end() && xiter->key().degree()<=r.degree()) {
        r.expansion().append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y.end() && yiter->key().degree()<=r.degree()) {
        r.expansion().append(yiter->key(),-yiter->data());
        ++yiter;
    }
    //std::cerr<<"x="<<x<<" y="<<y<<" x-y="<<r<<"\n";
    return r;
}

template<class X>
Differential<X> operator*(const Differential<X>& x, const Differential<X>& y)
{
    typedef typename Differential<X>::ConstIterator ConstIterator;
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<" y="<<y);
    Differential<X> r(x.argument_size(),std::min(x._degree,y._degree));
    MultiIndex a(x.argument_size());
    X c(0);
    for(ConstIterator xiter=x._expansion.begin(); xiter!=x._expansion.end(); ++xiter) {
        if(xiter->key().degree()>r.degree()) { break; }
        for(ConstIterator yiter=y._expansion.begin(); yiter!=y._expansion.end(); ++yiter) {
            if(xiter->key().degree()+yiter->key().degree()>r.degree()) { break; }
            a=xiter->key()+yiter->key();
            c=static_cast<const X&>(xiter->data())*static_cast<const X&>(yiter->data());
            r._expansion.append(a,c);
        }
    }
    r.cleanup();
    //std::cerr<<"x="<<x<<" y="<<y<<" x*y="<<r<<"\n";
    return r;
}

template<class X>
Differential<X> operator/(const Differential<X>& x, const Differential<X>& y)
{
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<" y="<<y);
    return x*rec(y);
}

template<class X>
Differential<X>
min(const Differential<X>& x1, const Differential<X>& x2)
{
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"min(Differential<X> x1, Differential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()<x2.value() ? x1 : x2;
}


template<class X>
Differential<X>
max(const Differential<X>& x1,const Differential<X>& x2)
{
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"max(Differential<X> x1, Differential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()>x2.value() ? x1 : x2;
}

template<class X>
Differential<X>
abs(const Differential<X>& x)
{
    // FIXME: Maybe need different code for validated and approximate paradigms
    if(decide(x.value()==0)) {
        ARIADNE_THROW(std::runtime_error,"abs(Differential<X> x)","x[0]==0");
    }
    return decide(x.value()>0) ? pos(x) : neg(x);
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
    return -x;
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
Differential<X> pow(const Differential<X>& x, Int n)
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
    ARIADNE_ASSERT_MSG(y.argument_size()==x.size(), "y="<<y<<" x="<<x);
    //std::cerr << "y=" << y << std::endl;
    //std::cerr << "x=" << x << std::endl;
    DegreeType d=y.degree();
    SizeType ms=x.size();
    ARIADNE_ASSERT(d>=1);

    Y zero = x.zero_element();
    Y one = zero; one+=1;

    // Use inefficient brute-force approach with lots of storage...
    Array< Array< Y > > val(ms, Array< Y >(d+1));
    for(SizeType j=0; j!=ms; ++j) {
        val[j][0]=one;
        val[j][1]=x[j];
        for(SizeType k=2; k<=d; ++k) {
            val[j][k]=val[j][k-1]*x[j];
        }
    }

    Y r(zero);
    for(typename Differential<X>::ConstIterator iter=y.begin();
        iter!=y.end(); ++iter)
        {
            const MultiIndex& j=iter->key();
            const X& c=iter->data();
            Y t=one;
            for(SizeType k=0; k!=ms; ++k) {
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
    SizeType as=y.argument_size();
    DegreeType d=std::min(x.degree(),y.degree());

    Differential<X> w=y;
    if(w.begin()->key().degree()==0) { w.begin()->data()=0; }
    Differential<X> r(as,d);
    r[MultiIndex(as)]=x[d];
    for(SizeType n=1; n<=d; ++n) {
        r=r*w;
        r+=x[d-n];
    }
    return r;
}


template<class X>
Vector<X>
gradient(const Differential<X>& x)
{
    Vector<X> r(x.argument_size());
    for(SizeType j=0; j!=x.argument_size(); ++j) {
        r[j]=x.gradient(j);
    }
    return r;
}


template<class X>
Differential<X> derivative(const Differential<X>& x, SizeType i)
{
    if(x.degree()==0) { return Differential<X>(x.argument_size(),0u); }
    Differential<X> r(x.argument_size(), x.degree()-1);
    MultiIndex a(x.argument_size());
    for(typename Differential<X>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter) {
        a=iter->key();
        SizeType n=a[i];
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
Differential<X> antiderivative(const Differential<X>& x, SizeType i)
{
    Differential<X> r(x.argument_size(), x.degree()+1);
    MultiIndex a(x.argument_size());
    MultiIndex ra=MultiIndex(x.argument_size());
    MultiIndex ai=MultiIndex::unit(x.argument_size(),i);
    for(typename Differential<X>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter) {
        a=iter->key();
        const X& xc=x[a];
        ++a[i];
        SizeType n=a[i];
        X& rc=r[a];
        rc=xc/n;
    }
    return r;
}


template<class X>
OutputStream& operator<<(OutputStream& os, const Differential<X>& x)
{
    Expansion<X> e=x.expansion();
    //e.graded_sort();
    os << "SD("<<x.argument_size()<<","<<x.degree()<<"){";
    for(typename Expansion<X>::ConstIterator iter=e.begin(); iter!=e.end(); ++iter) {
        if(iter!=e.begin()) { os << ","; } os << " ";
        for(SizeType i=0; i!=e.argument_size(); ++i) {
            if(i!=0) { os << ","; }
            os << SizeType(iter->key()[i]);
        }
        os << ":" << X(iter->data());
    }
    return os << " }";
}


//! Embed the arguments in a space of dimension \a size, starting at position \a start.
template<class X>
Differential<X>
embed(SizeType before_size, const Differential<X>& x,
      SizeType after_size)
{
    return Differential<X>(x.degree(),x.expansion().embed(before_size,after_size));
}





template<class X> template<class XX>
Vector<Differential<X>>::Vector(SizeType rs, SizeType as, DegreeType d, const XX* ptr)
    : _chars(as,d), _ary(rs,Differential<X>(as,d))
{
    for(SizeType i=0; i!=rs; ++i) {
        for(MultiIndex j(as); j.degree()<=d; ++j) {
            if(*ptr!=0) { (*this)[i][j]=*ptr; } ++ptr; } }
}

template<class X>
Vector<Differential<X>>::Vector(SizeType rs, SizeType as, DegreeType d,
                                  const Vector<X>& v, const Matrix<X>& A)
    :  _chars(as,d), _ary(rs,Differential<X>(as,d))
{
    ARIADNE_ASSERT(rs==v.size());
    ARIADNE_ASSERT(rs==A.row_size());
    ARIADNE_ASSERT(as==A.column_size());
    for(SizeType i=0; i!=this->result_size(); ++i) {
        (*this)[i]=v[i];
        for(SizeType j=0; j!=this->argument_size(); ++j) {
            const X& x=A[i][j];
            if(x!=0) { (*this)[i][j]=x; }
        }
    }
}


/*! \brief Compose the series of \a x and the series of \a y, assuming that \a x is centred at the value of \a y. The value of \a y is therefore unused by this method. */
template<class X, class Y>
Vector<Y>
evaluate(const Vector<Differential<X>>& x,
         const Vector<Y>& y)
{
    Vector<Y> r(x.result_size());
    for(SizeType i=0; i!=r.size(); ++i) {
        r[i]=evaluate(x[i],y);
    }
    return r;
}

/*! \brief Compose the series of \a x and the series of \a y, assuming that \a x is centred at the value of \a y. The value of \a y is therefore unused by this method. */
template<class X>
Differential<X>
compose(const Differential<X>& x,
        const Vector<Differential<X>>& y)
{
    Vector<Differential<X>>& ync=const_cast< Vector<Differential<X>>&>(y);
    Vector<X> yv(y.size());
    for(SizeType i=0; i!=ync.result_size(); ++i) { yv[i]=ync[i].value(); ync[i].set_value(0); }
    Differential<X> r=evaluate(x,ync);
    for(SizeType i=0; i!=ync.result_size(); ++i) { ync[i].set_value(yv[i]); }
    return r;
}


/*! \brief Compose the series of \a x and the series of \a y, assuming that \a x is centred at the value of \a y. The value of \a y is therefore unused by this method. */
template<class X>
Vector<Differential<X>>
compose(const Vector<Differential<X>>& x,
        const Vector<Differential<X>>& y)
{
    ARIADNE_ASSERT(x.degree()==y.degree());
    //std::cerr<<"compose(DV x, DV y)\n x="<<x<<"\n y="<<y<<std::endl;
    Vector<Differential<X>>& ync=const_cast< Vector<Differential<X>>&>(y);
    Vector<X> yv(y.size());
    for(SizeType i=0; i!=ync.result_size(); ++i) { yv[i]=ync[i].value(); ync[i].set_value(0); }
    Vector<Differential<X>> r(x.size(),y.argument_size(),y.degree());
    for(SizeType i=0; i!=x.result_size(); ++i) { r[i]=evaluate(x[i],y); }
    for(SizeType i=0; i!=ync.result_size(); ++i) { ync[i].set_value(yv[i]); }
    return r;
}


template<class X>
Vector<Differential<X>>
lie_derivative(const Vector<Differential<X> >& df, const Vector<Differential<X> >& dg)
{
    Vector<Differential<X>> r(df.result_size(),df.argument_size(),df.degree()-1);
    Differential<X> t(df.argument_size(), df.degree()-1);
    MultiIndex a; X c;
    for(SizeType i=0; i!=df.result_size(); ++i) {
        Expansion<X> const& dfi_expansion = df[i].expansion();
        Expansion<X>& t_expansion = t.expansion();
        for(SizeType j=0; j!=df.argument_size(); ++j) {
            for(typename Expansion<X>::ConstIterator iter=dfi_expansion.begin(); iter!=dfi_expansion.end(); ++iter) {
                if(iter->key()[j]!=0) {
                    a=iter->key();
                    c=iter->data()*SizeType(a[j]);
                    a[j]-=1;
                    t_expansion.append(a,c);
                }
            }
            r[i]+=t*dg[j];
            t.clear();
        }
    }
    return r;
}

template<class X>
Vector<Differential<X>>
antiderivative(const Vector<Differential<X>>& x, SizeType k) {
    Vector<Differential<X>> r(x.size(), Differential<X>(x.argument_size(),x.degree()+1));
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(x[i],k); }
    return r;
}

template<class X>
Vector<X>
value(const Vector<Differential<X>>& x)
{
    Vector<X> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].value();
    }
    return r;
}

template<class X>
Matrix<X>
jacobian(const Vector<Differential<X>>& x)
{
    if(x.size()==0) { return Matrix<X>(); }
    for(SizeType i=1; i!=x.size(); ++i) {
        ARIADNE_ASSERT(x[i].argument_size()==x[0].argument_size());
    }

    Matrix<X> r(x.size(),x[0].argument_size());
    for(SizeType i=0; i!=x.size(); ++i) {
        for(SizeType j=0; j!=x[0].argument_size(); ++j) {
            r[i][j]=x[i].gradient(j);

        }
    }
    return r;
}


} //namespace Ariadne
