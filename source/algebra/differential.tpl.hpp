/***************************************************************************
 *            differential.tpl.hpp
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

#include <map>

#include "utility/macros.hpp"
#include "utility/array.hpp"
#include "numeric/float.decl.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/series.hpp"
#include "algebra/expansion.hpp"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Series;

template<class X> class Expansion;
template<class X> class Differential;




template<class X> UnivariateDifferential<X>::UnivariateDifferential()
    : _ary(1u,X(0)) { }

template<class X> UnivariateDifferential<X>::UnivariateDifferential(DegreeType d)
    : _ary(d+1,X(0)) { }

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

template<class X> OutputStream& UnivariateDifferential<X>::write(OutputStream& os) const {
    os << this->_ary[0];
    for(DegreeType i=1; i<=this->degree(); ++i) {
        os <<" ";
        if(decide(this->_ary[i]>=X(0))) { os << "+"; }
        os << this->_ary[i] << "*dx^" << (uint)i;
    }
    return os;
}




//template<class X> Differential<X>::Differential() : _expansion(0), _degree(0) { }

template<class X> Differential<X>::Differential(SizeType as, DegreeType deg, X const& z) : _expansion(as,z), _degree(deg) { }

template<class X> Differential<X>::Differential(const Map<MultiIndex,X>& map, DegreeType deg) : _expansion(0u) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class X>
Differential<X>::Differential(SizeType as, DegreeType deg,
                              InitializerList< Pair<InitializerList<DegreeType>,X> > lst)
    : _expansion(Expansion<X>(lst)), _degree(deg)
{
    this->cleanup();
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
    X one(v); one=1;
    r._expansion.append(a,one);
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
    X zero=x.zero_element();
    Vector<Differential<X>> result(x.size(),Differential<X>(x.size(),deg,zero));
    X& one=zero; one=1;
    MultiIndex a(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        result[i]._expansion.append(a,x[i]);
        a[i]=1; result[i]._expansion.append(a,one); a[i]=0;
    }
    return result;
}

//template<class X> UnivariateDifferential<X> Differential<X>::identity(DegreeType deg, const X& x) {
//    return UnivariateDifferential<X>::variable(deg,x);
//}

template<class X> Differential<X> Differential<X>::identity(DegreeType deg, const X& x) {
    return Differential<X>::variable(1u,deg,x,0u);
}

template<class X> Vector<Differential<X>> Differential<X>::identity(DegreeType deg, const Vector<X>& x) {
    return Differential<X>::variables(deg,x);
}

template<class X> EqualityType<X> Differential<X>::operator==(const Differential<X>& other) const {
    Differential<X> const& self=*this;
    if(self.argument_size()!=other.argument_size()) { return false; }
    EqualityType<X> result=true;
    IndexComparisonType less;
    ConstIterator self_iter=self.begin(); ConstIterator other_iter=other.begin();
    while(self_iter!=self.end() && other_iter!=other.end()) {
        if(self_iter->key()==other_iter->key()) {
            result = result && (self_iter->data()==other_iter->data());
            ++self_iter; ++other_iter;
        } else if (less(self_iter->key(), other_iter->key())) {
            result = result && (self_iter->data()==0);
            ++self_iter;
        } else { // (self_iter->key() > other_iter->key())
            result = result && (other_iter->data()==0);
            ++other_iter;
        }
    }
    while(self_iter!=self.end()) {
        result = result && (self_iter->data()==0);
        ++self_iter;
    }
    while(other_iter!=other.end()) {
        result = result && (other_iter->data()==0);
        ++other_iter;
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

template<class X> X Differential<X>::zero_coefficient() const {
    return this->_expansion.zero_coefficient();
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

template<class X> Covector<X> Differential<X>::gradient() const {
    Covector<X> g(this->argument_size());
    for(SizeType j=0; j!=g.size(); ++j) { g[j]=this->gradient(j); }
    return g;
}

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
    if(iter==this->_expansion.end()) { return this->_expansion.zero_coefficient(); }
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

template<class X> Differential<X> Differential<X>::create_constant(NumericType const& c) const {
    return Differential<X>::constant(this->argument_size(),this->degree(),c);
}

template<class X> Differential<X> Differential<X>::create_constant(Int c) const {
    X xc=this->value(); xc=c; return Differential<X>::constant(this->argument_size(),this->degree(),xc);
}


template<class X> Void Differential<X>::clear() {
    this->_expansion.clear();
}

template<class X> Void Differential<X>::cleanup() {
    this->_expansion.sort();
    this->_expansion.combine_terms();
    this->_expansion.remove_zeros();
}

template<class X> Void Differential<X>::check() const {
    for(auto iter=this->begin(); iter!=this->end(); ++iter) {
        ARIADNE_ASSERT_MSG(iter->key().degree()<=this->degree(), *this);
        auto next = iter; ++next;
        ARIADNE_ASSERT_MSG(graded_less(iter->key(),next->key()),"ErrorTag in ordering Differential "<<this->expansion());
    }
}


template<class X>
Differential<X> AlgebraOperations<Differential<X>>::_pos(Differential<X> x)
{
    for(auto iter=x.begin(); iter!=x.end(); ++iter) {
        X& xa=iter->data();
        xa=+xa;
    }
    return std::move(x);
}

template<class X>
Differential<X> AlgebraOperations<Differential<X>>::_neg(Differential<X> x)
{
    for(auto iter=x.begin(); iter!=x.end(); ++iter) {
        X& xa=iter->data();
        xa=-xa;
    }
    return std::move(x);
}

template<class X>
Differential<X> AlgebraOperations<Differential<X>>::_add(Differential<X> x, const X& c)
{
    MultiIndex a(x.argument_size());
    if(x.expansion().empty()) {
        x.expansion().append(a,c);
    } else if(x.begin()->key()!=a) {
        x.expansion().prepend(a,c);
    } else {
        x.begin()->data()+=c;
    }
    return std::move(x);
}

template<class X>
Differential<X> AlgebraOperations<Differential<X>>::_mul(Differential<X> x, const X& c)
{
    if(decide(c==static_cast<X>(0))) {
        x.clear();
    } else {
        for(auto iter=x.begin(); iter!=x.end(); ++iter) {
            static_cast<X&>(iter->data())*=c;
        }
    }
    return std::move(x);
}



template<class X>
Differential<X> AlgebraOperations<Differential<X>>::_add(const Differential<X>& x, const Differential<X>& y)
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
Differential<X> AlgebraOperations<Differential<X>>::_sub(const Differential<X>& x, const Differential<X>& y)
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
Differential<X> AlgebraOperations<Differential<X>>::_mul(const Differential<X>& x, const Differential<X>& y)
{
    typedef typename Differential<X>::ConstIterator ConstIterator;
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<" y="<<y);
    Differential<X> r(x.argument_size(),std::min(x.degree(),y.degree()));
    MultiIndex a(x.argument_size());
    X c(0);
    for(ConstIterator xiter=x.expansion().begin(); xiter!=x.expansion().end(); ++xiter) {
        if(xiter->key().degree()>r.degree()) { break; }
        for(ConstIterator yiter=y.expansion().begin(); yiter!=y.expansion().end(); ++yiter) {
            if(xiter->key().degree()+yiter->key().degree()>r.degree()) { break; }
            a=xiter->key()+yiter->key();
            c=static_cast<const X&>(xiter->data())*static_cast<const X&>(yiter->data());
            r.expansion().append(a,c);
        }
    }
    r.cleanup();
    //std::cerr<<"x="<<x<<" y="<<y<<" x*y="<<r<<"\n";
    return r;
}

template<class X>
Differential<X> AlgebraOperations<Differential<X>>::_div(const Differential<X>& x, const Differential<X>& y)
{
    return x * rec(y);
}

template<class X>
Differential<X> AlgebraOperations<Differential<X>>::_rec(const Differential<X>& x)
{
    return GradedAlgebraOperations<Differential<X>>::_rec(x);
}


template<class X> Differential<X> AlgebraOperations<Differential<X>>::_min(const Differential<X>& x1, const Differential<X>& x2) {
    // FIXME: Maybe need different code for validated and approximate paradigms
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(decide(x1.value()==x2.value())) {
        ARIADNE_THROW(std::runtime_error,"min(Differential<X> x1, Differential<X> x2)","x1[0]==x2[0]");
    }
    return decide(x1.value()<x2.value()) ? x1 : x2;
}


template<class X> Differential<X> AlgebraOperations<Differential<X>>::_max(const Differential<X>& x1,const Differential<X>& x2) {
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(decide(x1.value()==x2.value())) {
        ARIADNE_THROW(std::runtime_error,"max(Differential<X> x1, Differential<X> x2)","x1[0]==x2[0]");
    }
    return decide(x1.value()>x2.value()) ? x1 : x2;
}

template<class X> Differential<X> AlgebraOperations<Differential<X>>::_abs(const Differential<X>& x) {
    // FIXME: Maybe need different code for validated and approximate paradigms
    if(decide(x.value()==X(0))) {
        ARIADNE_THROW(std::runtime_error,"abs(Differential<X> x)","x[0]==0");
    }
    return decide(x.value()>X(0)) ? pos(x) : neg(x);
}




template<class X> Differential<X> _evaluate(const Differential<X>& x, const Vector<Differential<X>>& a)
{
    ARIADNE_ASSERT_MSG(x.argument_size()==a.size(), "x="<<x<<" a="<<a);
    DegreeType d=x.degree();
    SizeType ms=a.size();
    ARIADNE_ASSERT(d>=1);

    Differential<X> zero = a.zero_element();
    Differential<X> one = zero.create_constant(X(1));

    // Use inefficient brute-force approach with lots of storage...
    Array< Array< Differential<X> > > val(ms, Array< Differential<X> >(d+1,zero));
    for(SizeType j=0; j!=ms; ++j) {
        val[j][0]=one;
        val[j][1]=a[j];
        for(SizeType k=2; k<=d; ++k) {
            val[j][k]=val[j][k-1]*a[j];
        }
    }

    Differential<X> r(zero);
    for(auto iter=x.begin(); iter!=x.end(); ++iter)
    {
        const MultiIndex& j=iter->key();
        const X& c=iter->data();
        Differential<X> t=one;
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
Differential<X> Differential<X>::_compose(const UnivariateDifferential<X>& x, const Differential<X>& y)
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
Differential<X> Differential<X>::_compose(const Differential<X>& x, const Vector<Differential<X>>& y)
{
    Vector<Differential<X>>& ync=const_cast< Vector<Differential<X>>&>(y);
    Vector<X> yv(y.size());
    X zero(0);
    for(SizeType i=0; i!=ync.result_size(); ++i) { yv[i]=ync[i].value(); ync[i].set_value(zero); }
    Differential<X> r=_evaluate(x,ync);
    for(SizeType i=0; i!=ync.result_size(); ++i) { ync[i].set_value(yv[i]); }
    return r;
}


template<class X>
Vector<Differential<X>> Differential<X>::_compose(const Vector<Differential<X>>& x, const Vector<Differential<X>>& y)
{
    ARIADNE_ASSERT(x.degree()==y.degree());
    //std::cerr<<"compose(DV x, DV y)\n x="<<x<<"\n y="<<y<<std::endl;
    Vector<Differential<X>>& ync=const_cast< Vector<Differential<X>>&>(y);
    Vector<X> yv(y.size());
    X zero(0);
    for(SizeType i=0; i!=ync.result_size(); ++i) { yv[i]=ync[i].value(); ync[i].set_value(zero); }
    Vector<Differential<X>> r(x.size(),y.argument_size(),y.degree());
    for(SizeType i=0; i!=x.result_size(); ++i) { r[i]=_evaluate(x[i],y); }
    for(SizeType i=0; i!=ync.result_size(); ++i) { ync[i].set_value(yv[i]); }
    return r;
}


template<class X>
Differential<X> Differential<X>::_derivative(const Differential<X>& x, SizeType i)
{
    if(x.degree()==0) { return Differential<X>(x.argument_size(),0u); }
    Differential<X> r(x.argument_size(), x.degree()-1);
    MultiIndex a(x.argument_size());
    for(typename Differential<X>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter) {
        a=iter->key();
        Nat n=a[i];
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
Differential<X> Differential<X>::_antiderivative(const Differential<X>& x, SizeType i)
{
    Differential<X> r(x.argument_size(), x.degree()+1);
    MultiIndex a(x.argument_size());
    MultiIndex ra=MultiIndex(x.argument_size());
    MultiIndex ai=MultiIndex::unit(x.argument_size(),i);
    for(typename Differential<X>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter) {
        a=iter->key();
        const X& xc=x[a];
        ++a[i];
        Nat n=a[i];
        X& rc=r[a];
        rc=xc/n;
    }
    return r;
}


template<class X>
OutputStream& Differential<X>::_write(OutputStream& os) const
{
    Differential<X> const& x=*this;
    Expansion<X> e=x.expansion();
    //e.graded_sort();
    os << "SD("<<x.argument_size()<<","<<(uint)x.degree()<<"){";
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





template<class X> Vector<Differential<X>>::Vector(SizeType rs, SizeType as, DegreeType d, X const& z)
    : _chars(as,d), _ary(rs,Differential<X>(as,d,z)) {
}

template<class X> Vector<Differential<X>>::Vector(SizeType rs, const Differential<X>& sd)
    : _chars(sd), _ary(rs,sd) {
}

template<class X> Vector<Differential<X>>::Vector(SizeType rs, const Differential<X>* p)
    : _chars(), _ary(p,p+rs)
{
    ARIADNE_ASSERT(rs>0); _chars=DifferentialCharacteristics<X>(p[0]);
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
            if(decide(x!=0)) { (*this)[i][j]=x; }
        }
    }
}


template<class X> Vector<X> Vector<Differential<X>>::value() const {
    Vector<Differential<X>> const& x=*this;
    Vector<X> r(x.result_size(),x.zero_coefficient());
    for(SizeType i=0; i!=r.size(); ++i) {
        r[i]=x[i].value();
    }
    return r;
}

template<class X> Matrix<X> Vector<Differential<X>>::jacobian() const {
    Vector<Differential<X>> const& x=*this;
    Matrix<X> r(x.size(),x.argument_size(),x.zero_coefficient());
    for(SizeType i=0; i!=x.size(); ++i) {
        for(SizeType j=0; j!=x[0].argument_size(); ++j) {
            r[i][j]=x[i].gradient(j);

        }
    }
    return r;
}

template<class X> Void Vector<Differential<X>>::set_value(const Vector<X>& c) {
    ARIADNE_ASSERT(this->result_size()==c.size());
    for(SizeType i=0; i!=c.size(); ++i) {
        (*this)[i].set_value(c[i]);
    }
}

template<class X> Vector<Differential<X>> Vector<Differential<X>>::constant(SizeType rs, SizeType as, DegreeType d, const Vector<X>& c) {
    ARIADNE_ASSERT(c.size()==rs);
    Vector< Differential<X> > result(rs,as,d);
    for(SizeType i=0; i!=rs; ++i) { result[i]=c[i]; }
    return result;
}

template<class X> Vector<Differential<X>> Vector<Differential<X>>::variable(SizeType rs, SizeType as, DegreeType d, const Vector<X>& x) {
    ARIADNE_ASSERT(x.size()==rs);
    Vector< Differential<X> > result(rs,as,d);
    for(SizeType i=0; i!=rs; ++i) { result[i]=x[i]; result[i][i]=X(1.0); }
    return result;
}

template<class X> Vector<Differential<X>> Vector<Differential<X>>::affine(SizeType rs, SizeType as, DegreeType d, const Vector<X>& b, const Matrix<X>& A) {
    ARIADNE_ASSERT(b.size()==rs);
    ARIADNE_ASSERT(A.row_size()==rs);
    ARIADNE_ASSERT(A.column_size()==as);
    Vector< Differential<X> > result(rs,as,d);
    for(SizeType i=0; i!=rs; ++i) {
        result[i]=b[i];
        for(SizeType j=0; j!=as; ++j) {
            result[i][j]=A[i][j];
        }
    }
    return result;
}



template<class X> Vector<Differential<X>> Vector<Differential<X>>::_compose(Vector<Differential<X>> const& x, Vector<Differential<X>> const& y) {
    Vector<Differential<X>> r(x.size(),y.zero_element());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=compose(x[i],y); }
    return r;
}

DegreeType min(DegreeType d1, DegreeType d2) { return std::min(d1,d2); }

template<class X>
Vector<Differential<X>>
Vector<Differential<X>>::_derivative(const Vector<Differential<X>>& x, SizeType k) {
    Vector<Differential<X>> r(x.size(), Differential<X>(x.argument_size(),min(x.degree(),(DegreeType)1u)-1u));
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=derivative(x[i],k); }
    return r;
}

template<class X>
Vector<Differential<X>>
Vector<Differential<X>>::_antiderivative(const Vector<Differential<X>>& x, SizeType k) {
    Vector<Differential<X>> r(x.size(), Differential<X>(x.argument_size(),x.degree()+1));
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(x[i],k); }
    return r;
}




template<class X>
Vector<Differential<X>>
Vector<Differential<X>>::_lie_derivative(const Vector<Differential<X> >& df, const Vector<Differential<X> >& dg)
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




} //namespace Ariadne
