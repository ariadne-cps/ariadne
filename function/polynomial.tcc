/***************************************************************************
 *            polynomial.tcc
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

#include "algebra/evaluate.h"

namespace Ariadne {

template<class X> inline auto operator==(X const& x1, Int n2) -> decltype(x1==X(n2)) { return x1==X(n2); }
template<class X> inline auto operator!=(X const& x1, Int n2) -> decltype(x1!=X(n2)) { return x1!=X(n2); }
template<class X> inline auto operator> (X const& x1, Int n2) -> decltype(x1> X(n2)) { return x1> X(n2); }

template<class X> inline bool is_null(X const& x) { return decide(x==0); }
template<class X> inline bool is_unit(X const& x) { return decide(x==1); }
template<class X> inline bool is_positive(X const& x) { return decide(x>=0); }


template<class X> Polynomial<X>::Polynomial(SizeType as)
    : _expansion(as)
{
}

template<class X>
Polynomial<X>::Polynomial(InitializerList<PairType<InitializerList<Int>,X>> lst)
    : _expansion(lst)
{
    this->cleanup();
}

template<class X>
Polynomial<X>::Polynomial(SizeType as, DegreeType deg, InitializerList<X> lst)
    : _expansion(as,deg,lst)
{
    this->cleanup();
}


template<class X> Polynomial<X> Polynomial<X>::create_zero() const {
    return Polynomial<X>(this->argument_size());
}


template<class X> Polynomial<X> Polynomial<X>::constant(SizeType as, const X& c) {
    Polynomial<X> r(as); r[MultiIndex::zero(as)]=c; return r;
}

template<class X> Polynomial<X> Polynomial<X>::variable(SizeType as, SizeType j) {
    ARIADNE_ASSERT(j<as); Polynomial<X> r(as); r[MultiIndex::unit(as,j)]=1; return r;
}

template<class X> Polynomial<X> Polynomial<X>::coordinate(SizeType as, SizeType j) {
    ARIADNE_ASSERT(j<as); Polynomial<X> r(as); r[MultiIndex::unit(as,j)]=1; return r;
}


template<class X> Vector<Polynomial<X>> Polynomial<X>::variables(SizeType as) {
    Vector<Polynomial<X>> r(as); for(SizeType i=0; i!=as; ++i) { r[i]=variable(as,i); } return r;
}


template<class X> Polynomial<X>& Polynomial<X>::operator=(const X& x) {
    this->_expansion.clear();
    this->_expansion.append(MultiIndex(this->argument_size()),x);
    return *this;
}



template<class X> SizeType Polynomial<X>::argument_size() const { return this->_expansion.argument_size(); }

template<class X> SizeType Polynomial<X>::number_of_nonzeros() const { return this->_expansion.number_of_nonzeros(); }

template<class X> SizeType Polynomial<X>::degree() const { return this->_expansion.degree(); }

template<class X> const X& Polynomial<X>::value() const { return this->_expansion[MultiIndex::zero(this->argument_size())]; }

template<class X> X& Polynomial<X>::operator[](const MultiIndex& a) { return this->_expansion.at(a); }

template<class X> const X& Polynomial<X>::operator[](const MultiIndex& a) const { return this->_expansion.get(a); }

template<class X> const Expansion<X>& Polynomial<X>::expansion() const { return this->_expansion; }

template<class X> Expansion<X>& Polynomial<X>::expansion() { return this->_expansion; }



template<class X> typename Polynomial<X>::Iterator Polynomial<X>::begin() { return this->_expansion.begin(); }

template<class X> typename Polynomial<X>::Iterator Polynomial<X>::end() { return this->_expansion.end(); }

template<class X> typename Polynomial<X>::Iterator Polynomial<X>::find(const MultiIndex& a) { return this->_expansion.find(a); }

template<class X> typename Polynomial<X>::ConstIterator Polynomial<X>::begin() const { return this->_expansion.begin(); }

template<class X> typename Polynomial<X>::ConstIterator Polynomial<X>::end() const { return this->_expansion.end(); }

template<class X> typename Polynomial<X>::ConstIterator Polynomial<X>::find(const MultiIndex& a) const { return this->_expansion.find(a); }



template<class X> Void Polynomial<X>::_append(const MultiIndex& a, const X& c) { this->_expansion.append(a,c); }

template<class X> Void Polynomial<X>::insert(const MultiIndex& a, const X& c) { this->_expansion.insert(a,c); }

template<class X> Void Polynomial<X>::reserve(SizeType n) { this->_expansion.reserve(n); }

template<class X> Void Polynomial<X>::erase(Iterator iter) { this->_expansion.erase(iter); }

template<class X> Void Polynomial<X>::clear() { this->_expansion.clear(); }



template<class X>
Polynomial<X>&
Polynomial<X>::differentiate(SizeType j) {
    for(typename Polynomial<X>::Iterator iter=this->begin(); iter!=this->end(); ++iter) {
        MultiIndex& a=iter->key();
        X& c=iter->data();
        c*=static_cast<Nat>(a[j]);
        if(a[j]!=0u) { ++a[j]; }
    }
    return *this;
}

template<class X>
Polynomial<X>&
Polynomial<X>::antidifferentiate(SizeType j) {
    for(typename Polynomial<X>::Iterator iter=this->begin(); iter!=this->end(); ++iter) {
        MultiIndex& a=iter->key();
        X& c=iter->data();
        ++a[j];
        c/=static_cast<Nat>(a[j]);
    }
    return *this;
}

template<class X>
Polynomial<X>& Polynomial<X>::truncate(DegreeType d) {
    Polynomial<X> r(this->argument_size());
    for(typename Polynomial<X>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->key().degree()<=d && iter->data()!=X(0)) {
            r._append(iter->key(),iter->data());
        }
    }
    this->_expansion.swap(r._expansion);
    return *this;
}

template<class FwdIter, class Op>
FwdIter unique_key(FwdIter first, FwdIter last, Op op) {
    FwdIter curr=first;
    FwdIter next=curr;
    while(next!=last) {
        if(curr!=next) { *curr=*next; }
        ++next;
        while(next!=last && curr->key()==next->key()) {
            if(curr->key()==next->key()) {
                curr->data()=op(curr->data(),next->data());
                ++next;
            }
        }
        // Removes zero entries; the code below is preferred to the case "curr->data()!=0" for Tribool results
        if(curr->data()==0) { }
        else { ++curr; }
    }
    return curr;
}

template<class X>
Void Polynomial<X>::cleanup()
{
    Polynomial<X>* self=const_cast<Polynomial<X>*>(this);
    self->_expansion.reverse_lexicographic_sort();
    Iterator new_end=unique_key(self->_expansion.begin(), self->_expansion.end(), std::plus<X>());
    self->_expansion.resize(new_end-self->_expansion.begin());
}

template<class X>
Void Polynomial<X>::check() const
{
    this->_expansion.check();
}




template<class X> Polynomial<X> Polynomial<X>::_neg(const Polynomial<X>& p) {
    Polynomial<X> r(p.argument_size());
    for(auto iter=p.begin(); iter!=p.end(); ++iter) {
        r[iter->key()]-=iter->data();
    }
    return r;
}


template<class X> Polynomial<X> Polynomial<X>::_add(const Polynomial<X>& p, const X& c) {
    Polynomial<X> r(p);
    r[MultiIndex(p.argument_size())]+=c;
    return r;
}

template<class X> Polynomial<X> Polynomial<X>::_mul(const Polynomial<X>& p, const X& c) {
    if(is_null(c)) { return Polynomial<X>(p.argument_size()); }
    Polynomial<X> r(p);
    for(auto iter=r.begin(); iter!=r.end(); ++iter) {
        iter->data()*=c;
    }
    return r;
}

template<class X> Polynomial<X> Polynomial<X>::_add(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    ComparisonType less;
    Polynomial<X> r(p1.argument_size());
    auto iter1=p1.begin(); auto iter2=p2.begin();
    while (iter1!=p1.end() && iter2!=p2.end()) {
        if (iter1->key()==iter2->key()) {
            r._expansion.append(iter1->key(),iter1->data()+iter2->data());
           ++iter1; ++iter2;
         } else if (less(iter1->key(),iter2->key())) {
            r._expansion.append(iter1->key(),iter1->data());
            ++iter1;
        } else { //  (greater(iter1->key(),iter2->key()))
            r._expansion.append(iter2->key(),iter2->data());
            ++iter2;
        }
    }
    while (iter1!=p1.end()) {
        r._expansion.append(iter1->key(),iter1->data());
        ++iter1;
    }
    while (iter2!=p2.end()) {
            r._expansion.append(iter2->key(),iter2->data());
            ++iter2;
    }
    return std::move(r);
}

template<class X> Polynomial<X> Polynomial<X>::_sub(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    ComparisonType less;
    Polynomial<X> r(p1.argument_size());
    auto iter1=p1.begin(); auto iter2=p2.begin();
    while (iter1!=p1.end() && iter2!=p2.end()) {
        if (iter1->key()==iter2->key()) {
            r._expansion.append(iter1->key(),iter1->data()-iter2->data());
           ++iter1; ++iter2;
         } else if (less(iter1->key(),iter2->key())) {
            r._expansion.append(iter1->key(),iter1->data());
            ++iter1;
        } else { //  (greater(iter1->key(),iter2->key()))
            r._expansion.append(iter2->key(),-iter2->data());
            ++iter2;
        }
    }
    while (iter1!=p1.end()) {
        r._expansion.append(iter1->key(),iter1->data());
        ++iter1;
    }
    while (iter2!=p2.end()) {
            r._expansion.append(iter2->key(),-iter2->data());
            ++iter2;
    }
    return std::move(r);
}

template<class X> Polynomial<X> Polynomial<X>::_mul(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    Polynomial<X> r(p1.argument_size());
    for(auto iter1=p1.begin(); iter1!=p1.end(); ++iter1) {
        for(auto iter2=p2.begin(); iter2!=p2.end(); ++iter2) {
            MultiIndex a=iter1->key()+iter2->key();
            r[a]+=iter1->data()*iter2->data();
        }
    }
    return r;
}

template<class X> Polynomial<X>& Polynomial<X>::_imul(Polynomial<X>& p, const Monomial<X>& m) {
    if(is_null(m.data())) { p.clear(); }
    for(auto iter=p.begin(); iter!=p.end(); ++iter) {
        iter->key()+=m.key();
        iter->data()*=m.data();
    }
    return p;
}



template<class X> X Polynomial<X>::_evaluate(const Polynomial<X>& p, const Vector<X>& x) {
    return horner_evaluate(p._expansion,x);
}

template<class X> Polynomial<X> Polynomial<X>::_compose(const Polynomial<X>& p, const Vector<Polynomial<X>>& q) {
    return evaluate(p,q);
}




template<class X>
Polynomial<X>
Polynomial<X>::_partial_evaluate(const Polynomial<X>& x, SizeType k, const X& c)
{
    Polynomial<X> r(x.argument_size()-1);
    MultiIndex ra(r.argument_size());
    if(is_null(c)) {
        for(typename Polynomial<X>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            MultiIndex::IndexType xak=xa[k];
            if(xak==0) {
                const X& xv=xiter->data();
                for(Nat i=0; i!=k; ++i) { ra[i]=xa[i]; }
                for(Nat i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
                r.expansion().append(ra,xv);
            }
        }
    } else if(is_unit(c)) {
        Polynomial<X> s(x.argument_size()-1);
        Array< Polynomial<X> > p(x.degree()+1,Polynomial<X>(x.argument_size()-1));

        for(typename Polynomial<X>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            const X& xv=xiter->data();
            MultiIndex::IndexType xak=xa[k];
            for(Nat i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(Nat i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
            assert(ra.degree()+xak==xa.degree());
            p[xak].expansion().append(ra,xv);
        }

        r=p[0];
        for(Nat i=1; i!=p.size(); ++i) {
            r+=p[i];
        }
    } else {
        Polynomial<X> s(x.argument_size()-1);
        Array< Polynomial<X> > p(x.degree()+1,Polynomial<X>(x.argument_size()-1));

        Array<X> cpowers(x.degree()+1);
        cpowers[0]=static_cast<X>(1); cpowers[1]=c;
        if(x.degree()>=2) { cpowers[2]=sqr(c); }
        for(Nat j=3; j<=x.degree(); ++j) {
            cpowers[j]=cpowers[j-2]*cpowers[2];
        }

        for(typename Polynomial<X>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            const X& xv=xiter->data();
            MultiIndex::IndexType xak=xa[k];
            for(Nat i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(Nat i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
            assert(ra.degree()+xak==xa.degree());
            p[xak].expansion().append(ra,xv);
        }
        for(Nat i=1; i!=p.size(); ++i) {
            p[i]*=cpowers[i];
        }

        r=p[0];
        for(Nat i=1; i!=p.size(); ++i) {
            r+=p[i];
        }
    }

    return r;
}



/*
template<class X>
OutputStream& operator<<(OutputStream& os, const Polynomial<X>& p) {
    if(p.begin()==p.end()) {
        return os << "{"<<MultiIndex::zero(p.argument_size())<<":0.0}"; }

    os << "{";
    for(typename Polynomial<X>::ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
        os << (iter==p.begin() ? "" : ",");
        for(SizeType i=0; i!=iter->key().size(); ++i) {
            os << (i==0?" ":",") << Int(iter->key()[i]); }
        os << ":" << iter->data(); }
    return os << " }";
}
*/

template<class X>
OutputStream& Polynomial<X>::_write(OutputStream& os) const {
    List<String> argument_names;
    for(SizeType i=0; i!=this->argument_size(); ++i) {
        StringStream ss;
        ss << "x" << i;
        argument_names.append(ss.str());
    }
    return this->_write(os,argument_names);
}

template<class X>
OutputStream& Polynomial<X>::_write(OutputStream& os, List<String> const& argument_names) const {
    const Polynomial<X>& p=*this;
    const std::vector<String>& n=argument_names;
    Bool first_term=true;
    Bool identically_zero=true;
    for(typename Polynomial<X>::ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
        MultiIndex a=iter->key();
        X v=iter->data();
        if(decide(v!=0)) {
            identically_zero=false;
            Bool first_factor=true;
            if(decide(v>0) && !first_term) { os << "+"; }
            first_term=false;
            if(v==1) { } else if (v==-1) { os << '-'; }
            else { os << v; first_factor=false; }
            for(Nat j=0; j!=a.size(); ++j) {
                if(a[j]!=0) {
                    if(first_factor) { first_factor=false; } else { os <<"*"; }
                    os<<n[j]; if(a[j]!=1) { os<<"^"<<Int(a[j]); } }
            }
            if(first_factor) { os << v; }
        }
    }
    if(identically_zero) { os << "0"; }
    return os;
}




} // namespace Ariadne
