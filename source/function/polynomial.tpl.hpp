/***************************************************************************
 *            polynomial.tpl.hpp
 *
 *  Copyright 2008-17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "../algebra/evaluate.tpl.hpp"
#include "../algebra/expansion.inl.hpp"

namespace Ariadne {

template<class X> inline bool is_null(X const& x) { return decide(x==0); }
template<class X> inline bool is_unit(X const& x) { return decide(x==1); }
template<class X> inline bool is_positive(X const& x) { return decide(x>=0); }


template<class X> MultivariatePolynomial<X>::MultivariatePolynomial(SizeType as)
    : _expansion(as)
{
}

template<class X>
MultivariatePolynomial<X>::MultivariatePolynomial(InitializerList<Pair<InitializerList<DegreeType>,X>> lst)
    : _expansion(lst)
{
    this->cleanup();
}


template<class X> MultivariatePolynomial<X> MultivariatePolynomial<X>::create_zero() const {
    return MultivariatePolynomial<X>(this->argument_size());
}


template<class X> MultivariatePolynomial<X> MultivariatePolynomial<X>::constant(SizeType as, const X& c) {
    MultivariatePolynomial<X> r(as); r[MultiIndex::zero(as)]=c; return r;
}

template<class X> MultivariatePolynomial<X> MultivariatePolynomial<X>::variable(SizeType as, SizeType j) {
    ARIADNE_ASSERT(j<as); MultivariatePolynomial<X> r(as); r[MultiIndex::unit(as,j)]=1; return r;
}

template<class X> MultivariatePolynomial<X> MultivariatePolynomial<X>::coordinate(SizeType as, SizeType j) {
    ARIADNE_ASSERT(j<as); MultivariatePolynomial<X> r(as); r[MultiIndex::unit(as,j)]=1; return r;
}


template<class X> Vector<MultivariatePolynomial<X>> MultivariatePolynomial<X>::variables(SizeType as) {
    Vector<MultivariatePolynomial<X>> r(as); for(SizeType i=0; i!=as; ++i) { r[i]=variable(as,i); } return r;
}

template<class X> Vector<MultivariatePolynomial<X>> MultivariatePolynomial<X>::coordinates(SizeType as) {
    Vector<MultivariatePolynomial<X>> r(as); for(SizeType i=0; i!=as; ++i) { r[i]=coordinate(as,i); } return r;
}


template<class X> MultivariatePolynomial<X>& MultivariatePolynomial<X>::operator=(const X& x) {
    this->_expansion.clear();
    this->_expansion.append(MultiIndex(this->argument_size()),x);
    return *this;
}



template<class X> SizeType MultivariatePolynomial<X>::argument_size() const { return this->_expansion.argument_size(); }

template<class X> SizeType MultivariatePolynomial<X>::number_of_terms() const { return this->_expansion.number_of_terms(); }

template<class X> DegreeType MultivariatePolynomial<X>::degree() const {
    DegreeType deg=0u; for(auto iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        deg=std::max(deg,iter->index().degree());
    }
    return deg;
}

template<class X> const X& MultivariatePolynomial<X>::value() const { return this->_expansion[MultiIndex::zero(this->argument_size())]; }

template<class X> X& MultivariatePolynomial<X>::operator[](const MultiIndex& a) { return this->_expansion.at(a); }

template<class X> const X& MultivariatePolynomial<X>::operator[](const MultiIndex& a) const { return this->_expansion.get(a); }

template<class X> const Expansion<MultiIndex,X>& MultivariatePolynomial<X>::expansion() const { return this->_expansion; }

template<class X> Expansion<MultiIndex,X>& MultivariatePolynomial<X>::expansion() { return this->_expansion; }



template<class X> typename MultivariatePolynomial<X>::Iterator MultivariatePolynomial<X>::begin() { return this->_expansion.begin(); }

template<class X> typename MultivariatePolynomial<X>::Iterator MultivariatePolynomial<X>::end() { return this->_expansion.end(); }

template<class X> typename MultivariatePolynomial<X>::Iterator MultivariatePolynomial<X>::find(const MultiIndex& a) { return this->_expansion.find(a); }

template<class X> typename MultivariatePolynomial<X>::ConstIterator MultivariatePolynomial<X>::begin() const { return this->_expansion.begin(); }

template<class X> typename MultivariatePolynomial<X>::ConstIterator MultivariatePolynomial<X>::end() const { return this->_expansion.end(); }

template<class X> typename MultivariatePolynomial<X>::ConstIterator MultivariatePolynomial<X>::find(const MultiIndex& a) const { return this->_expansion.find(a); }



template<class X> Void MultivariatePolynomial<X>::_append(const MultiIndex& a, const X& c) { this->_expansion.append(a,c); }

template<class X> Void MultivariatePolynomial<X>::insert(const MultiIndex& a, const X& c) { this->_expansion.insert(a,c); }

template<class X> Void MultivariatePolynomial<X>::reserve(SizeType n) { this->_expansion.reserve(n); }

template<class X> Void MultivariatePolynomial<X>::erase(Iterator iter) { this->_expansion.erase(iter); }

template<class X> Void MultivariatePolynomial<X>::clear() { this->_expansion.clear(); }



template<class X>
MultivariatePolynomial<X>&
MultivariatePolynomial<X>::differentiate(SizeType j) {
    for(typename MultivariatePolynomial<X>::Iterator iter=this->begin(); iter!=this->end(); ++iter) {
        IndexReference a=iter->index();
        CoefficientReference c=iter->coefficient();
        c*=static_cast<Nat>(a[j]);
        if(a[j]!=0u) { ++a[j]; }
    }
    return *this;
}

template<class X>
MultivariatePolynomial<X>&
MultivariatePolynomial<X>::antidifferentiate(SizeType j) {
    for(typename MultivariatePolynomial<X>::Iterator iter=this->begin(); iter!=this->end(); ++iter) {
        IndexReference a=iter->index();
        CoefficientReference c=iter->coefficient();
        ++a[j];
        c/=static_cast<Nat>(a[j]);
    }
    return *this;
}

template<class X>
MultivariatePolynomial<X>& MultivariatePolynomial<X>::truncate(DegreeType d) {
    MultivariatePolynomial<X> r(this->argument_size());
    for(typename MultivariatePolynomial<X>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->index().degree()<=d && decide(iter->coefficient()!=X(0))) {
            r._append(iter->index(),iter->coefficient());
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
        while(next!=last && curr->index()==next->index()) {
            if(curr->index()==next->index()) {
                curr->coefficient()=op(curr->coefficient(),next->coefficient());
                ++next;
            }
        }
        // Removes zero entries; the code below is preferred to the case "curr->coefficient()!=0" for ValidatedKleenean results
        if(definitely(curr->coefficient()==0)) { }
        else { ++curr; }
    }
    return curr;
}

template<class X>
Void MultivariatePolynomial<X>::cleanup()
{
    MultivariatePolynomial<X>* self=const_cast<MultivariatePolynomial<X>*>(this);
    self->_expansion.index_sort(IndexComparisonType());
    Iterator new_end=unique_key(self->_expansion.begin(), self->_expansion.end(), std::plus<X>());
    self->_expansion.resize(static_cast<SizeType>(new_end-self->_expansion.begin()));
}

template<class X>
Void MultivariatePolynomial<X>::check() const
{
    this->_expansion.check();
}




template<class X> MultivariatePolynomial<X> AlgebraOperations<MultivariatePolynomial<X>>::apply(Pos, const MultivariatePolynomial<X>& p) {
    MultivariatePolynomial<X> r(p.argument_size());
    for(auto iter=p.begin(); iter!=p.end(); ++iter) {
        r[iter->index()]=+iter->coefficient();
    }
    return r;
}

template<class X> MultivariatePolynomial<X> AlgebraOperations<MultivariatePolynomial<X>>::apply(Neg, const MultivariatePolynomial<X>& p) {
    MultivariatePolynomial<X> r(p.argument_size());
    for(auto iter=p.begin(); iter!=p.end(); ++iter) {
        r[iter->index()]=-iter->coefficient();
    }
    return r;
}


template<class X> MultivariatePolynomial<X> AlgebraOperations<MultivariatePolynomial<X>>::apply(Add, MultivariatePolynomial<X> p, const X& c) {
    p[MultiIndex(p.argument_size())]+=c;
    return std::move(p);
}

template<class X> MultivariatePolynomial<X>& AlgebraOperations<MultivariatePolynomial<X>>::iapply(Add, MultivariatePolynomial<X>& p, const X& c) {
    p[MultiIndex(p.argument_size())]+=c;
    return p;
}

template<class X> MultivariatePolynomial<X> AlgebraOperations<MultivariatePolynomial<X>>::apply(Mul, MultivariatePolynomial<X> p, const X& c) {
    if(is_null(c)) {
        p.expansion().clear();
    } else {
        for(auto iter=p.begin(); iter!=p.end(); ++iter) {
            iter->coefficient()*=c;
        }
    }
    return std::move(p);
}

template<class X> MultivariatePolynomial<X>& AlgebraOperations<MultivariatePolynomial<X>>::iapply(Mul, MultivariatePolynomial<X>& p, const X& c) {
    if(is_null(c)) {
        p.expansion().clear();
    } else {
        for(auto iter=p.begin(); iter!=p.end(); ++iter) {
            iter->coefficient()*=c;
        }
    }
    return p;
}

template<class X> MultivariatePolynomial<X> AlgebraOperations<MultivariatePolynomial<X>>::apply(Add, const MultivariatePolynomial<X>& p1, const MultivariatePolynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    typename MultivariatePolynomial<X>::IndexComparisonType less;
    MultivariatePolynomial<X> r(p1.argument_size());
    auto iter1=p1.begin(); auto iter2=p2.begin();
    while (iter1!=p1.end() && iter2!=p2.end()) {
        if (iter1->index()==iter2->index()) {
            r._expansion.append(iter1->index(),iter1->coefficient()+iter2->coefficient());
           ++iter1; ++iter2;
         } else if (less(iter1->index(),iter2->index())) {
            r._expansion.append(iter1->index(),iter1->coefficient());
            ++iter1;
        } else { //  (greater(iter1->index(),iter2->index()))
            r._expansion.append(iter2->index(),iter2->coefficient());
            ++iter2;
        }
    }
    while (iter1!=p1.end()) {
        r._expansion.append(iter1->index(),iter1->coefficient());
        ++iter1;
    }
    while (iter2!=p2.end()) {
            r._expansion.append(iter2->index(),iter2->coefficient());
            ++iter2;
    }
    return std::move(r);
}

template<class X> MultivariatePolynomial<X> AlgebraOperations<MultivariatePolynomial<X>>::apply(Sub, const MultivariatePolynomial<X>& p1, const MultivariatePolynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    typename MultivariatePolynomial<X>::IndexComparisonType less;
    MultivariatePolynomial<X> r(p1.argument_size());
    auto iter1=p1.begin(); auto iter2=p2.begin();
    while (iter1!=p1.end() && iter2!=p2.end()) {
        if (iter1->index()==iter2->index()) {
            r._expansion.append(iter1->index(),iter1->coefficient()-iter2->coefficient());
           ++iter1; ++iter2;
         } else if (less(iter1->index(),iter2->index())) {
            r._expansion.append(iter1->index(),iter1->coefficient());
            ++iter1;
        } else { //  (greater(iter1->index(),iter2->index()))
            r._expansion.append(iter2->index(),-iter2->coefficient());
            ++iter2;
        }
    }
    while (iter1!=p1.end()) {
        r._expansion.append(iter1->index(),iter1->coefficient());
        ++iter1;
    }
    while (iter2!=p2.end()) {
            r._expansion.append(iter2->index(),-iter2->coefficient());
            ++iter2;
    }
    return std::move(r);
}

template<class X> MultivariatePolynomial<X> AlgebraOperations<MultivariatePolynomial<X>>::apply(Mul, const MultivariatePolynomial<X>& p1, const MultivariatePolynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    MultivariatePolynomial<X> r(p1.argument_size());
    for(auto iter1=p1.begin(); iter1!=p1.end(); ++iter1) {
        for(auto iter2=p2.begin(); iter2!=p2.end(); ++iter2) {
            MultiIndex a=iter1->index()+iter2->index();
            r[a]+=iter1->coefficient()*iter2->coefficient();
        }
    }
    return r;
}

template<class X> MultivariatePolynomial<X> AlgebraOperations<MultivariatePolynomial<X>>::apply(Mul, MultivariatePolynomial<X> p, const MultivariateMonomial<X>& m) {
    if(is_null(m.coefficient())) { p.clear(); }
    for(auto iter=p.begin(); iter!=p.end(); ++iter) {
        iter->index()+=m.index();
        iter->coefficient()*=m.coefficient();
    }
    return std::move(p);
}

template<class X> MultivariatePolynomial<X>& AlgebraOperations<MultivariatePolynomial<X>>::iapply(Mul, MultivariatePolynomial<X>& p, const MultivariateMonomial<X>& m) {
    if(is_null(m.coefficient())) { p.clear(); }
    for(auto iter=p.begin(); iter!=p.end(); ++iter) {
        iter->index()+=m.index();
        iter->coefficient()*=m.coefficient();
    }
    return p;
}



template<class X> X MultivariatePolynomial<X>::_evaluate(const MultivariatePolynomial<X>& p, const Vector<X>& x) {
    return horner_evaluate(p._expansion,x);
}

template<class X> MultivariatePolynomial<X> MultivariatePolynomial<X>::_compose(const MultivariatePolynomial<X>& p, const Vector<MultivariatePolynomial<X>>& q) {
    return evaluate(p,q);
}




template<class X>
MultivariatePolynomial<X>
MultivariatePolynomial<X>::_partial_evaluate(const MultivariatePolynomial<X>& x, SizeType k, const X& c)
{
    MultivariatePolynomial<X> r(x.argument_size()-1u);
    MultiIndex ra(r.argument_size());
    if(is_null(c)) {
        for(typename MultivariatePolynomial<X>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            IndexConstReference xa=xiter->index();
            MultiIndex::IndexType xak=xa[k];
            if(xak==0) {
                CoefficientConstReference xv=xiter->coefficient();
                for(Nat i=0; i!=k; ++i) { ra[i]=xa[i]; }
                for(Nat i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1u]; }
                r.expansion().append(ra,xv);
            }
        }
    } else if(is_unit(c)) {
        MultivariatePolynomial<X> s(x.argument_size()-1u);
        Array< MultivariatePolynomial<X> > p(x.degree()+1u,MultivariatePolynomial<X>(x.argument_size()-1u));

        for(typename MultivariatePolynomial<X>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            IndexConstReference xa=xiter->index();
            CoefficientConstReference xv=xiter->coefficient();
            MultiIndex::IndexType xak=xa[k];
            for(Nat i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(Nat i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1u]; }
            assert(ra.degree()+xak==xa.degree());
            p[xak].expansion().append(ra,xv);
        }

        r=p[0];
        for(Nat i=1; i!=p.size(); ++i) {
            r+=p[i];
        }
    } else {
        MultivariatePolynomial<X> s(x.argument_size()-1u);
        Array< MultivariatePolynomial<X> > p(x.degree()+1u,MultivariatePolynomial<X>(x.argument_size()-1u));

        Array<X> cpowers(x.degree()+1u);
        cpowers[0]=static_cast<X>(1); cpowers[1]=c;
        if(x.degree()>=2) { cpowers[2]=sqr(c); }
        for(Nat j=3; j<=x.degree(); ++j) {
            cpowers[j]=cpowers[j-2]*cpowers[2];
        }

        for(typename MultivariatePolynomial<X>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            IndexConstReference xa=xiter->index();
            CoefficientConstReference xv=xiter->coefficient();
            MultiIndex::IndexType xak=xa[k];
            for(Nat i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(Nat i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1u]; }
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
OutputStream& operator<<(OutputStream& os, const MultivariatePolynomial<X>& p) {
    if(p.begin()==p.end()) {
        return os << "{"<<MultiIndex::zero(p.argument_size())<<":0.0}"; }

    os << "{";
    for(typename MultivariatePolynomial<X>::ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
        os << (iter==p.begin() ? "" : ",");
        for(SizeType i=0; i!=iter->index().size(); ++i) {
            os << (i==0?" ":",") << Int(iter->index()[i]); }
        os << ":" << iter->coefficient(); }
    return os << " }";
}
*/

template<class X>
OutputStream& MultivariatePolynomial<X>::_write(OutputStream& os) const {
    List<String> argument_names;
    for(SizeType i=0; i!=this->argument_size(); ++i) {
        StringStream ss;
        ss << "x" << i;
        argument_names.append(ss.str());
    }
    return this->_write(os,argument_names);
}

template<class X>
OutputStream& MultivariatePolynomial<X>::_write(OutputStream& os, List<String> const& argument_names) const {
    const MultivariatePolynomial<X>& p=*this;
    const std::vector<String>& n=argument_names;
    Bool first_term=true;
    Bool identically_zero=true;
    for(typename MultivariatePolynomial<X>::ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
        MultiIndex a=iter->index();
        X v=iter->coefficient();
        if(decide(v!=0)) {
            identically_zero=false;
            Bool first_factor=true;
            if(decide(v>0) && !first_term) { os << "+"; }
            first_term=false;
            if(decide(v==1)) { } else if (decide(v==-1)) { os << '-'; }
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
