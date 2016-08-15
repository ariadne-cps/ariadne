/***************************************************************************
 *            algebra/expansion.tcc
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

/*! \file algebra/expansion.tcc
 *  \brief
 */

#include <algorithm>
#include <exception>
#include <stdexcept>

#include "expansion.h"
#include "numeric/logical.h"

namespace Ariadne {

inline SizeType word_size(SizeType as) { return (1u+as)/sizeof(int)+1; }


template<class X> Expansion<X>::~Expansion()
{
    delete[] _indices;
    _indices=nullptr;
    delete[] _coefficients;
    _coefficients=nullptr;
}

template<class X> Expansion<X>::Expansion(SizeType as)
    : Expansion<X>(as,X()) {
}

template<class X> Expansion<X>::Expansion(SizeType as, PrecisionType pr, SizeType cap)
    : Expansion<X>(as,X(pr),cap) {
}

template<class X> Expansion<X>::Expansion(SizeType as, X const& z, SizeType cap)
    : _zero_coefficient(z), _capacity(cap), _size(0u), _argument_size(as)
    , _indices(new DegreeType[_capacity*(_argument_size+1)]), _coefficients(new CoefficientType[_capacity])
{ }

template<class X> Expansion<X>::Expansion(InitializerList<Pair<InitializerList<DegreeType>,X>> lst)
    : Expansion(lst.size()==0 ? 0u : lst.begin()->first.size(), nul(lst.begin()->second))
{
    MultiIndex a;
    X x;
    for(auto iter=lst.begin();
        iter!=lst.end(); ++iter)
    {
        a=iter->first;
        x=iter->second;
        if(decide(x!=0)) { this->append(a,x); }
    }
}



/*
// Call std::memcpy(dest,src,count) with same arguments as std::copy(src_begin,src_end,dest)
namespace std {
    template<> inline void copy(const unsigned char* b, const unsigned char* e, unsigned char* t) { memcpy(t,b,e-b); }
}
*/

template<class X> Expansion<X>::Expansion(const Expansion<X>& e)
    : _capacity(e.capacity()), _size(e.size()), _argument_size(e.argument_size())
    , _indices(new DegreeType[_capacity*(_argument_size+1)]), _coefficients(new CoefficientType[_capacity])
{
    std::copy(e._indices,e._indices+_size*(_argument_size+1),_indices);
    std::copy(e._coefficients,e._coefficients+_size,_coefficients);
}

template<class X> Expansion<X>& Expansion<X>::operator=(const Expansion<X>& e)
{
    if(*this!=e) {
        // Perform memory reallocation if necessary
        if(this->_capacity<e._size) {
            SizeType new_capacity=e._capacity;
            while(new_capacity/2>=e._size) { new_capacity/=2; }
            DegreeType* new_indices=new DegreeType[new_capacity*(e._argument_size+1)];
            CoefficientType* new_coefficients=new CoefficientType[new_capacity];
            _capacity=new_capacity;
            delete[] _indices;
            _indices=new_indices;
            delete[] _coefficients;
            _coefficients=new_coefficients;
        } else if(this->_argument_size!=e._argument_size) {
            DegreeType* new_indices=new DegreeType[this->_capacity*(e._argument_size+1)];
            delete[] _indices;
            _indices=new_indices;
        }
        _size=e._size;
        _argument_size=e._argument_size;
        std::copy(e._indices,e._indices+_size*(_argument_size+1),_indices);
        std::copy(e._coefficients,e._coefficients+_size,_coefficients);
    }
    return *this;
}

template<class X> Expansion<X>::Expansion(Expansion<X>&& e)
    : _capacity(e.capacity()), _size(e.size()), _argument_size(e.argument_size())
    , _indices(e._indices), _coefficients(e._coefficients)
{
    e._indices=nullptr;
    e._coefficients=nullptr;
}

template<class X> Expansion<X>& Expansion<X>::operator=(Expansion<X>&& e)
{
    if(this!=&e) {
        delete[] _indices;
        delete[] _coefficients;
        _capacity=e._capacity;
        _size=e._size;
        _argument_size=e._argument_size;
        _indices=e._indices;
        _coefficients=e._coefficients;
        e._indices=nullptr;
        e._coefficients=nullptr;
    }
    return *this;
}

template<class X> Void Expansion<X>::swap(Expansion<X>& other) {
    std::swap(this->_coefficients,other._coefficients);
    std::swap(this->_indices,other._indices);
    std::swap(this->_size,other._size);
    std::swap(this->_capacity,other._capacity);
    std::swap(this->_argument_size,other._argument_size);
}

template<class X> SizeType Expansion<X>::number_of_terms() const {
    return this->size();
}

template<class X> SizeType Expansion<X>::number_of_nonzeros() const {
    return this->size();
}

template<class X> Bool Expansion<X>::empty() const {
    return this->_size==0;
}

template<class X> SizeType Expansion<X>::size() const {
    return this->_size;
}

template<class X> SizeType Expansion<X>::argument_size() const {
    return this->_argument_size;
}

template<class X> X const& Expansion<X>::zero_coefficient() const {
    return this->_zero_coefficient;
}

template<class X> Void Expansion<X>::reserve(SizeType new_capacity) {
    if(this->_capacity < new_capacity) {
        DegreeType* new_indices = new DegreeType[new_capacity*(this->_argument_size+1)];
        CoefficientType* new_coefficients = new CoefficientType[new_capacity];
        std::copy(this->_indices, this->_indices+this->_size*(this->_argument_size+1), new_indices);
        std::copy(this->_coefficients, this->_coefficients+this->_size, new_coefficients);
        this->_capacity=new_capacity;
        delete[] this->_indices;
        this->_indices=new_indices;
        delete[] this->_coefficients;
        this->_coefficients=new_coefficients;
    }
}

template<class X> Void Expansion<X>::resize(SizeType new_size) {
    if(new_size<this->size()) {
        this->_size=new_size;
    } else {
        if(this->_capacity < new_size) {
            this->reserve(new_size);
        }
        MultiIndex a(this->argument_size());
        X c=_zero_coefficient;
        for (SizeType i=this->_size; i!=new_size; ++i) {
            this->append(a,c);
        }
    }
}

template<class X> SizeType Expansion<X>::capacity() const {
    return this->_capacity;
}

template<class X> Void Expansion<X>::clear() {
    this->_size=0;
}

template<class X> Void Expansion<X>::remove_zeros() {
    this->resize(std::remove_if(this->begin(),this->end(),CoefficientIsZero())-this->begin());
}

template<class X, EnableIf<IsSame<SumType<X>,X>> =dummy> Void combine_terms(Expansion<X>& e) {
    auto begin=e.begin();
    auto end=e.end();
    auto curr=begin;
    auto adv=begin;
    while (adv!=end) {
        curr->index()=adv->index();
        curr->coefficient()=adv->coefficient();
        ++adv;
        while (adv!=end && adv->index()==curr->index()) {
            curr->coefficient() += adv->coefficient();
            ++adv;
        }
        ++curr;
    }
    e.resize(curr-begin);
}

template<class X, DisableIf<IsSame<SumType<X>,X>> =dummy> Void combine_terms(Expansion<X>& e) {
    ARIADNE_ASSERT_MSG(false, "Cannot combine terms of an expansion if the coefficients do not support inplace addition.");
}

template<class X> Void Expansion<X>::combine_terms() {
    Ariadne::combine_terms(*this);
}

template<class X> Void Expansion<X>::check() const {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class X> ExpansionValueReference<X> Expansion<X>::operator[](const MultiIndex& a) {
    return ExpansionValueReference<X>(*this,a);
}

template<class X> const X& Expansion<X>::operator[](const MultiIndex& a) const {
    return this->get(a);
}

template<class X> X& Expansion<X>::at(const MultiIndex& a) {
    auto iter=this->find(a);
    if(iter==this->end()) { this->append(a,this->_zero_coefficient); iter=this->end()-1; }
    return iter->coefficient();
}

template<class X> Void Expansion<X>::set(const MultiIndex& a, const X& c) {
    auto iter=this->find(a);
    if(iter==this->end()) { this->append(a,c); }
    else { iter->coefficient() = c; }
}

template<class X> const X& Expansion<X>::get(const MultiIndex& a) const {
    auto iter=this->find(a);
    if(iter==this->end()) { return this->_zero_coefficient; }
    else { return iter->coefficient(); }
}

template<class X> Bool Expansion<X>::operator==(const Expansion<X>& other) const {
    auto iter1=this->begin();
    auto iter2=other.begin();
    auto end1=this->end();
    auto end2=other.end();

    while(true) {
        if(iter1!=end1 && iter2!=end2) {
            if(!decide(*iter1 == *iter2)) { return false; }
            ++iter1; ++iter2;
        } else if (iter1==end1 && iter2==end2) {
            return true;
        } else {
            return false;
        }
    }
}

template<class X> Bool Expansion<X>::operator!=(const Expansion<X>& other) const {
    return !(*this==other);
}

template<class X> auto Expansion<X>::insert(Iterator pos, const MultiIndex& a, const X& c) -> Iterator {
    if(this->size()==this->capacity()) {
        SizeType where=pos-this->begin();
        this->append(a,c);
        pos=this->begin()+where;
    } else {
        this->append(a,c);
    }
    auto curr=this->end();
    auto prev=curr-1;
    while(prev!=pos) {
        --curr;
        --prev;
        *curr=*prev;
    }
    pos->index()=a;
    pos->coefficient()=c;
    return pos;
}

template<class X> auto Expansion<X>::erase(Iterator pos) -> Iterator {
    auto curr=pos;
    auto next=curr;
    ++next;
    while(next!=this->end()) {
        *curr=*next;
        ++curr;
        ++next;
    }
    --this->_size;
    return pos;
}

template<class X> Void Expansion<X>::prepend(const MultiIndex& a, const X& c) {
    this->append(a,c);
    auto curr=this->end();
    auto prev=curr-1;
    while(prev!=this->begin()) {
        --curr;
        --prev;
        *curr=*prev;
    }
    --curr;
    curr->index()=a;
    curr->coefficient()=c;
}

template<class X> Void Expansion<X>::append(const MultiIndex& a, const X& c) {
    if(this->number_of_terms()==this->capacity()) {
        this->reserve(2*this->capacity());
    }
    DegreeType* _p=this->_indices+(this->_size)*(this->_argument_size+1);
    for(SizeType j=0; j!=this->_argument_size; ++j) { _p[j]=a[j]; }
    _p[this->_argument_size]=a.degree();
    _coefficients[this->_size]=c;
    ++this->_size;
}

template<class X> Void Expansion<X>::append_sum(const MultiIndex& a1, const MultiIndex& a2, const X& c) {
    if(this->number_of_terms()==this->capacity()) {
        this->reserve(2*this->capacity());
    }
    DegreeType* _p=this->_indices+(this->_size)*(this->_argument_size+1);
    for(SizeType j=0; j<=this->_argument_size; ++j) { _p[j]=a1[j]+a2[j]; }
    _coefficients[this->_size]=c;
    ++this->_size;
}

template<class X> ExpansionIterator<X> Expansion<X>::begin() {
    return ExpansionIterator<X>(_argument_size,_indices,_coefficients);
}

template<class X> ExpansionIterator<X> Expansion<X>::end() {
    return ExpansionIterator<X>(_argument_size,_indices+_size*(_argument_size+1u),_coefficients+_size);
}

template<class X> ExpansionConstIterator<X> Expansion<X>::begin() const {
    return ExpansionConstIterator<X>(_argument_size,const_cast<DegreeType*>(_indices),_coefficients);
}

template<class X> ExpansionConstIterator<X> Expansion<X>::end() const {
    return ExpansionConstIterator<X>(_argument_size,const_cast<DegreeType*>(_indices)+_size*(_argument_size+1u),_coefficients+_size);
}

template<class X> ExpansionIterator<X> Expansion<X>::find(const MultiIndex& a) {
    ExpansionIterator<X> iter=this->begin();
    ExpansionIterator<X> end=this->end();
    while(iter!=end && iter->index()!=a) {
        ++iter;
    }
    return iter;
}

template<class X> ExpansionConstIterator<X> Expansion<X>::find(const MultiIndex& a) const {
    ExpansionConstIterator<X> iter=this->begin();
    ExpansionConstIterator<X> end=this->end();
    while(iter!=end && iter->index()!=a) {
        ++iter;
    }
    return iter;
}

template<class CMP> struct IndexComparison {
    CMP cmp;
    template<class M1, class M2> auto operator()(M1 const& m1, M2 const& m2) -> decltype(cmp(m1.index(),m2.index())) {
        return cmp(m1.index(),m2.index()); }
};

template<class X> Void Expansion<X>::index_sort(ReverseLexicographicLess) {
    std::sort(this->begin(),this->end(),IndexComparison<ReverseLexicographicLess>());
}

template<class X> Void Expansion<X>::index_sort(GradedLess) {
    std::sort(this->begin(),this->end(),IndexComparison<GradedLess>());
}

template<class X> Void Expansion<X>::sort(ReverseLexicographicIndexLess) {
    std::sort(this->begin(),this->end(),IndexComparison<ReverseLexicographicLess>());
}

template<class X> Void Expansion<X>::sort(GradedIndexLess) {
    std::sort(this->begin(),this->end(),IndexComparison<GradedLess>());
}


template<class X>
Expansion<X> Expansion<X>::_embed(SizeType before_size, SizeType after_size) const
{
    const Expansion<X>& x=*this;
    SizeType old_size=x.argument_size();
    SizeType new_size=before_size+old_size+after_size;
    Expansion<X> r(new_size, this->_zero_coefficient);
    MultiIndex old_index(old_size);
    MultiIndex new_index(new_size);
    for(typename Expansion<X>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter) {
        old_index=iter->key();
        for(Nat j=0; j!=old_size; ++j) {
            Nat aj=old_index[j];
            new_index[j+before_size]=aj;
        }
        r.append(new_index,iter->data());
    }
    return r;
}

template<class X> OutputStream& Expansion<X>::write(OutputStream& os) const {
    os << "\nExpansion<X>{"<<this->number_of_terms()<<"/"<<this->capacity()<<",\n";
    for(auto iter=this->begin(); iter!=this->end(); ++iter) {
        os << "  "<<iter->index()<<":"<<iter->coefficient()<<",\n";
    }
    os << "}\n";
    return os;
}

template<class X>
OutputStream& Expansion<X>::write(OutputStream& os, const Array<String>& variable_names) const
{
    const Expansion<X>& p=*this;
    ARIADNE_ASSERT(p.argument_size()==variable_names.size());
    if(p.size()==0) {
        os << "0";
    } else {
        bool first_term=true;
        for(auto iter=p.begin(); iter!=p.end(); ++iter) {
            MultiIndex a=iter->index();
            X v=iter->coefficient();
            os << " ";
            if(decide(v>=0) && !first_term) { os<<"+"; }
            first_term=false;
            bool first_factor=true;
            if(decide(v<0)) { os<<"-"; }
            if(possibly(abs(v)!=1) || a.degree()==0) { os<<abs(v); first_factor=false; }
            for(SizeType j=0; j!=a.size(); ++j) {
                if(a[j]!=0) {
                    if(first_factor) { first_factor=false; } else { os <<"*"; }
                    os<<variable_names[j]; if(a[j]!=1) { os<<"^"<<int(a[j]); }
                }
            }
        }
    }
    return os;
}

template<class X, class CMP> SortedExpansion<X,CMP>::SortedExpansion(Expansion<X> e)
    : Expansion<X>(std::move(e))
{
    this->sort();
}

/*
template<class X, class CMP> auto SortedExpansion<X,CMP>::find(const MultiIndex& a) -> Iterator {
    ExpansionValue<X> term(a,Expansion<X>::_zero_coefficient);
    return  std::lower_bound(this->begin(),this->end(),a,CMP());
}

template<class X, class CMP> auto SortedExpansion<X,CMP>::find(const MultiIndex& a) const -> ConstIterator {
    ExpansionValue<X> term(a,Expansion<X>::_zero_coefficient);
    return std::lower_bound(this->begin(),this->end(),a,CMP());
}
*/
template<class X, class CMP> auto SortedExpansion<X,CMP>::at(const MultiIndex& a) -> CoefficientType& {
    ExpansionValue<X> term(a,Expansion<X>::_zero_coefficient);
    auto iter=std::lower_bound(this->begin(),this->end(),term,CMP());
    if (iter->index()!=a) {
        iter=this->Expansion<X>::insert(iter,a,X(0));
    }
    return iter->coefficient();
}

template<class X, class CMP> Void SortedExpansion<X,CMP>::insert(const MultiIndex& a, const X& c) {
    ExpansionValue<X> term(a,Expansion<X>::_zero_coefficient);
    auto iter=std::lower_bound(this->begin(),this->end(),term,CMP());
    if (iter->index()!=a) {
        iter=this->Expansion<X>::insert(iter,a,c);
    } else {
        ARIADNE_THROW(std::runtime_error,"Expansion<X>::set(const MultiIndex& a, const X& c):\n    e="<<*this,
                      " Index "<<a<<" already has a coefficient.");
    }
}

template<class X, class CMP> void SortedExpansion<X,CMP>::set(const MultiIndex& a, CoefficientType const& c) {
    this->at(a)=c;
}

template<class X, class CMP> void SortedExpansion<X,CMP>::sort() {
    std::sort(this->begin(),this->end(),CMP());
}


} // namespace Ariadne
