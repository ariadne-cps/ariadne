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

template<class X> X Expansion<X>::_zero = static_cast<X>(0);

template<class X> Expansion<X>::~Expansion() { }

template<class X> Expansion<X>::Expansion() : Expansion(0u) { }

template<class X> Expansion<X>::Expansion(SizeType as) : _argument_size(as) { }

template<class X> Expansion<X>::Expansion(SizeType as, DegreeType deg, InitializerList<X> lst)
    : _argument_size(as)
{
    MultiIndex a(as); X x;
    auto iter = lst.begin();
    while(a.degree()<=deg) {
        x=*iter;
        if(x!=X(0)) { this->append(a,x); }
        ++a;
        ++iter;
    }
}


template<class X> Expansion<X>::Expansion(SizeType as, InitializerList<PairType<InitializerList<Int>,X>> lst)
    : _argument_size(as)
{
    MultiIndex a;
    X x;
    for(auto iter=lst.begin();
        iter!=lst.end(); ++iter)
    {
        a=iter->first;
        x=iter->second;
        if(x!=X(0)) { this->append(a,x); }
    }
}

template<class X> Expansion<X>::Expansion(InitializerList<PairType<InitializerList<Int>,X>> lst)
    : Expansion(lst.size()==0?0u:lst.begin()->first.size(),lst)
{ }


template<class X> Void Expansion<X>::swap(Expansion<X>& other) {
    std::swap(this->_argument_size,other._argument_size);
    std::swap(this->_coefficients,other._coefficients);
}


template<class X> Bool Expansion<X>::operator==(const Expansion<X>& other) const {
    if(this->argument_size()!=other.argument_size()) { return false; }
    if(this->number_of_nonzeros()!=other.number_of_nonzeros()) { return false; }
    ConstIterator self_iter=this->begin(); ConstIterator other_iter=other.begin();
    while(self_iter!=this->end()) {
        if(self_iter->key()!=other_iter->key() || self_iter->data()!=other_iter->data()) { return false; }
        ++self_iter; ++other_iter; }
    return true;
}

template<class X> Bool Expansion<X>::operator!=(const Expansion<X>& other) const { return !this->operator==(other); }

template<class X> SizeType Expansion<X>::argument_size() const { return this->_argument_size; }

template<class X> SizeType Expansion<X>::number_of_nonzeros() const { return _coefficients.size()/_element_size(); }

template<class X> DegreeType Expansion<X>::degree() const {
   DegreeType deg=0u; for(ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        deg=std::max(deg,iter->key().degree());
    }
    return deg;
}

template<class X> Bool Expansion<X>::empty() const { return this->_coefficients.empty(); }
template<class X> SizeType Expansion<X>::size() const { return this->_coefficients.size()/_element_size(); }
template<class X> Void Expansion<X>::reserve(SizeType nnz) { this->_coefficients.reserve(nnz*_element_size()); }
template<class X> Void Expansion<X>::resize(SizeType nnz) { this->_coefficients.resize(nnz*_element_size()); }



template<class X> Void Expansion<X>::append(const MultiIndex& a, const RealType& c) {
    this->_append(a,c);
}
template<class X> Void Expansion<X>::prepend(const MultiIndex& a, const RealType& c) {
    this->_prepend(a,c);
}
template<class X> Void Expansion<X>::append_sum(const MultiIndex& a1, const MultiIndex& a2, const RealType& c) {
    this->_append(a1,a2,c);
}

template<class X> const X& Expansion<X>::operator[](const MultiIndex& a) const {
    ConstIterator iter=this->find(a);
    if(iter==this->end()) { return _zero; }
    else { return iter->data(); }
}

template<class X> typename Expansion<X>::Iterator Expansion<X>::begin() { return Iterator(_argument_size,_index_size(),_begin_ptr()); }
template<class X> typename Expansion<X>::Iterator Expansion<X>::end() { return Iterator(_argument_size,_index_size(),_end_ptr()); }
template<class X> typename Expansion<X>::Iterator Expansion<X>::find(const MultiIndex& a) {
    Iterator iter=this->end(); while(iter!=this->begin()) { --iter; if(iter->key()==a) { return iter; } } return this->end(); }

template<class X> typename Expansion<X>::ConstIterator Expansion<X>::begin() const { return ConstIterator(_argument_size,_index_size(),const_cast<WordType*>(_begin_ptr())); }
template<class X> typename Expansion<X>::ConstIterator Expansion<X>::end() const { return ConstIterator(_argument_size,_index_size(),const_cast<WordType*>(_end_ptr())); }
template<class X> typename Expansion<X>::ConstIterator Expansion<X>::find(const MultiIndex& a) const {
    ConstIterator iter=this->end(); while(iter!=this->begin()) { --iter; if(iter->key()==a) { return iter; } } return this->end();
}

template<class X> Void Expansion<X>::erase(Iterator iter) { iter->data()=static_cast<X>(0u); }
template<class X> Void Expansion<X>::clear() { _coefficients.clear(); }

template<class X> Void Expansion<X>::remove_zeros() {
    this->resize(std::remove_if(this->begin(),this->end(),DataIsZero())-this->begin()); }

template<class X> Void Expansion<X>::combine_terms() {
    Iterator curr=this->begin(); ConstIterator adv=this->begin(); ConstIterator end=this->end();
    while(adv!=end) { curr->key()=adv->key(); curr->data()=adv->data(); ++adv;
        while(adv!=end && curr->key()==adv->key()) { curr->data()+=adv->data(); ++adv; } ++curr; }
    this->resize(curr-this->begin()); }

template<> Void Expansion<ExactFloat>::combine_terms() { ARIADNE_NOT_IMPLEMENTED; }
template<> Void Expansion<ExactInterval>::combine_terms() { ARIADNE_NOT_IMPLEMENTED; }

template<class X> Void Expansion<X>::graded_sort() {
    this->_sort(GradedKeyLess());
}

template<class X> Void Expansion<X>::lexicographic_sort() {
    this->_sort(LexicographicKeyLess());
}

template<class X> Void Expansion<X>::reverse_lexicographic_sort() {
    this->_sort(ReverseLexicographicKeyLess());
}

template<class X> Void Expansion<X>::check() const { }

template<class X> typename Expansion<X>::SizeType Expansion<X>::_vector_size() const {
    return _coefficients.size(); }
template<class X> typename Expansion<X>::SizeType Expansion<X>::_index_size() const {
    return (_argument_size+sizeof_word)/sizeof_word; }

template<class X> typename Expansion<X>::SizeType Expansion<X>::_element_size() const {
    return (_argument_size+sizeof_word)/sizeof_word+sizeof_data/sizeof_word; }

template<class X> const typename Expansion<X>::WordType* Expansion<X>::_begin_ptr() const { return _coefficients.begin().operator->(); }
template<class X> const typename Expansion<X>::WordType* Expansion<X>::_end_ptr() const { return _coefficients.end().operator->(); }
template<class X> typename Expansion<X>::WordType* Expansion<X>::_begin_ptr() { return _coefficients.begin().operator->(); }
template<class X> typename Expansion<X>::WordType* Expansion<X>::_end_ptr() { return _coefficients.end().operator->(); }

template<class X> typename Expansion<X>::Iterator Expansion<X>::_insert(Iterator p, const MultiIndex& a, const RealType& x) {
    //std::cerr<<"_insert "<<*this<<" "<<p._ptr()<<" "<<a<<" "<<x<<std::endl;
    if(_coefficients.size()+_element_size()>_coefficients.capacity()) {
        DifferenceType i=p-begin();
        _coefficients.resize(_coefficients.size()+_element_size());
        p=begin()+i;
    } else {
        _coefficients.resize(_coefficients.size()+_element_size());
    }
    return _allocated_insert(p,a,x);
}

template<class X> typename Expansion<X>::Iterator Expansion<X>::_allocated_insert(Iterator p, const MultiIndex& a, const RealType& x) {
    //std::cerr<<"_allocated_insert "<<*this<<" "<<p<<" "<<p-begin()<<" "<<a<<" "<<x<<std::endl;
    Iterator curr=this->end()-1; Iterator prev=curr;
    while(curr!=p) { --prev; *curr=*prev; curr=prev; }
    curr->key()=a; curr->data()=x; return p;
}

template<class X> Void Expansion<X>::_prepend(const MultiIndex& a, const RealType& x) {
    //std::cerr<<"_prepend "<<*this<<" "<<a<<" "<<x<<std::endl;
    _coefficients.resize(_coefficients.size()+_element_size());
    _allocated_insert(begin(),a,x);
}

template<class X> Void Expansion<X>::_append(const MultiIndex& a, const RealType& x) {
    //std::cerr<<"_append "<<*this<<" "<<a<<" "<<x<<"... "<<std::flush;
    _coefficients.resize(_coefficients.size()+_element_size());
    WordType* vp=_end_ptr()-_element_size(); const WordType* ap=a.word_begin();
    for(SizeType j=0; j!=_index_size(); ++j) { vp[j]=ap[j]; }
    DataType* xp=reinterpret_cast<DataType*>(this->_end_ptr())-1; *xp=x;
    //std::cerr<<"done"<<std::endl;
}

template<class X> Void Expansion<X>::_append(const MultiIndex&  a1, const MultiIndex&  a2, const RealType& x) {
    //std::cerr<<"_append "<<*this<<" "<<a1<<" "<<a2<<" "<<x<<std::endl;
    _coefficients.resize(_coefficients.size()+_element_size());
    WordType* vp=_end_ptr()-_element_size();
    const WordType* ap1=a1.word_begin(); const WordType* ap2=a2.word_begin();
    for(SizeType j=0; j!=_index_size(); ++j) { vp[j]=ap1[j]+ap2[j]; }
    DataType* xp=reinterpret_cast<DataType*>(this->_end_ptr())-1; *xp=x;
}




template<class X>
Expansion<X> Expansion<X>::_embed(SizeType before_size, SizeType after_size) const
{
    const Expansion<X>& x=*this;
    Nat old_size=x.argument_size();
    Nat new_size=before_size+old_size+after_size;
    Expansion<X> r(new_size);
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

template<class X> Expansion<RawFloat>& Expansion<X>::raw() { return reinterpret_cast<Expansion<RawFloat>&>(*this); }

template<class X> Expansion<RawFloat>const& Expansion<X>::raw() const { return reinterpret_cast<Expansion<RawFloat>const&>(*this); }


template<class X>
OutputStream& Expansion<X>::write(OutputStream& os, const Array<std::string>& variable_names) const
{
    ARIADNE_ASSERT(this->argument_size()==variable_names.size());
    const Expansion<X>& p=*this;
    if(p.size()==0) {
        os << "0";
    } else {
        Bool first_term=true;
        for(ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
            MultiIndex a=iter->key();
            X v=iter->data();
            if(decide(v>=X(0)) && !first_term) { os<<"+"; }
            first_term=false;
            Bool first_factor=true;
            if(decide(v<X(0))) { os<<"-"; }
            if(abs(v)!=X(1) || a.degree()==0) { os<<abs(v); first_factor=false; }
            for(Nat j=0; j!=a.size(); ++j) {
                if(a[j]!=0) {
                    if(first_factor) { first_factor=false; } else { os <<"*"; }
                    os<<variable_names[j]; if(a[j]!=1) { os<<"^"<<Int(a[j]); }
                }
            }
        }
    }
    return os;
}

template<class X>
OutputStream& Expansion<X>::write_polynomial(OutputStream& os) const {
    const Expansion<X>& p=*this;
    Array<std::string> variable_names(p.argument_size());
    for(Nat j=0; j!=p.argument_size(); ++j) {
        StringStream sstr;
        sstr << 'x' << j;
        variable_names[j]=sstr.str();
    }
    return p.write(os,variable_names);
}

template<class X>
OutputStream& Expansion<X>::write_map(OutputStream& os) const {
    os << "E[" << this->argument_size() << ";" << this->number_of_nonzeros() << "]";
    os << "{";
    for(auto iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter!=this->begin()) { os << ", "; }
        os << iter->key() << ":" << iter->data();
    }
    os << "}";
    return os;
}

template<class X>
OutputStream& Expansion<X>::write(OutputStream& os) const {
    return this->write_map(os);
}





template<class T>
inline Vector< Expansion<MidpointType<T>> > midpoint(const Vector< Expansion<T> >& pse) {
    Vector< Expansion<MidpointType<T>> > r(pse.size(),Expansion<MidpointType<T>>());
    for(Nat i=0; i!=pse.size(); ++i) {
        r[i]=midpoint(pse[i]); }
    return r;
}




template<class X> template<class CMP> Void Expansion<X>::_sort(const CMP& cmp) {
    std::sort(this->begin(),this->end(),cmp);
}

template<class X> template<class CMP> typename Expansion<X>::CoefficientType& Expansion<X>::_at(const MultiIndex& a,const CMP& cmp) {
    Iterator p=std::lower_bound(this->begin(),this->end(),a,cmp);
    if(p!=this->end() && p->key()==a) { return p->data(); }
    else { p=this->_insert(p,a,X(0)); return p->data(); }
}

template<class X> template<class CMP> typename Expansion<X>::Iterator Expansion<X>::_insert(const MultiIndex& a, const CoefficientType& x,const CMP& cmp) {
    //std::cerr<<"_insert "<<*this<<" "<<a<<" "<<x<<std::endl;
    _coefficients.resize(_coefficients.size()+_element_size());
    Iterator p=std::lower_bound(this->begin(),this->end()-1,a,cmp);
    return _allocated_insert(p,a,x);
}



template<class X, class CMP> SortedExpansion<X,CMP>::SortedExpansion(Expansion<X> e)
    : Expansion<X>(std::move(e))
{
    this->sort();
}

template<class X, class CMP> Void SortedExpansion<X,CMP>::insert(const MultiIndex& a, const CoefficientType& c) {
    this->Expansion<X>::_insert(a,c,CMP());
}

template<class X, class CMP> typename SortedExpansion<X,CMP>::CoefficientType& SortedExpansion<X,CMP>::at(const MultiIndex& a) {
    return this->Expansion<X>::_at(a,CMP());
}

template<class X, class CMP> Void SortedExpansion<X,CMP>::set(const MultiIndex& a, const CoefficientType& c) {
    this->at(a)=c;
}

template<class X, class CMP> typename SortedExpansion<X,CMP>::CoefficientType const& SortedExpansion<X,CMP>::get(const MultiIndex& a) const {
    return const_cast<SortedExpansion<X,CMP>*>(this)->at(a);
}

template<class X, class CMP> typename Expansion<X>::Iterator SortedExpansion<X,CMP>::find(const MultiIndex& a) {
    Iterator iter=std::lower_bound(this->begin(),this->end(),a,CMP()); if(iter!=this->end() && iter->key()!=a) { iter=this->end(); } return iter;
}

template<class X, class CMP> typename Expansion<X>::ConstIterator SortedExpansion<X,CMP>::find(const MultiIndex& a) const {
    Iterator iter=std::lower_bound(this->begin(),this->end(),a,CMP()); if(iter!=this->end() && iter->key()!=a) { iter=this->end(); } return iter;
}

template<class X, class CMP> Void SortedExpansion<X,CMP>::sort() {
    this->Expansion<X>::_sort(CMP()); }




}

