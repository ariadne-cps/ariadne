/***************************************************************************
 *            multi_index.inl.hpp
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

namespace Ariadne {

inline MultiIndexValueReference& MultiIndexValueReference::operator--() {
    if(_p[_i]==0) {
        ARIADNE_THROW(std::runtime_error,"--MultiIndex[i]"," decrementing zero value at "<<_i<<" in "<<reinterpret_cast<const MultiIndex&>(*this)); }
    --_p[_n]; --_p[_i]; return *this;
}


inline MultiIndex::~MultiIndex()
{
    _deallocate(_p);
}

inline MultiIndex::MultiIndex()
    : MultiIndexData(0,_allocate(0))
{
    _p[0]=0;
}

inline MultiIndex::MultiIndex(SizeType n)
    : MultiIndexData(n,_allocate(n))
{
    std::fill(this->begin(),this->end(),0u); _p[n]=0;
}

inline MultiIndex::MultiIndex(SizeType n, const DegreeType* ary)
    : MultiIndexData(n,_allocate(n))
{
    DegreeType& _d=*(_p+_n); _d=0u;
    for(SizeType i=0; i!=n; ++i) { _p[i]=ary[i]; _d+=ary[i]; }
}

inline MultiIndex::MultiIndex(InitializerList<DegreeType> lst)
    : MultiIndexData(lst.size(),_allocate(lst.size()))
{
    DegreeType* _c=_p;
    DegreeType& _d=*(_p+_n); _d=0u;
    auto iter=lst.begin();
    while(iter!=lst.end()) {
        *_c=*iter;
        _d+=*_c;
        ++_c; ++iter;
    }
}

inline MultiIndex::MultiIndex(const MultiIndex& a)
    : MultiIndexData(a._n,_allocate(a._n))
{
    this->assign(a);
}

inline MultiIndex& MultiIndex::operator=(const MultiIndex& a) {
    if(this!=&a) { this->resize(a.size()); this->assign(a); }
    return *this;
}

inline Void MultiIndex::assign(const MultiIndex& a) {
    for(SizeType j=0; j!=this->size()+1; ++j) { this->_p[j]=a._p[j]; }
    //std::copy(a.word_begin(),a.word_end(),this->word_begin());
}

inline MultiIndex MultiIndex::zero(SizeType n)
{
    return MultiIndex(n);
}

inline MultiIndex MultiIndex::unit(SizeType n, SizeType i)
{
    MultiIndex r(n);
    r._p[i]=1u;
    r._p[r._n]=1u;
    return r;
}

inline MultiIndex MultiIndex::first(SizeType n, DegreeType d)
{
    MultiIndex r(n);
    r._p[0]=d;
    r._p[r._n]=d;
    return r;
}

inline Void MultiIndex::clear()
{
    std::fill(this->begin(),this->end(),0u);
    this->_p[_n]=0u;
}

inline Void MultiIndex::resize(SizeType n) {
    if(this->_n!=n) { _deallocate(this->_p); this->_n=n; this->_p=_allocate(this->_n); }
}

inline SizeType MultiIndex::size() const {
    return this->_n;
}

inline DegreeType MultiIndex::degree() const {
    return this->_p[_n];
}

inline SizeType MultiIndex::number_of_variables() const {
    return this->_n;
}

inline DegreeType MultiIndex::get(SizeType i) const {
    assert(i<this->size()); return this->_p[i];
}

inline Void MultiIndex::set(SizeType i, DegreeType k) {
    assert(i<this->size());
    DegreeType& ai=this->_p[i];
    DegreeType& d=this->_p[this->_n];
    d+=k; d-=ai; ai=k;
}

inline DegreeType const& MultiIndex::operator[](SizeType i) const {
    if(i>=this->size()) { std::cerr<<*this<<"["<<i<<"]"; }
    assert(i<this->size()); return this->_p[i];
}

inline MultiIndexValueReference MultiIndex::operator[](SizeType i) {
    return MultiIndexValueReference(this->_n,this->_p,i);
}

inline Void MultiIndex::increment(SizeType i) {
    ++this->_p[i]; ++this->_p[_n];
}

inline Void MultiIndex::decrement(SizeType i) {
    ARIADNE_ASSERT(this->_p[i]>0u);
    --this->_p[i]; --this->_p[_n];
}


inline Bool operator==(const MultiIndex& a1, const MultiIndex& a2) {
    if(a1.size()!=a2.size()) { return false; }
    for(MultiIndex::SizeType i=0; i!=a1.size(); ++i) {
        if(a1.at(i)!=a2.at(i)) { return false; } }
    return true;
}

inline Bool operator!=(const MultiIndex& a1, const MultiIndex& a2) {
    return !(a1==a2);
}



inline
MultiIndex& MultiIndex::operator+=(const MultiIndex& a) {
    for(SizeType i=0; i!=this->size()+1; ++i) {
        this->_p[i]+=a._p[i];
    }
    return *this;
}

inline
MultiIndex& MultiIndex::operator-=(const MultiIndex& a) {
    for(SizeType i=0; i!=this->size()+1; ++i) {
        ARIADNE_ASSERT(this->_p[i]>=a._p[i]);
        this->_p[i]-=a._p[i];
    }
    return *this;
}

inline
MultiIndex& MultiIndex::operator*=(const IndexType& s) {
    for(SizeType i=0; i!=this->size()+1; ++i) {
        this->_p[i]*=s;
    }
    return *this;
}

inline
Void iadd(MultiIndex& r, const MultiIndex& a1, const MultiIndex& a2) {
    for(SizeType j=0; j!=r.size()+1; ++j) {
        r._p[j]=a1._p[j]+a2._p[j];
    }
}

inline
MultiIndex operator+(MultiIndex a1, const MultiIndex& a2) {
    a1+=a2; return a1;
}

inline
MultiIndex operator-(MultiIndex a1, const MultiIndex& a2) {
    a1-=a2; return a1;
}

inline
MultiIndex operator*(MultiIndex a, DegreeType s) {
    a*=s; return a;
}

inline
MultiIndex operator*(DegreeType s, MultiIndex a) {
    a*=s; return a;
}

inline
Void swap(MultiIndex& a1, MultiIndex& a2) {
    ARIADNE_ASSERT(a1._n==a2._n); DegreeType t;
    for(SizeType i=0; i!=a1.size()+1; ++i) { t=a1._p[i]; a1._p[i]=a2._p[i]; a2._p[i]=t; }
    return;
    // FIXME: The code below works, but only if the MultiIndex objects are not casts of references to objects in lists
    std::swap(a1._n,a2._n);
    std::swap(a1._p,a2._p);
}

inline
unsigned int
number(const MultiIndex& i)
{
    unsigned int result=fac(i.degree());
    for(unsigned int k=0; k!=i.size(); ++k) {
        result/=fac(i[k]);
    }
    return result;
}

inline
unsigned int
fac(const MultiIndex& i)
{
    unsigned int result=1;
    for(unsigned int k=0; k!=i.size(); ++k) {
        result*=fac(i[k]);
    }
    return result;
}

inline
unsigned int
bin(const MultiIndex& n, const MultiIndex& k)
{
    assert(n.size()==k.size());
    unsigned int result=1;
    for(unsigned int i=0; i!=n.size(); ++i) {
        result*=bin(n[i],k[i]);
    }
    return result;
}


inline MultiIndexValueReference MultiIndexData::operator[](SizeType i) {
    return static_cast<MultiIndex&>(*this)[i];
}




//! \brief \brief A bound on a MultiIndex object, allowing different groups of
//!variables to have different maximum degrees.

inline MultiIndexBound::MultiIndexBound(SizeType as, SizeType d)
    : _groups(as), _max_degrees(1u)
{
    for(SizeType i=0; i!=as; ++i) {
        _groups[i]=0u;
    }
    _max_degrees[0]=d;
}


inline MultiIndexBound::MultiIndexBound(const MultiIndex& a)
    : _groups(a.size()), _max_degrees(a.size())
{
    for(SizeType i=0; i!=a.size(); ++i) {
        _groups[i]=i;
        _max_degrees[i]=a[i];
    }
}


inline Bool operator<=(const MultiIndex& a, const MultiIndexBound& b) {
    typedef MultiIndex::SizeType SizeType;
    Array<SizeType> degrees(b._max_degrees.size());
    for(SizeType j=0; j!=a.size(); ++j) {
        degrees[b._groups[j]]+=a[j];
    }
    for(SizeType k=0; k!=degrees.size(); ++k) {
        if(degrees[k]>b._max_degrees[k]) {
            return false;
        }
    }
    return true;
}


} // namespace Ariadne
