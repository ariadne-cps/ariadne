
/***************************************************************************
 *            multi_index.h
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

/*! \brief \file multi_index.h
 *  \brief An index specifying the degree of differentiation.
 */

#ifndef ARIADNE_MULTI_INDEX_H
#define ARIADNE_MULTI_INDEX_H

#include <cassert>
#include <cstdarg>
#include <cstring>
#include <iostream>
#include <vector>

#include <boost/iterator/iterator_facade.hpp>

#include "macros.h"
#include "array.h"
#include "numeric.h"

namespace Ariadne {

uint32_t fac(uint8_t);
uint32_t bin(uint8_t,uint8_t);

typedef unsigned char uchar;
typedef double Float;

class MultiIndexReference;
class MultiIndexValueReference;

class MultiIndexListIterator;

class MultiIndexBound;

template<class T> class Reference;
template<> class Reference<MultiIndexReference>;
template<> class Reference<const MultiIndexReference>;

//struct MultiIndexData { unsigned int _n; unsigned int* _p; };

//A reference to a MultiIndex
class MultiIndexReference {
    friend class MultiIndexListIterator;
    friend class MultiIndexListConstIterator;
  public:
    typedef unsigned int size_type;
    typedef unsigned char byte_type;
    typedef unsigned int word_type;
    typedef char raw_data_type;

    typedef byte_type index_type;
    typedef MultiIndexValueReference reference;
    typedef const index_type& const_reference;
    //typedef unsigned int word_type;
    //typedef unsigned long long int word_type;
  protected:
  public:
    //! \brief Constructors.
    MultiIndexReference(size_type n, index_type* p);
  public:
    //! \brief Assign values.
    MultiIndexReference& operator=(const MultiIndexReference& a);
    //! \brief Resize to hold n variables.
    void resize(size_type n);
    //! \brief Assigns values from another index. Precondition: the size of \a a must equal the current size.
    void assign(const MultiIndexReference& a);
    //! \brief Set all values to zero.
    void clear();
    //! \brief The number of variables.
    size_type size() const;
    //! \brief The degree of the multi-index, equal to the sum of the number of occurrences of the variables.
    index_type degree() const;
     //! \brief The number of variables.
    size_type number_of_variables() const;
    //! \brief The number of occurrences of the \a i th variable.
    index_type get(size_type i) const;
    //! \brief Set the number of occurrences of the \a i th variable to \a n.
    void set(size_type i, index_type n);
    //! \brief The number of occurrences of the \a i th variable.
    const_reference operator[](size_type i) const;
    //! \brief The number of occurrences of the \a i th variable.
    reference operator[](size_type i);
    //! \brief Increment the value of the \a ith element
    void increment(size_type i);
    //! \brief Decrement the value of the \a ith element
    void decrement(size_type i);

    //! \brief Equality operator.
    friend bool operator==(const MultiIndexReference& a1, const MultiIndexReference& a2); // inline
    //! \brief Inequality operator.
    friend bool operator!=(const MultiIndexReference& a1, const MultiIndexReference& a2); // inline
    //! \brief Comparison operator.
    friend bool operator<(const MultiIndexReference& a1, const MultiIndexReference& a2); // inline

    //! \brief Increment. No post-increment operator as we sometimes pass MultiIndexReference by reference.
    MultiIndexReference& operator++(); // inline
    // No post-increment operator as we sometimes pass MultiIndexReference by reference.Post increment.
    // MultiIndexReference operator++(int);
    //! \brief Inplace sum.
    MultiIndexReference& operator+=(const MultiIndexReference& a); // inline
    //! \brief Inplace difference.
    MultiIndexReference& operator-=(const MultiIndexReference& a); // inline
    //! \brief Inplace scalar product.
    MultiIndexReference& operator*=(const index_type& a); // inline
    //! \brief Sum.
    friend MultiIndexReference operator+(const MultiIndexReference& a1, const MultiIndexReference& a2); // inline
    //! \brief Difference.
    friend MultiIndexReference operator-(const MultiIndexReference& a1, const MultiIndexReference& a2); // inline
    //! \brief Scalar product.
    friend MultiIndexReference operator*(const MultiIndexReference& a, index_type s); // inline
    //! \brief Scalar product.
    friend MultiIndexReference operator*(index_type s, const MultiIndexReference& a); // inline

    //! \brief The position of the element in the Array of tensor values.
    unsigned int position() const;
    //! \brief The product of the factorials of the indices.
    unsigned int factorial() const;
    //! \brief The number of ordered index arrays with each element occurring the number of times specified by the multi index.
    unsigned int number() const;

    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream&, const MultiIndexReference&);
  public:
    //index_type& at(size_type i) { return _p[i]; }
    const index_type& at(size_type i) const { return _p[i]; }
    index_type& at(size_type i) { return _p[i]; }
    index_type* begin() { return _p; }
    const index_type* begin() const { return _p; }
    index_type* end() { return _p+(_n+1u); }
    const index_type* end() const { return _p+(_n+1u); }
  public:
    //size_type memory_size() { return sizeof(word_type)*(1+_n/sizeof(word_type)); }
    size_type memory_size() { return _n+1; }
    void allocate() { _p=new index_type[memory_size()]; }
    void deallocate() { delete[] _p; }
  private:
    size_type _n;
    index_type* _p;
};


class MultiIndex : public MultiIndexReference {
  public:
    //! \brief Destructor.
    ~MultiIndex();
    //! \brief Construct a multi index with no coefficients.
    explicit MultiIndex();
    //! \brief Construct a multi index of degree \a 0 with \a nv variables.
    explicit MultiIndex(size_type nv);
    //! \brief Construct a multi index with \a nv variables from the Array \a ary.
    explicit MultiIndex(size_type nv, const int* ary);
    //! \brief Construct a multi index with \a nv variables from variable arguments.
    explicit MultiIndex(size_type nv, int a1, ...);

    //! \brief Copy constructor.
    MultiIndex(const MultiIndex& a);
    //! \brief Copy from a reference.
    MultiIndex(const MultiIndexReference& a);
    //! \brief Copy assignment operator.
    MultiIndex& operator=(const MultiIndexReference& a);

    //! \brief Construct the zero multi index with \a nv variables.
    static MultiIndex zero(size_type nv);
    //! \brief Construct the unit multi index in variable \a j with \a nv variables.
    static MultiIndex unit(size_type nv, size_type j);
    //! \brief Construct the first multi index of degree \a d with \a nv variables.
    static MultiIndex first(size_type nv, index_type d);
};


class MultiIndexValueReference {
    typedef MultiIndexReference::size_type size_type;
    typedef MultiIndexReference::index_type index_type;
    typedef MultiIndexReference::word_type word_type;
  private:
    size_type _n; index_type* _p; size_type _i;
  public:
    MultiIndexValueReference(size_type n, index_type* p, size_type i) : _n(n), _p(p), _i(i) { }
    operator const index_type& () { return _p[_i]; }
    MultiIndexValueReference& operator=(const index_type& d) { _p[_n]+=(d-_p[_i]); _p[_i]=d; return *this; }
    MultiIndexValueReference& operator++() { ++_p[_n]; ++_p[_i]; return *this; }
    MultiIndexValueReference& operator--();
    MultiIndexValueReference& operator+=(int k) { _p[_n]+=k; _p[_i]+=k; return *this; }
    MultiIndexValueReference& operator-=(int k) { _p[_n]-=k; _p[_i]-=k; return *this; }
};

inline MultiIndexValueReference& MultiIndexValueReference::operator--() {
    if(_p[_i]==0) {
        ARIADNE_THROW(std::runtime_error,"--MultiIndexReference[i]"," decrementing zero value at "<<_i<<" in "<<MultiIndexReference(_n,_p)); }
    --_p[_n]; --_p[_i]; return *this;
}


inline MultiIndex::~MultiIndex()
{
    deallocate();
}

inline MultiIndex::MultiIndex()
    : MultiIndexReference(0,0)
{
    allocate();
    std::fill(begin(),end(),0);
}

inline MultiIndex::MultiIndex(size_type n)
    : MultiIndexReference(n,0)
{
    allocate();
    std::fill(this->begin(),this->end(),0);
}

inline MultiIndex::MultiIndex(size_type n, const int* ary)
    : MultiIndexReference(n,0)
{
    allocate();
    index_type* p=begin();
    p[n]=0; for(size_type i=0; i!=n; ++i) { p[i]=ary[i]; p[n]+=ary[i]; }
}

inline MultiIndex::MultiIndex(size_type n, int a0, ...)
    : MultiIndexReference(n,0)
{
    ARIADNE_ASSERT(n>0);
    allocate();
    index_type* p=begin();
    p[0]=a0; p[n]=a0;
    va_list args; va_start(args,a0);
    for(size_type i=1; i!=n; ++i) {
        p[i]=va_arg(args,int);
        p[n]+=p[i];
    }
    va_end(args);
}

inline MultiIndex::MultiIndex(const MultiIndex& a)
    : MultiIndexReference(a.size(),0)
{
    this->allocate();
    this->assign(a);
}

inline MultiIndex::MultiIndex(const MultiIndexReference& a)
    : MultiIndexReference(a.size(),0)
{
    this->allocate();
    this->assign(a);
}

inline MultiIndex& MultiIndex::operator=(const MultiIndexReference& a) {
    if(this!=&a) { this->resize(a.size()); this->assign(a); }
    return *this;
}

inline MultiIndex MultiIndex::zero(size_type n)
{
    return MultiIndex(n);
}

inline MultiIndex MultiIndex::unit(size_type n, size_type i)
{
    MultiIndex result(n);
    ++result[i];
    return result;
}

inline MultiIndex MultiIndex::first(size_type n, index_type d)
{
    MultiIndex result(n);
    result[0]=d;
    return result;
}


inline MultiIndexReference::MultiIndexReference(size_type n, index_type* p)
    : _n(n), _p(p) { }

inline MultiIndexReference& MultiIndexReference::operator=(const MultiIndexReference& a) {
    assert(this->size()==a.size()); this->assign(a); return *this;
}

inline void MultiIndexReference::assign(const MultiIndexReference& a) {
    memcpy(this->begin(),a.begin(),this->memory_size());
}

inline void MultiIndexReference::clear()
{
    memset(this->begin(),this->size(),0);
    //std::fill(this->begin(),this->end(),0);
}

inline void MultiIndexReference::resize(size_type n) {
    if(this->_n!=n) { deallocate(); this->_n=n; allocate(); }
}

inline MultiIndexReference::size_type MultiIndexReference::size() const {
    return this->_n;
}

inline MultiIndexReference::index_type MultiIndexReference::degree() const {
    return this->_p[this->_n];
}

inline MultiIndexReference::size_type MultiIndexReference::number_of_variables() const {
    return this->_n;
}

inline MultiIndexReference::index_type MultiIndexReference::get(size_type i) const {
    assert(i<this->size()); return this->_p[i];
}

inline void MultiIndexReference::set(size_type i, index_type k) {
    assert(i<this->size());
    index_type& ai=this->_p[i];
    index_type& d=this->_p[this->_n];
    d+=k; d-=ai; ai=k;
}

inline MultiIndexReference::index_type const& MultiIndexReference::operator[](size_type i) const {
    assert(i<this->size()); return this->_p[i];
}

inline MultiIndexValueReference MultiIndexReference::operator[](size_type i) {
    return MultiIndexValueReference(this->_n,this->_p,i);
}

inline void MultiIndexReference::increment(size_type i) {
    ++this->_p[i]; ++this->_p[this->_n];
}

inline void MultiIndexReference::decrement(size_type i) {
    ARIADNE_ASSERT(this->_p[i]>0u);
    --this->_p[i]; --this->_p[this->_n];
}


inline bool operator==(const MultiIndexReference& a1, const MultiIndexReference& a2) {
    if(a1.size()!=a2.size()) { return false; }
    //return memcmp(a1.begin(),a2.begin(),a1.size());
    for(MultiIndexReference::size_type i=0; i!=a1.size(); ++i) {
        if(a1.at(i)!=a2.at(i)) { return false; } }
    return true;
    //for(MultiIndexReference::size_type j=0; j!=a1.word_size(); ++j) {
    //    if(a1.word_at(j)!=a2.word_at(j)) { return false; } }
    //return true;
}

inline bool operator!=(const MultiIndexReference& a1, const MultiIndexReference& a2) {
    return !(a1==a2);
}

inline bool operator<(const MultiIndexReference& a1, const MultiIndexReference& a2) {
    if(a1.size()!=a2.size()) {
        throw std::runtime_error("operator<(MultiIndexReference,MultiIndexReference): number of variables must match");
    }

    if(a1.degree()!=a2.degree()) {
        return a1.degree()<a2.degree();
    } else {
        for(MultiIndexReference::size_type i=0; i!=a1.size(); ++i) {
            if(a1[i]!=a2[i]) {
                return a1[i]>a2[i];
            }
        }
        return false;
    }
}


inline
unsigned int MultiIndexReference::number() const
{
    unsigned int result=fac(this->degree());
    for(unsigned int k=0; k!=this->size(); ++k) {
        result/=fac((*this)[k]);
    }
    return result;
}

inline
unsigned int MultiIndexReference::factorial() const
{
    unsigned int result=1;
    for(unsigned int k=0; k!=this->size(); ++k) {
        result*=fac((*this)[k]);
    }
    return result;
}

inline
unsigned int MultiIndexReference::position() const
{
    unsigned int deg=this->degree()-1;
    unsigned int nvar=this->size();
    unsigned int result=bin(deg+nvar,nvar);
    for(unsigned int k=0; k!=this->size()-1; ++k) {
        --nvar;
        deg-=(*this)[k];
        result+=bin(deg+nvar,nvar);
    }
    return result;
}


inline
MultiIndexReference& MultiIndexReference::operator++()
{
    //std::cerr<<"MultiIndexReference::operator++() with *this="<<*this<<" "<<std::flush;
    assert(_n>0);

    size_type const n=this->_n;
    index_type* const p=this->_p;

    if(n==1) {
        ++p[0];
        ++p[n];
        return *this;
    }
    if(p[n-2]!=0) {
        --p[n-2];
        ++p[n-1];
        return *this;
    } else {
        index_type li=p[n-1];
        p[n-1]=0;
        for(size_type k=n-1; k!=0; --k) {
            if(p[k-1]!=0) {
                --p[k-1];
                p[k]=li+1;
                return *this;
            }
        }
        p[0]=li+1;
        ++p[n];
    }
    return *this;
}


inline
MultiIndexReference& MultiIndexReference::operator+=(const MultiIndexReference& a) {
    //for(size_type j=0; j!=this->word_size(); ++j) {
    //    this->word_at(j)+=a.word_at(j);
    //}
    //return *this;
    for(size_type i=0; i!=this->size()+1; ++i) {
        this->_p[i]+=a._p[i];
    }
    return *this;
}

inline
MultiIndexReference& MultiIndexReference::operator-=(const MultiIndexReference& a) {
    for(size_type i=0; i!=this->size()+1; ++i) {
        ARIADNE_ASSERT(this->_p[i]>=a._p[i]);
        this->_p[i]-=a._p[i];
    }
    return *this;
}

inline
MultiIndexReference& MultiIndexReference::operator*=(const index_type& s) {
    //for(size_type j=0; j!=this->word_size(); ++j) {
    //    this->word_at(j)*=s;
    //}
    //return *this;
    for(size_type i=0; i!=this->size()+1; ++i) {
        this->_p[i]*=s;
    }
    return *this;
}

inline
void iadd(MultiIndexReference& r, const MultiIndexReference& a1, const MultiIndexReference& a2) {
    for(MultiIndexReference::size_type j=0; j!=r.size(); ++j) {
        r.at(j)=a1.at(j)+a2.at(j);
    }
}

inline
MultiIndexReference operator+(const MultiIndexReference& a1, const MultiIndexReference& a2) {
    MultiIndexReference r(a1); r+=a2; return r;
}

inline
MultiIndexReference operator-(const MultiIndexReference& a1, const MultiIndexReference& a2) {
    MultiIndexReference r(a1); r-=a2; return r;
}

inline
MultiIndexReference operator*(const MultiIndexReference& a, MultiIndexReference::index_type s) {
    MultiIndexReference r(a); r*=s; return r;
}

inline
MultiIndexReference operator*(MultiIndexReference::index_type s, const MultiIndexReference& a) {
    MultiIndexReference r(a); r*=s; return r;
}


inline
unsigned int
number(const MultiIndexReference& i)
{
    unsigned int result=fac(i.degree());
    for(unsigned int k=0; k!=i.size(); ++k) {
        result/=fac(i[k]);
    }
    return result;
}

inline
unsigned int
fac(const MultiIndexReference& i)
{
    unsigned int result=1;
    for(unsigned int k=0; k!=i.size(); ++k) {
        result*=fac(i[k]);
    }
    return result;
}

inline
unsigned int
bin(const MultiIndexReference& n, const MultiIndexReference& k)
{
    assert(n.size()==k.size());
    unsigned int result=1;
    for(unsigned int i=0; i!=n.size(); ++i) {
        result*=bin(n[i],k[i]);
    }
    return result;
}

inline
std::ostream& operator<<(std::ostream& os, const MultiIndexReference& a) {
    //os<<"<"<<static_cast<const void*>(a.begin())<<">";
    if(a.size()==0) { os << '('; }
    for(MultiIndex::size_type i=0; i!=a.size(); ++i) { os << (i==0?'(':',') << int(a[i]); }
    //os<<";"<<int(a.degree());
    return os << ')';
}






//! \brief \brief A bound on a MultiIndexReference object, allowing different groups of
//!variables to have different maximum degrees.

class MultiIndexBound {
  public:
    typedef MultiIndexReference::size_type size_type;
    MultiIndexBound(size_type as, size_type d);
    MultiIndexBound(const MultiIndexReference& a);
    MultiIndexBound(size_type as, size_type ng, size_type g1s, size_type g1d, ...);
    size_type size() const { return _groups.size(); }
    friend bool operator<=(const MultiIndexReference& a, const MultiIndexBound& b);
  private:
    Array<size_type> _groups;
    Array<size_type> _max_degrees;
};

inline MultiIndexBound::MultiIndexBound(size_type as, size_type d)
    : _groups(as), _max_degrees(1u)
{
    for(size_type i=0; i!=as; ++i) {
        _groups[i]=0u;
    }
    _max_degrees[0]=d;
}


inline MultiIndexBound::MultiIndexBound(const MultiIndexReference& a)
    : _groups(a.size()), _max_degrees(a.size())
{
    for(size_type i=0; i!=a.size(); ++i) {
        _groups[i]=i;
        _max_degrees[i]=a[i];
    }
}


inline MultiIndexBound::MultiIndexBound(size_type as, size_type ng, size_type g1s, size_type g1d, ...)
    : _groups(as), _max_degrees(ng)
{
    va_list args; va_start(args,g1d);
    size_type k=0;
    for( ; k!=g1s; ++k) {
        _groups[k]=0;
    }
    _max_degrees[0]=g1d;
    for(size_type i=1; i!=ng; ++i) {
        size_type nk=k+va_arg(args,int);
        for( ; k!=nk; ++k) {
            _groups[k]=i;
        }
        _max_degrees[i]=va_arg(args,int);
    }
    va_end(args);
    ARIADNE_ASSERT(k==as);
}

inline bool operator<=(const MultiIndexReference& a, const MultiIndexBound& b) {
    typedef MultiIndexReference::size_type size_type;
    Array<size_type> degrees(b._max_degrees.size());
    for(size_type j=0; j!=a.size(); ++j) {
        degrees[b._groups[j]]+=a[j];
    }
    for(size_type k=0; k!=degrees.size(); ++k) {
        if(degrees[k]>b._max_degrees[k]) {
            return false;
        }
    }
    return true;
}



class MultiIndexListIterator
    : public boost::iterator_facade<MultiIndexListIterator,MultiIndex,boost::random_access_traversal_tag,MultiIndexReference>
{
    MultiIndexReference _ref;
  public:
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::index_type index_type;

    MultiIndexListIterator() : _ref(0u,0) { }
    MultiIndexListIterator(size_type as, index_type* p) : _ref(as,p) { }
    MultiIndexListIterator(const MultiIndexListIterator& iter) : _ref(iter._ref) { }
    MultiIndexListIterator& operator=(const MultiIndexListIterator& iter) { _ref._n=iter._ref._n; _ref._p=iter._ref._p; return *this; }

    bool equal(const MultiIndexListIterator& other) const { return this->_ref._p==other._ref._p; }
    MultiIndexReference dereference() const { return _ref; }
    void increment() { advance(1); }
    void decrement() { advance(-1); }
    void advance(int m) { _ref._p+=m*(_ref._n+1); }

    friend std::ostream& operator<<(std::ostream& os, const MultiIndexListIterator& iter) {
        return os << "<" << &iter._ref[0] << ">" << iter._ref << "\n";
    }
};

class MultiIndexListConstIterator
    : public boost::iterator_facade<MultiIndexListConstIterator,MultiIndex,boost::random_access_traversal_tag,const MultiIndexReference>
{
    MultiIndexReference _ref;
  public:
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::index_type index_type;

    MultiIndexListConstIterator() : _ref(0u,0) { }
    MultiIndexListConstIterator(size_type as, const index_type* p) : _ref(as,const_cast<index_type*>(p)) { }
    MultiIndexListConstIterator(MultiIndexListIterator iter) : _ref(*iter) { }
    MultiIndexListConstIterator(const MultiIndexListConstIterator& iter) : _ref(iter._ref) { }
    MultiIndexListConstIterator& operator=(const MultiIndexListConstIterator& iter) { _ref._n=iter._ref._n; _ref._p=iter._ref._p; return *this; }
    bool equal(const MultiIndexListConstIterator& other) const { return this->_ref._p==other._ref._p; }
    const MultiIndexReference dereference() const { return _ref; }
    void increment() { advance(1); }
    void decrement() { advance(-1); }
    void advance(int m) { _ref._p+=m*(_ref._n+1); }

    friend std::ostream& operator<<(std::ostream& os, const MultiIndexListConstIterator& iter) {
        return os << "<" << &iter._ref[0] << ">" << iter._ref << "\n";
    }
};


class MultiIndexList
{
  public:
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::index_type index_type;

    typedef MultiIndex value_type;
    typedef MultiIndexReference reference;
    typedef const MultiIndexReference const_reference;
    typedef MultiIndexListIterator iterator;
    typedef MultiIndexListConstIterator const_iterator;

    MultiIndexList() : _as(0u), _d() { }
    MultiIndexList(size_type as) : _as(as), _d() { }
    bool operator==(const MultiIndexList& other) const { return this==&other || (this->_as==other._as && this->_d==other._d); }
    bool operator!=(const MultiIndexList& other) const { return !(*this==other); }
    size_type element_size() const { return _as; }
    size_type size() const { return _d.size()/(_as+1); }
    size_type capacity() const { return _d.capacity()/(_as+1); }
    void clear() { _d.clear(); }
    void resize(size_type n) { _d.resize(n*(_as+1)); }
    void reserve(size_type c) { _d.reserve(c*(_as+1)); };
    void append(const MultiIndexReference& a) { this->resize(this->size()+1); (*this)[this->size()-1]=a; }

    reference operator[](size_type i) { return MultiIndexReference(_as,_begin_pointer()+i*(_as+1)); }
    const_reference operator[](size_type i) const { return MultiIndexReference(_as,const_cast<index_type*>(_begin_pointer())+i*(_as+1)); }
    reference front() { return (*this)[0]; }
    const_reference front() const { return (*this)[0]; }
    reference back() { return (*this)[this->size()-1]; }
    const_reference back() const { return (*this)[this->size()-1]; }

    iterator begin() { return iterator(_as,_begin_pointer()); }
    iterator end() { return iterator(_as,_begin_pointer()+size()*(_as+1u)); }
    const_iterator begin() const { return const_iterator(_as,_begin_pointer()); }
    const_iterator end() const { return const_iterator(_as,_begin_pointer()+size()*(_as+1u)); }
  private:
  public:
    size_type _data_size() const { return _d.size(); }
    index_type* _begin_pointer() { return &_d[0]; }
    const index_type* _begin_pointer() const { return &_d[0]; }
    const std::vector<index_type>& _data() const { return _d; }
  private:
    size_type _as; std::vector<index_type> _d;
};

inline std::ostream& operator<<(std::ostream& os, const MultiIndexList& lst) {
    os << "MultiIndexList";
    os << "<" << lst.size() << "," << lst.capacity() <<  "," << lst.element_size() << "," << lst._data_size() << ">";
    for(MultiIndexList::size_type i=0; i!=lst.size(); ++i) { os << (i==0?"[":",") << lst[i]; }
    return os << "]";
}

} // namespace Ariadne

#endif /* ARIADNE_MULTI_INDEX_H */
