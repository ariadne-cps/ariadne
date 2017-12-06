/***************************************************************************
 *            algebra/expansion.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file algebra/expansion.hpp
 *  \brief
 */



#ifndef ARIADNE_EXPANSION_HPP
#define ARIADNE_EXPANSION_HPP

#include "multi_index.hpp"

#include "utility/typedefs.hpp"
#include "utility/iterator.hpp"
#include "utility/macros.hpp"

namespace Ariadne {

/************ Expansion ******************************************************/

template<class I, class X> class ExpansionValue;
template<class I, class X> class ExpansionReference;
template<class I, class X> class ExpansionConstReference;
template<class I, class X> class ExpansionIterator;
template<class I, class X> class ExpansionConstIterator;
template<class I, class X> class ExpansionValueReference;

template<class T> struct ReferenceTypedef { typedef T& Type; };
template<class T> using ReferenceType = typename ReferenceTypedef<T>::Type;
template<class T> using ConstReferenceType = typename ReferenceTypedef<const T>::Type;

struct GradedIndexLess {
    template<class M1, class M2> bool operator()(M1 const& m1, M2 const& m2) const {
        return graded_less(m1.index(),m2.index());
    }
    template<class M> bool operator()(M const& m, typename M::IndexType const& a) const {
        return graded_less(m.index(),a);
    }
};

struct LexicographicIndexLess {
    template<class M1, class M2> bool operator()(M1 const& m1, M2 const& m2) const {
        return lexicographic_less(m1.index(),m2.index());
    }
    template<class M> bool operator()(M const& m, typename M::IndexType const& a) const {
        return lexicographic_less(m.index(),a);
    }
};

struct ReverseLexicographicIndexLess {
    template<class M1, class M2> bool operator()(M1 const& m1, M2 const& m2) const {
        return reverse_lexicographic_less(m1.index(),m2.index());
    }
    template<class M> bool operator()(M const& m, typename M::IndexType const& a) const {
        return reverse_lexicographic_less(m.index(),a);
    }
};

struct CoefficientLess {
    template<class M1, class M2> bool operator()(M1 const& m1, M2 const& m2) const {
        return m1.coefficient() < m2.coefficient();
    }
};

struct CoefficientIsZero {
    template<class M> bool operator()(M const& m) const {
        return decide(m.coefficient() == 0);
    }
};

template<class I, class X> class Expansion;

template<class X> class Expansion<MultiIndex,X> {
    using I=MultiIndex;
    static const SizeType DEFAULT_CAPACITY=16;
  public:
    X _zero_coefficient;
    SizeType _capacity;
    SizeType _size;
    SizeType _argument_size;
    DegreeType* _indices;
    X* _coefficients;
  public:
    typedef ExpansionIterator<I,X> Iterator;
    typedef ExpansionConstIterator<I,X> ConstIterator;
    typedef MultiIndex IndexType;
    typedef X CoefficientType;
    typedef ExpansionValue<I,X> ValueType;
    typedef ExpansionReference<I,X> Reference;
    typedef ExpansionConstReference<I,X> ConstReference;

    ~Expansion();
    explicit Expansion(SizeType as);
    explicit Expansion(SizeType as, X const& z, SizeType cap=DEFAULT_CAPACITY);
    Expansion(InitializerList<Pair<InitializerList<DegreeType>,X>> lst);
    template<class PR, EnableIf<IsConstructible<X,PR>> =dummy>
        explicit Expansion(SizeType as, PR pr, SizeType cap=DEFAULT_CAPACITY);
    template<class PR, EnableIf<IsConstructible<X,Dbl,PR>> =dummy>
        Expansion(InitializerList<Pair<InitializerList<DegreeType>,Dbl>> lst, PR prs);
    template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>> =dummy>
        explicit Expansion(Expansion<I,Y> const&, PRS... prs);
    Expansion(const Expansion<I,X>&);
    Expansion<I,X>& operator=(const Expansion<I,X>&);
    Expansion(Expansion<I,X>&&);
    Expansion<I,X>& operator=(Expansion<I,X>&&);
    Void swap(Expansion<I,X>&);

    Bool operator==(const Expansion<I,X>& other) const;
    Bool operator!=(const Expansion<I,X>& other) const;

    Bool empty() const;

    SizeType argument_size() const;
    SizeType number_of_terms() const;
    SizeType number_of_nonzeros() const; // DEPRECATED
    SizeType size() const;
    SizeType capacity() const;
    Void resize(SizeType sz);
    Void reserve(SizeType cap);

    CoefficientType const& zero_coefficient() const;
    const CoefficientType& operator[](const IndexType& a) const;
    ExpansionValueReference<I,X> operator[](const IndexType& a);

    Void set(const IndexType& a, const CoefficientType& c);
    CoefficientType& at(const IndexType& a);
    const CoefficientType& get(const IndexType& a) const;

    Iterator begin();
    Iterator end();
    ConstIterator begin() const;
    ConstIterator end() const;

    Void prepend(const IndexType& a, const CoefficientType& x);
    Void append(const IndexType& a, const CoefficientType& x);
    Void append_sum(const IndexType& a1, const IndexType& a2, const CoefficientType& x);

    Iterator find(const IndexType& a);
    ConstIterator find(const IndexType& a) const;

    Iterator insert(Iterator pos, const IndexType& a, const CoefficientType& x);
    Iterator erase(Iterator pos);

    Void clear();
    Void remove_zeros();
    Void combine_terms();
    Void check() const;

    //template<class CMP> Void sort(CMP cmp);
    Void index_sort(ReverseLexicographicLess cmp);
    Void index_sort(GradedLess cmp);
    Void sort(ReverseLexicographicIndexLess cmp);
    Void sort(GradedIndexLess cmp);

    Void reverse_lexicographic_sort();
    Void graded_sort();

    OutputStream& write(OutputStream& os) const;
    OutputStream& write(OutputStream& os, Array<String> const& vars) const;
  public:
    friend OutputStream& operator<<(OutputStream& os, Expansion<IndexType,CoefficientType> const& self) { return self.write(os); }
    friend Expansion<IndexType,CoefficientType> embed(SizeType as1, Expansion<IndexType,CoefficientType> const& e2, SizeType as3) {
        return e2._embed(as1,as3); }
  private:
    Expansion<IndexType,CoefficientType> _embed(SizeType as1, SizeType as3) const;
};

template<class I, class X> inline Bool same(const Expansion<I,X>& e1, const Expansion<I,X>& e2) {
    return e1==e2; }

template<class I, class X> inline OutputStream& operator<<(OutputStream& os, const Expansion<I,X>& e) {
    e.write(os); return os; }

template<class I, class X, class CMP> class SortedExpansion : public Expansion<I,X> {
public:
    typedef I IndexType;
    typedef X CoefficientType;
    typedef typename Expansion<I,X>::Iterator Iterator;
    typedef typename Expansion<I,X>::ConstIterator ConstIterator;
  public:
    using Expansion<I,X>::Expansion;
    SortedExpansion(Expansion<I,X> e);
    Void sort();
    Void insert(const IndexType& a, const CoefficientType& c);
    Void set(const IndexType& a, const CoefficientType& c);
    CoefficientType& at(const IndexType& a);
    CoefficientType const& get(const IndexType& a) const;
//    Iterator find(const IndexType& a);
//    ConstIterator find(const IndexType& a) const;
};


template<class I, class X> class ExpansionValue;

template<class X> class ExpansionValue<MultiIndex,X>
{
    using I = MultiIndex;
    I _a; X _c;
  public:
    typedef I IndexType;
    typedef X CoefficientType;
  public:
    ExpansionValue(SizeType as, const DegreeType* ip, const CoefficientType* cp) : _a(as,ip), _c(*cp) { }
  public:
    ExpansionValue(IndexType const& a, CoefficientType const& c) : _a(a), _c(c) { }
    IndexType& index() { return _a; }
    const IndexType& index() const { return _a; }
    CoefficientType& coefficient() { return _c; }
    const CoefficientType& coefficient() const { return _c; }
    friend OutputStream& operator<<(OutputStream& os, const ExpansionValue<I,X>& m) {
        return os << m.index()<<":" << m.coefficient(); }
};


template<class I, class X> class ExpansionConstReference;

template<class I, class X> class ExpansionReference;
template<class I, class X> Void swap(ExpansionReference<I,X>, ExpansionReference<I,X>);

template<class X> class ExpansionReference<MultiIndex,X>
{
    typedef MultiIndex I;
  public:
    typedef X CoefficientType;
  private:
    MultiIndexReference _a; CoefficientType* _cp;
  public:
    ExpansionReference(SizeType as, DegreeType* ip, CoefficientType* cp) : _a(as,ip), _cp(cp) { }
    ExpansionReference(MultiIndexReference& a, CoefficientType& c) : _a(a), _cp(&c) { }
    ReferenceType<MultiIndex> index() { return static_cast<MultiIndex&>(static_cast<MultiIndexData&>(_a)); }
    //MultiIndexReference& index() { return _a; }
    ConstReferenceType<MultiIndex> index() const { return _a; }
    ReferenceType<CoefficientType> coefficient() { return *_cp; }
    ConstReferenceType<CoefficientType> coefficient() const { return *_cp; }
    operator ExpansionValue<I,X>() const { return ExpansionValue<I,X>(_a,*_cp); }
    ExpansionReference<I,X>& operator=(const ExpansionValue<I,X>&);
    ExpansionReference<I,X>& operator=(const ExpansionReference<I,X>&);
    ExpansionReference<I,X>& operator=(const ExpansionConstReference<I,X>&);
};
template<class I, class X> Void swap(ExpansionReference<I,X>, ExpansionReference<I,X>);

template<class I, class X> class ExpansionConstReference;

template<class X> class ExpansionConstReference<MultiIndex,X>
{
    typedef MultiIndex I;
  public:
    typedef X CoefficientType;
  private:
    const MultiIndexReference _a; const CoefficientType* _cp;
  public:
    ExpansionConstReference(SizeType as, const DegreeType* ip, const CoefficientType* cp) : _a(as,const_cast<DegreeType*>(ip)), _cp(cp) { }
    ExpansionConstReference(const MultiIndexReference& a, const CoefficientType& c) : _a(a), _cp(&c) { }
    ConstReferenceType<MultiIndex> index() const { return _a; }
    ConstReferenceType<CoefficientType> coefficient() const { return *_cp; }
    ExpansionConstReference(const ExpansionValue<I,X>& other) : ExpansionConstReference(other._a._n,other._a._ip,other._cp) { }
    ExpansionConstReference(const ExpansionReference<I,X>& other) : ExpansionConstReference(other.index().size(),other.index().begin(),other._cp) { }
    operator ExpansionValue<I,X>() const { return ExpansionValue<I,X>(_a,*_cp); }
};

template<class I, class X> Void swap(ExpansionReference<I,X> m1, ExpansionReference<I,X> m2) {
    swap(m1.index(),m2.index());
    std::swap(m1.coefficient(),m2.coefficient());
}

template<class X> ExpansionReference<MultiIndex,X>& ExpansionReference<MultiIndex,X>::operator=(const ExpansionValue<MultiIndex,X>& other) {
    this->index()=other.index(); this->coefficient()=other.coefficient(); return *this; }

template<class X> ExpansionReference<MultiIndex,X>& ExpansionReference<MultiIndex,X>::operator=(const ExpansionReference<MultiIndex,X>& other) {
    this->index()=other.index(); this->coefficient()=other.coefficient(); return *this; }

template<class X> ExpansionReference<MultiIndex,X>& ExpansionReference<MultiIndex,X>::operator=(const ExpansionConstReference<MultiIndex,X>& other) {
    this->index()=other.index(); this->coefficient()=other.coefficient(); return *this; }

template<class I, class X> decltype(auto) operator==(const ExpansionConstReference<I,X>& m1, const ExpansionConstReference<I,X>& m2) {
    return m1.index()==m2.index() && (m1.coefficient()==m2.coefficient());
}

template<class I, class X> inline decltype(auto) operator!=(const ExpansionConstReference<I,X>& m1, const ExpansionConstReference<I,X>& m2) {
    return !(m1==m2);
}

template<class I, class X> OutputStream& operator<<(OutputStream& os, const ExpansionReference<I,X>& m) {
    return os << ExpansionValue<I,X>(m); }

template<class I, class X> OutputStream& operator<<(OutputStream& os, const ExpansionConstReference<I,X>& m) {
    return os << ExpansionValue<I,X>(m); }


template<class I, class X> class ExpansionIterator;

template<class X> class ExpansionIterator<MultiIndex,X>
    : public IteratorFacade<ExpansionIterator<MultiIndex,X>, ExpansionValue<MultiIndex,X>, RandomAccessTraversalTag, ExpansionReference<MultiIndex,X>>
{
    typedef MultiIndex I;
    typedef X CoefficientType;
    SizeType _as; DegreeType* _ip; CoefficientType* _cp;
    friend class ExpansionConstIterator<I,X>;
  public:
    ExpansionIterator(SizeType as, DegreeType* ip, CoefficientType* cp) : _as(as), _ip(ip), _cp(cp) { }
  public:
    template<class XX> Bool equal(const ExpansionIterator<I,XX>& other) const { return this->_cp==other._cp; }
    template<class XX> Bool equal(const ExpansionConstIterator<I,XX>& other) const { return this->_cp==other._cp; }
    Void advance(PointerDifferenceType k) { this->_ip+=k*(_as+1u); _cp+=k; }
    template<class XX> PointerDifferenceType distance_to(const ExpansionIterator<I,XX>& other) const { return other._cp - this->_cp; }
    ExpansionReference<I,X>& dereference() const { return reinterpret_cast<ExpansionReference<I,X>&>(const_cast<ExpansionIterator<I,X>&>(*this)); }
    ExpansionReference<I,X>* operator->() const { return reinterpret_cast<ExpansionReference<I,X>*>(const_cast<ExpansionIterator<I,X>*>(this)); }
    Void write(OutputStream& os) const { os << "{" << (void*)_ip << ":" << MultiIndex(_as,_ip) << ", " << _cp << ":" << *_cp << "}"; }
};
template<class I, class X> OutputStream& operator<<(OutputStream& os, const ExpansionIterator<I,X>& e) {
    e.write(os); return os; }

template<class I,class X> class ExpansionConstIterator;

template<class X> class ExpansionConstIterator<MultiIndex,X>
    : public IteratorFacade<ExpansionConstIterator<MultiIndex,X>, ExpansionValue<MultiIndex,X>, RandomAccessTraversalTag, ExpansionConstReference<MultiIndex,X>>
{
    typedef MultiIndex I;
    typedef X CoefficientType;
    friend class ExpansionIterator<I,X>;
    SizeType _as; DegreeType* _ip; CoefficientType* _cp;
  public:
    ExpansionConstIterator(SizeType as, DegreeType* ip, CoefficientType* cp) : _as(as), _ip(ip), _cp(cp) { }
    ExpansionConstIterator(const ExpansionIterator<I,X>& other) : ExpansionConstIterator(other._as,other._ip,other._cp) { }
  public:
    template<class XX> Bool equal(const ExpansionIterator<I,XX>& other) const { return this->_cp==other._cp; }
    template<class XX> Bool equal(const ExpansionConstIterator<I,XX>& other) const { return this->_cp==other._cp; }
    Void advance(PointerDifferenceType k) { this->_ip+=k*(_as+1u); _cp+=k; }
    template<class XX> PointerDifferenceType distance_to(const ExpansionConstIterator<I,XX>& other) const { return other._cp - this->_cp; }
    ExpansionConstReference<I,X> dereference() const { return ExpansionConstReference<I,X>(_as,_ip,_cp); }
    ExpansionConstReference<I,X>* operator->() const { return reinterpret_cast<ExpansionConstReference<I,X>*>(const_cast<ExpansionConstIterator<I,X>*>(this)); }
    Void write(OutputStream& os) const { os << "{" << (void*)_ip << "," << _cp << "}"; }
};
template<class I, class X> OutputStream& operator<<(OutputStream& os, const ExpansionConstIterator<I,X>& e) {
    e.write(os); return os; }

template<class I, class X> class ExpansionValueReference {
    Expansion<I,X>& _e; I const& _a;
  public:
    ExpansionValueReference(Expansion<I,X>& e, const I& a) : _e(e), _a(a) { }
    operator const X& () const { return _e.get(_a); }
    //operator X& () { return _e.at(_a); }
    ExpansionValueReference<I,X>& operator=(const X& x) { _e.at(_a)=x; return *this; }
    ExpansionValueReference<I,X>& operator+=(X const& c) { _e.at(_a)+=c; return *this; }
    ExpansionValueReference<I,X>& operator-=(X const& c) { _e.at(_a)-=c; return *this; }
    ExpansionValueReference<I,X>& operator*=(X const& c) { _e.at(_a)*=c; return *this; }
    ExpansionValueReference<I,X>& operator/=(X const& c) { _e.at(_a)/=c; return *this; }
};


template<class X> template<class PR, EnableIf<IsConstructible<X,PR>>>
Expansion<MultiIndex,X>::Expansion(SizeType as, PR pr, SizeType cap)
    : Expansion(as,X(pr),cap)
{
}

template<class X> template<class PR, EnableIf<IsConstructible<X,Dbl,PR>>>
Expansion<MultiIndex,X>::Expansion(InitializerList<Pair<InitializerList<DegreeType>,Dbl>> lst, PR pr)
    : Expansion( ((ARIADNE_PRECONDITION(lst.size()!=0)),lst.begin()->first.size()), X(pr) )
{
    MultiIndex a;
    X x;
    for(auto iter=lst.begin();
        iter!=lst.end(); ++iter)
    {
        a=iter->first;
        x=X(iter->second,pr);
        if(decide(x!=0)) { this->append(a,x); }
    }
}

template<class X> template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>>>
Expansion<MultiIndex,X>::Expansion(Expansion<MultiIndex,Y> const& other, PRS... prs)
    : Expansion(other.argument_size(),X(prs...))
{
    MultiIndex a;
    X x;
    for(auto iter=other.begin();
        iter!=other.end(); ++iter)
    {
        a=iter->index();
        x=X(iter->coefficient(),prs...);
        if(decide(x!=0)) { this->append(a,x); }
    }
}


} // namespace Ariadne

#endif
