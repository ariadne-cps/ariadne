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

template<class X> class ExpansionValue;
template<class X> class ExpansionReference;
template<class X> class ExpansionConstReference;
template<class X> class ExpansionIterator;
template<class X> class ExpansionConstIterator;
template<class X> class ExpansionValueReference;

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


template<class X> class Expansion {
    static const SizeType DEFAULT_CAPACITY=16;
  public:
    X _zero_coefficient;
    SizeType _capacity;
    SizeType _size;
    SizeType _argument_size;
    DegreeType* _indices;
    X* _coefficients;
  public:
    typedef ExpansionIterator<X> Iterator;
    typedef ExpansionConstIterator<X> ConstIterator;
    typedef MultiIndex MultiIndexType;
    typedef X CoefficientType;
    typedef ExpansionValue<X> ValueType;
    typedef ExpansionReference<X> Reference;
    typedef ExpansionConstReference<X> ConstReference;

    ~Expansion();
    explicit Expansion(SizeType as);
    explicit Expansion(SizeType as, X const& z, SizeType cap=DEFAULT_CAPACITY);
    Expansion(InitializerList<Pair<InitializerList<DegreeType>,X>> lst);
    template<class PR, EnableIf<IsConstructible<X,PR>> =dummy>
        explicit Expansion(SizeType as, PR pr, SizeType cap=DEFAULT_CAPACITY);
    template<class PR, EnableIf<IsConstructible<X,Dbl,PR>> =dummy>
        Expansion(InitializerList<Pair<InitializerList<DegreeType>,Dbl>> lst, PR prs);
    template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>> =dummy>
        explicit Expansion(Expansion<Y> const&, PRS... prs);
    Expansion(const Expansion<X>&);
    Expansion<X>& operator=(const Expansion<X>&);
    Expansion(Expansion<X>&&);
    Expansion<X>& operator=(Expansion<X>&&);
    Void swap(Expansion<X>&);

    Bool operator==(const Expansion<X>& other) const;
    Bool operator!=(const Expansion<X>& other) const;

    Bool empty() const;

    SizeType argument_size() const;
    SizeType number_of_terms() const;
    SizeType number_of_nonzeros() const; // DEPRECATED
    SizeType size() const;
    SizeType capacity() const;
    Void resize(SizeType sz);
    Void reserve(SizeType cap);

    CoefficientType const& zero_coefficient() const;
    const CoefficientType& operator[](const MultiIndex& a) const;
    ExpansionValueReference<X> operator[](const MultiIndex& a);

    Void set(const MultiIndex& a, const X& c);
    CoefficientType& at(const MultiIndex& a);
    const CoefficientType& get(const MultiIndex& a) const;

    Iterator begin();
    Iterator end();
    ConstIterator begin() const;
    ConstIterator end() const;

    Void prepend(const MultiIndex& a, const X& x);
    Void append(const MultiIndex& a, const X& x);
    Void append_sum(const MultiIndex& a1, const MultiIndex& a2, const X& x);

    Iterator find(const MultiIndex& a);
    ConstIterator find(const MultiIndex& a) const;

    Iterator insert(Iterator pos, const MultiIndex& a, const X& x);
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
    friend OutputStream& operator<<(OutputStream& os, Expansion<X> const& self) { return self.write(os); }
    friend Expansion<X> embed(SizeType as1, Expansion<X> const& e2, SizeType as3) {
        return e2._embed(as1,as3); }
  private:
    Expansion<X> _embed(SizeType as1, SizeType as3) const;
};

template<class X> inline Bool same(const Expansion<X>& e1, const Expansion<X>& e2) {
    return e1==e2; }

template<class X> inline OutputStream& operator<<(OutputStream& os, const Expansion<X>& e) {
    e.write(os); return os; }

template<class X, class CMP> class SortedExpansion : public Expansion<X> {
public:
    typedef X CoefficientType;
    typedef typename Expansion<X>::Iterator Iterator;
    typedef typename Expansion<X>::ConstIterator ConstIterator;
  public:
    using Expansion<X>::Expansion;
    SortedExpansion(Expansion<X> e);
    Void sort();
    Void insert(const MultiIndex& a, const X& c);
    Void set(const MultiIndex& a, const X& c);
    CoefficientType& at(const MultiIndex& a);
    CoefficientType const& get(const MultiIndex& a) const;
//    Iterator find(const MultiIndex& a);
//    ConstIterator find(const MultiIndex& a) const;
};



template<class X> class ExpansionValue
{
    typedef X CoefficientType;
    MultiIndex _a; CoefficientType _c;
  public:
    ExpansionValue(SizeType as, const DegreeType* ip, const CoefficientType* cp) : _a(as,ip), _c(*cp) { }
  public:
    ExpansionValue(MultiIndex const& a, CoefficientType const& c) : _a(a), _c(c) { }
    MultiIndex& index() { return _a; }
    const MultiIndex& index() const { return _a; }
    CoefficientType& coefficient() { return _c; }
    const CoefficientType& coefficient() const { return _c; }
        MultiIndex& key() { return _a; }
        const MultiIndex& key() const { return _a; }
        CoefficientType& data() { return _c; }
        const CoefficientType& data() const { return _c; }
};
template<class X> OutputStream& operator<<(OutputStream& os, const ExpansionValue<X>& m) {
    return os << m.index()<<":" << m.coefficient(); }

template<class X> class ExpansionConstReference;

template<class X> class ExpansionReference;
template<class X> Void swap(ExpansionReference<X>, ExpansionReference<X>);

template<class X> class ExpansionReference
{
  public:
    typedef X CoefficientType; typedef X& CoefficientReference;
  private:
    MultiIndexReference _a; CoefficientType* _cp;
  public:
    ExpansionReference(SizeType as, DegreeType* ip, CoefficientType* cp) : _a(as,ip), _cp(cp) { }
    ExpansionReference(MultiIndexReference& a, CoefficientType& c) : _a(a), _cp(&c) { }
    MultiIndexReference& index() { return _a; }
    const MultiIndex& index() const { return _a; }
    CoefficientReference coefficient() { return *_cp; }
    const CoefficientType& coefficient() const { return *_cp; }
        MultiIndexReference& key() { return _a; }
        const MultiIndex& key() const { return _a; }
        CoefficientReference data() { return *_cp; }
        const CoefficientType& data() const { return *_cp; }
    operator ExpansionValue<X>() const { return ExpansionValue<X>(_a,*_cp); }
    ExpansionReference<X>& operator=(const ExpansionValue<X>&);
    ExpansionReference<X>& operator=(const ExpansionReference<X>&);
    ExpansionReference<X>& operator=(const ExpansionConstReference<X>&);
};
template<class X> Void swap(ExpansionReference<X>, ExpansionReference<X>);

template<class X> class ExpansionConstReference
{
  public:
    typedef X CoefficientType;
  private:
    const MultiIndexReference _a; const CoefficientType* _cp;
  public:
    ExpansionConstReference(SizeType as, const DegreeType* ip, const CoefficientType* cp) : _a(as,const_cast<DegreeType*>(ip)), _cp(cp) { }
    ExpansionConstReference(const MultiIndexReference& a, const CoefficientType& c) : _a(a), _cp(&c) { }
    const MultiIndex& index() const { return _a; }
    const CoefficientType& coefficient() const { return *_cp; }
        const MultiIndex& key() const { return _a; }
        const CoefficientType& data() const { return *_cp; }
    ExpansionConstReference(const ExpansionValue<X>& other) : ExpansionConstReference(other._a._n,other._a._ip,other._cp) { }
    ExpansionConstReference(const ExpansionReference<X>& other) : ExpansionConstReference(other.index().size(),other.index().begin(),other._cp) { }
    operator ExpansionValue<X>() const { return ExpansionValue<X>(_a,*_cp); }
};

template<class X> Void swap(ExpansionReference<X> m1, ExpansionReference<X> m2) {
    swap(m1.index(),m2.index());
    std::swap(m1.coefficient(),m2.coefficient());
}

template<class X> ExpansionReference<X>& ExpansionReference<X>::operator=(const ExpansionValue<X>& other) {
    this->index()=other.index(); this->coefficient()=other.coefficient(); return *this; }

template<class X> ExpansionReference<X>& ExpansionReference<X>::operator=(const ExpansionReference<X>& other) {
    this->index()=other.index(); this->coefficient()=other.coefficient(); return *this; }

template<class X> ExpansionReference<X>& ExpansionReference<X>::operator=(const ExpansionConstReference<X>& other) {
    this->index()=other.index(); this->coefficient()=other.coefficient(); return *this; }

template<class X> decltype(auto) operator==(const ExpansionConstReference<X>& m1, const ExpansionConstReference<X>& m2) {
    return m1.index()==m2.index() && (m1.coefficient()==m2.coefficient());
}

template<class X> inline decltype(auto) operator!=(const ExpansionConstReference<X>& m1, const ExpansionConstReference<X>& m2) {
    return !(m1==m2);
}

template<class X> OutputStream& operator<<(OutputStream& os, const ExpansionReference<X>& m) {
    return os << ExpansionValue<X>(m); }

template<class X> OutputStream& operator<<(OutputStream& os, const ExpansionConstReference<X>& m) {
    return os << ExpansionValue<X>(m); }


template<class X> class ExpansionIterator
    : public IteratorFacade<ExpansionIterator<X>, ExpansionValue<X>, RandomAccessTraversalTag, ExpansionReference<X>>
{
    typedef X CoefficientType;
    SizeType _as; DegreeType* _ip; CoefficientType* _cp;
    friend class ExpansionConstIterator<X>;
  public:
    ExpansionIterator(SizeType as, DegreeType* ip, CoefficientType* cp) : _as(as), _ip(ip), _cp(cp) { }
  public:
    template<class XX> Bool equal(const ExpansionIterator<XX>& other) const { return this->_cp==other._cp; }
    template<class XX> Bool equal(const ExpansionConstIterator<XX>& other) const { return this->_cp==other._cp; }
    Void advance(PointerDifferenceType k) { this->_ip+=k*(_as+1u); _cp+=k; }
    template<class XX> PointerDifferenceType distance_to(const ExpansionIterator<XX>& other) const { return other._cp - this->_cp; }
    ExpansionReference<X>& dereference() { return reinterpret_cast<ExpansionReference<X>&>(*this); }
    ExpansionReference<X>* operator->() { return reinterpret_cast<ExpansionReference<X>*>(this); }
    Void write(OutputStream& os) const { os << "{" << (void*)_ip << ":" << MultiIndex(_as,_ip) << ", " << _cp << ":" << *_cp << "}"; }
};
template<class X> OutputStream& operator<<(OutputStream& os, const ExpansionIterator<X>& e) {
    e.write(os); return os; }

template<class X> class ExpansionConstIterator
    : public IteratorFacade<ExpansionConstIterator<X>, ExpansionValue<X>, RandomAccessTraversalTag, ExpansionConstReference<X>>
{
    typedef X CoefficientType;
    friend class ExpansionIterator<X>;
    SizeType _as; DegreeType* _ip; CoefficientType* _cp;
  public:
    ExpansionConstIterator(SizeType as, DegreeType* ip, CoefficientType* cp) : _as(as), _ip(ip), _cp(cp) { }
    ExpansionConstIterator(const ExpansionIterator<X>& other) : ExpansionConstIterator(other._as,other._ip,other._cp) { }
  public:
    template<class XX> Bool equal(const ExpansionIterator<XX>& other) const { return this->_cp==other._cp; }
    template<class XX> Bool equal(const ExpansionConstIterator<XX>& other) const { return this->_cp==other._cp; }
    Void advance(PointerDifferenceType k) { this->_ip+=k*(_as+1u); _cp+=k; }
    template<class XX> PointerDifferenceType distance_to(const ExpansionConstIterator<XX>& other) const { return other._cp - this->_cp; }
    ExpansionConstReference<X> dereference() { return ExpansionConstReference<X>(_as,_ip,_cp); }
    ExpansionConstReference<X>* operator->() { return reinterpret_cast<ExpansionConstReference<X>*>(this); }
    Void write(OutputStream& os) const { os << "{" << (void*)_ip << "," << _cp << "}"; }
};
template<class X> OutputStream& operator<<(OutputStream& os, const ExpansionConstIterator<X>& e) {
    e.write(os); return os; }

template<class X> class ExpansionValueReference {
    Expansion<X>& _e; MultiIndex const& _a;
  public:
    ExpansionValueReference(Expansion<X>& e, const MultiIndex& a) : _e(e), _a(a) { }
    operator const X& () const { return _e.get(_a); }
    //operator X& () { return _e.at(_a); }
    ExpansionValueReference<X>& operator=(const X& x) { _e.at(_a)=x; return *this; }
    ExpansionValueReference<X>& operator+=(X const& c) { _e.at(_a)+=c; return *this; }
    ExpansionValueReference<X>& operator-=(X const& c) { _e.at(_a)-=c; return *this; }
    ExpansionValueReference<X>& operator*=(X const& c) { _e.at(_a)*=c; return *this; }
    ExpansionValueReference<X>& operator/=(X const& c) { _e.at(_a)/=c; return *this; }
};


template<class X> template<class PR, EnableIf<IsConstructible<X,PR>>>
Expansion<X>::Expansion(SizeType as, PR pr, SizeType cap)
    : Expansion(as,X(pr),cap)
{
}

template<class X> template<class PR, EnableIf<IsConstructible<X,Dbl,PR>>>
Expansion<X>::Expansion(InitializerList<Pair<InitializerList<DegreeType>,Dbl>> lst, PR pr)
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
Expansion<X>::Expansion(Expansion<Y> const& other, PRS... prs)
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
