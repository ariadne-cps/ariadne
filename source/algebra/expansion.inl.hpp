/***************************************************************************
 *            algebra/expansion.inl.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

#ifndef ARIADNE_EXPANSION_INL_HPP
#define ARIADNE_EXPANSION_INL_HPP

#include "multi_index.hpp"
#include "expansion.hpp"

#include "../utility/typedefs.hpp"
#include "../utility/iterator.hpp"
#include "../utility/macros.hpp"

namespace Ariadne {

/************ Expansion ******************************************************/


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


template<class I, class X> class ExpansionValue
{
    I _a; X _c;
  public:
    typedef I IndexType;
    typedef X CoefficientType;
  public:
    ExpansionValue(const IndexType& a, const CoefficientType& c) : _a(a), _c(c) { }
    IndexType& index() { return _a; }
    const IndexType& index() const { return _a; }
    CoefficientType& coefficient() { return _c; }
    const CoefficientType& coefficient() const { return _c; }
    decltype(auto) operator==(ExpansionValue<I,X> const& other) const { return this->_a==other._a && this->_c == other._c; }
    decltype(auto) operator!=(ExpansionValue<I,X> const& other) const { return !(*this==other); }
    friend OutputStream& operator<<(OutputStream& os, const ExpansionValue<I,X>& m) {
        return os << m.index()<<":" << m.coefficient(); }
};


template<class I, class X> class ExpansionConstReference;

template<class I, class X> class ExpansionReference;
template<class I, class X> Void swap(ExpansionReference<I,X>, ExpansionReference<I,X>);

template<class I, class X> class ExpansionReference
{
  public:
    typedef I IndexType; typedef X CoefficientType;
    typedef typename Expansion<I,X>::IndexReference IndexReference;
    typedef typename Expansion<I,X>::CoefficientReference CoefficientReference;
    typedef typename Expansion<I,X>::IndexConstReference IndexConstReference;
    typedef typename Expansion<I,X>::CoefficientConstReference CoefficientConstReference;
  private:
    IndexReference _a; CoefficientReference _c;
  public:
    ExpansionReference(IndexReference a, CoefficientReference c) : _a(a), _c(c) { }
    IndexReference index() { return _a; }
    IndexConstReference index() const { return _a; }
    CoefficientReference coefficient() { return _c; }
    CoefficientConstReference coefficient() const { return _c; }
    operator ExpansionValue<I,X>() const { return ExpansionValue<I,X>(_a,_c); }
    ExpansionReference<I,X>& operator=(const ExpansionValue<I,X>&);
    ExpansionReference<I,X>& operator=(const ExpansionReference<I,X>&);
    ExpansionReference<I,X>& operator=(const ExpansionConstReference<I,X>&);
    friend Void swap(ExpansionReference<I,X> m1, ExpansionReference<I,X> m2) {
        swap(m1.index(),m2.index()); std::swap(m1.coefficient(),m2.coefficient()); }
};

template<class I, class X> class ExpansionConstReference
{
  public:
    typedef I IndexType; typedef X CoefficientType;
    typedef typename Expansion<I,X>::IndexConstReference IndexConstReference;
    typedef typename Expansion<I,X>::CoefficientConstReference CoefficientConstReference;
  private:
    IndexConstReference _a; CoefficientConstReference _c;
  public:
    ExpansionConstReference(IndexConstReference a, CoefficientConstReference c) : _a(a), _c(c) { }
    IndexConstReference index() const { return _a; }
    CoefficientConstReference coefficient() const { return _c; }
    ExpansionConstReference(const ExpansionValue<I,X>& other) : ExpansionConstReference(other.index(),other.coefficient()) { }
    ExpansionConstReference(const ExpansionReference<I,X>& other) : ExpansionConstReference(other.index(),other.coefficient()) { }
    operator ExpansionValue<I,X>() const { return ExpansionValue<I,X>(_a,_c); }
};



template<class I, class X> ExpansionReference<I,X>& ExpansionReference<I,X>::operator=(const ExpansionValue<I,X>& other) {
    this->index()=other.index(); this->coefficient()=other.coefficient(); return *this; }

template<class I, class X> ExpansionReference<I,X>& ExpansionReference<I,X>::operator=(const ExpansionReference<I,X>& other) {
    this->index()=other.index(); this->coefficient()=other.coefficient(); return *this; }

template<class I, class X> ExpansionReference<I,X>& ExpansionReference<I,X>::operator=(const ExpansionConstReference<I,X>& other) {
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


template<class I, class X> class ExpansionPointer
{
    ExpansionReference<I,X> _r;
  public:
    ExpansionPointer(typename UniformList<I>::Pointer ap, typename UniformList<X>::Pointer cp) : _r(*ap,*cp) { }
    ExpansionPointer(ExpansionReference<I,X>* p) : _r(*p) { }
    ExpansionReference<I,X>& operator*() { return _r; }
    ExpansionReference<I,X>* operator->() { return &_r; }
};

template<class I, class X> class ExpansionConstPointer
{
    ExpansionConstReference<I,X> _r;
  public:
    ExpansionConstPointer(typename UniformList<I>::ConstPointer ap, typename UniformList<X>::ConstPointer cp) : _r(*ap,*cp) { }
    ExpansionConstPointer(ExpansionConstReference<I,X>* p) : _r(*p) { }
    ExpansionConstReference<I,X>& operator*() { return _r; }
    ExpansionConstReference<I,X>* operator->() { return &_r; }
};


template<class I, class X> class ExpansionIterator
    : public IteratorFacade<ExpansionIterator<I,X>, ExpansionValue<I,X>, RandomAccessTraversalTag, ExpansionReference<I,X>>
{
    typedef I IndexType;
    typedef X CoefficientType;
    typedef typename UniformList<I>::Iterator IndexIterator;
    typedef typename UniformList<X>::Iterator CoefficientIterator;
    friend class ExpansionConstIterator<I,X>;
  private: public:
    IndexIterator _ap; CoefficientIterator _cp;
  public:
    ExpansionIterator(IndexIterator ap, CoefficientIterator cp) : _ap(ap), _cp(cp) { }
    template<class XX> Bool equal(const ExpansionIterator<I,XX>& other) const { return this->_cp==other._cp; }
    template<class XX> Bool equal(const ExpansionConstIterator<I,XX>& other) const { return this->_cp==other._cp; }
    Void advance(PointerDifferenceType k) { _ap+=k; _cp+=k; }
    template<class XX> PointerDifferenceType distance_to(const ExpansionIterator<I,XX>& other) const { return other._cp - this->_cp; }
    ExpansionReference<I,X> dereference() const {
        auto ar=*(this->_ap);
        return ExpansionReference<I,X>(*_ap,*_cp); }
    ExpansionPointer<I,X> operator->() { return ExpansionPointer<I,X>(_ap.operator->(),_cp.operator->()); }
    friend OutputStream& operator<<(OutputStream& os, const ExpansionIterator<I,X>& e) {
        return os << "{" << e._ap.operator->() << "," << e._cp.operator->() << "}"; }
};

template<class I, class X> class ExpansionConstIterator
    : public IteratorFacade<ExpansionConstIterator<I,X>, ExpansionValue<I,X>, RandomAccessTraversalTag, ExpansionConstReference<I,X>>
{
    typedef I IndexType;
    typedef X CoefficientType;
    typedef typename UniformList<I>::ConstIterator IndexConstIterator;
    typedef typename UniformList<X>::ConstIterator CoefficientConstIterator;
    friend class ExpansionIterator<I,X>;
    IndexConstIterator _ap; CoefficientConstIterator _cp;
  public:
    ExpansionConstIterator(IndexConstIterator ap, CoefficientConstIterator cp) : _ap(ap), _cp(cp) { }
    ExpansionConstIterator(const ExpansionIterator<I,X>& other) : ExpansionConstIterator(other._ap,other._cp) { }
    template<class XX> Bool equal(const ExpansionIterator<I,XX>& other) const { return this->_cp==other._cp; }
    template<class XX> Bool equal(const ExpansionConstIterator<I,XX>& other) const { return this->_cp==other._cp; }
    Void advance(PointerDifferenceType k) { _ap+=k; _cp+=k; }
    template<class XX> PointerDifferenceType distance_to(const ExpansionConstIterator<I,XX>& other) const { return other._cp - this->_cp; }
    ExpansionConstReference<I,X> dereference() const { return ExpansionConstReference<I,X>(*_ap,*_cp); }
    ExpansionConstPointer<I,X> operator->() { return ExpansionConstPointer<I,X>(_ap.operator->(),_cp.operator->()); }
    friend OutputStream& operator<<(OutputStream& os, const ExpansionConstIterator<I,X>& e) {
        return os << "{" << e._ap.operator->() << "," << e._cp.operator->() << "}"; }
};

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


template<class I, class X> template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>>>
Expansion<I,X>::Expansion(Expansion<I,Y> const& other, PRS... prs)
    : Expansion(other.argument_size(),X(prs...))
{
    MultiIndex a(other.argument_size());
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
