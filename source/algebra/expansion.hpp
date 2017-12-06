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

struct GradedIndexLess;
struct LexicographicIndexLess;
struct ReverseLexicographicIndexLess;
struct CoefficientLess;
struct CoefficientIsZero;


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
    Bool same_as(const Expansion<I,X>& e) const;

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
    friend OutputStream& operator<<(OutputStream& os, Expansion<I,X> const& self) { return self.write(os); }
    friend Bool same(Expansion<I,X> const& e1, Expansion<I,X> const& e2) { return e1.same_as(e2); }
    friend Expansion<IndexType,CoefficientType> embed(SizeType as1, Expansion<IndexType,CoefficientType> const& e2, SizeType as3) {
        return e2._embed(as1,as3); }
  private:
    Expansion<IndexType,CoefficientType> _embed(SizeType as1, SizeType as3) const;
};

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


} // namespace Ariadne

#endif
