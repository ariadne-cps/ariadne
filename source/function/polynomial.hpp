/***************************************************************************
 *            polynomial.hpp
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

/*! \file polynomial.hpp
 *  \brief Base class for polynomial rings.
 */

#ifndef ARIADNE_POLYNOMIAL_HPP
#define ARIADNE_POLYNOMIAL_HPP

#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "../algebra/multi_index.hpp"
#include "../algebra/expansion.hpp"
#include "../algebra/operations.hpp"
#include "../algebra/differential.hpp"


namespace Ariadne {

template<class T> class Array;
template<class X> class Algebra;

//! \brief A monomial with index \a I and coefficients of some type \a X.
template<class X>
class Monomial
    : public ExpansionValue<MultiIndex,X>
{
    typedef MultiIndex I;
  public:
    Monomial(const MultiIndex& a, const X& x) : ExpansionValue<I,X>(a,x) { }
    Monomial(const ExpansionValue<I,X>& v) : ExpansionValue<I,X>(v) { }
};

//! \ingroup FunctionModule
//! \brief A polynomial with coefficients of some type \a X.
template<class X>
class Polynomial
    : public DispatchAlgebraOperations<Polynomial<X>,X>
{
    template<class XX> friend class Polynomial;
    friend struct AlgebraOperations<Polynomial<X>,X>;
  public:
    typedef typename Expansion<MultiIndex,X>::ValueType ValueType;
    typedef typename Expansion<MultiIndex,X>::Reference Reference;
    typedef typename Expansion<MultiIndex,X>::ConstReference ConstReference;
    typedef typename Expansion<MultiIndex,X>::Iterator Iterator;
    typedef typename Expansion<MultiIndex,X>::ConstIterator ConstIterator;

    typedef typename Expansion<MultiIndex,X>::IndexReference IndexReference;
    typedef typename Expansion<MultiIndex,X>::IndexConstReference IndexConstReference;
    typedef typename Expansion<MultiIndex,X>::CoefficientReference CoefficientReference;
    typedef typename Expansion<MultiIndex,X>::CoefficientConstReference CoefficientConstReference;

    typedef typename X::Paradigm Paradigm;
    typedef typename X::NumericType NumericType;
    typedef Polynomial<X> SelfType;
    typedef ReverseLexicographicIndexLess ComparisonType;
    typedef ReverseLexicographicLess IndexComparisonType;
  public:
    //@{
    //! \name Constructors

    //! \brief The zero polynomial in \a as variables.
    explicit Polynomial(SizeType as=0u);
    //! \brief Copy/conversion constructor.
    template<class XX> Polynomial(const Polynomial<XX>& p);
    //! \brief Copy/conversion constructor.
    template<class XX> explicit Polynomial(const Expansion<MultiIndex,XX>& e);
    //! \brief A sparse polynomial with coefficients given by an initializer list of indices and coefficients.
    Polynomial(InitializerList<Pair<InitializerList<DegreeType>,X>> lst);
    //@}

    //! \brief Create the null polynomial in the same number of variables.
    Polynomial<X> create_zero() const;

    //! \brief Create a constant polynomial in \a as variables with value \a c.
    static Polynomial<X> constant(SizeType as, const X& c);
    //! \brief Create a polynomial in \a as variables which returns the value of the \a j<sup>th</sup> variable.
    static Polynomial<X> coordinate(SizeType as, SizeType j);
    static Polynomial<X> variable(SizeType as, SizeType j);
    //! \brief Create an Array of polynomials in \a as variables,
    //! the i<sup>th</sup> of  which returns the value of the i<sup>th</sup> variable.
    static Vector<Polynomial<X>> coordinates(SizeType as);
    static Vector<Polynomial<X>> variables(SizeType as);

    //! \brief Set equal to a constant.
    Polynomial<X>& operator=(const X& x);
    //@{
    //! \name Comparisons

    //! \brief Equality operator.
    template<class XX> Bool operator==(const Polynomial<XX>& p) const;
    //! \brief Inequality operator.
    template<class XX> Bool operator!=(const Polynomial<XX>& p) const;
    //@}

    //@{
    //! \name Data access

    //! \brief The number of variables in the argument of the polynomial.
    SizeType argument_size() const;
    //! \brief The number of structural nonzero terms.
    SizeType number_of_terms() const;
    //! \brief The order of the highest term.
    DegreeType degree() const;
    //! \brief The value of the polynomial at zero.
    const X& value() const;
    //! \brief A reference to the coefficient of the term in \f$x^{a_1}\cdots x^{a_n}\f$.
    X& operator[](const MultiIndex& a);
    //! \brief A constant referent to the coefficient of the term in \f$x^{a_1}\cdots x^{a_n}\f$.
    const X& operator[](const MultiIndex& a) const;
    //! \brief A constant reference to the raw data expansion.
    const Expansion<MultiIndex,X>& expansion() const;
    //! \brief A reference to the raw data expansion.
    Expansion<MultiIndex,X>& expansion();
    //@}

    //@{
    //! \name Iterators

    //! \brief An Iterator to the beginning of the list of terms.
    Iterator begin();
    //! \brief An Iterator to the end of the list of terms..
    Iterator end();
    //! \brief An Iterator to the term in \f$x^a\f$. Returns \c end() if there is no term in \a a.
    Iterator find(const MultiIndex& a);
    //! \brief A constant Iterator to the beginning of the list of terms.
    ConstIterator begin() const;
    //! \brief A constant Iterator to the end of the list of terms.
    ConstIterator end() const;
    //! \brief An Iterator to the term in \f$x^a\f$. Returns \c end() if there is no term in \a a.
    ConstIterator find(const MultiIndex& a) const;
    //@}


    //@{
    //! \name Modifying operations

    //! \brief Insert the term \f$c x^{a_1}\f$ into a sorted list of terms.
    Void insert(const MultiIndex& a, const X& c);
    //! \brief Reserve space for a total of \a n terms.
    Void reserve(SizeType n);
    //! \brief Remove the term pointed to by \a iter. May be expensive if the term is near the beginning of the list of terms.
    Void erase(Iterator iter);
    //! \brief Set the polynomial to zero.
    Void clear();
    //! \brief Remove all zero terms from the expansion, and order the expansion reverse lexicographically by term.
    Void cleanup();
    //@}

    //@{
    //! \name Evaluation

    //! Evaluate on a vector of algebra elements.
    template<class A> A operator() (Vector<A> const&) const;
    //@}

    //@{
    //! \name Modifying operators

    //! \brief Truncate to degree \a d.
    Polynomial<X>& truncate(DegreeType d);
    //! \brief Differentiate with respect to the \a j<sup>th</sup> variable.
    Polynomial<X>& differentiate(SizeType j);
    //! \brief Antidifferentiate (integrate) with respect to the \a j<sup>th</sup> variable.
    Polynomial<X>& antidifferentiate(SizeType j);
    //@}

    //@{
    //! \name Related operations
    friend Polynomial<X>& operator*=(Polynomial<X>& p, const Monomial<X>& m) { return Polynomial<X>::_imul(p,m); }

    template<class XX, class A> friend A evaluate(const Polynomial<XX>& p, const Vector<A>& v);
    template<class XX> friend Polynomial<XX> compose(const Polynomial<XX>& p, const Vector<Polynomial<XX>>& q);
    template<class XX> friend Polynomial<XX> derivative(Polynomial<XX> dx, SizeType k);
    template<class XX> friend Polynomial<XX> antiderivative(Polynomial<XX> dx, SizeType k);
    template<class XX> friend Polynomial<XX> truncate(Polynomial<XX> dx, DegreeType deg);
    //@}

    Void check() const;
    static Polynomial<X> _compose(const Polynomial<X>& p, const Vector<Polynomial<X>>& q);
    static X _evaluate(const Polynomial<X>& p, const Vector<X>& vx);
    static Algebra<X> _evaluate(const Polynomial<X>& p, const Vector<Algebra<X>>& va);
    static Polynomial<X> _partial_evaluate(const Polynomial<X>& p, SizeType k, const X& c);
    OutputStream& _write(OutputStream& os) const;
    OutputStream& _write(OutputStream& os, List<String> const& names) const;
  private:
    Void _append(const MultiIndex& a, const X& c);
    Iterator _unique_key();
  private:
    SortedExpansion<MultiIndex,X,ReverseLexicographicIndexLess> _expansion;
  private: // FIXME: Put these concrete-generic operations in proper place
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend Polynomial<X> operator+(Polynomial<X> p, const Y& c) {
            X xc=p.value(); xc=c; return p+xc; }
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend Polynomial<X> operator-(Polynomial<X> p, const Y& c) {
            X xc=p.value(); xc=c; return p-xc; }
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend Polynomial<X> operator*(Polynomial<X> p, const Y& c) {
            X xc=p.value(); xc=c; return p*xc; }
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend Polynomial<X> operator/(Polynomial<X> p, const Y& c) {
            X xc=p.value(); xc=c; return p/xc; }

};

template<class X> struct AlgebraOperations<Polynomial<X>> {
  public:
    static Polynomial<X> apply(Pos, const Polynomial<X>& p);
    static Polynomial<X> apply(Neg, const Polynomial<X>& p);
    static Polynomial<X> apply(Add, const Polynomial<X>& p1, const Polynomial<X>& p2);
    static Polynomial<X> apply(Sub, const Polynomial<X>& p1, const Polynomial<X>& p2);
    static Polynomial<X> apply(Mul, const Polynomial<X>& p1, const Polynomial<X>& p2);
    static Polynomial<X> apply(Add, Polynomial<X> p, const X& c);
    static Polynomial<X> apply(Mul, Polynomial<X> p, const X& c);
    static Polynomial<X> apply(Mul, Polynomial<X> p, const Monomial<X>& m);
    static Polynomial<X>& iapply(Add, Polynomial<X>& p, const X& c);
    static Polynomial<X>& iapply(Mul, Polynomial<X>& p, const X& c);
    static Polynomial<X>& iapply(Mul, Polynomial<X>& p, const Monomial<X>& m);

};


template<class X> template<class XX> Polynomial<X>::Polynomial(const Polynomial<XX>& p)
    : _expansion(p._expansion) { }

template<class X> template<class XX> Polynomial<X>::Polynomial(const Expansion<MultiIndex,XX>& e)
    : _expansion(e) { this->cleanup(); }

template<class X> template<class XX> Bool Polynomial<X>::operator==(const Polynomial<XX>& p) const {
    const_cast<Polynomial<X>*>(this)->cleanup();
    const_cast<Polynomial<XX>&>(p).cleanup();
    return this->_expansion==p._expansion;
}

template<class X> template<class XX> Bool Polynomial<X>::operator!=(const Polynomial<XX>& p) const {
    return !(*this==p);
}

template<class X> inline Polynomial<X> partial_evaluate(const Polynomial<X>& p, SizeType k, const X& c) {
    return Polynomial<X>::_partial_evaluate(p,k,c); }

template<class X> inline X evaluate(const Polynomial<X>& p, const Vector<X>& v) {
    return Polynomial<X>::_evaluate(p,v); }

template<class X> inline Polynomial<X> compose(const Polynomial<X>& p, const Vector<Polynomial<X>>& q) {
    return Polynomial<X>::_compose(p,q); }

template<class X> inline Vector<Polynomial<X>> compose(const Vector<Polynomial<X>>& p, const Vector<Polynomial<X>>& q) {
    return Polynomial<X>::_compose(p,q); }

template<class X, class A> inline A evaluate(const Polynomial<X>& p, const Vector<A>& v) {
    return horner_evaluate(p.expansion(),v); }

template<class X> inline Polynomial<X> derivative(Polynomial<X> p, SizeType k) {
    p.differentiate(k); return std::move(p); }

template<class X> inline Polynomial<X> antiderivative(Polynomial<X> p, SizeType k) {
    p.antidifferentiate(k); return std::move(p); }

template<class X> inline Polynomial<X> truncate(Polynomial<X> p, DegreeType deg) {
    p.truncate(deg); return std::move(p); }

template<class X> inline OutputStream& operator<<(OutputStream& os, const Polynomial<X>& p) {
    return p._write(os); }

template<class X> inline Bool compatible(const Polynomial<X>& x1, const Polynomial<X>& x2) {
    return x1.argument_size()==x2.argument_size(); }

template<class F> struct NamedArgumentRepresentation {
    const F& function; const List<String>& argument_names;
};

template<class F> inline NamedArgumentRepresentation<F> named_argument_repr(const F& function, const List<String>& argument_names) {
    NamedArgumentRepresentation<F> r={function,argument_names}; return r; }

template<class X> inline OutputStream& operator<<(OutputStream& os, const NamedArgumentRepresentation<Polynomial<X>>& repr) {
    return repr.function._write(os,repr.argument_names); }

template<class X> inline Polynomial<MidpointType<X>> midpoint(const Polynomial<X>& p) {
    return Polynomial<MidpointType<X>>(midpoint(p.expansion())); }


// Vectorised operations
template<class X, class A> Vector<A> evaluate(const Vector<Polynomial<X>>& p, const Vector<A>& v) {
    Vector<A> r(p.size(),v.zero_element());
    for(Nat i=0; i!=p.size(); ++i) { r[i]=evaluate(p[i],v); }
    return r;
}

template<class X> Vector<Polynomial<X>> derivative(const Vector<Polynomial<X>>& p, SizeType j) {
    Vector<Polynomial<X>> r(p.size());
    for(Nat i=0; i!=p.size(); ++i) { r[i]=derivative(p[i],j); }
    return r;
}

template<class X> Vector<Polynomial<X>> antiderivative(const Vector<Polynomial<X>>& p, SizeType j) {
    Vector<Polynomial<X>> r(p.size());
    for(Nat i=0; i!=p.size(); ++i) { r[i]=antiderivative(p[i],j); }
    return r;
}

template<class X> Vector<Polynomial<X>> truncate(const Vector<Polynomial<X>>& p, DegreeType d) {
    Vector<Polynomial<X>> r(p.size());
    for(Nat i=0; i!=p.size(); ++i) { r[i]=truncate(p[i],d); }
    return r;
}

template<class X> Vector<Polynomial<MidpointType<X>>> midpoint(const Vector<Polynomial<X>>& p) {
    Vector<Polynomial<MidpointType<X>>> r(p.size());
    for(Nat i=0; i!=p.size(); ++i) { r[i]=midpoint(p[i]); }
    return r;
}


} // namespace Ariadne

#endif /* ARIADNE_POLYNOMIAL_HPP */
