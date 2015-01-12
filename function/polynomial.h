/***************************************************************************
 *            polynomial.h
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

/*! \file polynomial.h
 *  \brief Base class for polynomial rings.
 */

#ifndef ARIADNE_POLYNOMIAL_H
#define ARIADNE_POLYNOMIAL_H

#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "algebra/multi_index.h"
#include "algebra/expansion.h"
#include "function/taylor_model.h"
#include "algebra/differential.h"
#include "algebra/evaluate.h"


namespace Ariadne {

template<class T> class Array;

//! \brief A monomial with coefficients of some type \a X.
template<class X>
class Monomial
    : public ExpansionValue<X>
{
    Monomial(const MultiIndex& a, const X& x) : ExpansionValue<X>(a,x) { }
    Monomial(const ExpansionValue<X>& v) : ExpansionValue<X>(v) { }
};

//! \ingroup FunctionModule
//! \brief A polynomial with coefficients of some type \a X.
template<class X>
class Polynomial
{
    template<class XX> friend class Polynomial;
  public:
    typedef typename Expansion<X>::ValueType ValueType;
    typedef typename Expansion<X>::Reference Reference;
    typedef typename Expansion<X>::ConstReference ConstReference;
    typedef typename Expansion<X>::Iterator Iterator;
    typedef typename Expansion<X>::ConstIterator ConstIterator;

    typedef typename X::NumericType NumericType;
    typedef Polynomial<X> SelfType;
    typedef ReverseLexicographicKeyLess ComparisonType;
  public:
    //@{
    //! \name Constructors

    //! \brief The zero polynomial in \a as variables.
    explicit Polynomial(SizeType as=0u);
    //! \brief Copy/conversion constructor.
    template<class XX> Polynomial(const Polynomial<XX>& p);
    //! \brief Copy/conversion constructor.
    template<class XX> explicit Polynomial(const Expansion<XX>& e);
    //! \brief A dense polynomial with coefficients given by an initializer list of doubles.
    explicit Polynomial(SizeType as, DegreeType deg, InitializerList<X> lst);
    //! \brief A sparse polynomial with coefficients given by an initializer list of indices and coefficients.
    Polynomial(InitializerList<PairType<InitializerList<Int>,X>> lst);
    //@}

    //! \brief Create the null polynomial in the same number of variables.
    Polynomial<X> create_zero() const;

    //! \brief Create a constant polynomial in \a as variables with value \a c.
    static Polynomial<X> constant(SizeType as, const X& c);
    //! \brief Create a polynomial in \a as variables which returns the value of the \a j<sup>th</sup> variable.
    static Polynomial<X> variable(SizeType as, SizeType j);
    static Polynomial<X> coordinate(SizeType as, SizeType j);
    //! \brief Create an Array of polynomials in \a as variables,
    //! the i<sup>th</sup> of  which returns the value of the i<sup>th</sup> variable.
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
    SizeType number_of_nonzeros() const;
    //! \brief The order of the highest term.
    SizeType degree() const;
    //! \brief The value of the polynomial at zero.
    const X& value() const;
    //! \brief A reference to the coefficient of the term in \f$x^{a_1}\cdots x^{a_n}\f$.
    X& operator[](const MultiIndex& a);
    //! \brief A constant referent to the coefficient of the term in \f$x^{a_1}\cdots x^{a_n}\f$.
    const X& operator[](const MultiIndex& a) const;
    //! \brief A constant reference to the raw data expansion.
    const Expansion<X>& expansion() const;
    //! \brief A reference to the raw data expansion.
    Expansion<X>& expansion();
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

    //! \brief Append the term \f$c x^{a_1}\f$ to the list of terms.
    Void append(const MultiIndex& a, const X& c);
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
    template<class XX, class A> friend A evaluate(const Polynomial<XX>& p, const Vector<A>& v);
    template<class XX> friend Polynomial<XX> compose(const Polynomial<XX>& p, const Vector<Polynomial<XX>>& q);
    template<class XX> friend Polynomial<XX> derivative(Polynomial<XX> dx, SizeType k);
    template<class XX> friend Polynomial<XX> antiderivative(Polynomial<XX> dx, SizeType k);
    template<class XX> friend Polynomial<XX> truncate(Polynomial<XX> dx, DegreeType deg);
    //@}

    Void check() const;
  public:
    static Polynomial<X> _neg(const Polynomial<X>& p);
    static Polynomial<X> _add(const Polynomial<X>& p1, const Polynomial<X>& p2);
    static Polynomial<X> _sub(const Polynomial<X>& p1, const Polynomial<X>& p2);
    static Polynomial<X> _mul(const Polynomial<X>& p1, const Polynomial<X>& p2);
    static Polynomial<X> _add(const Polynomial<X>& p, const X& c);
    static Polynomial<X> _mul(const Polynomial<X>& p, const X& c);
    static Polynomial<X>& _imul(Polynomial<X>& p, const Monomial<X>& m);
    static Polynomial<X> _compose(const Polynomial<X>& p, const Vector<Polynomial<X>>& q);
    static X _evaluate(const Polynomial<X>& x, const Vector<X>& c);
    static Polynomial<X> _partial_evaluate(const Polynomial<X>& x, SizeType k, const X& c);
    OutputStream& _write(OutputStream& os) const;
    OutputStream& _write(OutputStream& os, List<String> const& names) const;
  private:
    Iterator _unique_key();
  private:
    SortedExpansion<X,ReverseLexicographicKeyLess> _expansion;
};




template<class X> template<class XX> Polynomial<X>::Polynomial(const Polynomial<XX>& p)
    : _expansion(p._expansion) { }

template<class X> template<class XX> Polynomial<X>::Polynomial(const Expansion<XX>& e)
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

template<class X> OutputStream& operator<<(OutputStream& os, const Polynomial<X>& p) {
    return p._write(os); }


template<class F> struct NamedArgumentRepresentation {
    const F& function; const List<String>& argument_names;
};

template<class F> NamedArgumentRepresentation<F> named_argument_repr(const F& function, const List<String>& argument_names) {
    NamedArgumentRepresentation<F> r={function,argument_names}; return r; }

template<class X> OutputStream& operator<<(OutputStream& os, const NamedArgumentRepresentation<Polynomial<X>>& repr) {
    return repr.function._write(os,repr.argument_names); }



template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p) {
    return Polynomial<X>::_pos(p); }
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p) {
    return Polynomial<X>::_neg(p); }
template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p, const typename Polynomial<X>::NumericType& c) {
    return Polynomial<X>::_add(p,c); }
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p, const typename Polynomial<X>::NumericType& c) {
    return Polynomial<X>::_add(p,neg(c)); }
template<class X> inline Polynomial<X> operator*(const Polynomial<X>& p, const typename Polynomial<X>::NumericType& c) {
    return Polynomial<X>::_mul(p,c); }
template<class X> inline Polynomial<X> operator/(const Polynomial<X>& p, const typename Polynomial<X>::NumericType& c) {
    return Polynomial<X>::_mul(p,rec(c)); }
template<class X> inline Polynomial<X> operator+(const typename Polynomial<X>::NumericType& c, const Polynomial<X>& p) {
    return Polynomial<X>::_add(p,c); }
template<class X> inline Polynomial<X> operator-(const typename Polynomial<X>::NumericType& c, const Polynomial<X>& p) {
    return Polynomial<X>::_add(neg(p),c); }
template<class X> inline Polynomial<X> operator*(const typename Polynomial<X>::NumericType& c, const Polynomial<X>& p) {
    return Polynomial<X>::_mul(p,c); }

template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p1, const Polynomial<X>& p2)
{ return Polynomial<X>::_add(p1,p2); }
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    return Polynomial<X>::_sub(p1,p2); }
template<class X> inline Polynomial<X> operator*(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    return Polynomial<X>::_mul(p1,p2); }

template<class X> inline Polynomial<X> sqr(const Polynomial<X>& p) {
    return p*p; }

template<class X> inline Polynomial<X> pow(const Polynomial<X>& p, Nat m) {
    Polynomial<X> r=Polynomial<X>::constant(p.argument_size(),1.0); Polynomial<X> q(p);
    while(m) { if(m%2) { r=r*q; } q=q*q; m/=2; } return r;
}


template<class X> inline Polynomial<X>& operator+=(Polynomial<X>& p, const typename Polynomial<X>::SelfType& q) {
    return p=p+q; }
template<class X> inline Polynomial<X>& operator-=(Polynomial<X>& p, const typename Polynomial<X>::SelfType& q) {
    return p=p-q; }

template<class X> inline Polynomial<X>& operator+=(Polynomial<X>& p, const typename Polynomial<X>::NumericType& c) {
    return p=p+c; }
template<class X> inline Polynomial<X>& operator-=(Polynomial<X>& p, const typename Polynomial<X>::NumericType& c) {
    return p=p-c; }
template<class X> inline Polynomial<X>& operator*=(Polynomial<X>& p, const typename Polynomial<X>::NumericType& c) {
    return p=p*c; }
template<class X> inline Polynomial<X>& operator/=(Polynomial<X>& p, const typename Polynomial<X>::NumericType& c) {
    return p=p/c; }

template<class X> inline Polynomial<X>& operator*=(Polynomial<X>& p, const Monomial<X>& m) {
    return Polynomial<X>::_imul(p,m); }



template<class X> Polynomial<MidpointType<X>> midpoint(const Polynomial<X>& p) {
    Polynomial<MidpointType<X>> r(p.argument_size());
    for(auto iter=p.begin(); iter!=p.end(); ++iter) {
        r.append(iter->key(),static_cast<MidpointType<X>>(midpoint(iter->data()))); }
    return r;
}


// Vectorised operations
template<class X, class A> Vector<A> evaluate(const Vector<Polynomial<X>>& p, const Vector<A>& v) {
    Vector<A> r(p.size(),v.zero_element());
    for(Nat i=0; i!=p.size(); ++i) { r[i]=evaluate(p[i],v); }
    return r;
}

template<class X> Vector<Polynomial<X>> derivative(const Vector<Polynomial<X>>& p, Nat j) {
    Vector<Polynomial<X>> r(p.size());
    for(Nat i=0; i!=p.size(); ++i) { r[i]=derivative(p[i],j); }
    return r;
}

template<class X> Vector<Polynomial<X>> antiderivative(const Vector<Polynomial<X>>& p, Nat j) {
    Vector<Polynomial<X>> r(p.size());
    for(Nat i=0; i!=p.size(); ++i) { r[i]=antiderivative(p[i],j); }
    return r;
}

template<class X> Vector<Polynomial<X>> truncate(const Vector<Polynomial<X>>& p, Nat d) {
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

#include "polynomial.tcc"

#endif /* ARIADNE_POLYNOMIAL_H */
