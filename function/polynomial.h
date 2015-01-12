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
    //! \name Modifying operators

    //! \brief Truncate to degree \a d.
    Polynomial<X>& truncate(DegreeType d);
    //! \brief Differentiate with respect to the \a j<sup>th</sup> variable.
    Polynomial<X>& differentiate(SizeType j);
    //! \brief Antidifferentiate (integrate) with respect to the \a j<sup>th</sup> variable.
    Polynomial<X>& antidifferentiate(SizeType j);
    //@}

    Void check() const;
  private:
    SortedExpansion<X,ReverseLexicographicKeyLess> _expansion;
};

} // namespace Ariadne

#include "polynomial.tcc"

#endif /* ARIADNE_POLYNOMIAL_H */
