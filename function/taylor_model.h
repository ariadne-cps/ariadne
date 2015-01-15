/***************************************************************************
 *            taylor_model.h
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

/*! \file taylor_model.h
 *  \brief Approximate functions on a bounded domain with a sparse representation.
 */

#ifndef ARIADNE_TAYLOR_MODEL_H
#define ARIADNE_TAYLOR_MODEL_H

#include <map>

#include "utility/macros.h"
#include "utility/declarations.h"
#include "utility/array.h"
#include "utility/pointer.h"
#include "algebra/vector.h"
#include "algebra/multi_index.h"
#include "algebra/expansion.h"
#include "algebra/sweeper.h"
#include "algebra/algebra_mixin.h"
#include "geometry/interval.h"

namespace Ariadne {

class UnitInterval;

template<class T1, class T2> struct Product;

template<class X> class TaylorModel;
typedef TaylorModel<ApproximateFloat> ApproximateTaylorModel;
typedef TaylorModel<ValidatedFloat> ValidatedTaylorModel;

template<class X> struct IsScalar< TaylorModel<X> > { static const Bool value = true; };
template<class X> struct IsAlgebra< TaylorModel<X> > { static const Bool value = true; };
template<class X> struct IsNormedAlgebra< TaylorModel<X> > { static const Bool value = true; };

class IntersectionException;

struct IntersectionException : public std::runtime_error {
    IntersectionException(const StringType& what) : std::runtime_error(what) { }
};




/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, ScalarTaylorFunction, VectorTaylorFunction, TaylorConstrainedImageSet.
 */
template<>
class TaylorModel<ValidatedFloat>
    : public NormedAlgebraMixin<TaylorModel<ValidatedFloat>,ValidatedNumber>
{
  public:
    typedef ExactFloat CoefficientType;
    typedef ErrorFloat ErrorType;
    typedef ErrorFloat NormType;
    typedef ReverseLexicographicKeyLess ComparisonType;
    typedef SortedExpansion<CoefficientType,ComparisonType> ExpansionType;

    typedef ExactInterval CodomainType;
    typedef ApproximateInterval RangeType;

    //! \brief The computational paradigm.
    typedef Validated Paradigm;

    //! \brief The type used for the coefficients.
    typedef ValidatedNumber NumericType;

    typedef ValidatedScalarFunction FunctionType;
    typedef ValidatedScalarFunction ScalarFunctionType;
    typedef ValidatedVectorFunction VectorFunctionType;

    //! \brief The type used to index the coefficients.
    typedef MultiIndex IndexType;
    //! \brief The type used for the coefficients.
    typedef CoefficientType ValueType;

    //! \brief An Iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::Iterator Iterator;
    //! \brief A constant Iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::ConstIterator ConstIterator;
  private:
    ExpansionType _expansion;
    ErrorType _error;
    mutable Sweeper _sweeper;
  public:
    //@{
    /*! \name Constructors and destructors. */
    //! \brief Default constructor.
    TaylorModel<ValidatedFloat>();
    //! \brief Construct a TaylorModel<ValidatedFloat> in \a as arguments with the given accuracy control.
    TaylorModel<ValidatedFloat>(SizeType as, Sweeper swp);
    //! \brief Construct from a map giving the expansion, a constant giving the error, and an accuracy parameter.
    TaylorModel<ValidatedFloat>(const Expansion<CoefficientType>& f, const ErrorType& e, Sweeper swp);
    TaylorModel<ValidatedFloat>(const Expansion<RawFloat>& f, const RawFloat& e, Sweeper swp);
    //! \brief Fast swap with another Taylor model.
    Void swap(TaylorModel<ValidatedFloat>& tm);
    //! \brief The zero element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<ValidatedFloat> create() const;
    //! \brief The zero element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<ValidatedFloat> create_zero() const;
    //! \brief The \a j<sup>th</sup> coordinate element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<ValidatedFloat> create_coordinate(SizeType j) const;
    //! \brief Set to zero.
    TaylorModel<ValidatedFloat> create_ball(ErrorType e) const;
    //! \brief Set to zero.
    Void clear();

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<ValidatedFloat>& operator=(const ValidatedNumber& c);
    template<class X, typename std::enable_if<std::is_same<X,RawFloat>::value,Int>::type=0>
        TaylorModel<ValidatedFloat>& operator=(const X& c) { return (*this)=static_cast<ValidatedNumber>(c); }
    template<class X, typename std::enable_if<std::is_same<X,double>::value,Int>::type=0>
        TaylorModel<ValidatedFloat>& operator=(const X& c) { return (*this)=static_cast<ValidatedNumber>(c); }
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct the zero quantity in \a as independent variables.
    static TaylorModel<ValidatedFloat> zero(SizeType as, Sweeper swp) {
        TaylorModel<ValidatedFloat> r(as,swp); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<ValidatedFloat> constant(SizeType as, double c, Sweeper swp) {
        return TaylorModel<ValidatedFloat>::constant(as,CoefficientType(c),swp); }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<ValidatedFloat> constant(SizeType as, const CoefficientType& c, Sweeper swp) {
        TaylorModel<ValidatedFloat> r(as,swp); r.set_value(c); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<ValidatedFloat> constant(SizeType as, const ValidatedNumber& c, Sweeper swp) {
        TaylorModel<ValidatedFloat> r(as,swp); r.set_value(1); r*=c; return r; }
    //! \brief Construct the quantity with expansion \f$x_j\f$ in \a as independent variables.
    static TaylorModel<ValidatedFloat> variable(SizeType as, SizeType j, Sweeper swp) {
        TaylorModel<ValidatedFloat> r(as,swp); r.set_gradient(j,1); return r; }
    //! \brief Construct the quantity which scales the unit interval into the domain \a dom.
    static TaylorModel<ValidatedFloat> scaling(SizeType as, SizeType j, const ExactInterval& dom, Sweeper swp) {
        TaylorModel<ValidatedFloat> r(as,swp); r.set_gradient(j,1); r.rescale(ExactInterval(-1,1),dom); return r; }
    //! \brief Construct the quantity which scales the codomain \a codom into the unit interval.
    static TaylorModel<ValidatedFloat> unscaling(SizeType as, SizeType j, const ExactInterval& codom, Sweeper swp) {
        TaylorModel<ValidatedFloat> r(as,swp); r.set_gradient(j,1); r.rescale(codom,ExactInterval(-1,+1)); return r; }
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a e
    static TaylorModel<ValidatedFloat> error(SizeType as, ErrorType e, Sweeper swp) {
        TaylorModel<ValidatedFloat> r(as,swp); r.set_error(e); return r; }

    //! \brief Return the vector of zero variables of size \a rs in \a as arguments.
    static Vector<TaylorModel<ValidatedFloat>> zeros(SizeType rs, SizeType as, Sweeper swp);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector<TaylorModel<ValidatedFloat>> constants(SizeType as, const Vector<ExactNumber>& c, Sweeper swp);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector<TaylorModel<ValidatedFloat>> constants(SizeType as, const Vector<ValidatedNumber>& c, Sweeper swp);
    //! \brief Return the vector of variables on the unit domain.
    static Vector<TaylorModel<ValidatedFloat>> variables(SizeType as, Sweeper swp);
    //! \brief Return the vector scaling the unit interval onto the domain \a d.
    static Vector<TaylorModel<ValidatedFloat>> scalings(const Vector<ExactInterval>& dom, Sweeper swp);
    //! \brief Return the vector scaling the unit interval onto the codomain \a cd.
    static Vector<TaylorModel<ValidatedFloat>> unscalings(const Vector<ExactInterval>& dom, Sweeper swp);
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    Bool operator==(const TaylorModel<ValidatedFloat>& sd) const {
        return this->_expansion==sd._expansion && this->_error == sd._error; }
    //! \brief Inequality operator.
    Bool operator!=(const TaylorModel<ValidatedFloat>& sd) const {
        return !(*this==sd); }
    //! \brief Comparison with another Taylor model.
    Tribool operator<(const TaylorModel<ValidatedFloat>& sd) const {
        return (sd-*this)>0; }
    //! \brief Comparison with another Taylor model.
    Tribool operator>(const TaylorModel<ValidatedFloat>& sd) const {
        return (*this-sd)>0; }

    //! \brief Comparison with a scalar.
    Tribool operator<(double c) const {
        return this->range()<c; }
    //! \brief Comparison with a scalar.
    Tribool operator>(double c) const {
        return this->range()>c; }
    //@}

    //@{
    /*! \name Data access */
    //! \brief The expansion.
    const ExpansionType& expansion() const { return this->_expansion; }
    //! \brief A reference to the expansion.
    ExpansionType& expansion() { return this->_expansion; }
    //! \brief The error of the expansion over the domain.
    const ErrorType& error() const { return this->_error; }
    //! \brief A reference to the error of the expansion over the domain.
    ErrorType& error() { return this->_error; }
    //! \brief The constant term in the expansion.
    const CoefficientType& value() const { return (*this)[MultiIndex::zero(this->argument_size())]; }
    //! \brief A reference to the constant term in the expansion.
    CoefficientType& value() { return (*this)[MultiIndex::zero(this->argument_size())]; }
    //! \brief The coefficient of the gradient term \f$df/dx_j\f$.
    const CoefficientType& gradient(SizeType j) const { return (*this)[MultiIndex::unit(this->argument_size(),j)]; }
    //! \brief A reference to the coefficient of the gradient term \f$df/dx_j\f$.
    CoefficientType& gradient(SizeType j) { return (*this)[MultiIndex::unit(this->argument_size(),j)]; }

    //! \brief The constant term in the expansion.
    CoefficientType average() const { return (*this)[MultiIndex::zero(this->argument_size())]; }
    //! \brief The radius of the smallest ball about a constant function containing the model.
    NormType radius() const;
    //! \brief An over-approximation to the supremum norm.
    NormType norm() const;
    //! \brief A value \c e such that analytic functions are evaluated to a tolerance of \c e. Equal to the sweep threshold.
    RawFloat tolerance() const;

    //! \brief Set the error of the expansion.
    Void set_error(const ErrorType& ne) {
        ARIADNE_ASSERT(ne>=0); this->_error=ne; }
    template<class E, EnableIf<IsSame<E,Float>> =dummy> Void set_error(const E& ne) { set_error(ErrorType(ne)); }
    template<class E, EnableIf<IsSame<E,double>> =dummy> Void set_error(const E& ne) { set_error(ErrorType(ne)); }
    //! \brief Set the constant term in the expansion.
    Void set_value(const CoefficientType& c) {
        this->_expansion.set(MultiIndex::zero(this->argument_size()),c); }
    //! \brief Set the coefficient of the term \f$df/dx_j\f$.
    Void set_gradient(SizeType j, const CoefficientType& c) {
        this->_expansion.set(MultiIndex::unit(this->argument_size(),j),c); }

    //! \brief The coefficient of the term in $x^a$.
    const CoefficientType& operator[](const MultiIndex& a) const { return this->_expansion[a]; }
    //! \brief A read/write reference to the coefficient of the term in $x^a$.
    CoefficientType& operator[](const MultiIndex& a) { return this->_expansion.at(a); }

    //! \brief An Iterator to the first term in the expansion.
    Iterator begin() { return this->_expansion.begin(); }
    //! \brief A constant Iterator to the first term in the expansion.
    ConstIterator begin() const { return this->_expansion.begin(); }
    //! \brief An Iterator to the end of the expansion.
    Iterator end() { return this->_expansion.end(); }
    //! \brief A constant Iterator to the end of the expansion.
    ConstIterator end() const { return this->_expansion.end(); }

    //! \brief The number of variables in the argument of the quantity.
    SizeType argument_size() const { return this->_expansion.argument_size(); }
    //! \brief The maximum degree of terms in the expansion.
    DegreeType degree() const;
    //! \brief The number of nonzero terms in the expansion.
    SizeType number_of_nonzeros() const { return this->_expansion.number_of_nonzeros(); }
    //@}

    //@{
    /*! \name Function evaluation. */
    //! \brief The domain of the quantity, always given by \f$[-1,1]^{\mathrm{as}}\f$.
    Vector<UnitInterval> domain() const;
    //! \brief The codomain of the quantity.
    ExactInterval codomain() const;
    //! \brief An over-approximation to the range of the quantity.
    UpperInterval range() const;
    //! \brief Compute the gradient of the expansion with respect to the \a jth variable over the domain.
    ExactInterval gradient_range(SizeType j) const;

    //! \brief Evaluate the quantity over the interval of points \a x.
    friend ValidatedNumber evaluate(const TaylorModel<ValidatedFloat>&, const Vector<ValidatedNumber>& x);
    //! \brief Evaluate the quantity over the interval of points \a x.
    friend TaylorModel<ValidatedFloat> compose(const TaylorModel<ValidatedFloat>&, const Vector<TaylorModel<ValidatedFloat>>& x);
    //! \brief Compose two models, where the second is scaled so that the codomain is a unit box.
    friend Vector<ValidatedNumber> evaluate(const Vector<TaylorModel<ValidatedFloat>>& f, const Vector<ValidatedNumber>& x);
    //! \brief Compose two models, where the second is scaled so that the codomain is a unit box.
    friend Vector<TaylorModel<ValidatedFloat>> compose(const Vector<TaylorModel<ValidatedFloat>>& f, const Vector<TaylorModel<ValidatedFloat>>& g);
    //@}

    //@{
    /*! \name Inplace modifications. */
    // TODO: Change these to return Void
    //! \brief Scale so that the old codomain maps into the new codomain.
    TaylorModel<ValidatedFloat>& rescale(const ExactInterval& old_codomain, const ExactInterval& new_codomain);
    //! \brief Restrict to a subdomain.
    TaylorModel<ValidatedFloat>& restrict(const Vector<ExactInterval>& new_domain);
    //! \brief Compute the antiderivative (in place).
    TaylorModel<ValidatedFloat>& antidifferentiate(SizeType k);
    //@}

    //@{
    /*! \name Set-based operations. */
    //! \brief Test if one model refines (is a subset of) another.
    friend Bool refines(const TaylorModel<ValidatedFloat>& tm1, const TaylorModel<ValidatedFloat>& tm2);
    //! \brief Test if one model is disjoint from (is incompatible with) another.
    friend Bool disjoint(const TaylorModel<ValidatedFloat>& tm1, const TaylorModel<ValidatedFloat>& tm2);
    //! \brief An over-approximation of the intersection of the sets of functions
    //! allowed by the two models. It is guaranteed that any function represented
    //! by both models is also represented by the result.
    friend TaylorModel<ValidatedFloat> intersection(const TaylorModel<ValidatedFloat>& tm1, const TaylorModel<ValidatedFloat>& tm2);
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorModel<ValidatedFloat>& sweep();
    //! \brief Remove all terms whose coefficient has magnitude less than \a eps.
    TaylorModel<ValidatedFloat>& sweep(const SweeperInterface& accuracy);
    //! \brief Set the error to zero.
    //! WARNING: This method does not preserve rigour of the model approximation.
    TaylorModel<ValidatedFloat>& clobber();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel<ValidatedFloat>& unique_sort();
    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \brief Specify a policy to use to remove low-impact terms.
    Void set_sweeper(Sweeper swp) { this->_sweeper=swp; }
    //! \brief A shared pointer to an object using for removing low-impact terms.
    Sweeper sweeper() const { return this->_sweeper; }
    //@}

    //@{
    /*! \name Inplace arithmetic operations. */
    //! \brief Add a constant numerical scalar \c r+=c .
    Void iadd(const ValidatedNumber& c);
    //! \brief Multiply by a numerical scalar \c r*=c .
    Void imul(const ValidatedNumber& c);
    //! \brief Scalar multiply and add \c r+=c*x .
    Void isma(const ValidatedNumber& c, const TaylorModel<ValidatedFloat>& x);
    //! \brief Fused multiply and add \c r+=x1*x2 .
    Void ifma(const TaylorModel<ValidatedFloat>& x1, const TaylorModel<ValidatedFloat>& x2);
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const TaylorModel<ValidatedFloat>& x);
    //@}

  public:
    OutputStream& str(const OutputStream&) const;
    OutputStream& repr(const OutputStream&) const;

};

// Rescale the vector \a x from the domain \a dom to the unit domain.
Vector<ValidatedNumber> unscale(const Vector<ValidatedNumber>& x, const Vector<ExactInterval>& dom);

//! \relates TaylorModel<ValidatedFloat> \brief The magnitude of the variable
ErrorType mag(const TaylorModel<ValidatedFloat>& tm);
//! \relates TaylorModel<ValidatedFloat> \brief Split the variable over two domains, subdividing along the independent variable j.
Pair< TaylorModel<ValidatedFloat>, TaylorModel<ValidatedFloat> > split(const TaylorModel<ValidatedFloat>& x, SizeType j);
//! \relates TaylorModel<ValidatedFloat>
//!\brief Split the variable, subdividing along the independent variable j
//! and taking the lower/middle/upper half depending on whether half is false, indeterminate or true.
TaylorModel<ValidatedFloat> split(const TaylorModel<ValidatedFloat>& x, SizeType j, Tribool half);

//! \relates TaylorModel<ValidatedFloat> \brief Scale the variable by post-composing with an affine map taking the unit interval to ivl.
TaylorModel<ValidatedFloat> scale(const TaylorModel<ValidatedFloat>& x, const ExactInterval& ivl);
//! \relates TaylorModel<ValidatedFloat> \brief Scale the variable by post-composing with an affine map taking the interval ivl to the unit interval
TaylorModel<ValidatedFloat> unscale(const TaylorModel<ValidatedFloat>& x, const ExactInterval& ivl);
//! \relates TaylorModel<ValidatedFloat> \brief Scale the variable by post-composing with an affine map taking the interval ivl1 to interval \a ivl2
TaylorModel<ValidatedFloat> rescale(const TaylorModel<ValidatedFloat>& x, const ExactInterval& ivl1, const ExactInterval& ivl2);

//! \relates TaylorModel<ValidatedFloat> \brief Evaluate an array of Taylor variables on a vector.
ValidatedNumber evaluate(const TaylorModel<ValidatedFloat>& x, const Vector<ValidatedNumber>& sy);
//! \relates TaylorModel<ValidatedFloat> \brief Substite \a c for the \a k th variable.
TaylorModel<ValidatedFloat> partial_evaluate(const TaylorModel<ValidatedFloat>& x, SizeType k, ValidatedNumber c);
//! \relates TaylorModel<ValidatedFloat>
//! Substitute the TaylorModel<ValidatedFloat> y in the  kth variable of \a x.
//! Precondition: x.argument_size()==y.argument_size()+1
TaylorModel<ValidatedFloat> substitute(const TaylorModel<ValidatedFloat>& x, SizeType k, const TaylorModel<ValidatedFloat>& y);

//! \relates TaylorModel<ValidatedFloat> \brief Embed the model in a space of higher dimension
TaylorModel<ValidatedFloat> embed(const TaylorModel<ValidatedFloat>& tm, SizeType d);
//! \relates TaylorModel<ValidatedFloat> \brief Embed the model in a space of higher dimension, placing the error in variable i.
TaylorModel<ValidatedFloat> embed_error(const TaylorModel<ValidatedFloat>& tm, SizeType d, SizeType i);
//! \relates TaylorModel<ValidatedFloat> \brief Embed the model in a space of higher dimension
TaylorModel<ValidatedFloat> embed(SizeType as, const TaylorModel<ValidatedFloat>& tv);

//! \relates TaylorModel<ValidatedFloat> \brief Test if a model refines another
Bool refines(const TaylorModel<ValidatedFloat>& tv1, const TaylorModel<ValidatedFloat>& tv2);
//! \relates TaylorModel<ValidatedFloat> \brief Test if a model is disjoint from
Bool disjoint(const TaylorModel<ValidatedFloat>& tv1, const TaylorModel<ValidatedFloat>& tv2);

//! \relates TaylorModel<ValidatedFloat> \brief Antidifferentiation operator
TaylorModel<ValidatedFloat> antiderivative(const TaylorModel<ValidatedFloat>& x, SizeType k);

//! \relates TaylorModel<ValidatedFloat> \brief Differentiation operator; discards error term
TaylorModel<ValidatedFloat> derivative(const TaylorModel<ValidatedFloat>& x, SizeType k);

//! \relates TaylorModel<ValidatedFloat> \brief Replace the variale x[k] with a*x[k]+b
TaylorModel<ValidatedFloat> preaffine(const TaylorModel<ValidatedFloat>&, SizeType k, const ValidatedNumber& a, const ValidatedNumber& b);
//! \relates TaylorModel<ValidatedFloat> \brief Restricts the range of the variable x[k] to the interval d.
//! \pre -1 <= d.lower() <= d.upper() <= 1 .
TaylorModel<ValidatedFloat> restrict(const TaylorModel<ValidatedFloat>&, SizeType k, const ExactInterval& d);

//! \brief Abstract away the given variables.
//! For example, the model tm(x0,x1,x2,x3) becomes tm'(x0,x1)=tm(x0,[-1,+1],x1,[-1,+1]) on discarding x1 and x3.
TaylorModel<ValidatedFloat>  discard(const TaylorModel<ValidatedFloat>&, const Array<SizeType>& variables);

TaylorModel<ValidatedFloat> recondition(const TaylorModel<ValidatedFloat>& tm, Array<SizeType>& discarded_variables,
                                  SizeType number_of_error_variables);

TaylorModel<ValidatedFloat> recondition(const TaylorModel<ValidatedFloat>& tm, Array<SizeType>& discarded_variables,
                                  SizeType number_of_error_variables, SizeType index_of_error);

//! \relates TaylorModel<ValidatedFloat>
//! An over-approximation to the intersection of two Taylor models.
//! Since the intersection cannot be represented exactly in the class of
//! TaylorModels, truncation errors as well as roundoff errors may be present.
//! In the absence of roundoff errors, the result is a subset of both arguments,
//! and is guaranteed to contain any function contained in both arguments.
TaylorModel<ValidatedFloat> intersection(const TaylorModel<ValidatedFloat>& x1, const TaylorModel<ValidatedFloat>& x2);

// Compose an Array of Taylor variables with another, assuming that y has been scaled to have unit codomain
TaylorModel<ValidatedFloat> compose(const TaylorModel<ValidatedFloat>& x, const Vector<TaylorModel<ValidatedFloat>>& y);

// Compose an Array of Taylor variables with another, after scaling by the interval vectors
TaylorModel<ValidatedFloat> compose(const TaylorModel<ValidatedFloat>& x, const Vector<ExactInterval>& bx, const Vector<TaylorModel<ValidatedFloat>>& y);

ErrorType norm(const TaylorModel<ValidatedFloat>& tm);
ErrorType norm(const Vector<TaylorModel<ValidatedFloat>>& tv);

TaylorModel<ValidatedFloat> max(const TaylorModel<ValidatedFloat>& x, const TaylorModel<ValidatedFloat>& y);
TaylorModel<ValidatedFloat> min(const TaylorModel<ValidatedFloat>& x, const TaylorModel<ValidatedFloat>& y);
TaylorModel<ValidatedFloat> abs(const TaylorModel<ValidatedFloat>& x);

TaylorModel<ValidatedFloat> sqrt(const TaylorModel<ValidatedFloat>& x);
TaylorModel<ValidatedFloat> rec(const TaylorModel<ValidatedFloat>& x);
TaylorModel<ValidatedFloat> exp(const TaylorModel<ValidatedFloat>& x);
TaylorModel<ValidatedFloat> log(const TaylorModel<ValidatedFloat>& x);
TaylorModel<ValidatedFloat> sin(const TaylorModel<ValidatedFloat>& x);
TaylorModel<ValidatedFloat> cos(const TaylorModel<ValidatedFloat>& x);
TaylorModel<ValidatedFloat> tan(const TaylorModel<ValidatedFloat>& x);


// Vector operations which can be evaluated componentwise
Bool refines(const Vector<TaylorModel<ValidatedFloat>>&,const Vector<TaylorModel<ValidatedFloat>>&);
Bool disjoint(const Vector<TaylorModel<ValidatedFloat>>&,const Vector<TaylorModel<ValidatedFloat>>&);
Pair< Vector<TaylorModel<ValidatedFloat>>, Vector<TaylorModel<ValidatedFloat>> > split(const Vector<TaylorModel<ValidatedFloat>>& x, SizeType j);
Vector<TaylorModel<ValidatedFloat>> split(const Vector<TaylorModel<ValidatedFloat>>& x, SizeType j, Bool half);
Vector<TaylorModel<ValidatedFloat>> unscale(const Vector<TaylorModel<ValidatedFloat>>& x, const Vector<ExactInterval>& bx);
Vector<TaylorModel<ValidatedFloat>> scale(const Vector<TaylorModel<ValidatedFloat>>& x, const Vector<ExactInterval>& bx);
Vector<ValidatedNumber> evaluate(const Vector<TaylorModel<ValidatedFloat>>& x, const Vector<ValidatedNumber>& sy);
Vector<TaylorModel<ValidatedFloat>> partial_evaluate(const Vector<TaylorModel<ValidatedFloat>>& x, SizeType k, ValidatedNumber sy);
Vector<TaylorModel<ValidatedFloat>> substitute(const Vector<TaylorModel<ValidatedFloat>>& x, SizeType k, const TaylorModel<ValidatedFloat>& y);
Vector<TaylorModel<ValidatedFloat>> antiderivative(const Vector<TaylorModel<ValidatedFloat>>& x, SizeType k);
Vector<TaylorModel<ValidatedFloat>> embed(const Vector<TaylorModel<ValidatedFloat>>& x, SizeType as);
Vector<TaylorModel<ValidatedFloat>> embed(SizeType as, const Vector<TaylorModel<ValidatedFloat>>& x);
Vector<ValidatedNumber> gradient(const TaylorModel<ValidatedFloat>& x, const Vector<ValidatedNumber>& y);
Matrix<ValidatedNumber> jacobian(const Vector<TaylorModel<ValidatedFloat>>& x, const Vector<ValidatedNumber>& y);
//Matrix<ExactInterval> jacobian(const Vector<TaylorModel<ValidatedFloat>>& x);
Bool refines(const Vector<TaylorModel<ValidatedFloat>>& x1, const Vector<TaylorModel<ValidatedFloat>>& x2);
Vector<TaylorModel<ValidatedFloat>> combine(const Vector<TaylorModel<ValidatedFloat>>& x1, const Vector<TaylorModel<ValidatedFloat>>& x2);
Vector<TaylorModel<ValidatedFloat>> combine(const Vector<TaylorModel<ValidatedFloat>>& x1, const TaylorModel<ValidatedFloat>& x2);
Vector<TaylorModel<ValidatedFloat>> combine(const TaylorModel<ValidatedFloat>& x1, const Vector<TaylorModel<ValidatedFloat>>& x2);
Vector<TaylorModel<ValidatedFloat>> combine(const TaylorModel<ValidatedFloat>& x1, const TaylorModel<ValidatedFloat>& x2);
Vector<TaylorModel<ValidatedFloat>> compose(const Vector<TaylorModel<ValidatedFloat>>& f, const Vector<TaylorModel<ValidatedFloat>>& g);
Vector<TaylorModel<ValidatedFloat>> compose(const Vector<TaylorModel<ValidatedFloat>>& f, const Vector<ExactInterval>& dom, const Vector<TaylorModel<ValidatedFloat>>& g);

TaylorModel<ValidatedFloat> unchecked_compose(const TaylorModel<ValidatedFloat>& x, const Vector<TaylorModel<ValidatedFloat>>& y);
Vector<TaylorModel<ValidatedFloat>> unchecked_compose(const Vector<TaylorModel<ValidatedFloat>>& x, const Vector<TaylorModel<ValidatedFloat>>& y);
Vector<TaylorModel<ValidatedFloat>> unchecked_compose(const Vector<TaylorModel<ValidatedFloat>>& x, const Vector<ExactInterval>& dom, const Vector<TaylorModel<ValidatedFloat>>& y);



/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, Polynomila, TaylorModel<ValidatedFloat><ExactInterval>.
 */
template<>
class TaylorModel<ApproximateFloat>
    : public NormedAlgebraMixin<TaylorModel<ApproximateFloat>,ApproximateNumber>
{
  public:
    typedef ApproximateFloat CoefficientType;
    typedef ReverseLexicographicKeyLess ComparisonType;
    typedef SortedExpansion<CoefficientType,ComparisonType> ExpansionType;

    typedef ExactInterval CodomainType;
    typedef ApproximateInterval RangeType;

    //! \brief The type used for the coefficients.
    typedef ApproximateNumber NumericType;
    //! \brief The type used to index the coefficients.
    typedef MultiIndex IndexType;
    //! \brief The type used for the coefficients.
    typedef CoefficientType ValueType;

    //! \brief An Iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::Iterator Iterator;
    //! \brief A constant Iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::ConstIterator ConstIterator;
  private:
    ExpansionType _expansion;
    mutable Sweeper _sweeper;
  private:
    static const CoefficientType _zero;

  public:
    //@{
    /*! \name Constructors and destructors. */
    //! \brief Construct a TaylorModel<ValidatedFloat> in \a as arguments.
    TaylorModel<ApproximateFloat>(SizeType as = 0u);
    TaylorModel<ApproximateFloat>(SizeType as, Sweeper swp);

    TaylorModel<ApproximateFloat> create() const { return TaylorModel<ApproximateFloat>(this->argument_size(),this->_sweeper); }
    TaylorModel<ApproximateFloat> create_zero() const { return TaylorModel<ApproximateFloat>(this->argument_size(),this->_sweeper); }
    TaylorModel<ApproximateFloat> create_ball(ApproximateErrorType r) const { return TaylorModel<ApproximateFloat>(this->argument_size(),this->_sweeper); }
    //! \brief Fast swap with another Taylor model.
    Void swap(TaylorModel<ApproximateFloat>& other) { this->_expansion.swap(other._expansion); }
    //! \brief Set to zero.
    Void clear() { this->_expansion.clear(); }
    //! \brief A zero element with the same parameters.
    TaylorModel<ApproximateFloat> null() const { return TaylorModel<ApproximateFloat>(this->argument_size()); }
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a built-in, keeping the same number of arguments.
    TaylorModel<ApproximateFloat>& operator=(double c) {
        this->_expansion.clear(); this->_expansion.append(MultiIndex(this->argument_size()),ApproximateFloat(c)); return *this; }
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorModel<ApproximateFloat>& operator=(const ApproximateNumber& c) {
        this->_expansion.clear(); this->_expansion.append(MultiIndex(this->argument_size()),c); return *this; }
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<ApproximateFloat>& operator=(const ValidatedNumber& c) { return (*this)=ApproximateNumber(c); }
    //! \brief Set equal to a real constant, keeping the same number of arguments.
    TaylorModel<ApproximateFloat>& operator=(const EffectiveNumber& c) { return (*this)=ApproximateNumber(c); }
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    Bool operator==(const TaylorModel<ApproximateFloat>& other) const {
        return this->_expansion==other._expansion; }
    //! \brief Inequality operator.
    Bool operator!=(const TaylorModel<ApproximateFloat>& other) const {
        return !(*this==other); }
    //@}

    //@{
    /*! \name Data access */
    //! \brief The number of variables in the argument of the quantity.
    SizeType argument_size() const { return this->_expansion.argument_size(); }
    //! \brief The coefficient of the term in $x^a$.
    const CoefficientType& operator[](const MultiIndex& a) const { return this->_expansion[a]; }
    //! \brief The error of the expansion over the domain.
    Void error() const { }
    //@}

    //@{
    /*! \name Function evaluation. */
    //! \brief The domain of the quantity.
    Vector<UnitInterval> domain() const;
    //! \brief A coarse over-approximation to the range of the quantity.
    ExactInterval codomain() const;
    //! \brief An over-approximation to the range of the quantity.
    ApproximateInterval range() const;
    //! \brief Compute the gradient of the expansion with respect to the \a jth variable over the domain.
    ApproximateInterval gradient_range(SizeType j) const;
    //@}

    //@{
    /*! \name Inplace modifications. */
    //! \brief Scale so that the old codomain maps into the unit interval.
    Void unscale(const ExactInterval& codomain);
    //! \brief Compute the antiderivative (in place).
    Void antidifferentiate(SizeType k);
    //! \brief Compute the derivative (in place).
    Void differentiate(SizeType k);
    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \brief Specify a policy to use to remove low-impact terms.
    Void set_sweeper(Sweeper swp) { this->_sweeper=swp; }
    //! \brief A shared pointer to an object using for removing low-impact terms.
    Sweeper sweeper() const { return this->_sweeper; }
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorModel<ApproximateFloat>& sweep();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel<ApproximateFloat>& unique_sort();
    //@}

  public:
    //@{
    /*! \name Standard algebra interface. */
    //! \brief An approximation to the norm of the function.
    virtual NormType norm() const;
    //! \brief An approximation to the average value of the function.
    virtual CoefficientType average() const;
    //! \brief The tolerance to which analytic functions should be computed.
    virtual RawFloat tolerance() const;
    //! \brief The radius of the ball containing the functions.
    virtual NormType radius() const;
    //! \brief Write to an output stream.
    virtual OutputStream& write(OutputStream&) const;
    //! \brief Inplace addition of a scalar constant.
    virtual Void iadd(const ApproximateNumber& c);
    //! \brief Inplace multiplication of a scalar constant.
    virtual Void imul(const ApproximateNumber& c);
    //! \brief Inplace addition of a scalar multiple of a Taylor model.
    virtual Void isma(const ApproximateNumber& c, const TaylorModel<ApproximateFloat>& x);
    //! \brief Inplace addition of a product of Taylor models.
    virtual Void ifma(const TaylorModel<ApproximateFloat>& x1, const TaylorModel<ApproximateFloat>& x2);


    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const TaylorModel<ApproximateFloat>& x);
    //@}

  public:
    OutputStream& str(OutputStream&) const;
    OutputStream& repr(OutputStream&) const;
};


inline OutputStream& operator<<(OutputStream& os, const TaylorModel<ApproximateFloat>& x) {
    x.str(os); return os; }

inline Vector<ExactInterval> codomain(const Vector< TaylorModel<ApproximateFloat> >& t) {
    Vector<ExactInterval> r(t.size()); for(SizeType i=0; i!=t.size(); ++i) { r[i]=t[i].codomain(); } return r; }



} // namespace Ariadne

#endif // ARIADNE_TAYLOR_MODEL_H
