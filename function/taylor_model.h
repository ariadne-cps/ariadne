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
#include "function/scaling.h"
#include "geometry/interval.h"

namespace Ariadne {

class UnitInterval;
enum class SplitPart : char { LOWER, MIDDLE, UPPER };

template<class T1, class T2> struct Product;

template<class P, class F> class TaylorModel;
typedef TaylorModel<Approximate,Float> ApproximateTaylorModel;
typedef TaylorModel<Validated,Float> ValidatedTaylorModel;

template<class P, class F> struct IsScalar< TaylorModel<P,F> > { static const Bool value = true; };
template<class P, class F> struct IsAlgebra< TaylorModel<P,F> > { static const Bool value = true; };
template<class P, class F> struct IsNormedAlgebra< TaylorModel<P,F> > { static const Bool value = true; };

class IntersectionException;

struct IntersectionException : public std::runtime_error {
    IntersectionException(const StringType& what) : std::runtime_error(what) { }
};




/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, ScalarTaylorFunction, VectorTaylorFunction, TaylorConstrainedImageSet.
 */
template<class F>
class TaylorModel<Validated,F>
    : public NormedAlgebraMixin<TaylorModel<Validated,F>,ValidatedNumber>
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
    TaylorModel<Validated,F>();
    //! \brief Construct a TaylorModel in \a as arguments with the given accuracy control.
    TaylorModel<Validated,F>(SizeType as, Sweeper swp);
    //! \brief Construct from a map giving the expansion, a constant giving the error, and an accuracy parameter.
    TaylorModel<Validated,F>(const Expansion<CoefficientType>& f, const ErrorType& e, Sweeper swp);
    TaylorModel<Validated,F>(const Expansion<RawFloat>& f, const RawFloat& e, Sweeper swp);
    //! \brief Fast swap with another Taylor model.
    Void swap(TaylorModel<Validated,F>& tm);
    //! \brief The zero element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<Validated,F> create() const;
    //! \brief The zero element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<Validated,F> create_zero() const;
    //! \brief A constant element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<Validated,F> create_constant(NumericType c) const;
    //! \brief The \a j<sup>th</sup> coordinate element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<Validated,F> create_coordinate(SizeType j) const;
    //! \brief Set to zero.
    TaylorModel<Validated,F> create_ball(ErrorType e) const;
    //! \brief Set to zero.
    Void clear();

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<Validated,F>& operator=(const ValidatedNumber& c);
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct the zero quantity in \a as independent variables.
    static TaylorModel<Validated,F> zero(SizeType as, Sweeper swp) {
        TaylorModel<Validated,F> r(as,swp); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<Validated,F> constant(SizeType as, const NumericType& c, Sweeper swp) {
        TaylorModel<Validated,F> r(as,swp); r.set_value(1); r*=c; return r; }
    //! \brief Construct the quantity with expansion \f$x_j\f$ in \a as independent variables.
    static TaylorModel<Validated,F> coordinate(SizeType as, SizeType j, Sweeper swp) {
        TaylorModel<Validated,F> r(as,swp); r.set_gradient(j,1); return r; }
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a e
    static TaylorModel<Validated,F> error(SizeType as, ErrorType e, Sweeper swp) {
        TaylorModel<Validated,F> r(as,swp); r.set_error(e); return r; }
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a e
    static TaylorModel<Validated,F> ball(SizeType as, ErrorType e, Sweeper swp) {
        TaylorModel<Validated,F> r(as,swp); r.set_error(e); return r; }

    //! \brief Construct the quantity which scales the interval \a codom onto the unit interval.
    static TaylorModel<Validated,F> scaling(SizeType as, SizeType j, const ExactInterval& codom, Sweeper swp);

    //! \brief Return the vector of zero variables of size \a rs in \a as arguments.
    static Vector<TaylorModel<Validated,F>> zeros(SizeType rs, SizeType as, Sweeper swp);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector<TaylorModel<Validated,F>> constants(SizeType as, const Vector<ValidatedNumber>& c, Sweeper swp);
    //! \brief Return the vector of variables on the unit domain.
    static Vector<TaylorModel<Validated,F>> coordinates(SizeType as, Sweeper swp);

    //! \brief Return the vector scaling the box \a codom onto the unit box.
    static Vector<TaylorModel<Validated,F>> scalings(const Vector<ExactInterval>& codom, Sweeper swp);
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    Bool operator==(const TaylorModel<Validated,F>& sd) const {
        return this->_expansion==sd._expansion && this->_error == sd._error; }
    //! \brief Inequality operator.
    Bool operator!=(const TaylorModel<Validated,F>& sd) const {
        return !(*this==sd); }
    //! \brief Comparison with another Taylor model.
    Tribool operator<(const TaylorModel<Validated,F>& sd) const {
        return (sd-*this)>0; }
    //! \brief Comparison with another Taylor model.
    Tribool operator>(const TaylorModel<Validated,F>& sd) const {
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
    //! \brief The coefficient of the gradient term \f$df/dx_j\f$.
    const CoefficientType& gradient_value(SizeType j) const { return (*this)[MultiIndex::unit(this->argument_size(),j)]; }

    //! \brief The constant term in the expansion.
    CoefficientType average() const;
    //! \brief The radius of the smallest ball about a constant function containing the model.
    NormType radius() const;
    //! \brief An over-approximation to the supremum norm.
    NormType norm() const;
    //! \brief A value \c e such that analytic functions are evaluated to a tolerance of \c e. Equal to the sweep threshold.
    RawFloat tolerance() const;

    //! \brief Set the error of the expansion.
    Void set_error(const ErrorType& ne) { ARIADNE_ASSERT(ne>=0); this->_error=ne; }
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
    Box<UnitInterval> domain() const;
    //! \brief The codomain of the quantity.
    ExactInterval codomain() const;
    //! \brief An over-approximation to the range of the quantity.
    UpperInterval range() const;
    //! \brief Compute the gradient of the expansion with respect to the \a jth variable over the domain.
    UpperInterval gradient_range(SizeType j) const;

    //! \brief Evaluate the quantity over the interval of points \a x.
    friend ValidatedNumber evaluate(const TaylorModel<Validated,F>&, const Vector<ValidatedNumber>& x);
    //! \brief Evaluate the gradient over the interval of points \a x.
    friend Covector<ValidatedNumber> gradient(const TaylorModel<Validated,F>&, const Vector<ValidatedNumber>& x);
    //! \brief Evaluate the quantity over the interval of points \a x.
    friend TaylorModel<Validated,F> compose(const TaylorModel<Validated,F>&, const Vector<TaylorModel<Validated,F>>& x);
    //! \brief Compose two models, where the second is scaled so that the codomain is a unit box.
    friend Vector<ValidatedNumber> evaluate(const Vector<TaylorModel<Validated,F>>& f, const Vector<ValidatedNumber>& x);
    //! \brief Compose two models, where the second is scaled so that the codomain is a unit box.
    friend Vector<TaylorModel<Validated,F>> compose(const Vector<TaylorModel<Validated,F>>& f, const Vector<TaylorModel<Validated,F>>& g);
    //@}

    //@{
    /*! \name Inplace modifications. */
    //! \brief Scales the model by a function mapping \a dom into the unit interval.
    Void unscale(ExactInterval dom);
    //! \brief Compute the antiderivative (in place).
    Void antidifferentiate(SizeType k);
    //! \brief Compute the weak derivative (in place).
    Void differentiate(SizeType k);
    //@}

    //@{
    /*! \name Validated paradigm-based operations. */
    //! \brief Test if one model refines (is a subset of) another.
    friend Bool refines(const TaylorModel<Validated,F>& tm1, const TaylorModel<Validated,F>& tm2);
    //! \brief Test if one model is disjoint from (is incompatible with) another.
    friend Bool inconsistent(const TaylorModel<Validated,F>& tm1, const TaylorModel<Validated,F>& tm2);
    //! \brief An over-approximation of the intersection of the sets of functions
    //! allowed by the two models. It is guaranteed that any function represented
    //! by both models is also represented by the result.
    friend TaylorModel<Validated,F> refinement(const TaylorModel<Validated,F>& tm1, const TaylorModel<Validated,F>& tm2);
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Remove all terms based on the \a swp conditions.
    TaylorModel<Validated,F>& sweep(const Sweeper& swp);

    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorModel<Validated,F>& sweep();
    //! \brief Sorts keys.
    TaylorModel<Validated,F>& sort();
    //! \brief Remove terms with the same keys. Assumes sorted.
    TaylorModel<Validated,F>& unique();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel<Validated,F>& cleanup();

    //! \brief Set the error to zero.
    //! WARNING: This method does not preserve rigour of the model approximation.
    TaylorModel<Validated,F>& clobber();

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
    Void isma(const ValidatedNumber& c, const TaylorModel<Validated,F>& x);
    //! \brief Fused multiply and add \c r+=x1*x2 .
    Void ifma(const TaylorModel<Validated,F>& x1, const TaylorModel<Validated,F>& x2);
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const TaylorModel<Validated,F>& x);
    //@}

  public:
    OutputStream& str(const OutputStream&) const;
    OutputStream& repr(const OutputStream&) const;
  public: // FIXME: Should be private
    Void _set_error(const RawFloat& ne) { ARIADNE_ASSERT(ne>=0); this->_error=ErrorType(ne); }
    Void _append(MultiIndex const& a, CoefficientType const& v) { this->_expansion.append(a,v); }
};


//! \relates TaylorModel<Validated,F> \brief Test if a model refines another
template<class F> Bool refines(const TaylorModel<Validated,F>& tv1, const TaylorModel<Validated,F>& tv2);
//! \relates TaylorModel<Validated,F> \brief Test if two models are possible inconstent. If \a true, then models may still be consistent.
template<class F> Bool inconsistent(const TaylorModel<Validated,F>& tv1, const TaylorModel<Validated,F>& tv2);
//! \relates TaylorModel<Validated,F>
//! An over-approximation to the common refinement of two Taylor models.
//! Since the intersection of represented sets of functions cannot be represented
//! exactly in the class of TaylorModels, truncation errors as well as roundoff errors
//! may be present. In the absence of roundoff errors, the result is a subset of both
//! arguments, and is guaranteed to contain any function contained in both arguments.
template<class F> TaylorModel<Validated,F> refinement(const TaylorModel<Validated,F>& x1, const TaylorModel<Validated,F>& x2);

//! \relates TaylorModel<Validated,F> \brief Antidifferentiation operator
template<class F> TaylorModel<Validated,F> antiderivative(const TaylorModel<Validated,F>& x, SizeType k);
//! \relates TaylorModel<Validated,F> \brief Differentiation operator; discards error term
template<class F> TaylorModel<Validated,F> weak_derivative(const TaylorModel<Validated,F>& x, SizeType k);
template<class F> TaylorModel<Validated,F> derivative(const TaylorModel<Validated,F>& x, SizeType k);
template<class F> Covector<ValidatedNumber> gradient(const TaylorModel<Validated,F>& x, const Vector<ValidatedNumber>& y);

//! \relates TaylorModel<Validated,F> \brief Evaluate an array of Taylor variables on a vector.
template<class F> ValidatedNumber evaluate(const TaylorModel<Validated,F>& x, const Vector<ValidatedNumber>& sy);
//! \relates TaylorModel<Validated,F> \brief Substite \a c for the \a k th variable.
template<class F> TaylorModel<Validated,F> partial_evaluate(const TaylorModel<Validated,F>& x, SizeType k, ValidatedNumber c);

//! \relates Rescale the vector \a x from the domain \a dom to the unit domain.
template<class F> TaylorModel<Validated,F> unscale(const TaylorModel<Validated,F>& tv, const ExactInterval& ivl);
//! \relates TaylorModel<Validated,F Compose a vector of Taylor models with another.
template<class F> TaylorModel<Validated,F> compose(const Unscaling& u, const TaylorModel<Validated,F>& tg);
//! \relates TaylorModel<Validated,F Compose a vector of Taylor models with another.
template<class F> TaylorModel<Validated,F> compose(const TaylorModel<Validated,F>& tf, const Vector<TaylorModel<Validated,F>>& tg);
//! \relates TaylorModel<Validated,F> \brief Scale the variable by post-composing with an affine map taking the interval ivl to the unit interval
template<class F> TaylorModel<Validated,F> compose(const TaylorModel<Validated,F>& tf, const VectorUnscaling& u);
template<class F> Vector<TaylorModel<Validated,F>> compose(VectorUnscaling const& u, const Vector<TaylorModel<Validated,F>>& tf);

//! \relates TaylorModel<Validated,F>
//!\brief Split the variable, subdividing along the independent variable j, taking the lower/middle/upper half.
template<class F> TaylorModel<Validated,F> split(const TaylorModel<Validated,F>& x, SizeType j, SplitPart part);
//! \relates TaylorModel<Validated,F> \brief Embed the model in a space of higher dimension
template<class F> TaylorModel<Validated,F> embed(SizeType as1, const TaylorModel<Validated,F>& tm2, SizeType as3);

//! \relates TaylorModel<Validated,F> \brief Embed the model in a space of higher dimension, placing the error in the final variable.
template<class F> TaylorModel<Validated,F> embed_error(const TaylorModel<Validated,F>& tm);
//! \relates TaylorModel<Validated,F> \brief Abstract away the given variables.
//! For example, the model tm(x0,x1,x2,x3) becomes tm'(x0,x1)=tm(x0,[-1,+1],x1,[-1,+1]) on discarding x1 and x3.
template<class F> TaylorModel<Validated,F>  discard_variables(const TaylorModel<Validated,F>&, const Array<SizeType>& variables);

//! \relates TaylorModel<Validated,F> \brief The magnitude of the model. Returns an over-approximation to the supremum norm.
template<class F> ErrorType mag(const TaylorModel<Validated,F>& tm);
//! \relates TaylorModel<Validated,F> \brief An over-approximation to the supremum norm.
template<class F> ErrorType norm(const TaylorModel<Validated,F>& tm);

template<class F> TaylorModel<Validated,F> max(const TaylorModel<Validated,F>& x, const TaylorModel<Validated,F>& y);
template<class F> TaylorModel<Validated,F> min(const TaylorModel<Validated,F>& x, const TaylorModel<Validated,F>& y);
template<class F> TaylorModel<Validated,F> abs(const TaylorModel<Validated,F>& x);




/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, Polynomila, TaylorModel<Validated,F><ExactInterval>.
 */
template<class F>
class TaylorModel<Approximate,F>
    : public NormedAlgebraMixin<TaylorModel<Approximate,F>,ApproximateNumber>
{
  public:
    typedef ApproximateFloat CoefficientType;
    typedef ApproximateErrorType ErrorType;
    typedef ReverseLexicographicKeyLess ComparisonType;
    typedef SortedExpansion<CoefficientType,ComparisonType> ExpansionType;

    typedef ExactInterval CodomainType;
    typedef ApproximateInterval RangeType;
    typedef ApproximateFloat NormType;

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
    //! \brief Construct a TaylorModel<Validated,F> in \a as arguments.
    TaylorModel<Approximate,F>(SizeType as = 0u);
    TaylorModel<Approximate,F>(SizeType as, Sweeper swp);

    TaylorModel<Approximate,F> create() const { return TaylorModel<Approximate,F>(this->argument_size(),this->_sweeper); }
    TaylorModel<Approximate,F> create_zero() const { return TaylorModel<Approximate,F>(this->argument_size(),this->_sweeper); }
    TaylorModel<Approximate,F> create_constant(NumericType) const;
    TaylorModel<Approximate,F> create_variable(SizeType i) const;
    TaylorModel<Approximate,F> create_ball(ErrorType r) const;

    //! \brief Fast swap with another Taylor model.
    Void swap(TaylorModel<Approximate,F>& other) { this->_expansion.swap(other._expansion); }
    //! \brief Set to zero.
    Void clear() { this->_expansion.clear(); }
    //! \brief A zero element with the same parameters.
    TaylorModel<Approximate,F> null() const { return TaylorModel<Approximate,F>(this->argument_size()); }
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a built-in, keeping the same number of arguments.
    TaylorModel<Approximate,F>& operator=(double c) {
        this->_expansion.clear(); this->_expansion.append(MultiIndex(this->argument_size()),CoefficientType(c)); return *this; }
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorModel<Approximate,F>& operator=(const ApproximateNumber& c) {
        this->_expansion.clear(); this->_expansion.append(MultiIndex(this->argument_size()),c); return *this; }
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<Approximate,F>& operator=(const ValidatedNumber& c) { return (*this)=ApproximateNumber(c); }
    //! \brief Set equal to a real constant, keeping the same number of arguments.
    TaylorModel<Approximate,F>& operator=(const EffectiveNumber& c) { return (*this)=ApproximateNumber(c); }
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    Bool operator==(const TaylorModel<Approximate,F>& other) const {
        return this->_expansion==other._expansion; }
    //! \brief Inequality operator.
    Bool operator!=(const TaylorModel<Approximate,F>& other) const {
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
    TaylorModel<Approximate,F>& sweep();
    //! \brief Combine terms with the same index; requires sorted.
    TaylorModel<Approximate,F>& unique();
    //! \brief Sort the terms in index order.
    TaylorModel<Approximate,F>& sort();
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
    virtual Void isma(const ApproximateNumber& c, const TaylorModel<Approximate,F>& x);
    //! \brief Inplace addition of a product of Taylor models.
    virtual Void ifma(const TaylorModel<Approximate,F>& x1, const TaylorModel<Approximate,F>& x2);


    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const TaylorModel<Approximate,F>& x);
    //@}

  public:
    OutputStream& str(OutputStream&) const;
    OutputStream& repr(OutputStream&) const;
};


template<class F> inline OutputStream& operator<<(OutputStream& os, const TaylorModel<Approximate,F>& x) {
    x.str(os); return os; }

template<class P, class F> inline Box<ExactInterval> codomain(const Vector<TaylorModel<P,F> >& vtm) {
    Box<ExactInterval> r(vtm.size()); for(SizeType i=0; i!=vtm.size(); ++i) { r[i]=vtm[i].codomain(); } return r; }

template<class F> TaylorModel<Approximate,F> unscale(const TaylorModel<Approximate,F>& tv, const ExactInterval& ivl);




// Vector operations which can be evaluated componentwise
template<class F> Bool refines(const Vector<TaylorModel<Validated,F>>&, const Vector<TaylorModel<Validated,F>>&);
template<class F> Bool inconsistent(const Vector<TaylorModel<Validated,F>>&, const Vector<TaylorModel<Validated,F>>&);
template<class F> Vector<TaylorModel<Validated,F>> refinement(const Vector<TaylorModel<Validated,F>>&, const Vector<TaylorModel<Validated,F>>&);

template<class F> Vector<TaylorModel<Validated,F>> split(const Vector<TaylorModel<Validated,F>>& x, SizeType k, SplitPart part);

template<class F> Vector<ValidatedNumber> evaluate(const Vector<TaylorModel<Validated,F>>& x, const Vector<ValidatedNumber>& sy);
template<class F> Vector<TaylorModel<Validated,F>> partial_evaluate(const Vector<TaylorModel<Validated,F>>& x, SizeType k, ValidatedNumber const& sy);
template<class F> Vector<TaylorModel<Validated,F>> compose(const VectorUnscaling& u, const Vector<TaylorModel<Validated,F>>& x);
template<class F> Vector<TaylorModel<Validated,F>> compose(const VectorScaling& u, const Vector<TaylorModel<Validated,F>>& x);
template<class F> Vector<TaylorModel<Validated,F>> compose(const Vector<TaylorModel<Validated,F>>& f, const Vector<TaylorModel<Validated,F>>& g);
template<class F> Vector<TaylorModel<Validated,F>> compose(const Vector<TaylorModel<Validated,F>>& f, const Vector<ExactInterval>& dom, const Vector<TaylorModel<Validated,F>>& g);

template<class F> Vector<TaylorModel<Validated,F>> antiderivative(const Vector<TaylorModel<Validated,F>>& x, SizeType k);
template<class F> Vector<TaylorModel<Validated,F>> weak_derivative(const Vector<TaylorModel<Validated,F>>& x, SizeType k);

template<class F> Matrix<ValidatedNumber> jacobian(const Vector<TaylorModel<Validated,F>>& x, const Vector<ValidatedNumber>& y);
template<class F> Matrix<ValidatedNumber> jacobian(const Vector<TaylorModel<Validated,F>>& x, const Vector<ValidatedNumber>& y, Array<SizeType>& p);
template<class F> Matrix<ExactFloat> jacobian_value(const Vector<TaylorModel<Validated,F>>& x);
template<class F> Matrix<ExactFloat> jacobian_value(const Vector<TaylorModel<Validated,F>>& x, const Array<SizeType>& p);
template<class F> Matrix<UpperInterval> jacobian_range(const Vector<TaylorModel<Validated,F>>& x);
template<class F> Matrix<UpperInterval> jacobian_range(const Vector<TaylorModel<Validated,F>>& x, const Array<SizeType>& p);


template<class F> Vector<TaylorModel<Validated,F>> embed(SizeType as1, const Vector<TaylorModel<Validated,F>>& tm2, SizeType as3);
template<class F> Vector<TaylorModel<Validated,F>> combine(const Vector<TaylorModel<Validated,F>>& x1, const Vector<TaylorModel<Validated,F>>& x2);
template<class F> Vector<TaylorModel<Validated,F>> combine(const Vector<TaylorModel<Validated,F>>& x1, const TaylorModel<Validated,F>& x2);
template<class F> Vector<TaylorModel<Validated,F>> combine(const TaylorModel<Validated,F>& x1, const Vector<TaylorModel<Validated,F>>& x2);
template<class F> Vector<TaylorModel<Validated,F>> combine(const TaylorModel<Validated,F>& x1, const TaylorModel<Validated,F>& x2);

} // namespace Ariadne

#endif // ARIADNE_TAYLOR_MODEL_H
