/***************************************************************************
 *            taylor_model.h
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

/*! \file taylor_model.h
 *  \brief Approximate functions on a is_bounded domain with a sparse representation.
 */

#ifndef ARIADNE_TAYLOR_MODEL_H
#define ARIADNE_TAYLOR_MODEL_H

#include <map>

#include "utility/macros.h"
#include "utility/declarations.h"
#include "utility/array.h"
#include "utility/pointer.h"
#include "algebra/vector.h"
#include "algebra/covector.h"
#include "algebra/multi_index.h"
#include "algebra/expansion.h"
#include "algebra/sweeper.h"
#include "algebra/algebra_mixin.h"
#include "algebra/algebra_operations.h"
#include "function/scaling.h"
#include "geometry/interval.h"

namespace Ariadne {

class UnitInterval;
enum class SplitPart : char;

template<class T1, class T2> struct Product;

template<class P, class F> class TaylorModel;
typedef TaylorModel<Approximate,Float64> ApproximateTaylorModel;
typedef TaylorModel<Validated,Float64> ValidatedTaylorModel;

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
    typedef ExactFloat64 CoefficientType;
    typedef ErrorFloat64 ErrorType;
    typedef ErrorFloat64 NormType;
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
    TaylorModel<Validated,F>(const Expansion<RawFloat64>& f, const RawFloat64& e, Sweeper swp);
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
    TaylorModel<Validated,F>& operator =(const ValidatedNumber& c);
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
    friend Bool same(const TaylorModel<Validated,F>& tm1, const TaylorModel<Validated,F>& tm2) {
        return same(tm1._expansion, tm2._expansion) && same(tm1._error, tm2._error); }
    Bool operator==(const TaylorModel<Validated,F>& sd) const {
        return same(*this,sd); }
    //! \brief Inequality operator.
    Bool operator!=(const TaylorModel<Validated,F>& sd) const {
        return !same(*this,sd); }
    Kleenean operator<(const TaylorModel<Validated,F>& sd) const {
        return (sd-*this)>0; }
    //! \brief Comparison with another Taylor model.
    Kleenean operator>(const TaylorModel<Validated,F>& sd) const {
        return (*this-sd)>0; }

    //! \brief Comparison with a scalar.
    Kleenean operator<(Int c) const {
        return this->range()<c; }
    //! \brief Comparison with a scalar.
    Kleenean operator>(Int c) const {
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
    Covector<CoefficientType> gradient_value() const { Covector<CoefficientType> r(this->argument_size());
        for(SizeType j=0; j!=this->argument_size(); ++j) { r[j]=(*this)[MultiIndex::unit(this->argument_size(),j)]; } return r; }

    //! \brief The constant term in the expansion.
    CoefficientType average() const;
    //! \brief The radius of the smallest ball about a constant function containing the model.
    NormType radius() const;
    //! \brief An over-approximation to the supremum norm.
    NormType norm() const;
    //! \brief An over-approximation to the supremum norm.
    friend typename TaylorModel<Validated,F>::NormType norm(const TaylorModel<Validated,F>& tm) { return tm.norm(); }
    friend typename TaylorModel<Validated,F>::NormType mag(const TaylorModel<Validated,F>& tm) { return tm.norm(); }


    //! \brief A value \c e such that analytic functions are evaluated to a tolerance of \c e. Equal to the sweep threshold.
    RawFloat64 tolerance() const;


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
    Covector<UpperInterval> gradient_range() const;

    //! \brief Evaluate the quantity over the interval of points \a x.
    friend ApproximateNumber evaluate(const TaylorModel<Validated,F>& f, const Vector<ApproximateNumber>& x) {
        return TaylorModel<Validated,F>::_evaluate(f,x); }
    friend ValidatedNumber evaluate(const TaylorModel<Validated,F>& f, const Vector<ValidatedNumber>& x) {
        return TaylorModel<Validated,F>::_evaluate(f,x); }
    //! \brief Evaluate the gradient over the interval of points \a x.
    friend Covector<ValidatedNumber> gradient(const TaylorModel<Validated,F>& f, const Vector<ValidatedNumber>& x) {
        return TaylorModel<Validated,F>::_gradient(f,x); }
    //! \brief Substite \a c for the \a k th variable.
    friend TaylorModel<Validated,F> partial_evaluate(const TaylorModel<Validated,F>& x, SizeType k, ValidatedNumber c) {
        return TaylorModel<Validated,F>::_partial_evaluate(x,k,c); }
    //! \relates TaylorModel<Validated,F Compose a vector of Taylor models with another.
    friend TaylorModel<Validated,F> compose(const Unscaling& uf, const TaylorModel<Validated,F>& tg) {
        return TaylorModel<Validated,F>::_compose(uf,tg); }
    //! \brief Evaluate the quantity over the interval of points \a x.
    friend TaylorModel<Validated,F> compose(const TaylorModel<Validated,F>& f, const Vector<TaylorModel<Validated,F>>& g) {
        return TaylorModel<Validated,F>::_compose(f,g); }
    //! \brief Scale the variable by post-composing with an affine map taking the interval ivl to the unit interval
    friend TaylorModel<Validated,F> compose(const TaylorModel<Validated,F>& tf, const VectorUnscaling& u);

    //@}

    //@{
    /*! \name Inplace modifications. */
    //! \brief Scales the model by a function mapping \a dom into the unit interval.
    Void unscale(ExactInterval const& dom);
    //! \brief Compute the antiderivative (in place).
    Void antidifferentiate(SizeType k);
    //! \brief Compute the weak derivative (in place).
    Void differentiate(SizeType k);

    friend TaylorModel<Validated,F> unscale(TaylorModel<Validated,F> tm, ExactInterval const& dom) {
        tm.unscale(dom); return std::move(tm); }
    friend TaylorModel<Validated,F> antiderivative(TaylorModel<Validated,F> tm, SizeType k) {
        tm.antidifferentiate(k); return std::move(tm); }
    friend TaylorModel<Validated,F> derivative(TaylorModel<Validated,F> tm, SizeType k) {
        tm.differentiate(k); return std::move(tm); }
    //@}

    //@{
    /*! \name Operations on the domain. */
    //!\brief Split the Taylor model \a tm, subdividing along the independent variable \a k, taking the lower/middle/upper \a part.
    friend TaylorModel<Validated,F> split(const TaylorModel<Validated,F>& tm, SizeType k, SplitPart part) {
        return TaylorModel<Validated,F>::_split(tm,k,part); }
    //! \relates TaylorModel<Validated,F> \brief Embed the model in a space of higher dimension
    friend TaylorModel<Validated,F> embed(SizeType as1, const TaylorModel<Validated,F>& tm2, SizeType as3) {
        return TaylorModel<Validated,F>::_embed(as1,tm2,as3); }
    //@}

    //@{
    /*! \name Validated paradigm-based operations. */
    //! \brief Test if one model is disjoint from (is incompatible with) another.
    friend Bool consistent(const TaylorModel<Validated,F>& tm1, const TaylorModel<Validated,F>& tm2) {
        return TaylorModel<Validated,F>::_consistent(tm1,tm2); }
    //! \brief Test if one model is disjoint from (is incompatible with) another.
    friend Bool inconsistent(const TaylorModel<Validated,F>& tm1, const TaylorModel<Validated,F>& tm2) {
        return TaylorModel<Validated,F>::_inconsistent(tm1,tm2); }
    //! \brief Test if one model refines (is a subset of) another.
    friend Bool refines(const TaylorModel<Validated,F>& tm1, const TaylorModel<Validated,F>& tm2) {
        return TaylorModel<Validated,F>::_refines(tm1,tm2); }
    //! \brief An over-approximation to the common refinement of two Taylor models.
    //! Since the intersection of represented sets of functions cannot be represented
    //! exactly in the class of TaylorModels, truncation errors as well as roundoff errors
    //! may be present. In the absence of roundoff errors, the result is a subset of both
    //! arguments, and is guaranteed to contain any function contained in both arguments.
    friend TaylorModel<Validated,F> refinement(const TaylorModel<Validated,F>& tm1, const TaylorModel<Validated,F>& tm2) {
        return TaylorModel<Validated,F>::_refinement(tm1,tm2); }
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \relates TaylorModel<Validated,F> \brief Embed the model in a space of higher dimension, placing the error in the final variable.
    friend TaylorModel<Validated,F> embed_error(const TaylorModel<Validated,F>& tm) {
        return TaylorModel<Validated,F>::_embed_error(tm); }
    //! \relates TaylorModel<Validated,F> \brief Abstract away the given variables.
    //! For example, the model tm(x0,x1,x2,x3) becomes tm'(x0,x1)=tm(x0,[-1,+1],x1,[-1,+1]) on discarding x1 and x3.
    friend TaylorModel<Validated,F>  discard_variables(const TaylorModel<Validated,F>& tm, const Array<SizeType>& variables) {
        return TaylorModel<Validated,F>::_discard_variables(tm,variables); }

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

    friend TaylorModel<Validated,F> operator-(TaylorModel<Validated,F> x) { x.imul(-1); return std::move(x); }
    friend TaylorModel<Validated,F>& operator+=(TaylorModel<Validated,F>& x, NumericType const& c) { x.iadd(c); return x; }
    friend TaylorModel<Validated,F>& operator*=(TaylorModel<Validated,F>& x, NumericType const& c) { x.imul(c); return x; }
    friend TaylorModel<Validated,F> operator+(TaylorModel<Validated,F> const& x1, TaylorModel<Validated,F> const& x2) {
        TaylorModel<Validated,F> r=x1; r.isma(+1,x2); return std::move(r); }
    friend TaylorModel<Validated,F> operator-(TaylorModel<Validated,F> const& x1, TaylorModel<Validated,F> const& x2) {
        TaylorModel<Validated,F> r=x1; r.isma(-1,x2); return std::move(r); }
    friend TaylorModel<Validated,F> operator*(TaylorModel<Validated,F> const& x1, TaylorModel<Validated,F> const& x2) {
        TaylorModel<Validated,F> r(x1.argument_size(),x1.sweeper()); r.ifma(x1,x2); return std::move(r); }
    //@}

    //@{
    /*! \name Order operators. */
    //! \brief The pointwise maximum.
    template<class FF> friend TaylorModel<Validated,FF> max(const TaylorModel<Validated,FF>& x, const TaylorModel<Validated,FF>& y);
    //! \brief The pointwise minimum.
    template<class FF> friend TaylorModel<Validated,FF> min(const TaylorModel<Validated,FF>& x, const TaylorModel<Validated,FF>& y);
    //! \brief The pointwise absolute value.
    //! \details If the range of \a x definitely does not include 0, returns +x or -x. Otherwise, uses a uniform polynomial approximation to abs.
    template<class FF> friend TaylorModel<Validated,FF> abs(const TaylorModel<Validated,FF>& x);
    //@}
    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const TaylorModel<Validated,F>& x) { return x.str(os); }
    //@}

  public:
    OutputStream& str(OutputStream&) const;
    OutputStream& repr(OutputStream&) const;
  public: // FIXME: Should be private
    Void _set_error(const RawFloat64& ne) { ARIADNE_ASSERT(ne>=0); this->_error=ErrorType(ne); }
    Void _append(MultiIndex const& a, CoefficientType const& v) { this->_expansion.append(a,v); }
    static Bool _consistent(const TaylorModel<Validated,F>& tm1, const TaylorModel<Validated,F>& tm2);
    static Bool _inconsistent(const TaylorModel<Validated,F>& tm1, const TaylorModel<Validated,F>& tm2);
    static Bool _refines(const TaylorModel<Validated,F>& tm1, const TaylorModel<Validated,F>& tm2);
    static TaylorModel<Validated,F> _refinement(const TaylorModel<Validated,F>& tm1, const TaylorModel<Validated,F>& tm2);
    static TaylorModel<Validated,F> _antiderivative(const TaylorModel<Validated,F>& tm, SizeType k);
    static TaylorModel<Validated,F> _weak_derivative(const TaylorModel<Validated,F>& tm, SizeType k);
    static TaylorModel<Validated,F> _embed_error(const TaylorModel<Validated,F>& tm);
    static TaylorModel<Validated,F>  _discard_variables(const TaylorModel<Validated,F>&, const Array<SizeType>& variables);
    static TaylorModel<Validated,F> _split(const TaylorModel<Validated,F>& tm, SizeType k, SplitPart part);
    static TaylorModel<Validated,F> _embed(SizeType as1, const TaylorModel<Validated,F>& tm2, SizeType as3);
    static TaylorModel<Validated,F> _compose(TaylorModel<Validated,F> const& tf, const Vector<TaylorModel<Validated,F>>& tg);
    static TaylorModel<Validated,F> _compose(Unscaling const& uf, TaylorModel<Validated,F> const& tg);
    static TaylorModel<Validated,F> _compose(TaylorModel<Validated,F> const& tf, VectorUnscaling const& u, const Vector<TaylorModel<Validated,F>>& tg);
    static TaylorModel<Validated,F> _partial_evaluate(const TaylorModel<Validated,F>& x, SizeType k, ValidatedNumber c);
    static Covector<ValidatedNumber> _gradient(const TaylorModel<Validated,F>& x, Vector<ValidatedNumber> const& v);
    static ValidatedNumber _evaluate(const TaylorModel<Validated,F>& x, Vector<ValidatedNumber> const& v);
    static ApproximateNumber _evaluate(const TaylorModel<Validated,F>& x, Vector<ApproximateNumber> const& v);
};

Covector<ValidatedNumber> gradient(const TaylorModel<Validated,Float64>& x, const Vector<ValidatedNumber>& v);



/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, Polynomila, TaylorModel<Validated,F><ExactInterval>.
 */
template<class F>
class TaylorModel<Approximate,F>
    : public NormedAlgebraMixin<TaylorModel<Approximate,F>,ApproximateNumber>
{
  public:
    typedef ApproximateFloat64 CoefficientType;
    typedef ApproximateErrorType ErrorType;
    typedef ReverseLexicographicKeyLess ComparisonType;
    typedef SortedExpansion<CoefficientType,ComparisonType> ExpansionType;

    typedef ExactInterval CodomainType;
    typedef ApproximateInterval RangeType;
    typedef ApproximateFloat64 NormType;

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
    Void unscale(const ExactInterval& dom);
    //! \brief Compute the antiderivative (in place).
    Void antidifferentiate(SizeType k);
    //! \brief Compute the derivative (in place).
    Void differentiate(SizeType k);

    friend TaylorModel<Approximate,F> unscale(TaylorModel<Approximate,F> tm, const ExactInterval& dom) {
        tm.unscale(dom); return std::move(tm); }

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
    virtual RawFloat64 tolerance() const;
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

    template<class FF> friend TaylorModel<Approximate,FF> max(const TaylorModel<Approximate,FF>& x, const TaylorModel<Approximate,FF>& y);
    template<class FF> friend TaylorModel<Approximate,FF> min(const TaylorModel<Approximate,FF>& x, const TaylorModel<Approximate,FF>& y);
    template<class FF> friend TaylorModel<Approximate,FF> abs(const TaylorModel<Approximate,FF>& x);

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const TaylorModel<Approximate,F>& x) {
        return x.str(os); }
    //@}

  public:
    OutputStream& str(OutputStream&) const;
    OutputStream& repr(OutputStream&) const;
};



template<class F> Bool inconsistent(const Vector<TaylorModel<Validated,F>>& tm1, const Vector<TaylorModel<Validated,F>>& tm2) {
    ARIADNE_ASSERT(tm1.size()==tm2.size());
    for(SizeType i=0; i!=tm1.size(); ++i) { if(inconsistent(tm1[i],tm2[i])) { return true; } }
    return false;
}

template<class F> Bool refines(const Vector<TaylorModel<Validated,F>>& tm1, const Vector<TaylorModel<Validated,F>>& tm2) {
    ARIADNE_ASSERT(tm1.size()==tm2.size());
    for(SizeType i=0; i!=tm1.size(); ++i) { if(not refines(tm1[i],tm2[i])) { return false; } }
    return true;
}

template<class F> Vector<TaylorModel<Validated,F>> refinement(const Vector<TaylorModel<Validated,F>>& tm1, const Vector<TaylorModel<Validated,F>>& tm2) {
    ARIADNE_ASSERT(tm1.size()==tm2.size());
    Vector<TaylorModel<Validated,F>> r(tm1.size());
    for(SizeType i=0; i!=tm1.size(); ++i) { r[i]=refinement(tm1[i],tm2[i]); }
    return std::move(r);
}





template<class F> Vector<TaylorModel<Validated,F>> partial_evaluate(const Vector<TaylorModel<Validated,F>>& tf, SizeType k, const ValidatedNumber& c) {
    Vector<TaylorModel<Validated,F>> r(tf.size(),ValidatedTaylorModel::zero(tf.zero_element().argument_size()-1,tf.zero_element().sweeper()));
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=partial_evaluate(tf[i],k,c); }
    return std::move(r);
}

template<class F> Vector<ValidatedNumber> evaluate(const Vector<TaylorModel<Validated,F>>& tf, const Vector<ValidatedNumber>& x) {
    Vector<ValidatedNumber> r(tf.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=evaluate(tf[i],x); }
    return std::move(r);
}

template<class F> Vector<TaylorModel<Validated,F>> compose(const Vector<TaylorModel<Validated,F>>& tf, const Vector<TaylorModel<Validated,F>>& tg) {
    Vector<TaylorModel<Validated,F>> r(tf.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=compose(tf[i],tg); }
    return std::move(r);
}

template<class P, class F> Vector<TaylorModel<P,F>> unscale(const Vector<TaylorModel<P,F>>& x, const Vector<ExactInterval>& dom) {
    Vector<TaylorModel<P,F>> r(x.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=unscale(x[i],dom[i]); }
    return std::move(r);
}

template<class P, class F> Vector<TaylorModel<P,F>> compose(const VectorUnscaling& u, const Vector<TaylorModel<P,F>>& g) {
    ARIADNE_ASSERT_MSG(u.size()==g.size(),"u="<<u.domain()<<", g="<<g);
    Vector<TaylorModel<P,F>> r(u.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=compose(u[i],g[i]); }
    return std::move(r);
}

template<class P, class F> Vector<TaylorModel<P,F>> compose(const Vector<TaylorModel<P,F>>& x,
                                                            const VectorUnscaling& u,
                                                            const Vector<TaylorModel<P,F>>& y) {
    return compose(x,compose(u,y));
}

template<class P, class F> Vector<TaylorModel<P,F>> antiderivative(const Vector<TaylorModel<P,F>>& x, SizeType k) {
    Vector<TaylorModel<P,F>> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=antiderivative(x[i],k); }
    return std::move(r);
}

template<class P, class F> Vector<TaylorModel<P,F>> derivative(const Vector<TaylorModel<P,F>>& x, SizeType k) {
    Vector<TaylorModel<P,F>> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=derivative(x[i],k); }
    return std::move(r);
}

template<class P, class F> Vector<TaylorModel<P,F>> combine(const Vector<TaylorModel<P,F>>& x1, const Vector<TaylorModel<P,F>>& x2) {
    return join(embed(0u,x1,x2.zero_element().argument_size()),embed(x1.zero_element().argument_size(),x2,0u));
}

template<class F> Vector<TaylorModel<Validated,F>> embed(SizeType as1, const Vector<TaylorModel<Validated,F>>& x2, SizeType as3) {
    Vector<TaylorModel<Validated,F>> r(x2.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=embed(as1,x2[i],as3); }
    return std::move(r);
}

template<class P, class F> Vector<TaylorModel<P,F>> split(const Vector<TaylorModel<P,F>>& x, SizeType j, SplitPart h) {
    Vector<TaylorModel<P,F>> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=split(x[i],j,h); }
    return std::move(r);
}

template<class P, class F> typename TaylorModel<P,F>::RangeType ranges(const Vector<TaylorModel<P,F>>& f) {
    Vector<typename TaylorModel<P,F>::RangeType> r(f.size()); for(SizeType i=0; i!=f.size(); ++i) { r[i]=f[i].range(); } return std::move(r);
}

template<class P, class F> Vector<typename TaylorModel<P,F>::ErrorType> errors(const Vector<TaylorModel<P,F>>& h) {
    Vector<typename TaylorModel<P,F>::ErrorType> e(h.size()); for(SizeType i=0; i!=h.size(); ++i) { e[i]=h[i].error(); } return e; }

template<class P, class F> Vector<typename TaylorModel<P,F>::NormType> norms(const Vector<TaylorModel<P,F>>& h) {
    Vector<NormType> r(h.size()); for(SizeType i=0; i!=h.size(); ++i) { r[i]=norm(h[i]); } return std::move(r); }

template<class P, class F> typename TaylorModel<P,F>::NormType norm(const Vector<TaylorModel<P,F>>& h) {
    typename TaylorModel<P,F>::NormType r=0u; for(SizeType i=0; i!=h.size(); ++i) { r=max(r,norm(h[i])); } return std::move(r);
}

template<class F> Matrix<ValidatedNumber> jacobian(const Vector<TaylorModel<Validated,F>>& x, const Vector<ValidatedNumber>& y);
template<class F> Matrix<ValidatedNumber> jacobian(const Vector<TaylorModel<Validated,F>>& x, const Vector<ValidatedNumber>& y, Array<SizeType>& p);
template<class F> Matrix<ExactFloat64> jacobian_value(const Vector<TaylorModel<Validated,F>>& x);
template<class F> Matrix<ExactFloat64> jacobian_value(const Vector<TaylorModel<Validated,F>>& x, const Array<SizeType>& p);
template<class F> Matrix<UpperInterval> jacobian_range(const Vector<TaylorModel<Validated,F>>& x);
template<class F> Matrix<UpperInterval> jacobian_range(const Vector<TaylorModel<Validated,F>>& x, const Array<SizeType>& p);




} // namespace Ariadne

#endif // ARIADNE_TAYLOR_MODEL_H
