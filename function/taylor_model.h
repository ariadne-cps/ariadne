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
typedef TaylorModel<ApproximateNumber> ApproximateTaylorModel;
typedef TaylorModel<ValidatedNumber> ValidatedTaylorModel;

template<class X> struct IsScalar< TaylorModel<X> > { static const bool value = true; };
template<class X> struct IsAlgebra< TaylorModel<X> > { static const bool value = true; };
template<class X> struct IsNormedAlgebra< TaylorModel<X> > { static const bool value = true; };

template<> struct Arithmetic< TaylorModel<ApproximateNumber>,ApproximateNumber > { typedef TaylorModel<ApproximateNumber> ResultType; };
template<> struct Arithmetic< ApproximateNumber,TaylorModel<ApproximateNumber> > { typedef TaylorModel<ApproximateNumber> ResultType; };
template<> struct Arithmetic< TaylorModel<ApproximateNumber>,TaylorModel<ApproximateNumber> > { typedef TaylorModel<ApproximateNumber> ResultType; };
template<> struct Arithmetic< ExactNumber,TaylorModel<ValidatedNumber> > { typedef TaylorModel<ValidatedNumber> ResultType; };
template<> struct Arithmetic< TaylorModel<ValidatedNumber>,ExactNumber > { typedef TaylorModel<ValidatedNumber> ResultType; };
template<> struct Arithmetic< TaylorModel<ValidatedNumber>,ValidatedNumber > { typedef TaylorModel<ValidatedNumber> ResultType; };
template<> struct Arithmetic< ValidatedNumber,TaylorModel<ValidatedNumber> > { typedef TaylorModel<ValidatedNumber> ResultType; };
template<> struct Arithmetic< TaylorModel<ValidatedNumber>,TaylorModel<ValidatedNumber> > { typedef TaylorModel<ValidatedNumber> ResultType; };

class IntersectionException;

struct IntersectionException : public std::runtime_error {
    IntersectionException(const StringType& what) : std::runtime_error(what) { }
};


ValidatedTaylorModel sqrt(const ValidatedTaylorModel& x);
ValidatedTaylorModel rec(const ValidatedTaylorModel& x);
ValidatedTaylorModel exp(const ValidatedTaylorModel& x);
ValidatedTaylorModel log(const ValidatedTaylorModel& x);
ValidatedTaylorModel sin(const ValidatedTaylorModel& x);
ValidatedTaylorModel cos(const ValidatedTaylorModel& x);
ValidatedTaylorModel tan(const ValidatedTaylorModel& x);


/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, ScalarTaylorFunction, VectorTaylorFunction, TaylorConstrainedImageSet.
 */
template<>
class TaylorModel<ValidatedNumber>
    : public NormedAlgebraMixin<TaylorModel<ValidatedNumber>,ValidatedNumber>
{
    friend class ScalarTaylorFunction;
    friend class VectorTaylorFunction;
    typedef Expansion<CoefficientType> ExpansionType;
    typedef ReverseLexicographicKeyLess ComparisonType;
  public:
    typedef ValidatedNumber NumericType;
  private:
    ExpansionType _expansion;
    ErrorType _error;
    mutable Sweeper _sweeper;
  public:
    //! \brief The type used for the coefficients.
    typedef CoefficientType ScalarType;
    //! \brief The type used to index the coefficients.
    typedef MultiIndex IndexType;
    //! \brief The type used for the coefficients.
    typedef CoefficientType ValueType;

    //! \brief An Iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::Iterator Iterator;
    //! \brief A constant Iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::ConstIterator ConstIterator;

    //@{
    /*! \name Constructors and destructors. */
    //! \brief Default constructor.
    TaylorModel<ValidatedNumber>();
    //! \brief Construct a TaylorModel<ValidatedNumber> in \a as arguments with the given accuracy control.
    TaylorModel<ValidatedNumber>(uint as, Sweeper swp);
    //! \brief Construct from a map giving the expansion, a constant giving the error, and an accuracy parameter.
    TaylorModel<ValidatedNumber>(const Expansion<CoefficientType>& f, const ErrorType& e, Sweeper swp);
    TaylorModel<ValidatedNumber>(const Expansion<RawFloat>& f, const RawFloat& e, Sweeper swp);
    //! \brief Fast swap with another Taylor model.
    void swap(TaylorModel<ValidatedNumber>& tm);
    //! \brief The zero element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<ValidatedNumber> create() const;
    //! \brief The zero element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<ValidatedNumber> create_zero() const;
    //! \brief The \a j<sup>th</sup> coordinate element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<ValidatedNumber> create_coordinate(uint j) const;
    //! \brief Set to zero.
    TaylorModel<ValidatedNumber> create_ball(ErrorType e) const;
    //! \brief Set to zero.
    void clear();

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<ValidatedNumber>& operator=(const ValidatedNumber& c);
    template<class X, typename std::enable_if<std::is_same<X,RawFloat>::value,int>::type=0>
        TaylorModel<ValidatedNumber>& operator=(const X& c) { return (*this)=static_cast<ValidatedNumber>(c); }
    template<class X, typename std::enable_if<std::is_same<X,double>::value,int>::type=0>
        TaylorModel<ValidatedNumber>& operator=(const X& c) { return (*this)=static_cast<ValidatedNumber>(c); }
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct the zero quantity in \a as independent variables.
    static TaylorModel<ValidatedNumber> zero(uint as, Sweeper swp) {
        TaylorModel<ValidatedNumber> r(as,swp); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<ValidatedNumber> constant(uint as, double c, Sweeper swp) {
        return TaylorModel<ValidatedNumber>::constant(as,CoefficientType(c),swp); }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<ValidatedNumber> constant(uint as, const CoefficientType& c, Sweeper swp) {
        TaylorModel<ValidatedNumber> r(as,swp); r.set_value(c); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<ValidatedNumber> constant(uint as, const ValidatedNumber& c, Sweeper swp) {
        TaylorModel<ValidatedNumber> r(as,swp); r.set_value(1); r*=c; return r; }
    //! \brief Construct the quantity with expansion \f$x_j\f$ in \a as independent variables.
    static TaylorModel<ValidatedNumber> variable(uint as, uint j, Sweeper swp) {
        TaylorModel<ValidatedNumber> r(as,swp); r.set_gradient(j,1); return r; }
    //! \brief Construct the quantity which scales the unit interval into the domain \a dom.
    static TaylorModel<ValidatedNumber> scaling(uint as, uint j, const ExactInterval& dom, Sweeper swp) {
        TaylorModel<ValidatedNumber> r(as,swp); r.set_gradient(j,1); r.rescale(ExactInterval(-1,1),dom); return r; }
    //! \brief Construct the quantity which scales the codomain \a codom into the unit interval.
    static TaylorModel<ValidatedNumber> unscaling(uint as, uint j, const ExactInterval& codom, Sweeper swp) {
        TaylorModel<ValidatedNumber> r(as,swp); r.set_gradient(j,1); r.rescale(codom,ExactInterval(-1,+1)); return r; }
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a e
    static TaylorModel<ValidatedNumber> error(uint as, ErrorType e, Sweeper swp) {
        TaylorModel<ValidatedNumber> r(as,swp); r.set_error(e); return r; }

    //! \brief Return the vector of zero variables of size \a rs in \a as arguments.
    static Vector< TaylorModel<ValidatedNumber> > zeros(uint rs, uint as, Sweeper swp);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector< TaylorModel<ValidatedNumber> > constants(uint as, const Vector<ExactNumber>& c, Sweeper swp);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector< TaylorModel<ValidatedNumber> > constants(uint as, const Vector<ValidatedNumber>& c, Sweeper swp);
    //! \brief Return the vector of variables on the unit domain.
    static Vector< TaylorModel<ValidatedNumber> > variables(uint as, Sweeper swp);
    //! \brief Return the vector scaling the unit interval onto the domain \a d.
    static Vector< TaylorModel<ValidatedNumber> > scalings(const Vector<ExactInterval>& dom, Sweeper swp);
    //! \brief Return the vector scaling the unit interval onto the codomain \a cd.
    static Vector< TaylorModel<ValidatedNumber> > unscalings(const Vector<ExactInterval>& dom, Sweeper swp);
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    bool operator==(const TaylorModel<ValidatedNumber>& sd) const {
        return this->_expansion==sd._expansion && this->_error == sd._error; }
    //! \brief Inequality operator.
    bool operator!=(const TaylorModel<ValidatedNumber>& sd) const {
        return !(*this==sd); }
    //! \brief Comparison with another Taylor model.
    Tribool operator<(const TaylorModel<ValidatedNumber>& sd) const {
        return (sd-*this)>0; }
    //! \brief Comparison with another Taylor model.
    Tribool operator>(const TaylorModel<ValidatedNumber>& sd) const {
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
    const CoefficientType& gradient(uint j) const { return (*this)[MultiIndex::unit(this->argument_size(),j)]; }
    //! \brief A reference to the coefficient of the gradient term \f$df/dx_j\f$.
    CoefficientType& gradient(uint j) { return (*this)[MultiIndex::unit(this->argument_size(),j)]; }

    //! \brief The constant term in the expansion.
    CoefficientType average() const { return (*this)[MultiIndex::zero(this->argument_size())]; }
    //! \brief The radius of the smallest ball containing the model.
    ErrorType radius() const;
    //! \brief An over-approximation to the supremum norm.
    NormType norm() const;
    //! \brief A value \c e such that analytic functions are evaluated to a tolerance of \c e. Equal to the sweep threshold.
    RawFloat tolerance() const;

    //! \brief Set the error of the expansion.
    void set_error(const ErrorType& ne) {
        ARIADNE_ASSERT(ne>=0); this->_error=ne; }
    template<class E, EnableIf<IsSame<E,Float>> =dummy> void set_error(const E& ne) { set_error(ErrorType(ne)); }
    template<class E, EnableIf<IsSame<E,double>> =dummy> void set_error(const E& ne) { set_error(ErrorType(ne)); }
    //! \brief Set the constant term in the expansion.
    void set_value(const CoefficientType& c) {
        this->_expansion.set(MultiIndex::zero(this->argument_size()),c,ReverseLexicographicKeyLess()); }
    //! \brief Set the coefficient of the term \f$df/dx_j\f$.
    void set_gradient(uint j, const CoefficientType& c) {
        this->_expansion.set(MultiIndex::unit(this->argument_size(),j),c,ReverseLexicographicKeyLess()); }

    //! \brief The coefficient of the term in $x^a$.
    const CoefficientType& operator[](const MultiIndex& a) const { return this->_expansion[a]; }
    //! \brief A read/write reference to the coefficient of the term in $x^a$.
    CoefficientType& operator[](const MultiIndex& a) { return this->_expansion.at(a,ReverseLexicographicKeyLess()); }

    //! \brief An Iterator to the first term in the expansion.
    Iterator begin() { return this->_expansion.begin(); }
    //! \brief A constant Iterator to the first term in the expansion.
    ConstIterator begin() const { return this->_expansion.begin(); }
    //! \brief An Iterator to the end of the expansion.
    Iterator end() { return this->_expansion.end(); }
    //! \brief A constant Iterator to the end of the expansion.
    ConstIterator end() const { return this->_expansion.end(); }

    //! \brief The number of variables in the argument of the quantity.
    uint argument_size() const { return this->_expansion.argument_size(); }
    //! \brief The maximum degree of terms in the expansion.
    uint degree() const;
    //! \brief The number of nonzero terms in the expansion.
    uint number_of_nonzeros() const { return this->_expansion.number_of_nonzeros(); }
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
    ExactInterval gradient_range(uint j) const;

    //! \brief Evaluate the quantity over the interval of points \a x.
    friend ValidatedNumber evaluate(const TaylorModel<ValidatedNumber>&, const Vector<ValidatedNumber>& x);
    //! \brief Evaluate the quantity over the interval of points \a x.
    friend TaylorModel<ValidatedNumber> compose(const TaylorModel<ValidatedNumber>&, const Vector< TaylorModel<ValidatedNumber> >& x);
    //! \brief Compose two models, where the second is scaled so that the codomain is a unit box.
    friend Vector<ValidatedNumber> evaluate(const Vector< TaylorModel<ValidatedNumber> >& f, const Vector<ValidatedNumber>& x);
    //! \brief Compose two models, where the second is scaled so that the codomain is a unit box.
    friend Vector< TaylorModel<ValidatedNumber> > compose(const Vector< TaylorModel<ValidatedNumber> >& f, const Vector< TaylorModel<ValidatedNumber> >& g);
    //@}

    //@{
    /*! \name Inplace modifications. */
    // TODO: Change these to return void
    //! \brief Scale so that the old codomain maps into the new codomain.
    TaylorModel<ValidatedNumber>& rescale(const ExactInterval& old_codomain, const ExactInterval& new_codomain);
    //! \brief Restrict to a subdomain.
    TaylorModel<ValidatedNumber>& restrict(const Vector<ExactInterval>& new_domain);
    //! \brief Compute the antiderivative (in place).
    TaylorModel<ValidatedNumber>& antidifferentiate(uint k);
    //@}

    //@{
    /*! \name Set-based operations. */
    //! \brief Test if one model refines (is a subset of) another.
    friend bool refines(const TaylorModel<ValidatedNumber>& tm1, const TaylorModel<ValidatedNumber>& tm2);
    //! \brief Test if one model is disjoint from (is incompatible with) another.
    friend bool disjoint(const TaylorModel<ValidatedNumber>& tm1, const TaylorModel<ValidatedNumber>& tm2);
    //! \brief An over-approximation of the intersection of the sets of functions
    //! allowed by the two models. It is guaranteed that any function represented
    //! by both models is also represented by the result.
    friend TaylorModel<ValidatedNumber> intersection(const TaylorModel<ValidatedNumber>& tm1, const TaylorModel<ValidatedNumber>& tm2);
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorModel<ValidatedNumber>& sweep();
    //! \brief Remove all terms whose coefficient has magnitude less than \a eps.
    TaylorModel<ValidatedNumber>& sweep(const SweeperInterface& accuracy);
    //! \brief Set the error to zero.
    //! WARNING: This method does not preserve rigour of the model approximation.
    TaylorModel<ValidatedNumber>& clobber();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel<ValidatedNumber>& unique_sort();
    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \brief Specify a policy to use to remove low-impact terms.
    void set_sweeper(Sweeper swp) { this->_sweeper=swp; }
    //! \brief A shared pointer to an object using for removing low-impact terms.
    Sweeper sweeper() const { return this->_sweeper; }
    //@}

    //@{
    /*! \name Inplace arithmetic operations. */
    //! \brief Add a constant numerical scalar \c r+=c .
    void iadd(const ValidatedNumber& c);
    //! \brief Multiply by a numerical scalar \c r*=c .
    void imul(const ValidatedNumber& c);
    //! \brief Scalar multiply and add \c r+=c*x .
    void isma(const ValidatedNumber& c, const TaylorModel<ValidatedNumber>& x);
    //! \brief Fused multiply and add \c r+=x1*x2 .
    void ifma(const TaylorModel<ValidatedNumber>& x1, const TaylorModel<ValidatedNumber>& x2);
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const TaylorModel<ValidatedNumber>& x);
    //@}

  public:
    OutputStream& str(const OutputStream&) const;
    OutputStream& repr(const OutputStream&) const;

};

// Rescale the vector \a x from the domain \a dom to the unit domain.
Vector<ValidatedNumber> unscale(const Vector<ValidatedNumber>& x, const Vector<ExactInterval>& dom);

//! \relates TaylorModel<ValidatedNumber> \brief The magnitude of the variable
ErrorType mag(const TaylorModel<ValidatedNumber>& tm);
//! \relates TaylorModel<ValidatedNumber> \brief Split the variable over two domains, subdividing along the independent variable j.
Pair< TaylorModel<ValidatedNumber>, TaylorModel<ValidatedNumber> > split(const TaylorModel<ValidatedNumber>& x, uint j);
//! \relates TaylorModel<ValidatedNumber>
//!\brief Split the variable, subdividing along the independent variable j
//! and taking the lower/middle/upper half depending on whether half is false, indeterminate or true.
TaylorModel<ValidatedNumber> split(const TaylorModel<ValidatedNumber>& x, uint j, Tribool half);

//! \relates TaylorModel<ValidatedNumber> \brief Scale the variable by post-composing with an affine map taking the unit interval to ivl.
TaylorModel<ValidatedNumber> scale(const TaylorModel<ValidatedNumber>& x, const ExactInterval& ivl);
//! \relates TaylorModel<ValidatedNumber> \brief Scale the variable by post-composing with an affine map taking the interval ivl to the unit interval
TaylorModel<ValidatedNumber> unscale(const TaylorModel<ValidatedNumber>& x, const ExactInterval& ivl);
//! \relates TaylorModel<ValidatedNumber> \brief Scale the variable by post-composing with an affine map taking the interval ivl1 to interval \a ivl2
TaylorModel<ValidatedNumber> rescale(const TaylorModel<ValidatedNumber>& x, const ExactInterval& ivl1, const ExactInterval& ivl2);

//! \relates TaylorModel<ValidatedNumber> \brief Evaluate an array of Taylor variables on a vector.
ValidatedNumber evaluate(const TaylorModel<ValidatedNumber>& x, const Vector<ValidatedNumber>& sy);
//! \relates TaylorModel<ValidatedNumber> \brief Substite \a c for the \a k th variable.
TaylorModel<ValidatedNumber> partial_evaluate(const TaylorModel<ValidatedNumber>& x, uint k, ValidatedNumber c);
//! \relates TaylorModel<ValidatedNumber>
//! Substitute the TaylorModel<ValidatedNumber> y in the  kth variable of \a x.
//! Precondition: x.argument_size()==y.argument_size()+1
TaylorModel<ValidatedNumber> substitute(const TaylorModel<ValidatedNumber>& x, uint k, const TaylorModel<ValidatedNumber>& y);

//! \relates TaylorModel<ValidatedNumber> \brief Embed the model in a space of higher dimension
TaylorModel<ValidatedNumber> embed(const TaylorModel<ValidatedNumber>& tm, uint d);
//! \relates TaylorModel<ValidatedNumber> \brief Embed the model in a space of higher dimension, placing the error in variable i.
TaylorModel<ValidatedNumber> embed_error(const TaylorModel<ValidatedNumber>& tm, uint d, uint i);
//! \relates TaylorModel<ValidatedNumber> \brief Embed the model in a space of higher dimension
TaylorModel<ValidatedNumber> embed(uint as, const TaylorModel<ValidatedNumber>& tv);

//! \relates TaylorModel<ValidatedNumber> \brief Test if a model refines another
bool refines(const TaylorModel<ValidatedNumber>& tv1, const TaylorModel<ValidatedNumber>& tv2);
//! \relates TaylorModel<ValidatedNumber> \brief Test if a model is disjoint from
bool disjoint(const TaylorModel<ValidatedNumber>& tv1, const TaylorModel<ValidatedNumber>& tv2);

//! \relates TaylorModel<ValidatedNumber> \brief Antidifferentiation operator
TaylorModel<ValidatedNumber> antiderivative(const TaylorModel<ValidatedNumber>& x, uint k);

//! \relates TaylorModel<ValidatedNumber> \brief Differentiation operator; discards error term
TaylorModel<ValidatedNumber> derivative(const TaylorModel<ValidatedNumber>& x, uint k);

//! \relates TaylorModel<ValidatedNumber> \brief Replace the variale x[k] with a*x[k]+b
TaylorModel<ValidatedNumber> preaffine(const TaylorModel<ValidatedNumber>&, uint k, const ValidatedNumber& a, const ValidatedNumber& b);
//! \relates TaylorModel<ValidatedNumber> \brief Restricts the range of the variable x[k] to the interval d.
//! \pre -1 <= d.lower() <= d.upper() <= 1 .
TaylorModel<ValidatedNumber> restrict(const TaylorModel<ValidatedNumber>&, uint k, const ExactInterval& d);

//! \brief Abstract away the given variables.
//! For example, the model tm(x0,x1,x2,x3) becomes tm'(x0,x1)=tm(x0,[-1,+1],x1,[-1,+1]) on discarding x1 and x3.
TaylorModel<ValidatedNumber>  discard(const TaylorModel<ValidatedNumber>&, const Array<uint>& variables);

TaylorModel<ValidatedNumber> recondition(const TaylorModel<ValidatedNumber>& tm, Array<uint>& discarded_variables,
                                  uint number_of_error_variables);

TaylorModel<ValidatedNumber> recondition(const TaylorModel<ValidatedNumber>& tm, Array<uint>& discarded_variables,
                                  uint number_of_error_variables, uint index_of_error);

//! \relates TaylorModel<ValidatedNumber>
//! An over-approximation to the intersection of two Taylor models.
//! Since the intersection cannot be represented exactly in the class of
//! TaylorModels, truncation errors as well as roundoff errors may be present.
//! In the absence of roundoff errors, the result is a subset of both arguments,
//! and is guaranteed to contain any function contained in both arguments.
TaylorModel<ValidatedNumber> intersection(const TaylorModel<ValidatedNumber>& x1, const TaylorModel<ValidatedNumber>& x2);

// Compose an Array of Taylor variables with another, assuming that y has been scaled to have unit codomain
TaylorModel<ValidatedNumber> compose(const TaylorModel<ValidatedNumber>& x, const Vector< TaylorModel<ValidatedNumber> >& y);

// Compose an Array of Taylor variables with another, after scaling by the interval vectors
TaylorModel<ValidatedNumber> compose(const TaylorModel<ValidatedNumber>& x, const Vector<ExactInterval>& bx, const Vector< TaylorModel<ValidatedNumber> >& y);

ErrorType norm(const TaylorModel<ValidatedNumber>& tm);
ErrorType norm(const Vector< TaylorModel<ValidatedNumber> >& tv);

TaylorModel<ValidatedNumber> max(const TaylorModel<ValidatedNumber>& x, const TaylorModel<ValidatedNumber>& y);
TaylorModel<ValidatedNumber> min(const TaylorModel<ValidatedNumber>& x, const TaylorModel<ValidatedNumber>& y);
TaylorModel<ValidatedNumber> abs(const TaylorModel<ValidatedNumber>& x);



// Vector operations which can be evaluated componentwise
bool refines(const Vector< TaylorModel<ValidatedNumber> >&,const Vector< TaylorModel<ValidatedNumber> >&);
bool disjoint(const Vector< TaylorModel<ValidatedNumber> >&,const Vector< TaylorModel<ValidatedNumber> >&);
Pair< Vector< TaylorModel<ValidatedNumber> >, Vector< TaylorModel<ValidatedNumber> > > split(const Vector< TaylorModel<ValidatedNumber> >& x, uint j);
Vector< TaylorModel<ValidatedNumber> > split(const Vector< TaylorModel<ValidatedNumber> >& x, uint j, bool half);
Vector< TaylorModel<ValidatedNumber> > unscale(const Vector< TaylorModel<ValidatedNumber> >& x, const Vector<ExactInterval>& bx);
Vector< TaylorModel<ValidatedNumber> > scale(const Vector< TaylorModel<ValidatedNumber> >& x, const Vector<ExactInterval>& bx);
Vector<ValidatedNumber> evaluate(const Vector< TaylorModel<ValidatedNumber> >& x, const Vector<ValidatedNumber>& sy);
Vector< TaylorModel<ValidatedNumber> > partial_evaluate(const Vector< TaylorModel<ValidatedNumber> >& x, uint k, ValidatedNumber sy);
Vector< TaylorModel<ValidatedNumber> > substitute(const Vector< TaylorModel<ValidatedNumber> >& x, uint k, const TaylorModel<ValidatedNumber>& y);
Vector< TaylorModel<ValidatedNumber> > antiderivative(const Vector< TaylorModel<ValidatedNumber> >& x, uint k);
Vector< TaylorModel<ValidatedNumber> > embed(const Vector< TaylorModel<ValidatedNumber> >& x, uint as);
Vector< TaylorModel<ValidatedNumber> > embed(uint as, const Vector< TaylorModel<ValidatedNumber> >& x);
Matrix<ValidatedNumber> jacobian(const Vector< TaylorModel<ValidatedNumber> >& x, const Vector<ValidatedNumber>& y);
//Matrix<ExactInterval> jacobian(const Vector< TaylorModel<ValidatedNumber> >& x);
bool refines(const Vector< TaylorModel<ValidatedNumber> >& x1, const Vector< TaylorModel<ValidatedNumber> >& x2);
Vector< TaylorModel<ValidatedNumber> > combine(const Vector< TaylorModel<ValidatedNumber> >& x1, const Vector< TaylorModel<ValidatedNumber> >& x2);
Vector< TaylorModel<ValidatedNumber> > combine(const Vector< TaylorModel<ValidatedNumber> >& x1, const TaylorModel<ValidatedNumber>& x2);
Vector< TaylorModel<ValidatedNumber> > combine(const TaylorModel<ValidatedNumber>& x1, const Vector< TaylorModel<ValidatedNumber> >& x2);
Vector< TaylorModel<ValidatedNumber> > combine(const TaylorModel<ValidatedNumber>& x1, const TaylorModel<ValidatedNumber>& x2);
Vector< TaylorModel<ValidatedNumber> > compose(const Vector< TaylorModel<ValidatedNumber> >& f, const Vector< TaylorModel<ValidatedNumber> >& g);
Vector< TaylorModel<ValidatedNumber> > compose(const Vector< TaylorModel<ValidatedNumber> >& f, const Vector<ExactInterval>& dom, const Vector< TaylorModel<ValidatedNumber> >& g);

TaylorModel<ValidatedNumber> unchecked_compose(const TaylorModel<ValidatedNumber>& x, const Vector< TaylorModel<ValidatedNumber> >& y);
Vector< TaylorModel<ValidatedNumber> > unchecked_compose(const Vector< TaylorModel<ValidatedNumber> >& x, const Vector< TaylorModel<ValidatedNumber> >& y);
Vector< TaylorModel<ValidatedNumber> > unchecked_compose(const Vector< TaylorModel<ValidatedNumber> >& x, const Vector<ExactInterval>& dom, const Vector< TaylorModel<ValidatedNumber> >& y);



/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, Polynomila, TaylorModel<ValidatedNumber><ExactInterval>.
 */
template<>
class TaylorModel<ApproximateNumber>
    : public NormedAlgebraMixin<TaylorModel<ApproximateNumber>,ApproximateNumber>
{
    typedef Expansion<ApproximateCoefficientType> ExpansionType;
  private:
    ExpansionType _expansion;
    mutable Sweeper _sweeper;
  private:
    static const ApproximateCoefficientType _zero;

  public:
    //! \brief The type used for the coefficients.
    typedef ApproximateNumber NumericType;
    //! \brief The type used to index the coefficients.
    typedef MultiIndex IndexType;
    //! \brief The type used for the coefficients.
    typedef ApproximateCoefficientType ValueType;

    //! \brief An Iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::Iterator Iterator;
    //! \brief A constant Iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::ConstIterator ConstIterator;

  public:
    //@{
    /*! \name Constructors and destructors. */
    //! \brief Construct a ValidatedTaylorModel in \a as arguments.
    TaylorModel<ApproximateNumber>(uint as = 0u);
    TaylorModel<ApproximateNumber>(uint as, Sweeper swp);

    TaylorModel<ApproximateNumber> create() const { return TaylorModel<ApproximateNumber>(this->argument_size(),this->_sweeper); }
    TaylorModel<ApproximateNumber> create_zero() const { return TaylorModel<ApproximateNumber>(this->argument_size(),this->_sweeper); }
    TaylorModel<ApproximateNumber> create_ball(ApproximateErrorType r) const { return TaylorModel<ApproximateNumber>(this->argument_size(),this->_sweeper); }
    //! \brief Fast swap with another Taylor model.
    void swap(TaylorModel<ApproximateNumber>& other) { this->_expansion.swap(other._expansion); }
    //! \brief Set to zero.
    void clear() { this->_expansion.clear(); }
    //! \brief A zero element with the same parameters.
    TaylorModel<ApproximateNumber> null() const { return TaylorModel<ApproximateNumber>(this->argument_size()); }
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a built-in, keeping the same number of arguments.
    TaylorModel<ApproximateNumber>& operator=(double c) {
        this->_expansion.clear(); this->_expansion.append(MultiIndex(this->argument_size()),ApproximateCoefficientType(c)); return *this; }
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorModel<ApproximateNumber>& operator=(const ApproximateNumber& c) {
        this->_expansion.clear(); this->_expansion.append(MultiIndex(this->argument_size()),c); return *this; }
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<ApproximateNumber>& operator=(const ValidatedNumber& c) { return (*this)=ApproximateNumber(c); }
    //! \brief Set equal to a real constant, keeping the same number of arguments.
    TaylorModel<ApproximateNumber>& operator=(const EffectiveNumber& c) { return (*this)=ApproximateNumber(c); }
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    bool operator==(const TaylorModel<ApproximateNumber>& other) const {
        return this->_expansion==other._expansion; }
    //! \brief Inequality operator.
    bool operator!=(const TaylorModel<ApproximateNumber>& other) const {
        return !(*this==other); }
    //@}

    //@{
    /*! \name Data access */
    //! \brief The number of variables in the argument of the quantity.
    uint argument_size() const { return this->_expansion.argument_size(); }
    //! \brief The coefficient of the term in $x^a$.
    const ApproximateCoefficientType& operator[](const MultiIndex& a) const { return this->_expansion[a]; }
    //! \brief The error of the expansion over the domain.
    const ApproximateErrorType& error() const { return _zero; }
    //@}

    //@{
    /*! \name Function evaluation. */
    //! \brief The domain of the quantity.
    Vector<UnitInterval> domain() const;
    //! \brief A coarse over-approximation to the range of the quantity.
    ExactInterval codomain() const;
    //! \brief An over-approximation to the range of the quantity.
    UpperInterval range() const;
    //! \brief Compute the gradient of the expansion with respect to the \a jth variable over the domain.
    UpperInterval gradient_range(uint j) const;
    //@}

    //@{
    /*! \name Inplace modifications. */
    //! \brief Scale so that the old codomain maps into the unit interval.
    void unscale(const ExactInterval& codomain);
    //! \brief Compute the antiderivative (in place).
    void antidifferentiate(uint k);
    //! \brief Compute the derivative (in place).
    void differentiate(uint k);
    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \brief Specify a policy to use to remove low-impact terms.
    void set_sweeper(Sweeper swp) { this->_sweeper=swp; }
    //! \brief A shared pointer to an object using for removing low-impact terms.
    Sweeper sweeper() const { return this->_sweeper; }
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorModel<ApproximateNumber>& sweep();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel<ApproximateNumber>& unique_sort();
    //@}

  public:
    //@{
    /*! \name Standard algebra interface. */
    //! \brief An approximation to the norm of the function.
    virtual ErrorType norm() const;
    //! \brief An approximation to the average value of the function.
    virtual CoefficientType average() const;
    //! \brief The tolerance to which analytic functions should be computed.
    virtual RawFloat tolerance() const;
    //! \brief The radius of the ball containing the functions.
    virtual ErrorType radius() const;
    //! \brief Write to an output stream.
    virtual OutputStream& write(OutputStream&) const;
    //! \brief Inplace addition of a scalar constant.
    virtual void iadd(const ApproximateNumber& c);
    //! \brief Inplace multiplication of a scalar constant.
    virtual void imul(const ApproximateNumber& c);
    //! \brief Inplace addition of a scalar multiple of a Taylor model.
    virtual void isma(const ApproximateNumber& c, const TaylorModel<ApproximateNumber>& x);
    //! \brief Inplace addition of a product of Taylor models.
    virtual void ifma(const TaylorModel<ApproximateNumber>& x1, const TaylorModel<ApproximateNumber>& x2);


    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const TaylorModel<ApproximateNumber>& x);
    //@}

  public:
    OutputStream& str(OutputStream&) const;
    OutputStream& repr(OutputStream&) const;
};


inline OutputStream& operator<<(OutputStream& os, const TaylorModel<ApproximateNumber>& x) {
    x.str(os); return os; }

inline Vector<ExactInterval> codomain(const Vector< TaylorModel<ApproximateNumber> >& t) {
    Vector<ExactInterval> r(t.size()); for(uint i=0; i!=t.size(); ++i) { r[i]=t[i].codomain(); } return r; }



} // namespace Ariadne

#endif // ARIADNE_TAYLOR_MODEL_H
