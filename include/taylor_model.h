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

#include "macros.h"
#include "declarations.h"
#include "array.h"
#include "pointer.h"
#include "vector.h"
#include "multi_index.h"
#include "expansion.h"
#include "sweeper.h"
#include "algebra_mixin.h"

namespace Ariadne {

class UnitInterval;

template<class T1, class T2> struct Product;

template<class X> class TaylorModel;
typedef TaylorModel<ApproximateTag> ApproximateTaylorModel;
typedef TaylorModel<ValidatedTag> ValidatedTaylorModel;

template<class X> struct IsScalar< TaylorModel<X> > { static const bool value = true; };
template<class X> struct IsAlgebra< TaylorModel<X> > { static const bool value = true; };
template<class X> struct IsNormedAlgebra< TaylorModel<X> > { static const bool value = true; };

template<> struct Arithmetic< TaylorModel<ApproximateTag>,ApproximateNumberType > { typedef TaylorModel<ApproximateTag> ResultType; };
template<> struct Arithmetic< ApproximateNumberType,TaylorModel<ApproximateTag> > { typedef TaylorModel<ApproximateTag> ResultType; };
template<> struct Arithmetic< TaylorModel<ApproximateTag>,TaylorModel<ApproximateTag> > { typedef TaylorModel<ApproximateTag> ResultType; };
template<> struct Arithmetic< ExactNumberType,TaylorModel<ValidatedTag> > { typedef TaylorModel<ValidatedTag> ResultType; };
template<> struct Arithmetic< TaylorModel<ValidatedTag>,ExactNumberType > { typedef TaylorModel<ValidatedTag> ResultType; };
template<> struct Arithmetic< TaylorModel<ValidatedTag>,ValidatedNumberType > { typedef TaylorModel<ValidatedTag> ResultType; };
template<> struct Arithmetic< ValidatedNumberType,TaylorModel<ValidatedTag> > { typedef TaylorModel<ValidatedTag> ResultType; };
template<> struct Arithmetic< TaylorModel<ValidatedTag>,TaylorModel<ValidatedTag> > { typedef TaylorModel<ValidatedTag> ResultType; };

class IntersectionException;

struct IntersectionException : public std::runtime_error {
    IntersectionException(const std::string& what) : std::runtime_error(what) { }
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
class TaylorModel<ValidatedTag>
    : public NormedAlgebraMixin<TaylorModel<ValidatedTag>,ValidatedNumberType>
{
    friend class ScalarTaylorFunction;
    friend class VectorTaylorFunction;
    typedef Expansion<CoefficientType> ExpansionType;
    typedef ReverseLexicographicKeyLess ComparisonType;
  public:
    typedef ValidatedNumberType NumericType;
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

    //! \brief An iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::iterator iterator;
    //! \brief A constant iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::const_iterator const_iterator;

    //@{
    /*! \name Constructors and destructors. */
    //! \brief Default constructor.
    TaylorModel<ValidatedTag>();
    //! \brief Construct a TaylorModel<ValidatedTag> in \a as arguments with the given accuracy control.
    TaylorModel<ValidatedTag>(uint as, Sweeper swp);
    //! \brief Construct from a map giving the expansion, a constant giving the error, and an accuracy parameter.
    TaylorModel<ValidatedTag>(const Expansion<CoefficientType>& f, const CoefficientType& e, Sweeper swp);
    //! \brief Fast swap with another Taylor model.
    void swap(TaylorModel<ValidatedTag>& tm);
    //! \brief The zero element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<ValidatedTag> create() const;
    //! \brief The \a j<sup>th</sup> coordinate element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<ValidatedTag> create_coordinate(uint j) const;
    //! \brief Set to zero.
    TaylorModel<ValidatedTag> create_ball(ErrorType e) const;
    //! \brief Set to zero.
    void clear();

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<ValidatedTag>& operator=(const ValidatedNumberType& c);
    template<class X, typename std::enable_if<std::is_same<X,RawFloatType>::value,int>::type=0>
        TaylorModel<ValidatedTag>& operator=(const X& c) { return (*this)=static_cast<ValidatedNumberType>(c); }
    template<class X, typename std::enable_if<std::is_same<X,double>::value,int>::type=0>
        TaylorModel<ValidatedTag>& operator=(const X& c) { return (*this)=static_cast<ValidatedNumberType>(c); }
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct the zero quantity in \a as independent variables.
    static TaylorModel<ValidatedTag> zero(uint as, Sweeper swp) {
        TaylorModel<ValidatedTag> r(as,swp); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<ValidatedTag> constant(uint as, double c, Sweeper swp) {
        return TaylorModel<ValidatedTag>::constant(as,CoefficientType(c),swp); }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<ValidatedTag> constant(uint as, const CoefficientType& c, Sweeper swp) {
        TaylorModel<ValidatedTag> r(as,swp); r.set_value(c); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<ValidatedTag> constant(uint as, const ValidatedNumberType& c, Sweeper swp) {
        TaylorModel<ValidatedTag> r(as,swp); r.set_value(1); r*=c; return r; }
    //! \brief Construct the quantity with expansion \f$x_j\f$ in \a as independent variables.
    static TaylorModel<ValidatedTag> variable(uint as, uint j, Sweeper swp) {
        TaylorModel<ValidatedTag> r(as,swp); r.set_gradient(j,1); return r; }
    //! \brief Construct the quantity which scales the unit interval into the domain \a dom.
    static TaylorModel<ValidatedTag> scaling(uint as, uint j, const Interval& dom, Sweeper swp) {
        TaylorModel<ValidatedTag> r(as,swp); r.set_gradient(j,1); r.rescale(Interval(-1,1),dom); return r; }
    //! \brief Construct the quantity which scales the codomain \a codom into the unit interval.
    static TaylorModel<ValidatedTag> unscaling(uint as, uint j, const Interval& codom, Sweeper swp) {
        TaylorModel<ValidatedTag> r(as,swp); r.set_gradient(j,1); r.rescale(codom,Interval(-1,+1)); return r; }
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a e
    static TaylorModel<ValidatedTag> error(uint as, ErrorType e, Sweeper swp) {
        TaylorModel<ValidatedTag> r(as,swp); r.set_error(e); return r; }

    //! \brief Return the vector of zero variables of size \a rs in \a as arguments.
    static Vector< TaylorModel<ValidatedTag> > zeros(uint rs, uint as, Sweeper swp);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector< TaylorModel<ValidatedTag> > constants(uint as, const Vector<ExactNumberType>& c, Sweeper swp);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector< TaylorModel<ValidatedTag> > constants(uint as, const Vector<ValidatedNumberType>& c, Sweeper swp);
    //! \brief Return the vector of variables on the unit domain.
    static Vector< TaylorModel<ValidatedTag> > variables(uint as, Sweeper swp);
    //! \brief Return the vector scaling the unit interval onto the domain \a d.
    static Vector< TaylorModel<ValidatedTag> > scalings(const Vector<Interval>& dom, Sweeper swp);
    //! \brief Return the vector scaling the unit interval onto the codomain \a cd.
    static Vector< TaylorModel<ValidatedTag> > unscalings(const Vector<Interval>& dom, Sweeper swp);
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    bool operator==(const TaylorModel<ValidatedTag>& sd) const {
        return this->_expansion==sd._expansion && this->_error == sd._error; }
    //! \brief Inequality operator.
    bool operator!=(const TaylorModel<ValidatedTag>& sd) const {
        return !(*this==sd); }
    //! \brief Comparison with another Taylor model.
    tribool operator<(const TaylorModel<ValidatedTag>& sd) const {
        return (sd-*this)>0; }
    //! \brief Comparison with another Taylor model.
    tribool operator>(const TaylorModel<ValidatedTag>& sd) const {
        return (*this-sd)>0; }

    //! \brief Comparison with a scalar.
    tribool operator<(double c) const {
        return this->range()<c; }
    //! \brief Comparison with a scalar.
    tribool operator>(double c) const {
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
    RawFloatType tolerance() const;

    //! \brief Set the error of the expansion.
    void set_error(const ErrorType& ne) {
        ARIADNE_ASSERT(ne>=0); this->_error=ne; }
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

    //! \brief An iterator to the first term in the expansion.
    iterator begin() { return this->_expansion.begin(); }
    //! \brief A constant iterator to the first term in the expansion.
    const_iterator begin() const { return this->_expansion.begin(); }
    //! \brief An iterator to the end of the expansion.
    iterator end() { return this->_expansion.end(); }
    //! \brief A constant iterator to the end of the expansion.
    const_iterator end() const { return this->_expansion.end(); }

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
    Interval codomain() const;
    //! \brief An over-approximation to the range of the quantity.
    Interval range() const;
    //! \brief Compute the gradient of the expansion with respect to the \a jth variable over the domain.
    Interval gradient_range(uint j) const;

    //! \brief Evaluate the quantity over the interval of points \a x.
    friend ValidatedNumberType evaluate(const TaylorModel<ValidatedTag>&, const Vector<ValidatedNumberType>& x);
    //! \brief Evaluate the quantity over the interval of points \a x.
    friend TaylorModel<ValidatedTag> compose(const TaylorModel<ValidatedTag>&, const Vector< TaylorModel<ValidatedTag> >& x);
    //! \brief Compose two models, where the second is scaled so that the codomain is a unit box.
    friend Vector<ValidatedNumberType> evaluate(const Vector< TaylorModel<ValidatedTag> >& f, const Vector<ValidatedNumberType>& x);
    //! \brief Compose two models, where the second is scaled so that the codomain is a unit box.
    friend Vector< TaylorModel<ValidatedTag> > compose(const Vector< TaylorModel<ValidatedTag> >& f, const Vector< TaylorModel<ValidatedTag> >& g);
    //@}

    //@{
    /*! \name Inplace modifications. */
    // TODO: Change these to return void
    //! \brief Scale so that the old codomain maps into the new codomain.
    TaylorModel<ValidatedTag>& rescale(const Interval& old_codomain, const Interval& new_codomain);
    //! \brief Restrict to a subdomain.
    TaylorModel<ValidatedTag>& restrict(const Vector<Interval>& new_domain);
    //! \brief Compute the antiderivative (in place).
    TaylorModel<ValidatedTag>& antidifferentiate(uint k);
    //@}

    //@{
    /*! \name Set-based operations. */
    //! \brief Test if one model refines (is a subset of) another.
    friend bool refines(const TaylorModel<ValidatedTag>& tm1, const TaylorModel<ValidatedTag>& tm2);
    //! \brief Test if one model is disjoint from (is incompatible with) another.
    friend bool disjoint(const TaylorModel<ValidatedTag>& tm1, const TaylorModel<ValidatedTag>& tm2);
    //! \brief An over-approximation of the intersection of the sets of functions
    //! allowed by the two models. It is guaranteed that any function represented
    //! by both models is also represented by the result.
    friend TaylorModel<ValidatedTag> intersection(const TaylorModel<ValidatedTag>& tm1, const TaylorModel<ValidatedTag>& tm2);
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorModel<ValidatedTag>& sweep();
    //! \brief Remove all terms whose coefficient has magnitude less than \a eps.
    TaylorModel<ValidatedTag>& sweep(const SweeperInterface& accuracy);
    //! \brief Set the error to zero.
    //! WARNING: This method does not preserve rigour of the model approximation.
    TaylorModel<ValidatedTag>& clobber();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel<ValidatedTag>& unique_sort();
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
    void iadd(const ValidatedNumberType& c);
    //! \brief Multiply by a numerical scalar \c r*=c .
    void imul(const ValidatedNumberType& c);
    //! \brief Scalar multiply and add \c r+=c*x .
    void isma(const ValidatedNumberType& c, const TaylorModel<ValidatedTag>& x);
    //! \brief Fused multiply and add \c r+=x1*x2 .
    void ifma(const TaylorModel<ValidatedTag>& x1, const TaylorModel<ValidatedTag>& x2);
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const TaylorModel<ValidatedTag>& x);
    //@}

  public:
    std::ostream& str(const std::ostream&) const;
    std::ostream& repr(const std::ostream&) const;

};

// Rescale the vector \a x from the domain \a dom to the unit domain.
Vector<ValidatedNumberType> unscale(const Vector<ValidatedNumberType>& x, const Vector<Interval>& dom);

//! \relates TaylorModel<ValidatedTag> \brief The magnitude of the variable
ErrorType mag(const TaylorModel<ValidatedTag>& tm);
//! \relates TaylorModel<ValidatedTag> \brief Split the variable over two domains, subdividing along the independent variable j.
std::pair< TaylorModel<ValidatedTag>, TaylorModel<ValidatedTag> > split(const TaylorModel<ValidatedTag>& x, uint j);
//! \relates TaylorModel<ValidatedTag>
//!\brief Split the variable, subdividing along the independent variable j
//! and taking the lower/middle/upper half depending on whether half is false, indeterminate or true.
TaylorModel<ValidatedTag> split(const TaylorModel<ValidatedTag>& x, uint j, tribool half);

//! \relates TaylorModel<ValidatedTag> \brief Scale the variable by post-composing with an affine map taking the unit interval to ivl.
TaylorModel<ValidatedTag> scale(const TaylorModel<ValidatedTag>& x, const Interval& ivl);
//! \relates TaylorModel<ValidatedTag> \brief Scale the variable by post-composing with an affine map taking the interval ivl to the unit interval
TaylorModel<ValidatedTag> unscale(const TaylorModel<ValidatedTag>& x, const Interval& ivl);
//! \relates TaylorModel<ValidatedTag> \brief Scale the variable by post-composing with an affine map taking the interval ivl1 to interval \a ivl2
TaylorModel<ValidatedTag> rescale(const TaylorModel<ValidatedTag>& x, const Interval& ivl1, const Interval& ivl2);

//! \relates TaylorModel<ValidatedTag> \brief Evaluate an array of Taylor variables on a vector.
ValidatedNumberType evaluate(const TaylorModel<ValidatedTag>& x, const Vector<ValidatedNumberType>& sy);
//! \relates TaylorModel<ValidatedTag> \brief Substite \a c for the \a k th variable.
TaylorModel<ValidatedTag> partial_evaluate(const TaylorModel<ValidatedTag>& x, uint k, ValidatedNumberType c);
//! \relates TaylorModel<ValidatedTag>
//! Substitute the TaylorModel<ValidatedTag> y in the  kth variable of \a x.
//! Precondition: x.argument_size()==y.argument_size()+1
TaylorModel<ValidatedTag> substitute(const TaylorModel<ValidatedTag>& x, uint k, const TaylorModel<ValidatedTag>& y);

//! \relates TaylorModel<ValidatedTag> \brief Embed the model in a space of higher dimension
TaylorModel<ValidatedTag> embed(const TaylorModel<ValidatedTag>& tm, uint d);
//! \relates TaylorModel<ValidatedTag> \brief Embed the model in a space of higher dimension, placing the error in variable i.
TaylorModel<ValidatedTag> embed_error(const TaylorModel<ValidatedTag>& tm, uint d, uint i);
//! \relates TaylorModel<ValidatedTag> \brief Embed the model in a space of higher dimension
TaylorModel<ValidatedTag> embed(uint as, const TaylorModel<ValidatedTag>& tv);

//! \relates TaylorModel<ValidatedTag> \brief Test if a model refines another
bool refines(const TaylorModel<ValidatedTag>& tv1, const TaylorModel<ValidatedTag>& tv2);
//! \relates TaylorModel<ValidatedTag> \brief Test if a model is disjoint from
bool disjoint(const TaylorModel<ValidatedTag>& tv1, const TaylorModel<ValidatedTag>& tv2);

//! \relates TaylorModel<ValidatedTag> \brief Antidifferentiation operator
TaylorModel<ValidatedTag> antiderivative(const TaylorModel<ValidatedTag>& x, uint k);

//! \relates TaylorModel<ValidatedTag> \brief Differentiation operator; discards error term
TaylorModel<ValidatedTag> derivative(const TaylorModel<ValidatedTag>& x, uint k);

//! \relates TaylorModel<ValidatedTag> \brief Replace the variale x[k] with a*x[k]+b
TaylorModel<ValidatedTag> preaffine(const TaylorModel<ValidatedTag>&, uint k, const ValidatedNumberType& a, const ValidatedNumberType& b);
//! \relates TaylorModel<ValidatedTag> \brief Restricts the range of the variable x[k] to the interval d.
//! \pre -1 <= d.lower() <= d.upper() <= 1 .
TaylorModel<ValidatedTag> restrict(const TaylorModel<ValidatedTag>&, uint k, const Interval& d);

//! \brief Abstract away the given variables.
//! For example, the model tm(x0,x1,x2,x3) becomes tm'(x0,x1)=tm(x0,[-1,+1],x1,[-1,+1]) on discarding x1 and x3.
TaylorModel<ValidatedTag>  discard(const TaylorModel<ValidatedTag>&, const Array<uint>& variables);

TaylorModel<ValidatedTag> recondition(const TaylorModel<ValidatedTag>& tm, Array<uint>& discarded_variables,
                                  uint number_of_error_variables);

TaylorModel<ValidatedTag> recondition(const TaylorModel<ValidatedTag>& tm, Array<uint>& discarded_variables,
                                  uint number_of_error_variables, uint index_of_error);

//! \relates TaylorModel<ValidatedTag>
//! An over-approximation to the intersection of two Taylor models.
//! Since the intersection cannot be represented exactly in the class of
//! TaylorModels, truncation errors as well as roundoff errors may be present.
//! In the absence of roundoff errors, the result is a subset of both arguments,
//! and is guaranteed to contain any function contained in both arguments.
TaylorModel<ValidatedTag> intersection(const TaylorModel<ValidatedTag>& x1, const TaylorModel<ValidatedTag>& x2);

// Compose an Array of Taylor variables with another, assuming that y has been scaled to have unit codomain
TaylorModel<ValidatedTag> compose(const TaylorModel<ValidatedTag>& x, const Vector< TaylorModel<ValidatedTag> >& y);

// Compose an Array of Taylor variables with another, after scaling by the interval vectors
TaylorModel<ValidatedTag> compose(const TaylorModel<ValidatedTag>& x, const Vector<Interval>& bx, const Vector< TaylorModel<ValidatedTag> >& y);

ErrorType norm(const TaylorModel<ValidatedTag>& tm);
ErrorType norm(const Vector< TaylorModel<ValidatedTag> >& tv);

TaylorModel<ValidatedTag> max(const TaylorModel<ValidatedTag>& x, const TaylorModel<ValidatedTag>& y);
TaylorModel<ValidatedTag> min(const TaylorModel<ValidatedTag>& x, const TaylorModel<ValidatedTag>& y);
TaylorModel<ValidatedTag> abs(const TaylorModel<ValidatedTag>& x);



// Vector operations which can be evaluated componentwise
bool refines(const Vector< TaylorModel<ValidatedTag> >&,const Vector< TaylorModel<ValidatedTag> >&);
bool disjoint(const Vector< TaylorModel<ValidatedTag> >&,const Vector< TaylorModel<ValidatedTag> >&);
std::pair< Vector< TaylorModel<ValidatedTag> >, Vector< TaylorModel<ValidatedTag> > > split(const Vector< TaylorModel<ValidatedTag> >& x, uint j);
Vector< TaylorModel<ValidatedTag> > split(const Vector< TaylorModel<ValidatedTag> >& x, uint j, bool half);
Vector< TaylorModel<ValidatedTag> > unscale(const Vector< TaylorModel<ValidatedTag> >& x, const Vector<Interval>& bx);
Vector< TaylorModel<ValidatedTag> > scale(const Vector< TaylorModel<ValidatedTag> >& x, const Vector<Interval>& bx);
Vector<ValidatedNumberType> evaluate(const Vector< TaylorModel<ValidatedTag> >& x, const Vector<ValidatedNumberType>& sy);
Vector< TaylorModel<ValidatedTag> > partial_evaluate(const Vector< TaylorModel<ValidatedTag> >& x, uint k, ValidatedNumberType sy);
Vector< TaylorModel<ValidatedTag> > substitute(const Vector< TaylorModel<ValidatedTag> >& x, uint k, const TaylorModel<ValidatedTag>& y);
Vector< TaylorModel<ValidatedTag> > antiderivative(const Vector< TaylorModel<ValidatedTag> >& x, uint k);
Vector< TaylorModel<ValidatedTag> > embed(const Vector< TaylorModel<ValidatedTag> >& x, uint as);
Vector< TaylorModel<ValidatedTag> > embed(uint as, const Vector< TaylorModel<ValidatedTag> >& x);
Matrix<ValidatedNumberType> jacobian(const Vector< TaylorModel<ValidatedTag> >& x, const Vector<ValidatedNumberType>& y);
//Matrix<Interval> jacobian(const Vector< TaylorModel<ValidatedTag> >& x);
bool refines(const Vector< TaylorModel<ValidatedTag> >& x1, const Vector< TaylorModel<ValidatedTag> >& x2);
Vector< TaylorModel<ValidatedTag> > combine(const Vector< TaylorModel<ValidatedTag> >& x1, const Vector< TaylorModel<ValidatedTag> >& x2);
Vector< TaylorModel<ValidatedTag> > combine(const Vector< TaylorModel<ValidatedTag> >& x1, const TaylorModel<ValidatedTag>& x2);
Vector< TaylorModel<ValidatedTag> > combine(const TaylorModel<ValidatedTag>& x1, const Vector< TaylorModel<ValidatedTag> >& x2);
Vector< TaylorModel<ValidatedTag> > combine(const TaylorModel<ValidatedTag>& x1, const TaylorModel<ValidatedTag>& x2);
Vector< TaylorModel<ValidatedTag> > compose(const Vector< TaylorModel<ValidatedTag> >& f, const Vector< TaylorModel<ValidatedTag> >& g);
Vector< TaylorModel<ValidatedTag> > compose(const Vector< TaylorModel<ValidatedTag> >& f, const Vector<Interval>& dom, const Vector< TaylorModel<ValidatedTag> >& g);

TaylorModel<ValidatedTag> unchecked_compose(const TaylorModel<ValidatedTag>& x, const Vector< TaylorModel<ValidatedTag> >& y);
Vector< TaylorModel<ValidatedTag> > unchecked_compose(const Vector< TaylorModel<ValidatedTag> >& x, const Vector< TaylorModel<ValidatedTag> >& y);
Vector< TaylorModel<ValidatedTag> > unchecked_compose(const Vector< TaylorModel<ValidatedTag> >& x, const Vector<Interval>& dom, const Vector< TaylorModel<ValidatedTag> >& y);



/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, Polynomila, TaylorModel<ValidatedTag><Interval>.
 */
template<>
class TaylorModel<ApproximateTag>
    : public NormedAlgebraMixin<TaylorModel<ApproximateTag>,ApproximateNumberType>
{
    typedef Expansion<ApproximateCoefficientType> ExpansionType;
  private:
    ExpansionType _expansion;
    mutable Sweeper _sweeper;
  private:
    static const ApproximateCoefficientType _zero;

  public:
    //! \brief The type used for the coefficients.
    typedef ApproximateNumberType NumericType;
    //! \brief The type used to index the coefficients.
    typedef MultiIndex IndexType;
    //! \brief The type used for the coefficients.
    typedef ApproximateCoefficientType ValueType;

    //! \brief An iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::iterator iterator;
    //! \brief A constant iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::const_iterator const_iterator;

  public:
    //@{
    /*! \name Constructors and destructors. */
    //! \brief Construct a ValidatedTaylorModel in \a as arguments.
    TaylorModel<ApproximateTag>(uint as = 0u);
    TaylorModel<ApproximateTag>(uint as, Sweeper swp);

    TaylorModel<ApproximateTag> create() const { return TaylorModel<ApproximateTag>(this->argument_size(),this->_sweeper); }
    TaylorModel<ApproximateTag> create_ball(ApproximateErrorType r) const { return TaylorModel<ApproximateTag>(this->argument_size(),this->_sweeper); }
    //! \brief Fast swap with another Taylor model.
    void swap(TaylorModel<ApproximateTag>& other) { this->_expansion.swap(other._expansion); }
    //! \brief Set to zero.
    void clear() { this->_expansion.clear(); }
    //! \brief A zero element with the same parameters.
    TaylorModel<ApproximateTag> null() const { return TaylorModel<ApproximateTag>(this->argument_size()); }
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a built-in, keeping the same number of arguments.
    TaylorModel<ApproximateTag>& operator=(double c) {
        this->_expansion.clear(); this->_expansion.append(MultiIndex(this->argument_size()),ApproximateCoefficientType(c)); return *this; }
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorModel<ApproximateTag>& operator=(const ApproximateNumberType& c) {
        this->_expansion.clear(); this->_expansion.append(MultiIndex(this->argument_size()),c); return *this; }
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<ApproximateTag>& operator=(const ValidatedNumberType& c) { return (*this)=midpoint(c); }
    //! \brief Set equal to a real constant, keeping the same number of arguments.
    TaylorModel<ApproximateTag>& operator=(const EffectiveNumberType& c) { return (*this)=ApproximateNumberType(c); }
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    bool operator==(const TaylorModel<ApproximateTag>& other) const {
        return this->_expansion==other._expansion; }
    //! \brief Inequality operator.
    bool operator!=(const TaylorModel<ApproximateTag>& other) const {
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
    Interval codomain() const;
    //! \brief An over-approximation to the range of the quantity.
    Interval range() const;
    //! \brief Compute the gradient of the expansion with respect to the \a jth variable over the domain.
    Interval gradient_range(uint j) const;
    //@}

    //@{
    /*! \name Inplace modifications. */
    //! \brief Scale so that the old codomain maps into the unit interval.
    void unscale(const Interval& codomain);
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
    TaylorModel<ApproximateTag>& sweep();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel<ApproximateTag>& unique_sort();
    //@}

  public:
    //@{
    /*! \name Standard algebra interface. */
    //! \brief An approximation to the norm of the function.
    virtual ErrorType norm() const;
    //! \brief An approximation to the average value of the function.
    virtual CoefficientType average() const;
    //! \brief The tolerance to which analytic functions should be computed.
    virtual RawFloatType tolerance() const;
    //! \brief The radius of the ball containing the functions.
    virtual ErrorType radius() const;
    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream&) const;
    //! \brief Inplace addition of a scalar constant.
    virtual void iadd(const ApproximateNumberType& c);
    //! \brief Inplace multiplication of a scalar constant.
    virtual void imul(const ApproximateNumberType& c);
    //! \brief Inplace addition of a scalar multiple of a Taylor model.
    virtual void isma(const ApproximateNumberType& c, const TaylorModel<ApproximateTag>& x);
    //! \brief Inplace addition of a product of Taylor models.
    virtual void ifma(const TaylorModel<ApproximateTag>& x1, const TaylorModel<ApproximateTag>& x2);


    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const TaylorModel<ApproximateTag>& x);
    //@}

  public:
    std::ostream& str(std::ostream&) const;
    std::ostream& repr(std::ostream&) const;
};


inline std::ostream& operator<<(std::ostream& os, const TaylorModel<ApproximateTag>& x) {
    x.str(os); return os; }

inline Vector<Interval> codomain(const Vector< TaylorModel<ApproximateTag> >& t) {
    Vector<Interval> r(t.size()); for(uint i=0; i!=t.size(); ++i) { r[i]=t[i].codomain(); } return r; }



} // namespace Ariadne

#endif // ARIADNE_TAYLOR_MODEL_H
