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
#include "array.h"
#include "pointer.h"
#include "vector.h"
#include "multi_index.h"
#include "expansion.h"
#include "sweeper.h"
#include "algebra_mixin.h"

namespace Ariadne {

class Float;
class Interval;
class Real;

template<class T1, class T2> class Product;
template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Expansion;

template<class X> class ScalarFunction;
template<class X> class VectorFunction;

template<class X> class TaylorModel;
typedef TaylorModel<Float> FloatTaylorModel;
typedef TaylorModel<Interval> IntervalTaylorModel;

template<class X> struct IsAlgebra< TaylorModel<X> > { static const bool value = true; };
template<class X> struct IsNormedAlgebra< TaylorModel<X> > { static const bool value = true; };


template<> struct Arithmetic< TaylorModel<Float>,Float > { typedef TaylorModel<Float> ResultType; };
template<> struct Arithmetic< Float,TaylorModel<Float> > { typedef TaylorModel<Float> ResultType; };
template<> struct Arithmetic< TaylorModel<Float>,TaylorModel<Float> > { typedef TaylorModel<Float> ResultType; };
template<> struct Arithmetic< Float,TaylorModel<Interval> > { typedef TaylorModel<Interval> ResultType; };
template<> struct Arithmetic< TaylorModel<Interval>,Float > { typedef TaylorModel<Interval> ResultType; };
template<> struct Arithmetic< TaylorModel<Interval>,Interval > { typedef TaylorModel<Interval> ResultType; };
template<> struct Arithmetic< Interval,TaylorModel<Interval> > { typedef TaylorModel<Interval> ResultType; };
template<> struct Arithmetic< TaylorModel<Interval>,TaylorModel<Interval> > { typedef TaylorModel<Interval> ResultType; };

class IntersectionException;

struct IntersectionException : public std::runtime_error {
    IntersectionException(const std::string& what) : std::runtime_error(what) { }
};





/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, ScalarTaylorFunction, VectorTaylorFunction, TaylorConstrainedImageSet.
 */
template<>
class TaylorModel<Interval>
    : public NormedAlgebraMixin<TaylorModel<Interval>,Interval>
{
    friend class ScalarTaylorFunction;
    friend class VectorTaylorFunction;
    typedef Expansion<Float> ExpansionType;
    typedef ReverseLexicographicKeyLess ComparisonType;
  public:
    typedef Interval NumericType;
  private:
    ExpansionType _expansion;
    Float _error;
    mutable Sweeper _sweeper;
  public:
    //! \brief The type used for the coefficients.
    typedef Float ScalarType;
    //! \brief The type used to index the coefficients.
    typedef MultiIndex IndexType;
    //! \brief The type used for the coefficients.
    typedef Float ValueType;

    //! \brief An iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::iterator iterator;
    //! \brief A constant iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::const_iterator const_iterator;

    //@{
    /*! \name Constructors and destructors. */
    //! \brief Default constructor.
    TaylorModel<Interval>();
    //! \brief Construct a TaylorModel<Interval> in \a as arguments with the given accuracy control.
    TaylorModel<Interval>(uint as, Sweeper swp);
    //! \brief Construct from a map giving the expansion, a constant giving the error, and an accuracy parameter.
    TaylorModel<Interval>(const Expansion<Float>& f, const Float& e, Sweeper swp);
    //! \brief Fast swap with another Taylor model.
    void swap(TaylorModel<Interval>& tm);
    //! \brief The zero element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<Interval> create() const;
    //! \brief The \a j<sup>th</sup> coordinate element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<Interval> create_coordinate(uint j) const;
    //! \brief Set to zero.
    TaylorModel<Interval> create_ball(Float e) const;
    //! \brief Set to zero.
    void clear();

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorModel<Interval>& operator=(double c);
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorModel<Interval>& operator=(const Real& c);
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorModel<Interval>& operator=(const Float& c);
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<Interval>& operator=(const Interval& c);
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct the zero quantity in \a as independent variables.
    static TaylorModel<Interval> zero(uint as, Sweeper swp) {
        TaylorModel<Interval> r(as,swp); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<Interval> constant(uint as, double c, Sweeper swp) {
        return TaylorModel<Interval>::constant(as,Float(c),swp); }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<Interval> constant(uint as, const Float& c, Sweeper swp) {
        TaylorModel<Interval> r(as,swp); r.set_value(c); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<Interval> constant(uint as, const Interval& d, Sweeper swp) {
        TaylorModel<Interval> r(as,swp); r.set_value(1.0); r*=d; return r; }
    //! \brief Construct the quantity with expansion \f$x_j\f$ in \a as independent variables.
    static TaylorModel<Interval> variable(uint as, uint j, Sweeper swp) {
        TaylorModel<Interval> r(as,swp); r.set_gradient(j,1.0); return r; }
    //! \brief Construct the quantity which scales the unit interval into the domain \a d.
    static TaylorModel<Interval> scaling(uint as, uint j, const Interval& d, Sweeper swp) {
        TaylorModel<Interval> r(as,swp); r.set_gradient(j,1.0); r.rescale(Interval(-1,1),d); return r; }
    //! \brief Construct the quantity which scales the codomain \a cd into the unit interval.
    static TaylorModel<Interval> unscaling(uint as, uint j, const Interval& d, Sweeper swp) {
        TaylorModel<Interval> r(as,swp); r.set_gradient(j,1.0); r.rescale(d,Interval(-1,+1)); return r; }
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a e
    static TaylorModel<Interval> error(uint as, Float e, Sweeper swp) {
        TaylorModel<Interval> r(as,swp); r.set_error(e); return r; }

    //! \brief Return the vector of zero variables of size \a rs in \a as arguments.
    static Vector< TaylorModel<Interval> > zeros(uint rs, uint as, Sweeper swp);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector< TaylorModel<Interval> > constants(uint as, const Vector<Float>& c, Sweeper swp);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector< TaylorModel<Interval> > constants(uint as, const Vector<Interval>& c, Sweeper swp);
    //! \brief Return the vector of variables on the unit domain.
    static Vector< TaylorModel<Interval> > variables(uint as, Sweeper swp);
    //! \brief Return the vector scaling the unit interval onto the domain \a d.
    static Vector< TaylorModel<Interval> > scalings(const Vector<Interval>& d, Sweeper swp);
    //! \brief Return the vector scaling the unit interval onto the codomain \a cd.
    static Vector< TaylorModel<Interval> > unscalings(const Vector<Interval>& d, Sweeper swp);
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    bool operator==(const TaylorModel<Interval>& sd) const {
        return this->_expansion==sd._expansion && this->_error == sd._error; }
    //! \brief Inequality operator.
    bool operator!=(const TaylorModel<Interval>& sd) const {
        return !(*this==sd); }
    //! \brief Comparison with another Taylor model.
    tribool operator<(const TaylorModel<Interval>& sd) const {
        return (sd-*this)>0; }
    //! \brief Comparison with another Taylor model.
    tribool operator>(const TaylorModel<Interval>& sd) const {
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
    const Float& error() const { return this->_error; }
    //! \brief A reference to the error of the expansion over the domain.
    Float& error() { return this->_error; }
    //! \brief The constant term in the expansion.
    const Float& value() const { return (*this)[MultiIndex::zero(this->argument_size())]; }
    //! \brief A reference to the constant term in the expansion.
    Float& value() { return (*this)[MultiIndex::zero(this->argument_size())]; }
    //! \brief The coefficient of the gradient term \f$df/dx_j\f$.
    const Float& gradient(uint j) const { return (*this)[MultiIndex::unit(this->argument_size(),j)]; }
    //! \brief A reference to the coefficient of the gradient term \f$df/dx_j\f$.
    Float& gradient(uint j) { return (*this)[MultiIndex::unit(this->argument_size(),j)]; }

    //! \brief The constant term in the expansion.
    Float average() const { return (*this)[MultiIndex::zero(this->argument_size())]; }
    //! \brief The radius of the smallest ball containing the model.
    Float radius() const;
    //! \brief An over-approximation to the supremum norm.
    Float norm() const;
    //! \brief A value \c e such that analytic functions are evaluated to a tolerance of \c e. Equal to the sweep threshold.
    Float tolerance() const;

    //! \brief Set the error of the expansion.
    void set_error(const Float& ne) {
        ARIADNE_ASSERT(ne>=0); this->_error=ne; }
    //! \brief Set the constant term in the expansion.
    void set_value(const Float& c) {
        this->_expansion.set(MultiIndex::zero(this->argument_size()),c,ReverseLexicographicKeyLess()); }
    //! \brief Set the coefficient of the term \f$df/dx_j\f$.
    void set_gradient(uint j, const Float& c) {
        this->_expansion.set(MultiIndex::unit(this->argument_size(),j),c,ReverseLexicographicKeyLess()); }

    //! \brief The coefficient of the term in $x^a$.
    const Float& operator[](const MultiIndex& a) const { return this->_expansion[a]; }
    //! \brief A read/write reference to the coefficient of the term in $x^a$.
    Float& operator[](const MultiIndex& a) { return this->_expansion.at(a,ReverseLexicographicKeyLess()); }

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
    Vector<Interval> domain() const;
    //! \brief An over-approximation to the range of the quantity.
    Interval range() const;
    //! \brief Compute the gradient of the expansion with respect to the \a jth variable over the domain.
    Interval gradient_range(uint j) const;

    //! \brief Evaluate the quantity over the interval of points \a x.
    friend Interval evaluate(const TaylorModel<Interval>&, const Vector<Interval>& x);
    //! \brief Evaluate the quantity over the interval of points \a x.
    friend TaylorModel<Interval> compose(const TaylorModel<Interval>&, const Vector< TaylorModel<Interval> >& x);
    //! \brief Compose two models, where the second is scaled so that the codomain is a unit box.
    friend Vector<Interval> evaluate(const Vector< TaylorModel<Interval> >& f, const Vector<Interval>& x);
    //! \brief Compose two models, where the second is scaled so that the codomain is a unit box.
    friend Vector< TaylorModel<Interval> > compose(const Vector< TaylorModel<Interval> >& f, const Vector< TaylorModel<Interval> >& g);
    //@}

    //@{
    /*! \name Inplace modifications. */
    // TODO: Change these to return void
    //! \brief Scale so that the old codomain maps into the new codomain.
    TaylorModel<Interval>& rescale(const Interval& old_codomain, const Interval& new_codomain);
    //! \brief Restrict to a subdomain.
    TaylorModel<Interval>& restrict(const Vector<Interval>& new_domain);
    //! \brief Compute the antiderivative (in place).
    TaylorModel<Interval>& antidifferentiate(uint k);
    //@}

    //@{
    /*! \name Set-based operations. */
    //! \brief Test if one model refines (is a subset of) another.
    friend bool refines(const TaylorModel<Interval>& tm1, const TaylorModel<Interval>& tm2);
    //! \brief Test if one model is disjoint from (is incompatible with) another.
    friend bool disjoint(const TaylorModel<Interval>& tm1, const TaylorModel<Interval>& tm2);
    //! \brief An over-approximation of the intersection of the sets of functions
    //! allowed by the two models. It is guaranteed that any function represented
    //! by both models is also represented by the result.
    friend TaylorModel<Interval> intersection(const TaylorModel<Interval>& tm1, const TaylorModel<Interval>& tm2);
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorModel<Interval>& sweep();
    //! \brief Remove all terms whose coefficient has magnitude less than \a eps.
    TaylorModel<Interval>& sweep(const SweeperInterface& accuracy);
    //! \brief Set the error to zero.
    //! WARNING: This method does not preserve rigour of the model approximation.
    TaylorModel<Interval>& clobber();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel<Interval>& unique_sort();
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
    void iadd(const Interval& c);
    //! \brief Multiply by a numerical scalar \c r*=c .
    void imul(const Interval& c);
    //! \brief Scalar multiply and add \c r+=c*x .
    void isma(const Interval& c, const TaylorModel<Interval>& x);
    //! \brief Fused multiply and add \c r+=x1*x2 .
    void ifma(const TaylorModel<Interval>& x1, const TaylorModel<Interval>& x2);
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const TaylorModel<Interval>& x);
    //@}

  public:
    std::ostream& str(const std::ostream&) const;
    std::ostream& repr(const std::ostream&) const;

};

// Rescale the vector x from the domain d to the unit domain.
Vector<Interval> unscale(const Vector<Interval>& x, const Vector<Interval>& d);

//! \relates TaylorModel<Interval> \brief The magnitude of the variable
Float mag(const TaylorModel<Interval>& tm);
//! \relates TaylorModel<Interval> \brief Split the variable over two domains, subdividing along the independent variable j.
std::pair< TaylorModel<Interval>, TaylorModel<Interval> > split(const TaylorModel<Interval>& x, uint j);
//! \relates TaylorModel<Interval>
//!\brief Split the variable, subdividing along the independent variable j
//! and taking the lower/middle/upper half depending on whether half is false, indeterminate or true.
TaylorModel<Interval> split(const TaylorModel<Interval>& x, uint j, tribool half);

//! \relates TaylorModel<Interval> \brief Scale the variable by post-composing with an affine map taking the unit interval to ivl.
TaylorModel<Interval> scale(const TaylorModel<Interval>& x, const Interval& ivl);
//! \relates TaylorModel<Interval> \brief Scale the variable by post-composing with an affine map taking the interval ivl to the unit interval
TaylorModel<Interval> unscale(const TaylorModel<Interval>& x, const Interval& ivl);
//! \relates TaylorModel<Interval> \brief Scale the variable by post-composing with an affine map taking the interval ivl1 to interval \a ivl2
TaylorModel<Interval> rescale(const TaylorModel<Interval>& x, const Interval& ivl1, const Interval& ivl2);

//! \relates TaylorModel<Interval> \brief Evaluate an Array of Taylor variables on a vector.
Interval evaluate(const TaylorModel<Interval>& x, const Vector<Interval>& sy);
//! \relates TaylorModel<Interval> \brief Evaluate an Array of Taylor variables on a vector.
TaylorModel<Interval> partial_evaluate(const TaylorModel<Interval>& x, uint k, Float c);
//! \relates TaylorModel<Interval> \brief Evaluate an Array of Taylor variables on a vector.
TaylorModel<Interval> partial_evaluate(const TaylorModel<Interval>& x, uint k, Interval c);
//! \relates TaylorModel<Interval>
//! Substitute the TaylorModel<Interval> y in the  kth variable of \a x.
//! Precondition: x.argument_size()==y.argument_size()+1
TaylorModel<Interval> substitute(const TaylorModel<Interval>& x, uint k, const TaylorModel<Interval>& y);

//! \relates TaylorModel<Interval> \brief Embed the model in a space of higher dimension
TaylorModel<Interval> embed(const TaylorModel<Interval>& tm, uint d);
//! \relates TaylorModel<Interval> \brief Embed the model in a space of higher dimension, placing the error in variable i.
TaylorModel<Interval> embed_error(const TaylorModel<Interval>& tm, uint d, uint i);
//! \relates TaylorModel<Interval> \brief Embed the model in a space of higher dimension
TaylorModel<Interval> embed(uint as, const TaylorModel<Interval>& tv);

//! \relates TaylorModel<Interval> \brief Test if a model refines another
bool refines(const TaylorModel<Interval>& tv1, const TaylorModel<Interval>& tv2);
//! \relates TaylorModel<Interval> \brief Test if a model is disjoint from
bool disjoint(const TaylorModel<Interval>& tv1, const TaylorModel<Interval>& tv2);

//! \relates TaylorModel<Interval> \brief Antidifferentiation operator
TaylorModel<Interval> antiderivative(const TaylorModel<Interval>& x, uint k);

//! \relates TaylorModel<Interval> \brief Differentiation operator; discards error term
TaylorModel<Interval> derivative(const TaylorModel<Interval>& x, uint k);

//! \relates TaylorModel<Interval> \brief Replace the variale x[k] with a*x[k]+b
TaylorModel<Interval> preaffine(const TaylorModel<Interval>&, uint k, const Interval& a, const Interval& b);
//! \relates TaylorModel<Interval> \brief Restricts the range of the variable x[k] to the interval d.
//! \pre -1 <= d.lower() <= d.upper() <= 1 .
TaylorModel<Interval> restrict(const TaylorModel<Interval>&, uint k, const Interval& d);

//! \brief Abstract away the given variables.
//! For example, the model tm(x0,x1,x2,x3) becomes tm'(x0,x1)=tm(x0,[-1,+1],x1,[-1,+1]) on discarding x1 and x3.
TaylorModel<Interval>  discard(const TaylorModel<Interval>&, const Array<uint>& variables);

TaylorModel<Interval> recondition(const TaylorModel<Interval>& tm, Array<uint>& discarded_variables,
                                  uint number_of_error_variables);

TaylorModel<Interval> recondition(const TaylorModel<Interval>& tm, Array<uint>& discarded_variables,
                                  uint number_of_error_variables, uint index_of_error);

//! \relates TaylorModel<Interval>
//! An over-approximation to the intersection of two Taylor models.
//! Since the intersection cannot be represented exactly in the class of
//! TaylorModels, truncation errors as well as roundoff errors may be present.
//! In the absence of roundoff errors, the result is a subset of both arguments,
//! and is guaranteed to contain any function contained in both arguments.
TaylorModel<Interval> intersection(const TaylorModel<Interval>& x1, const TaylorModel<Interval>& x2);

// Compose an Array of Taylor variables with another, assuming that y has been scaled to have unit codomain
TaylorModel<Interval> compose(const TaylorModel<Interval>& x, const Vector< TaylorModel<Interval> >& y);

// Compose an Array of Taylor variables with another, after scaling by the interval vectors
TaylorModel<Interval> compose(const TaylorModel<Interval>& x, const Vector<Interval>& bx, const Vector< TaylorModel<Interval> >& y);

Float norm(const TaylorModel<Interval>& tm);
Float norm(const Vector< TaylorModel<Interval> >& tv);

TaylorModel<Interval> max(const TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
TaylorModel<Interval> min(const TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
TaylorModel<Interval> abs(const TaylorModel<Interval>& x);



// Vector operations which can be evaluated componentwise
bool refines(const Vector< TaylorModel<Interval> >&,const Vector< TaylorModel<Interval> >&);
bool disjoint(const Vector< TaylorModel<Interval> >&,const Vector< TaylorModel<Interval> >&);
std::pair< Vector< TaylorModel<Interval> >, Vector< TaylorModel<Interval> > > split(const Vector< TaylorModel<Interval> >& x, uint j);
Vector< TaylorModel<Interval> > split(const Vector< TaylorModel<Interval> >& x, uint j, bool half);
Vector< TaylorModel<Interval> > unscale(const Vector< TaylorModel<Interval> >& x, const Vector<Interval>& bx);
Vector< TaylorModel<Interval> > scale(const Vector< TaylorModel<Interval> >& x, const Vector<Interval>& bx);
Vector<Interval> evaluate(const Vector< TaylorModel<Interval> >& x, const Vector<Interval>& sy);
Vector< TaylorModel<Interval> > partial_evaluate(const Vector< TaylorModel<Interval> >& x, uint k, Float sy);
Vector< TaylorModel<Interval> > partial_evaluate(const Vector< TaylorModel<Interval> >& x, uint k, Interval sy);
Vector< TaylorModel<Interval> > substitute(const Vector< TaylorModel<Interval> >& x, uint k, const TaylorModel<Interval>& y);
Vector< TaylorModel<Interval> > antiderivative(const Vector< TaylorModel<Interval> >& x, uint k);
Vector< TaylorModel<Interval> > embed(const Vector< TaylorModel<Interval> >& x, uint as);
Vector< TaylorModel<Interval> > embed(uint as, const Vector< TaylorModel<Interval> >& x);
Matrix<Interval> jacobian(const Vector< TaylorModel<Interval> >& x, const Vector<Interval>& d);
//Matrix<Interval> jacobian(const Vector< TaylorModel<Interval> >& x);
bool refines(const Vector< TaylorModel<Interval> >& x1, const Vector< TaylorModel<Interval> >& x2);
Vector< TaylorModel<Interval> > combine(const Vector< TaylorModel<Interval> >& x1, const Vector< TaylorModel<Interval> >& x2);
Vector< TaylorModel<Interval> > combine(const Vector< TaylorModel<Interval> >& x1, const TaylorModel<Interval>& x2);
Vector< TaylorModel<Interval> > combine(const TaylorModel<Interval>& x1, const Vector< TaylorModel<Interval> >& x2);
Vector< TaylorModel<Interval> > combine(const TaylorModel<Interval>& x1, const TaylorModel<Interval>& x2);
Vector< TaylorModel<Interval> > compose(const Vector< TaylorModel<Interval> >& f, const Vector< TaylorModel<Interval> >& g);
Vector< TaylorModel<Interval> > compose(const Vector< TaylorModel<Interval> >& f, const Vector<Interval>& d, const Vector< TaylorModel<Interval> >& g);

TaylorModel<Interval> unchecked_compose(const TaylorModel<Interval>& x, const Vector< TaylorModel<Interval> >& y);
Vector< TaylorModel<Interval> > unchecked_compose(const Vector< TaylorModel<Interval> >& x, const Vector< TaylorModel<Interval> >& y);
Vector< TaylorModel<Interval> > unchecked_compose(const Vector< TaylorModel<Interval> >& x, const Vector<Interval>& d, const Vector< TaylorModel<Interval> >& y);

Vector< TaylorModel<Interval> > operator*(const Matrix<Float>& A, const Vector< TaylorModel<Interval> >& x);
Vector< TaylorModel<Interval> > operator*(const Matrix<Interval>& A, const Vector< TaylorModel<Interval> >& x);


/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, Polynomila, TaylorModel<Interval><Interval>.
 */
template<>
class TaylorModel<Float>
    : public NormedAlgebraMixin<TaylorModel<Float>,Float>
{
    typedef Expansion<Float> ExpansionType;
  private:
    ExpansionType _expansion;
    mutable Sweeper _sweeper;
  private:
    static const Float _zero;

  public:
    //! \brief The type used for the coefficients.
    typedef Float NumericType;
    //! \brief The type used to index the coefficients.
    typedef MultiIndex IndexType;
    //! \brief The type used for the coefficients.
    typedef Float ValueType;

    //! \brief An iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::iterator iterator;
    //! \brief A constant iterator through the (index,coefficient) pairs of the expansion.
    typedef ExpansionType::const_iterator const_iterator;

  public:
    //@{
    /*! \name Constructors and destructors. */
    //! \brief Construct a IntervalTaylorModel in \a as arguments.
    TaylorModel<Float>(uint as = 0u);
    TaylorModel<Float>(uint as, Sweeper swp);

    TaylorModel<Float> create() const { return TaylorModel<Float>(this->argument_size(),this->_sweeper); }
    TaylorModel<Float> create_ball(Float r) const { return TaylorModel<Float>(this->argument_size(),this->_sweeper); }
    //! \brief Fast swap with another Taylor model.
    void swap(TaylorModel<Float>& other) { this->_expansion.swap(other._expansion); }
    //! \brief Set to zero.
    void clear() { this->_expansion.clear(); }
    //! \brief A zero element with the same parameters.
    TaylorModel<Float> null() const { return TaylorModel<Float>(this->argument_size()); }
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a built-in, keeping the same number of arguments.
    TaylorModel<Float>& operator=(double c) {
        this->_expansion.clear(); this->_expansion.append(MultiIndex(this->argument_size()),Float(c)); return *this; }
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorModel<Float>& operator=(const Float& c) {
        this->_expansion.clear(); this->_expansion.append(MultiIndex(this->argument_size()),c); return *this; }
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<Float>& operator=(const Interval& c) { return (*this)=midpoint(c); }
    //! \brief Set equal to a real constant, keeping the same number of arguments.
    TaylorModel<Float>& operator=(const Real& c) { return (*this)=Float(c); }
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    bool operator==(const TaylorModel<Float>& other) const {
        return this->_expansion==other._expansion; }
    //! \brief Inequality operator.
    bool operator!=(const TaylorModel<Float>& other) const {
        return !(*this==other); }
    //@}

    //@{
    /*! \name Data access */
    //! \brief The number of variables in the argument of the quantity.
    uint argument_size() const { return this->_expansion.argument_size(); }
    //! \brief The coefficient of the term in $x^a$.
    const Float& operator[](const MultiIndex& a) const { return this->_expansion[a]; }
    //! \brief The error of the expansion over the domain.
    const Float& error() const { return _zero; }
    //@}

    //@{
    /*! \name Function evaluation. */
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
    TaylorModel<Float>& sweep();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel<Float>& unique_sort();
    //@}

  public:
    //@{
    /*! \name Standard algebra interface. */
    //! \brief An approximation to the norm of the function.
    virtual Float norm() const;
    //! \brief An approximation to the average value of the function.
    virtual Float average() const;
    //! \brief The tolerance to which analytic functions should be computed.
    virtual Float tolerance() const;
    //! \brief The radius of the ball containing the functions.
    virtual Float radius() const;
    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream&) const;
    //! \brief Inplace addition of a scalar constant.
    virtual void iadd(const Float& c);
    //! \brief Inplace multiplication of a scalar constant.
    virtual void imul(const Float& c);
    //! \brief Inplace addition of a scalar multiple of a Taylor model.
    virtual void isma(const Float& c, const TaylorModel<Float>& x);
    //! \brief Inplace addition of a product of Taylor models.
    virtual void ifma(const TaylorModel<Float>& x1, const TaylorModel<Float>& x2);


    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const TaylorModel<Float>& x);
    //@}

  public:
    std::ostream& str(std::ostream&) const;
    std::ostream& repr(std::ostream&) const;
};


inline std::ostream& operator<<(std::ostream& os, const TaylorModel<Float>& x) {
    x.str(os); return os; }

inline Vector<Interval> codomain(const Vector< TaylorModel<Float> >& t) {
    Vector<Interval> r(t.size()); for(uint i=0; i!=t.size(); ++i) { r[i]=t[i].codomain(); } return r; }

Vector< TaylorModel<Float> > operator*(const Matrix<Float>& A, const Vector< TaylorModel<Float> >& x);



} // namespace Ariadne

#endif // ARIADNE_TAYLOR_MODEL_H
