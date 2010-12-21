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

class IntersectionException;

struct IntersectionException : public std::runtime_error {
    IntersectionException(const std::string& what) : std::runtime_error(what) { }
};

class TaylorModelAccuracy
{
    friend class TaylorModel<Interval>;
    friend class TaylorModel<Float>;

    TaylorModelAccuracy();

    friend TaylorModelAccuracy max(const TaylorModelAccuracy& acc1, const TaylorModelAccuracy& acc2);
    friend TaylorModelAccuracy min(const TaylorModelAccuracy& acc1, const TaylorModelAccuracy& acc2);
    friend bool operator==(const TaylorModelAccuracy& acc1, const TaylorModelAccuracy& acc2);
    friend std::ostream& operator<<(std::ostream& os, const TaylorModelAccuracy& acc);
  public:
    TaylorModelAccuracy(double st, uint md);
  public:
    static void set_default_sweep_threshold(double dst) { ARIADNE_ASSERT(dst>=0.0); _default_sweep_threshold=dst; }
    static void set_default_maximum_degree(int dmd) { ARIADNE_ASSERT(dmd>=0); _default_maximum_degree=dmd; }
  public:
    bool discard(const Float& x) const;
    bool discard(const MultiIndex& a) const;
    bool discard(const MultiIndex& a, const Float& x) const;
  private:
    static double _default_sweep_threshold;
    static uint _default_maximum_degree;
  private:
    double _sweep_threshold;
    uint _maximum_degree;
};

/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, ScalarTaylorFunction, VectorTaylorFunction, TaylorConstrainedImageSet.
 */
template<>
class TaylorModel<Interval>
{
    friend class ScalarTaylorFunction;
    friend class VectorTaylorFunction;
    typedef Expansion<Float> ExpansionType;
    typedef ReverseLexicographicKeyLess ComparisonType;
  public:
    typedef TaylorModelAccuracy Accuracy;
    typedef Interval NumericType;
  private:
    ExpansionType _expansion;
    Float _error;
    mutable shared_ptr<Accuracy> _accuracy_ptr;
  private:
    static const double _zero;
    static double _default_sweep_threshold;
    static uint _default_maximum_degree;
  public:
    const Accuracy& accuracy() const { return *this->_accuracy_ptr; }
    shared_ptr<Accuracy> accuracy_ptr() const { return this->_accuracy_ptr; }
    void set_accuracy(shared_ptr<Accuracy> acc) { this->_accuracy_ptr=acc; }

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
    //! \brief Construct a TaylorModel<Interval> in \a as arguments.
    TaylorModel<Interval>(uint as);
    //! \brief Construct a TaylorModel<Interval> in \a as arguments with the given accuracy control.
    TaylorModel<Interval>(uint as, shared_ptr<Accuracy> acc);
    //! \brief Construct from a map giving the expansion, a constant giving the error, and an accuracy parameter.
    TaylorModel<Interval>(const std::map<MultiIndex,Float>& d, const Float& e);
    //! \brief Construct from a map giving the expansion, and a constant giving the error.
    TaylorModel<Interval>(const Expansion<Float>& f, const Float& e=0.0);
    //! \brief Construct from a map giving the expansion, a constant giving the error, and an accuracy parameter.
    TaylorModel<Interval>(const Expansion<Float>& f, const Float& e, shared_ptr<Accuracy> a);
    //! \brief Fast swap with another Taylor model.
    void swap(TaylorModel<Interval>& tm);
    //! \brief Set to zero.
    void clear();

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a builtin, keeping the same number of arguments.
    TaylorModel<Interval>& operator=(double c);
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorModel<Interval>& operator=(const Real& c);
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorModel<Interval>& operator=(const Float& c);
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<Interval>& operator=(const Interval& c);
    //! \brief Test if the quantity is a better approximation than \a t throughout the domain.
    bool refines(const TaylorModel<Interval>& t);
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct the zero quantity in \a as independent variables.
    static TaylorModel<Interval> zero(uint as) {
        TaylorModel<Interval> r(as); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<Interval> constant(uint as, double c) {
        TaylorModel<Interval> r(as); r.set_value(static_cast<Float>(c)); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<Interval> constant(uint as, const Float& c) {
        TaylorModel<Interval> r(as); r.set_value(c); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<Interval> constant(uint as, const Interval& d) {
        TaylorModel<Interval> r(as); r.set_value(1.0); r*=d; return r; }
    //! \brief Construct the quantity with expansion \f$x_j\f$ in \a as independent variables.
    static TaylorModel<Interval> variable(uint as, uint j) {
        TaylorModel<Interval> r(as); r.set_gradient(j,1.0); return r; }
    //! \brief Construct the quantity which scales the unit interval into the domain \a d.
    static TaylorModel<Interval> scaling(uint as, uint j, const Interval& d) {
        TaylorModel<Interval> r(as); r.set_gradient(j,1.0); r.rescale(Interval(-1,1),d); return r; }
    //! \brief Construct the quantity which scales the codomain \a cd into the unit interval.
    static TaylorModel<Interval> unscaling(uint as, uint j, const Interval& d) {
        TaylorModel<Interval> r(as); r.set_gradient(j,1.0); r.rescale(d,Interval(-1,+1)); return r; }
    //! \brief Construct the quantity which scales the interval \a cd onto the interval \a d.
    static TaylorModel<Interval> rescaling(uint as, uint j, const Interval& cd, const Interval& d) {
        TaylorModel<Interval> r(as); r.set_gradient(j,1.0); r.rescale(cd,d); return r; }
    //! \brief Construct the quantity \f$c+\sum g_jx_j\f$.
    static TaylorModel<Interval> affine(const Float& c, const Vector<Float>& g) {
        TaylorModel<Interval> r(g.size()); r.set_value(c);
        for(uint j=0; j!=g.size(); ++j) { r.set_gradient(j,g[j]); } return r; }
    //! \brief Construct the quantity \f$c+\sum g_jx_j \pm e\f$.
    static TaylorModel<Interval> affine(const Float& x, const Vector<Float>& g, const Float& e) {
        TaylorModel<Interval> r(g.size()); r.set_value(x); r.set_error(e);
        for(uint j=0; j!=g.size(); ++j) { r.set_gradient(j,g[j]); } return r; }

    //! \brief Return the vector of zero variables of size \a rs in \a as arguments.
    static Vector< TaylorModel<Interval> > zeros(uint rs, uint as);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector< TaylorModel<Interval> > constants(uint as, const Vector<Float>& c);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector< TaylorModel<Interval> > constants(uint as, const Vector<Interval>& c);
    //! \brief Return the vector of variables on the unit domain.
    static Vector< TaylorModel<Interval> > variables(uint as);
    //! \brief Return the vector scaling the unit interval onto the domain \a d.
    static Vector< TaylorModel<Interval> > scalings(const Vector<Interval>& d);
    //! \brief Return the vector scaling the unit interval onto the codomain \a cd.
    static Vector< TaylorModel<Interval> > unscalings(const Vector<Interval>& d);
    //! \brief Return the vector scaling the codomain \a cd onto the domain \a d.
    static Vector< TaylorModel<Interval> > rescalings(const Vector<Interval>& cd, const Vector<Interval>& d);
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

    //! \brief The coefficient of the term \f$df/dx_j\f$.
    const Float& operator[](uint j) const {
        return (*this)[MultiIndex::unit(this->argument_size(),j)]; }
    //! \brief A reference to the coefficient of the term \f$df/dx_j\f$.
    Float& operator[](uint j) {
        return (*this)[MultiIndex::unit(this->argument_size(),j)]; }

    //! \brief An iterator to the first term in the expansion.
    iterator begin() { return this->_expansion.begin(); }
    //! \brief A constant iterator to the first term in the expansion.
    const_iterator begin() const { return this->_expansion.begin(); }
    //! \brief An iterator to the end of the expansion.
    iterator end() { return this->_expansion.end(); }
    //! \brief A constant iterator to the end of the expansion.
    const_iterator end() const { return this->_expansion.end(); }
    //! \brief An iterator to the term with index \a a.
    iterator find(const MultiIndex& a) { return this->_expansion.find(a); }
    //! \brief A constant iterator to the term with index \a a.
    const_iterator find(const MultiIndex& a) const { return this->_expansion.find(a); }

    //! \brief The number of variables in the argument of the quantity.
    uint argument_size() const { return this->_expansion.argument_size(); }
    //! \brief The maximum degree of terms in the expansion.
    uint degree() const;
    //! \brief The number of nonzero terms in the expansion.
    uint number_of_nonzeros() const { return this->_expansion.number_of_nonzeros(); }
    //@}

    //@{
    /*! \name Function evaluation. */
    //! \brief An over-approximation to the supremum norm.
    Float norm() const;
    //! \brief The domain of the quantity, always given by \f$[-1,1]^{\mathrm{as}}\f$.
    Vector<Interval> domain() const;
    //! \brief An over-approximation to the range of the quantity.
    Interval range() const;
    //! \brief Compute the gradient of the expansion with respect to the \a jth variable over the domain.
    Interval gradient_range(uint j) const;
    //! \brief Evaluate the quantity at the point \a x.
    Interval evaluate(const Vector<Float>& x) const;
    //! \brief Evaluate the quantity over the interval of points \a x.
    Interval evaluate(const Vector<Interval>& x) const;
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
    //! \brief The supremum norm of the model, given by \f$|e|+\sum_\alpha |c_\alpha|\f$.
    friend Float norm(const TaylorModel<Interval>& tm);
    //@{
    /*! \name Vectoral function operators. */
    //! \brief Compose two models, where the second is scaled so that the codomain is a unit box.
    friend Vector< TaylorModel<Interval> > compose(const Vector< TaylorModel<Interval> >& f, const Vector< TaylorModel<Interval> >& g);
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Truncate to the default maximum degree of the quantity.
    TaylorModel<Interval>& truncate();
    //! \brief Truncate to degree \a deg.
    TaylorModel<Interval>& truncate(uint deg);
    //! \brief Truncate all terms with any coefficient higher than \a a.
    TaylorModel<Interval>& truncate(const MultiIndex& a);
    //! \brief Truncate all terms with any coefficient higher than those given by \a a.
    TaylorModel<Interval>& truncate(const MultiIndexBound& a);
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorModel<Interval>& sweep();
    //! \brief Remove all terms whose coefficient has magnitude less than \a eps.
    TaylorModel<Interval>& sweep(double eps);
    //! \brief Remove all terms whose degree is higher than \a deg or
    //! whose coefficient has magnitude less than \a eps.
    TaylorModel<Interval>& clean(const Accuracy& accuracy);
    //! \brief Remove all terms which have high degree or small magnitude.
    TaylorModel<Interval>& clean();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel<Interval>& unique_sort();
    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \brief Specify a bound on the terms which may be present in the expansion.
    void set_maximum_index(MultiIndexBound ma);
    //! \brief Specify the maximum degree \a md for terms which may be present in the expansion.
    void set_maximum_degree(uint md);
    //! \brief Specify the minimum absolute value \a me for coefficients of terms which may be present in the expansion.
    void set_sweep_threshold(double me);
    //! \brief The maximum index of terms which may be present in the expansion.
    //! Any term with index \f$a\not\leq a_{\max}\f$ will be assimilated into the error term when truncate() or clean() are called.
    MultiIndexBound maximum_index() const;
    //! \brief The maximum degree for terms which may be present in the expansion.
    //! Any term with degree \f$d>d_{\max}\f$ will be assimilated into the error term when truncate() or clean() are called.
    uint maximum_degree() const;
     //! \brief The minimum absolute value for coefficients of terms which may be present in the expansion.
    //! Any term with coefficient \f$c\f$ with \f$|c|<e_{\max}\f$ will be assimilated into the error term when sweep() or clean() are called.
   double sweep_threshold() const;
    //@}

    //@{
    /*! \name Arithmetic operations. */
    //! \brief Inplace addition of another variable.
    friend TaylorModel<Interval>& operator+=(TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
    //! \brief Inplace subtraction of another variable.
    friend TaylorModel<Interval>& operator-=(TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
    //! \brief Inplace multiplication of another variable. Not any more efficient than ordinary multiplication.
    friend TaylorModel<Interval>& operator*=(TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
    //! \brief Inplace division of another variable. Not any more efficient than ordinary division.
    friend TaylorModel<Interval>& operator/=(TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
    //! \brief Inplace addition of a product of two variables.
    friend TaylorModel<Interval>& operator+=(TaylorModel<Interval>& x, const Product< TaylorModel<Interval>, TaylorModel<Interval> >& y);
    //! \brief Inplace addition of a built-in floating-point constant.
    friend TaylorModel<Interval>& operator+=(TaylorModel<Interval>& x, double c);
    //! \brief Inplace subtraction of a built-in floating-point constant.
    friend TaylorModel<Interval>& operator-=(TaylorModel<Interval>& x, double c);
    //! \brief Inplace multiplication by a builting scalar.
    friend TaylorModel<Interval>& operator*=(TaylorModel<Interval>& x, double c);
    //! \brief Inplace division by an built-in scalar.
    friend TaylorModel<Interval>& operator/=(TaylorModel<Interval>& x, double c);
    //! \brief Inplace addition of an exact real number.
    friend TaylorModel<Interval>& operator+=(TaylorModel<Interval>& x, const Real& c);
    //! \brief Inplace subtraction of an exact real number.
    friend TaylorModel<Interval>& operator-=(TaylorModel<Interval>& x, const Real& c);
    //! \brief Inplace multiplication by an exact real number.
    friend TaylorModel<Interval>& operator*=(TaylorModel<Interval>& x, const Real& c);
    //! \brief Inplace division by an exact real number.
    friend TaylorModel<Interval>& operator/=(TaylorModel<Interval>& x, const Real& c);
    //! \brief Inplace addition of an exact floating-point constant.
    friend TaylorModel<Interval>& operator+=(TaylorModel<Interval>& x, const Float& c);
    //! \brief Inplace subtraction of an exact floating-point constant.
    friend TaylorModel<Interval>& operator-=(TaylorModel<Interval>& x, const Float& c);
    //! \brief Inplace multiplication by an exact scalar.
    friend TaylorModel<Interval>& operator*=(TaylorModel<Interval>& x, const Float& c);
    //! \brief Inplace division by an exact scalar.
    friend TaylorModel<Interval>& operator/=(TaylorModel<Interval>& x, const Float& c);
    //! \brief Inplace addition of an interval constant.
    friend TaylorModel<Interval>& operator+=(TaylorModel<Interval>& x, const Interval& c);
    //! \brief Inplace subtraction of an interval constant.
    friend TaylorModel<Interval>& operator-=(TaylorModel<Interval>& x, const Interval& c);
    //! \brief Inplace multiplication by an approximate scalar.
    friend TaylorModel<Interval>& operator*=(TaylorModel<Interval>& x, const Interval& c);
    //! \brief Inplace division by an approximate scalar.
    friend TaylorModel<Interval>& operator/=(TaylorModel<Interval>& x, const Interval& c);

    //! \brief Unary plus.
    friend TaylorModel<Interval> operator+(const TaylorModel<Interval>& x);
    //! \brief Unary minus.
    friend TaylorModel<Interval> operator-(const TaylorModel<Interval>& x);
    //! \brief Addition.
    friend TaylorModel<Interval> operator+(const TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
    //! \brief Subtraction.
    friend TaylorModel<Interval> operator-(const TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
    //! \brief Multiplication.
    friend TaylorModel<Interval> operator*(const TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
    //! \brief Division.
    friend TaylorModel<Interval> operator/(const TaylorModel<Interval>& x, const TaylorModel<Interval>& y);

    //! \brief Addition of a scalar.
    friend TaylorModel<Interval> operator+(const TaylorModel<Interval>& x, double c);
    //! \brief Subtraction of a scalar.
    friend TaylorModel<Interval> operator-(const TaylorModel<Interval>& x, double c);
    //! \brief Multiplication by a scalar.
    friend TaylorModel<Interval> operator*(const TaylorModel<Interval>& x, double c);
    //! \brief Division by a scalar.
    friend TaylorModel<Interval> operator/(const TaylorModel<Interval>& x, double c);
    //! \brief Addition of a scalar.
    friend TaylorModel<Interval> operator+(const TaylorModel<Interval>& x, const Real& c);
    //! \brief Addition of a scalar.
    friend TaylorModel<Interval> operator-(const TaylorModel<Interval>& x, const Real& c);
    //! \brief Addition of a scalar.
    friend TaylorModel<Interval> operator*(const TaylorModel<Interval>& x, const Real& c);
    //! \brief Addition of a scalar.
    friend TaylorModel<Interval> operator/(const TaylorModel<Interval>& x, const Real& c);
    //! \brief Addition of a scalar.
    friend TaylorModel<Interval> operator+(const TaylorModel<Interval>& x, const Float& c);
    //! \brief Subtraction of a scalar.
    friend TaylorModel<Interval> operator-(const TaylorModel<Interval>& x, const Float& c);
    //! \brief Multiplication by a scalar.
    friend TaylorModel<Interval> operator*(const TaylorModel<Interval>& x, const Float& c);
    //! \brief Division by a scalar.
    friend TaylorModel<Interval> operator/(const TaylorModel<Interval>& x, const Float& c);
    //! \brief Addition of a scalar.
    friend TaylorModel<Interval> operator+(const TaylorModel<Interval>& x, const Interval& c);
    //! \brief Subtraction of a scalar.
    friend TaylorModel<Interval> operator-(const TaylorModel<Interval>& x, const Interval& c);
    //! \brief Multiplication by a scalar.
    friend TaylorModel<Interval> operator*(const TaylorModel<Interval>& x, const Interval& c);
    //! \brief Division by a scalar.
    friend TaylorModel<Interval> operator/(const TaylorModel<Interval>& x, const Interval& c);
    //! \brief Addition of a built-in scalar.
    friend TaylorModel<Interval> operator+(double c, const TaylorModel<Interval>& x);
    //! \brief Subtraction from a built-in scalar.
    friend TaylorModel<Interval> operator-(double c, const TaylorModel<Interval>& x);
    //! \brief Multiplication by a built-in scalar.
    friend TaylorModel<Interval> operator*(double c, const TaylorModel<Interval>& x);
    //! \brief Division through a built-in scalar.
    friend TaylorModel<Interval> operator/(double c, const TaylorModel<Interval>& x);
    //! \brief Addition of a scalar.
    friend TaylorModel<Interval> operator+(const Real& c, const TaylorModel<Interval>& x);
    //! \brief Subtraction from a scalar.
    friend TaylorModel<Interval> operator-(const Real& c, const TaylorModel<Interval>& x);
    //! \brief Multiplication by a scalar.
    friend TaylorModel<Interval> operator*(const Real& c, const TaylorModel<Interval>& x);
    //! \brief Division through a scalar.
    friend TaylorModel<Interval> operator/(const Real& c, const TaylorModel<Interval>& x);
    //! \brief Addition of a scalar.
    friend TaylorModel<Interval> operator+(const Float& c, const TaylorModel<Interval>& x);
    //! \brief Subtraction from a scalar.
    friend TaylorModel<Interval> operator-(const Float& c, const TaylorModel<Interval>& x);
    //! \brief Multiplication by a scalar.
    friend TaylorModel<Interval> operator*(const Float& c, const TaylorModel<Interval>& x);
    //! \brief Division through a scalar.
    friend TaylorModel<Interval> operator/(const Float& c, const TaylorModel<Interval>& x);
    //! \brief Addition of a scalar.
    friend TaylorModel<Interval> operator+(const Interval& c, const TaylorModel<Interval>& x);
    //! \brief Subtraction from a scalar.
    friend TaylorModel<Interval> operator-(const Interval& c, const TaylorModel<Interval>& x);
    //! \brief Multiplication by a scalar.
    friend TaylorModel<Interval> operator*(const Interval& c, const TaylorModel<Interval>& x);
    //! \brief Division through a scalar.
    friend TaylorModel<Interval> operator/(const Interval& c, const TaylorModel<Interval>& x);
    //@}

    //@{
    /*! \name Algebraic and transcendental functions. */

    //! \brief Maximum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend TaylorModel<Interval> max(const TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
    //! \brief Minimum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend TaylorModel<Interval> min(const TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
    //! \brief Addition.
    friend TaylorModel<Interval> add(const TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
    //! \brief Multiplication.
    friend TaylorModel<Interval> mul(const TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
    //! \brief Absolute value. Throws an error if one variable is neither greater than
    //! zero over the entire domain, nor less than zero over the entire domain.
    friend TaylorModel<Interval> abs(const TaylorModel<Interval>& x);
    //! \brief Negation.
    friend TaylorModel<Interval> neg(const TaylorModel<Interval>& x);
    //! \brief Reciprocal.
    friend TaylorModel<Interval> rec(const TaylorModel<Interval>& x);
    //! \brief Square.
    friend TaylorModel<Interval> sqr(const TaylorModel<Interval>& x);
    //! \brief Power.
    friend TaylorModel<Interval> pow(const TaylorModel<Interval>& x, int n);
    //! \brief Square root.
    friend TaylorModel<Interval> sqrt(const TaylorModel<Interval>& x);
    //! \brief Natural exponent.
    friend TaylorModel<Interval> exp(const TaylorModel<Interval>& x);
    //! \brief Natural logarithm.
    friend TaylorModel<Interval> log(const TaylorModel<Interval>& x);
    //! \brief Sine.
    friend TaylorModel<Interval> sin(const TaylorModel<Interval>& x);
    //! \brief Cosine.
    friend TaylorModel<Interval> cos(const TaylorModel<Interval>& x);
    //! \brief Tangent.
    friend TaylorModel<Interval> tan(const TaylorModel<Interval>& x);
    //! \brief Inverse sine.
    friend TaylorModel<Interval> asin(const TaylorModel<Interval>& x);
    //! \brief Inverse cosine.
    friend TaylorModel<Interval> acos(const TaylorModel<Interval>& x);
    //! \brief Inverse tangent.
    friend TaylorModel<Interval> atan(const TaylorModel<Interval>& x);
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const TaylorModel<Interval>& x);
    //@}

  public:
    std::ostream& str(const std::ostream&) const;
    std::ostream& repr(const std::ostream&) const;

    TaylorModel<Interval>& clobber();
    TaylorModel<Interval>& clobber(uint o);
    TaylorModel<Interval>& clobber(uint so, uint to);
};

// Rescale the vector x from the domain d to the unit domain.
Vector<Interval> unscale(const Vector<Interval>& x, const Vector<Interval>& d);

//! \relates TaylorModel<Interval> \brief The magnitude of the variable
Float mag(const TaylorModel<Interval>& tm);
//! \relates TaylorModel<Interval> \brief Split the variable over two domains, subdividing along the independent variable j.
pair< TaylorModel<Interval>, TaylorModel<Interval> > split(const TaylorModel<Interval>& x, uint j);
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

//! \relates TaylorModel<Interval> \brief Evaluate an array of Taylor variables on a vector.
Interval evaluate(const TaylorModel<Interval>& x, const Vector<Interval>& sy);
//! \relates TaylorModel<Interval> \brief Evaluate an array of Taylor variables on a vector.
TaylorModel<Interval> partial_evaluate(const TaylorModel<Interval>& x, uint k, Float c);
//! \relates TaylorModel<Interval> \brief Evaluate an array of Taylor variables on a vector.
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
TaylorModel<Interval>  discard(const TaylorModel<Interval>&, const array<uint>& variables);

TaylorModel<Interval> recondition(const TaylorModel<Interval>& tm, array<uint>& discarded_variables,
                                  uint number_of_error_variables);

TaylorModel<Interval> recondition(const TaylorModel<Interval>& tm, array<uint>& discarded_variables,
                                  uint number_of_error_variables, uint index_of_error);

//! \relates TaylorModel<Interval>
//! An over-approximation to the intersection of two Taylor models.
//! Since the intersection cannot be represented exactly in the class of
//! TaylorModels, truncation errors as well as roundoff errors may be present.
//! In the absence of roundoff errors, the result is a subset of both arguments,
//! and is guaranteed to contain any function contained in both arguments.
TaylorModel<Interval> intersection(const TaylorModel<Interval>& x1, const TaylorModel<Interval>& x2);

// Compose an array of Taylor variables with another, assuming that y has been scaled to have unit codomain
TaylorModel<Interval> compose(const TaylorModel<Interval>& x, const Vector< TaylorModel<Interval> >& y);

// Compose an array of Taylor variables with another, after scaling by the interval vectors
TaylorModel<Interval> compose(const TaylorModel<Interval>& x, const Vector<Interval>& bx, const Vector< TaylorModel<Interval> >& y);

Float norm(const Vector< TaylorModel<Interval> >& tv);

TaylorModel<Interval> max(const TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
TaylorModel<Interval> min(const TaylorModel<Interval>& x, const TaylorModel<Interval>& y);
TaylorModel<Interval> abs(const TaylorModel<Interval>& x);
TaylorModel<Interval> neg(const TaylorModel<Interval>& x);
TaylorModel<Interval> rec(const TaylorModel<Interval>& x);
TaylorModel<Interval> sqr(const TaylorModel<Interval>& x);
TaylorModel<Interval> pow(const TaylorModel<Interval>& x, int n);
TaylorModel<Interval> sqrt(const TaylorModel<Interval>& x);
TaylorModel<Interval> exp(const TaylorModel<Interval>& x);
TaylorModel<Interval> log(const TaylorModel<Interval>& x);
TaylorModel<Interval> sin(const TaylorModel<Interval>& x);
TaylorModel<Interval> cos(const TaylorModel<Interval>& x);
TaylorModel<Interval> tan(const TaylorModel<Interval>& x);
TaylorModel<Interval> asin(const TaylorModel<Interval>& x);
TaylorModel<Interval> acos(const TaylorModel<Interval>& x);
TaylorModel<Interval> atan(const TaylorModel<Interval>& x);

std::ostream& operator<<(std::ostream&, const TaylorModel<Interval>::Accuracy&);



// Vector operations which can be evaluated componentwise
bool refines(const Vector< TaylorModel<Interval> >&,const Vector< TaylorModel<Interval> >&);
bool disjoint(const Vector< TaylorModel<Interval> >&,const Vector< TaylorModel<Interval> >&);
pair< Vector< TaylorModel<Interval> >, Vector< TaylorModel<Interval> > > split(const Vector< TaylorModel<Interval> >& x, uint j);
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




/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, Polynomila, TaylorModel<Interval><Interval>.
 */
template<>
class TaylorModel<Float>
{
    typedef TaylorModelAccuracy Accuracy;
    typedef Expansion<Float> ExpansionType;
  private:
    ExpansionType _expansion;
    mutable shared_ptr<Accuracy> _accuracy_ptr;
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
    TaylorModel<Float>(uint as, shared_ptr<Accuracy> acc);
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
    shared_ptr<Accuracy> accuracy_ptr() const { return this->_accuracy_ptr; }
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
    //! \brief Evaluate the quantity at the point \a x.
    Float evaluate(const Vector<Float>& x) const;
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
    /*! \name Simplification operations. */
    //! \brief Truncate to the default maximum degree of the quantity.
    TaylorModel<Float>& truncate();
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorModel<Float>& sweep();
    //! \brief Remove all terms which have high degree or small magnitude.
    TaylorModel<Float>& clean();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel<Float>& unique_sort();
    //@}

  public:
    //@{
    /*! \name Standard algebra interface. */
    //! \brief Create a dynamically-allocated model which is the zero constant.
    virtual TaylorModel<Float>* create() const;
    //! \brief Create a dynamically-allocated copy.
    virtual TaylorModel<Float>* clone() const;
    //! \brief An approximation to the norm of the function.
    virtual Float norm() const;
    //! \brief An approximation to the average value of the function.
    virtual Float average() const;
    //! \brief The tolerance to which analytic functions should be computed.
    virtual Float tolerance() const;
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
    //! \brief Addition of a scalar multiple of a Taylor model.
    TaylorModel<Float> sma(const Float& c, const TaylorModel<Float>& x) const;


    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const TaylorModel<Float>& x);
    //@}

  public:
    std::ostream& str(std::ostream&) const;
    std::ostream& repr(std::ostream&) const;
};

TaylorModel<Float> neg(const TaylorModel<Float>& t);
TaylorModel<Float> rec(const TaylorModel<Float>& t);
TaylorModel<Float> sqr(const TaylorModel<Float>& t);
TaylorModel<Float> pow(const TaylorModel<Float>& t, int n);
TaylorModel<Float> sqrt(const TaylorModel<Float>& t);
TaylorModel<Float> exp(const TaylorModel<Float>& t);
TaylorModel<Float> log(const TaylorModel<Float>& t);
TaylorModel<Float> sin(const TaylorModel<Float>& t);
TaylorModel<Float> cos(const TaylorModel<Float>& t);
TaylorModel<Float> tan(const TaylorModel<Float>& t);

inline TaylorModel<Float> operator+(const TaylorModel<Float>& t) {
    TaylorModel<Float> r=t; return r; }
inline TaylorModel<Float> operator-(const TaylorModel<Float>& t) {
    TaylorModel<Float> r=t; r.imul(-1); return r; }
inline TaylorModel<Float> operator+(const TaylorModel<Float>& t1, const TaylorModel<Float>& t2) {
    TaylorModel<Float> r(t1); r.isma(+1,t2); return r; }
inline TaylorModel<Float> operator-(const TaylorModel<Float>& t1, const TaylorModel<Float>& t2) {
    TaylorModel<Float> r(t1); r.isma(-1,t2); return r; }
inline TaylorModel<Float> operator*(const TaylorModel<Float>& t1, const TaylorModel<Float>& t2) {
    TaylorModel<Float> r=t1.null(); r.ifma(t1,t2); return r; }
inline TaylorModel<Float> operator/(const TaylorModel<Float>& t1, const TaylorModel<Float>& t2) {
    return t1*rec(t2); }
inline TaylorModel<Float>& operator+=(TaylorModel<Float>& t1, const TaylorModel<Float>& t2) {
    t1.isma(+1,t2); return t1; }
inline TaylorModel<Float>& operator-=(TaylorModel<Float>& t1, const TaylorModel<Float>& t2) {
    t1.isma(+1,t2); return t1; }
inline TaylorModel<Float> operator+(const TaylorModel<Float>& t, const Float& c) {
    TaylorModel<Float> r=t; r.iadd(c); return r; }
inline TaylorModel<Float> operator-(const TaylorModel<Float>& t, const Float& c) {
    TaylorModel<Float> r=t; r.iadd(-c); return r; }
inline TaylorModel<Float> operator*(const TaylorModel<Float>& t, const Float& c) {
    TaylorModel<Float> r=t; r.imul(c); return r; }
inline TaylorModel<Float> operator/(const TaylorModel<Float>& t, const Float& c) {
    TaylorModel<Float> r=t; r.imul(static_cast<Float>(1)/c); return r; }
inline TaylorModel<Float> operator+(const Float& c, const TaylorModel<Float>& t) {
    TaylorModel<Float> r=t; r.iadd(c); return r; }
inline TaylorModel<Float> operator-(const Float& c, const TaylorModel<Float>& t) {
    TaylorModel<Float> r=neg(t); r.iadd(c); return r; }
inline TaylorModel<Float> operator*(const Float& c, const TaylorModel<Float>& t) {
    TaylorModel<Float> r=t; r.imul(c); return r; }
inline TaylorModel<Float> operator/(const Float& c, const TaylorModel<Float>& t) {
    TaylorModel<Float> r=rec(t); r.imul(c); return r; }
inline TaylorModel<Float>& operator+=(TaylorModel<Float>& t, const Float& c) {
    t.iadd(c); return t; }
inline TaylorModel<Float>& operator-=(TaylorModel<Float>& t, const Float& c) {
    t.iadd(-c); return t; }
inline TaylorModel<Float>& operator*=(TaylorModel<Float>& t, const Float& c) {
    t.imul(c); return t; }
inline TaylorModel<Float>& operator/=(TaylorModel<Float>& t, const Float& c) {
    t.imul(static_cast<Float>(1)/c); return t; }


// FIXME: These should not need to be given explicitly
inline TaylorModel<Float>& operator+=(TaylorModel<Float>& t, const Real& c) { return t+=Float(c); }
inline TaylorModel<Float>& operator-=(TaylorModel<Float>& t, const Real& c) { return t-=Float(c); }
inline TaylorModel<Float>& operator*=(TaylorModel<Float>& t, const Real& c) { return t*=Float(c); }
inline TaylorModel<Float>& operator/=(TaylorModel<Float>& t, const Real& c) { return t/=Float(c); }
inline TaylorModel<Float> operator+(const TaylorModel<Float>& t, const Real& c) { return t+Float(c); }
inline TaylorModel<Float> operator-(const TaylorModel<Float>& t, const Real& c) { return t-Float(c); }
inline TaylorModel<Float> operator*(const TaylorModel<Float>& t, const Real& c) { return t*Float(c); }
inline TaylorModel<Float> operator/(const TaylorModel<Float>& t, const Real& c) { return t/Float(c); }
inline TaylorModel<Float> operator+(const Real& c, const TaylorModel<Float>& t) { return Float(c)+t; }
inline TaylorModel<Float> operator-(const Real& c, const TaylorModel<Float>& t) { return Float(c)-t; }
inline TaylorModel<Float> operator*(const Real& c, const TaylorModel<Float>& t) { return Float(c)*t; }
inline TaylorModel<Float> operator/(const Real& c, const TaylorModel<Float>& t) { return Float(c)/t; }

inline TaylorModel<Float>& operator+=(TaylorModel<Float>& t, double c) { return t+=Float(c); }
inline TaylorModel<Float>& operator-=(TaylorModel<Float>& t, double c) { return t-=Float(c); }
inline TaylorModel<Float>& operator*=(TaylorModel<Float>& t, double c) { return t*=Float(c); }
inline TaylorModel<Float>& operator/=(TaylorModel<Float>& t, double c) { return t/=Float(c); }
inline TaylorModel<Float> operator+(const TaylorModel<Float>& t, double c) { return t+Float(c); }
inline TaylorModel<Float> operator-(const TaylorModel<Float>& t, double c) { return t-Float(c); }
inline TaylorModel<Float> operator*(const TaylorModel<Float>& t, double c) { return t*Float(c); }
inline TaylorModel<Float> operator/(const TaylorModel<Float>& t, double c) { return t/Float(c); }
inline TaylorModel<Float> operator+(double c, const TaylorModel<Float>& t) { return Float(c)+t; }
inline TaylorModel<Float> operator-(double c, const TaylorModel<Float>& t) { return Float(c)-t; }
inline TaylorModel<Float> operator*(double c, const TaylorModel<Float>& t) { return Float(c)*t; }
inline TaylorModel<Float> operator/(double c, const TaylorModel<Float>& t) { return Float(c)/t; }

inline std::ostream& operator<<(std::ostream& os, const TaylorModel<Float>& x) {
    x.str(os); return os; }

inline Vector<Interval> codomain(const Vector< TaylorModel<Float> >& t) {
    Vector<Interval> r(t.size()); for(uint i=0; i!=t.size(); ++i) { r[i]=t[i].codomain(); } return r; }



} // namespace Ariadne

#endif // ARIADNE_TAYLOR_MODEL_H
