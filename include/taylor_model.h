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

template<class T1, class T2> class Product;
template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Expansion;

class ScalarFunctionInterface;
class FunctionInterface;
class TaylorModel;
class TaylorCalculus;

class IntersectionException;
class ImplicitFunctionException;
class FlowBoundsException;

// Rescale the vector x from the domain d to the unit domain.
Vector<Interval> unscale(const Vector<Interval>& x, const Vector<Interval>& d);

//! \relates TaylorModel \brief The magnitude of the variable
Float mag(const TaylorModel& tm);
//! \relates TaylorModel \brief Split the variable over two domains, subdividing along the independent variable j.
pair<TaylorModel,TaylorModel> split(const TaylorModel& x, uint j);
//! \relates TaylorModel
//!\brief Split the variable, subdividing along the independent variable j
//! and taking the lower/upper half depending on whether half is false or true
TaylorModel split(const TaylorModel& x, uint j, bool half);

//! \relates TaylorModel \brief Scale the variable by post-composing with an affine map taking the unit interval to ivl.
TaylorModel scale(const TaylorModel& x, const Interval& ivl);
//! \relates TaylorModel \brief Scale the variable by post-composing with an affine map taking the interval ivl to the unit interval
TaylorModel unscale(const TaylorModel& x, const Interval& ivl);
//! \relates TaylorModel \brief Scale the variable by post-composing with an affine map taking the interval ivl1 to interval \a ivl2
TaylorModel rescale(const TaylorModel& x, const Interval& ivl1, const Interval& ivl2);

//! \relates TaylorModel \brief Evaluate an array of Taylor variables on a vector.
Interval evaluate(const TaylorModel& x, const Vector<Interval>& sy);
//! \relates TaylorModel \brief Evaluate an array of Taylor variables on a vector.
TaylorModel partial_evaluate(const TaylorModel& x, uint k, Float c);
//! \relates TaylorModel \brief Evaluate an array of Taylor variables on a vector.
TaylorModel partial_evaluate(const TaylorModel& x, uint k, Interval c);
//! \relates TaylorModel
//! Substitute the TaylorModel y in the  kth variable of \a x.
//! Precondition: x.argument_size()==y.argument_size()+1
TaylorModel substitute(const TaylorModel& x, uint k, const TaylorModel& y);

//! \relates TaylorModel \brief Embed the model in a space of higher dimension
TaylorModel embed(const TaylorModel& tv, uint as);
//! \relates TaylorModel \brief Embed the model in a space of higher dimension
TaylorModel embed(uint as, const TaylorModel& tv);

//! \relates TaylorModel \brief Test if a model refines another
bool refines(const TaylorModel& tv1, const TaylorModel& tv2);
//! \relates TaylorModel \brief Test if a model is disjoint from
bool disjoint(const TaylorModel& tv1, const TaylorModel& tv2);

//! \relates TaylorModel \brief Antidifferentiation operator
TaylorModel antiderivative(const TaylorModel& x, uint k);

//! \relates TaylorModel \brief Differentiation operator; discards error term
TaylorModel derivative(const TaylorModel& x, uint k);

//! \relates TaylorModel \brief Replace the variale x[k] with a*x[k]+b
TaylorModel preaffine(const TaylorModel&, uint k, const Interval& a, const Interval& b);

//! \relates TaylorModel
//! An over-approximation to the intersection of two Taylor models.
//! Since the intersection cannot be represented exactly in the class of
//! TaylorModels, truncation errors as well as roundoff errors may be present.
//! In the absence of roundoff errors, the result is a subset of both arguments,
//! and is guaranteed to contain any function contained in both arguments.
TaylorModel intersection(const TaylorModel& x1, const TaylorModel& x2);

// Compose an array of Taylor variables with another, assuming that y has been scaled to have unit codomain
TaylorModel compose(const TaylorModel& x, const Vector<TaylorModel>& y);

// Compose an array of Taylor variables with another, after scaling by the interval vectors
TaylorModel compose(const TaylorModel& x, const Vector<Interval>& bx, const Vector<TaylorModel>& y);

// Compute the implicit function f(x,h(x))=0
TaylorModel implicit(const TaylorModel& f);
TaylorModel implicit_step(const TaylorModel& f, const TaylorModel& h);

// Compute the implicit function f(g(x),h(x))=0
TaylorModel implicit(const ScalarFunctionInterface& f, const Vector<TaylorModel>& g);

// Solve the equation f(x)=0
Interval solve(const TaylorModel& f);

// Vector operations which can be evaluated componentwise
bool refines(const Vector<TaylorModel>&,const Vector<TaylorModel>&);
bool disjoint(const Vector<TaylorModel>&,const Vector<TaylorModel>&);
pair< Vector<TaylorModel>, Vector<TaylorModel> > split(const Vector<TaylorModel>& x, uint j);
Vector<TaylorModel> split(const Vector<TaylorModel>& x, uint j, bool half);
Vector<TaylorModel> unscale(const Vector<TaylorModel>& x, const Vector<Interval>& bx);
Vector<TaylorModel> scale(const Vector<TaylorModel>& x, const Vector<Interval>& bx);
Vector<Interval> evaluate(const Vector<TaylorModel>& x, const Vector<Interval>& sy);
Vector<TaylorModel> partial_evaluate(const Vector<TaylorModel>& x, uint k, Float sy);
Vector<TaylorModel> partial_evaluate(const Vector<TaylorModel>& x, uint k, Interval sy);
Vector<TaylorModel> substitute(const Vector<TaylorModel>& x, uint k, const TaylorModel& y);
Vector<TaylorModel> antiderivative(const Vector<TaylorModel>& x, uint k);
Vector<TaylorModel> embed(const Vector<TaylorModel>& x, uint as);
Vector<TaylorModel> embed(uint as, const Vector<TaylorModel>& x);
Matrix<Interval> jacobian(const Vector<TaylorModel>& x, const Vector<Interval>& d);
//Matrix<Interval> jacobian(const Vector<TaylorModel>& x);
bool refines(const Vector<TaylorModel>& x1, const Vector<TaylorModel>& x2);
Vector<TaylorModel> combine(const Vector<TaylorModel>& x1, const Vector<TaylorModel>& x2);
Vector<TaylorModel> combine(const Vector<TaylorModel>& x1, const TaylorModel& x2);
Vector<TaylorModel> combine(const TaylorModel& x1, const Vector<TaylorModel>& x2);
Vector<TaylorModel> combine(const TaylorModel& x1, const TaylorModel& x2);
Vector<TaylorModel> compose(const Vector<TaylorModel>& f, const Vector<TaylorModel>& g);

//Vector operations which cannot be computed componentwise
Vector<Interval> solve(const Vector<TaylorModel>& f);
Vector<TaylorModel> implicit(const Vector<TaylorModel>& f);
Vector<TaylorModel> implicit_step(const Vector<TaylorModel>& f, const Vector<TaylorModel>& h);
Vector<TaylorModel> flow(const Vector<TaylorModel>& x, const Vector<TaylorModel>& y0, uint order);
Vector<TaylorModel> parameterised_flow(const Vector<TaylorModel>& x, const Vector<TaylorModel>& y0, uint order);

TaylorModel unchecked_compose(const TaylorModel& x, const Vector<TaylorModel>& y);
Vector<TaylorModel> unchecked_compose(const Vector<TaylorModel>& x, const Vector<TaylorModel>& y);
Vector<TaylorModel> unchecked_flow(const Vector<TaylorModel>& x, const Vector<TaylorModel>& y0, uint order);

Float norm(const Vector<TaylorModel>& tv);

/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, TaylorExpression, TaylorFunction, TaylorSet.
 */
class TaylorModel
{
    friend class TaylorExpression;
    typedef Expansion<Float> ExpansionType;
  public:
    class Accuracy;
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
    TaylorModel();
    //! \brief Construct a TaylorModel in \a as arguments.
    TaylorModel(uint as);
    //! \brief Construct a TaylorModel in \a as arguments with the given accuracy control.
    TaylorModel(uint as, shared_ptr<Accuracy> acc);
    //! \brief Construct from a map giving the expansion, a constant giving the error, and an accuracy parameter.
    TaylorModel(const std::map<MultiIndex,Float>& d, const Float& e);
    //! \brief Construct from a map giving the expansion, and a constant giving the error.
    TaylorModel(const Expansion<Float>& f, const Float& e=0.0);
    //! \brief Construct from a map giving the expansion, a constant giving the error, and an accuracy parameter.
    TaylorModel(const Expansion<Float>& f, const Float& e, shared_ptr<Accuracy> a);
    //! \brief Fast swap with another Taylor model.
    void swap(TaylorModel& tm);
    //! \brief Set to zero.
    void clear();

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorModel& operator=(const Float& c);
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel& operator=(const Interval& c);
    //! \brief Test if the quantity is a better approximation than \a t throughout the domain.
    bool refines(const TaylorModel& t);
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct the zero quantity in \a as independent variables.
    static TaylorModel zero(uint as) {
        TaylorModel r(as); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel constant(uint as, const Float& c) {
        TaylorModel r(as); r.set_value(c); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel constant(uint as, const Interval& d) {
        TaylorModel r(as); r.set_value(1.0); r*=d; return r; }
    //! \brief Construct the quantity with expansion \f$x_j\f$ in \a as independent variables.
    static TaylorModel variable(uint as, uint j) {
        TaylorModel r(as); r.set_gradient(j,1.0); return r; }
    //! \brief Construct the quantity which scales the unit interval into the domain \a d.
    static TaylorModel scaling(uint as, uint j, const Interval& d) {
        TaylorModel r(as); r.set_gradient(j,1.0); r.rescale(Interval(-1,1),d); return r; }
    //! \brief Construct the quantity which scales the codomain \a cd into the unit interval.
    static TaylorModel unscaling(uint as, uint j, const Interval& d) {
        TaylorModel r(as); r.set_gradient(j,1.0); r.rescale(d,Interval(-1,+1)); return r; }
    //! \brief Construct the quantity which scales the interval \a cd onto the interval \a d.
    static TaylorModel rescaling(uint as, uint j, const Interval& cd, const Interval& d) {
        TaylorModel r(as); r.set_gradient(j,1.0); r.rescale(cd,d); return r; }
    //! \brief Construct the quantity \f$c+\sum g_jx_j\f$.
    static TaylorModel affine(const Float& c, const Vector<Float>& g) {
        TaylorModel r(g.size()); r.set_value(c);
        for(uint j=0; j!=g.size(); ++j) { r.set_gradient(j,g[j]); } return r; }
    //! \brief Construct the quantity \f$c+\sum g_jx_j \pm e\f$.
    static TaylorModel affine(const Float& x, const Vector<Float>& g, const Float& e) {
        TaylorModel r(g.size()); r.set_value(x); r.set_error(e);
        for(uint j=0; j!=g.size(); ++j) { r.set_gradient(j,g[j]); } return r; }

    //! \brief Return the vector of zero variables of size \a rs in \a as arguments.
    static Vector<TaylorModel> zeros(uint rs, uint as);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector<TaylorModel> constants(uint as, const Vector<Float>& c);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector<TaylorModel> constants(uint as, const Vector<Interval>& c);
    //! \brief Return the vector of variables on the unit domain.
    static Vector<TaylorModel> variables(uint as);
    //! \brief Return the vector scaling the unit interval onto the domain \a d.
    static Vector<TaylorModel> scalings(const Vector<Interval>& d);
    //! \brief Return the vector scaling the unit interval onto the codomain \a cd.
    static Vector<TaylorModel> unscalings(const Vector<Interval>& d);
    //! \brief Return the vector scaling the codomain \a cd onto the domain \a d.
    static Vector<TaylorModel> rescalings(const Vector<Interval>& cd, const Vector<Interval>& d);
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    bool operator==(const TaylorModel& sd) const {
        return this->_expansion==sd._expansion && this->_error == sd._error; }
    //! \brief Inequality operator.
    bool operator!=(const TaylorModel& sd) const {
        return !(*this==sd); }
    //! \brief Comparison with another Taylor model.
    tribool operator<(const TaylorModel& sd) const {
        return (sd-*this)>0; }
    //! \brief Comparison with another Taylor model.
    tribool operator>(const TaylorModel& sd) const {
        return (*this-sd)>0; }

    //! \brief Comparison with a scalar.
    tribool operator<(const Float& c) const {
        return this->range()<c; }
    //! \brief Comparison with a scalar.
    tribool operator>(const Float& c) const {
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
        this->_expansion[MultiIndex::zero(this->argument_size())]=c; }
    //! \brief Set the coefficient of the term \f$df/dx_j\f$.
    void set_gradient(uint j, const Float& c) {
        this->_expansion[MultiIndex::unit(this->argument_size(),j)]=c; }

    //! \brief The coefficient of the term in $x^a$.
    const Float& operator[](const MultiIndex& a) const { return this->_expansion[a]; }
    //! \brief A read/write reference to the coefficient of the term in $x^a$.
    Float& operator[](const MultiIndex& a) { return this->_expansion[a]; }

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
    uint degree() const { if(this->_expansion.empty()) { return 0; } else { return (--this->_expansion.end())->key().degree(); } }
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
    //! \brief Evaluate the quantity at the point \a x.
    Interval evaluate(const Vector<Float>& x) const;
    //! \brief Evaluate the quantity over the interval of points \a x.
    Interval evaluate(const Vector<Interval>& x) const;
    //@}

    //@{
    /*! \name Inplace modifications. */
    //! \brief Scale so that the old codomain maps into the new codomain.
    TaylorModel& rescale(const Interval& old_codomain, const Interval& new_codomain);
    //! \brief Restrict to a subdomain.
    TaylorModel& restrict(const Vector<Interval>& new_domain);
    //! \brief Compute the antiderivative (in place).
    TaylorModel& antidifferentiate(uint k);
    //@}

    //@{
    /*! \name Set-based operations. */
    //! \brief Test if one model refines (is a subset of) another.
    friend bool refines(const TaylorModel& tm1, const TaylorModel& tm2);
    //! \brief Test if one model is disjoint from (is incompatible with) another.
    friend bool disjoint(const TaylorModel& tm1, const TaylorModel& tm2);
    //! \brief An over-approximation of the intersection of the sets of functions
    //! allowed by the two models. It is guaranteed that any function represented
    //! by both models is also represented by the result.
    friend TaylorModel intersection(const TaylorModel& tm1, const TaylorModel& tm2);
    //! \brief The supremum norm of the model, given by \f$|e|+\sum_\alpha |c_\alpha|\f$.
    friend Float norm(const TaylorModel& tm);
    //@{
    /*! \name Vectoral function operators. */
    //! \brief Solve the equation \f$f(x)=0\f$ in the unit box.
    friend Vector<Interval> solve(const Vector<TaylorModel>& f);
    //! \brief Compose two models, where the second is scaled so that the codomain is a unit box.
    friend Vector<TaylorModel> compose(const Vector<TaylorModel>& f, const Vector<TaylorModel>& g);
    //! \brief Compute the implicit function h satisfying f(x,h(x))=0.
    friend Vector<TaylorModel> implicit(const Vector<TaylorModel>& f);
    //! \brief Compute the flow of the vector field \a vf, starting from the box \a d,
    //! over the time interval \a h, using temporal order \a o.
    friend Vector<TaylorModel> flow(const Vector<TaylorModel>& vf,
                                    const Vector<Interval>& d, const Interval& h, uint o);
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Truncate to the default maximum degree of the quantity.
    TaylorModel& truncate();
    //! \brief Truncate to degree \a deg.
    TaylorModel& truncate(uint deg);
    //! \brief Truncate all terms with any coefficient higher than \a a.
    TaylorModel& truncate(const MultiIndex& a);
    //! \brief Truncate all terms with any coefficient higher than those given by \a a.
    TaylorModel& truncate(const MultiIndexBound& a);
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorModel& sweep();
    //! \brief Remove all terms whose coefficient has magnitude less than \a eps.
    TaylorModel& sweep(double eps);
    //! \brief Remove all terms whose degree is higher than \a deg or
    //! whose coefficient has magnitude less than \a eps.
    TaylorModel& clean(const Accuracy& accuracy);
    //! \brief Remove all terms which have high degree or small magnitude.
    TaylorModel& clean();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel& unique_sort();
    //@}

    //@{
    /*! \name Default accuracy parameters. */
    //! \brief Set the default maximum degree.
    static void set_default_maximum_degree(uint md);
    //! \brief Set the default sweep threshold.
    static void set_default_sweep_threshold(double me);
    //! \brief The default maximum degree.
    static uint default_maximum_degree();
    //! \brief The default sweep threshold.
    static double default_sweep_threshold();
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
    friend TaylorModel& operator+=(TaylorModel& x, const TaylorModel& y);
    //! \brief Inplace subtraction of another variable.
    friend TaylorModel& operator-=(TaylorModel& x, const TaylorModel& y);
    //! \brief Inplace multiplication of another variable. Not any more efficient than ordinary multiplication.
    friend TaylorModel& operator*=(TaylorModel& x, const TaylorModel& y);
    //! \brief Inplace division of another variable. Not any more efficient than ordinary division.
    friend TaylorModel& operator/=(TaylorModel& x, const TaylorModel& y);
    //! \brief Inplace addition of a product of two variables.
    friend TaylorModel& operator+=(TaylorModel& x, const Product<TaylorModel,TaylorModel>& y);
    //! \brief Inplace addition of an exact floating-point constant.
    friend TaylorModel& operator+=(TaylorModel& x, const Float& c);
    //! \brief Inplace addition of an interval constant.
    friend TaylorModel& operator+=(TaylorModel& x, const Interval& c);
    //! \brief Inplace subtraction of an exact floating-point constant.
    friend TaylorModel& operator-=(TaylorModel& x, const Float& c);
    //! \brief Inplace subtraction of an interval constant.
    friend TaylorModel& operator-=(TaylorModel& x, const Interval& c);
    //! \brief Inplace multiplication by an exact scalar.
    friend TaylorModel& operator*=(TaylorModel& x, const Float& c);
    //! \brief Inplace multiplication by an approximate scalar.
    friend TaylorModel& operator*=(TaylorModel& x, const Interval& c);
    //! \brief Inplace division by an exact scalar.
    friend TaylorModel& operator/=(TaylorModel& x, const Float& c);
    //! \brief Inplace division by an approximate scalar.
    friend TaylorModel& operator/=(TaylorModel& x, const Interval& c);

    //! \brief Unary plus.
    friend TaylorModel operator+(const TaylorModel& x);
    //! \brief Unary minus.
    friend TaylorModel operator-(const TaylorModel& x);
    //! \brief Addition.
    friend TaylorModel operator+(const TaylorModel& x, const TaylorModel& y);
    //! \brief Subtraction.
    friend TaylorModel operator-(const TaylorModel& x, const TaylorModel& y);
    //! \brief Multiplication.
    friend TaylorModel operator*(const TaylorModel& x, const TaylorModel& y);
    //! \brief Division.
    friend TaylorModel operator/(const TaylorModel& x, const TaylorModel& y);

    //! \brief Addition of a scakar.
    friend TaylorModel operator+(const TaylorModel& x, const Float& c);
    //! \brief Subtraction of a scakar.
    friend TaylorModel operator-(const TaylorModel& x, const Float& c);
    //! \brief Multiplication by a scakar.
    friend TaylorModel operator*(const TaylorModel& x, const Float& c);
    //! \brief Division by a scakar.
    friend TaylorModel operator/(const TaylorModel& x, const Float& c);
    //! \brief Addition of a scakar.
    friend TaylorModel operator+(const TaylorModel& x, const Interval& c);
    //! \brief Subtraction of a scakar.
    friend TaylorModel operator-(const TaylorModel& x, const Interval& c);
    //! \brief Multiplication by a scakar.
    friend TaylorModel operator*(const TaylorModel& x, const Interval& c);
    //! \brief Division by a scakar.
    friend TaylorModel operator/(const TaylorModel& x, const Interval& c);
    //! \brief Addition of a scakar.
    friend TaylorModel operator+(const Float& c, const TaylorModel& x);
    //! \brief Subtraction from a scakar.
    friend TaylorModel operator-(const Float& c, const TaylorModel& x);
    //! \brief Multiplication by a scakar.
    friend TaylorModel operator*(const Float& c, const TaylorModel& x);
    //! \brief Division through a scalar.
    friend TaylorModel operator/(const Float& c, const TaylorModel& x);
    //! \brief Addition of a scakar.
    friend TaylorModel operator+(const Interval& c, const TaylorModel& x);
    //! \brief Subtraction from a scakar.
    friend TaylorModel operator-(const Interval& c, const TaylorModel& x);
    //! \brief Multiplication by a scakar.
    friend TaylorModel operator*(const Interval& c, const TaylorModel& x);
    //! \brief Division through a scalar.
    friend TaylorModel operator/(const Interval& c, const TaylorModel& x);
    //@}

    //@{
    /*! \name Algebraic and transcendental functions. */

    //! \brief Maximum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend TaylorModel max(const TaylorModel& x, const TaylorModel& y);
    //! \brief Minimum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend TaylorModel min(const TaylorModel& x, const TaylorModel& y);
    //! \brief Addition.
    friend TaylorModel add(const TaylorModel& x, const TaylorModel& y);
    //! \brief Multiplication.
    friend TaylorModel mul(const TaylorModel& x, const TaylorModel& y);
    //! \brief Absolute value. Throws an error if one variable is neither greater than
    //! zero over the entire domain, nor less than zero over the entire domain.
    friend TaylorModel abs(const TaylorModel& x);
    //! \brief Negation.
    friend TaylorModel neg(const TaylorModel& x);
    //! \brief Reciprocal.
    friend TaylorModel rec(const TaylorModel& x);
    //! \brief Square.
    friend TaylorModel sqr(const TaylorModel& x);
    //! \brief Power.
    friend TaylorModel pow(const TaylorModel& x, int n);
    //! \brief Square root.
    friend TaylorModel sqrt(const TaylorModel& x);
    //! \brief Natural exponent.
    friend TaylorModel exp(const TaylorModel& x);
    //! \brief Natural logarithm.
    friend TaylorModel log(const TaylorModel& x);
    //! \brief Sine.
    friend TaylorModel sin(const TaylorModel& x);
    //! \brief Cosine.
    friend TaylorModel cos(const TaylorModel& x);
    //! \brief Tangent.
    friend TaylorModel tan(const TaylorModel& x);
    //! \brief Inverse sine.
    friend TaylorModel asin(const TaylorModel& x);
    //! \brief Inverse cosine.
    friend TaylorModel acos(const TaylorModel& x);
    //! \brief Inverse tangent.
    friend TaylorModel atan(const TaylorModel& x);
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const TaylorModel& x);
    //@}

  public:
    std::ostream& str(const std::ostream&) const;
    std::ostream& repr(const std::ostream&) const;

    TaylorModel& clobber();
    TaylorModel& clobber(uint o);
    TaylorModel& clobber(uint so, uint to);
};

TaylorModel max(const TaylorModel& x, const TaylorModel& y);
TaylorModel min(const TaylorModel& x, const TaylorModel& y);
TaylorModel abs(const TaylorModel& x);
TaylorModel neg(const TaylorModel& x);
TaylorModel rec(const TaylorModel& x);
TaylorModel sqr(const TaylorModel& x);
TaylorModel pow(const TaylorModel& x, int n);
TaylorModel sqrt(const TaylorModel& x);
TaylorModel exp(const TaylorModel& x);
TaylorModel log(const TaylorModel& x);
TaylorModel sin(const TaylorModel& x);
TaylorModel cos(const TaylorModel& x);
TaylorModel tan(const TaylorModel& x);
TaylorModel asin(const TaylorModel& x);
TaylorModel acos(const TaylorModel& x);
TaylorModel atan(const TaylorModel& x);

std::ostream& operator<<(std::ostream&, const TaylorModel::Accuracy&);


struct TaylorModel::Accuracy {
    friend class TaylorModel;
    friend class TaylorCalculus;

    Accuracy();
    Accuracy(double st, uint md);

    friend Accuracy max(const Accuracy& acc1, const Accuracy& acc2);
    friend Accuracy min(const Accuracy& acc1, const Accuracy& acc2);
    friend std::ostream& operator<<(std::ostream& os, const Accuracy& acc);
  public:
    bool discard(const Float& x) const;
    bool discard(const MultiIndex& a) const;
    bool discard(const MultiIndex& a, const Float& x) const;
  private:
    double _sweep_threshold;
    uint _maximum_degree;
};

struct IntersectionException : public std::runtime_error {
    IntersectionException(const std::string& what) : std::runtime_error(what) { }
};

struct ImplicitFunctionException : public std::runtime_error {
    ImplicitFunctionException(const std::string& what) : std::runtime_error(what) { }
};

struct FlowBoundsException : public std::runtime_error {
    FlowBoundsException(const std::string& what) : std::runtime_error(what) { }
};


} // namespace Ariadne

#endif // ARIADNE_TAYLOR_MODEL_H
