/***************************************************************************
 *            taylor_variable.h
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

/*! \file taylor_variable.h
 *  \brief Approximate functions on a bounded domain with a sparse representation.
 */

#ifndef ARIADNE_TAYLOR_VARIABLE_H
#define ARIADNE_TAYLOR_VARIABLE_H

#include <map>

#include "macros.h"
#include "array.h"
#include "vector.h"
#include "multi_index.h"
#include "expansion.h"

namespace Ariadne {

template<class T1, class T2> class Product;
template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Expansion;
template<class X> class Polynomial;

class TaylorVariable;
template<> class Vector<TaylorVariable>;

// Split the variable over two domains, subdividing along the independent variable j.
pair<TaylorVariable,TaylorVariable> split(const TaylorVariable& x, uint j);
pair< Vector<TaylorVariable>, Vector<TaylorVariable> > split(const Vector<TaylorVariable>& x, uint j);

// Scale the variabe by post-composing with an affine map taking the interval \a ivl to the unit interval
TaylorVariable unscale(const TaylorVariable& x, const Interval& ivl);
Vector<TaylorVariable> unscale(const Vector<TaylorVariable>& x, const Vector<Interval>& bx);

// Scale the variable by post-composing with an affine map taking the unit interval to \a ivl.
TaylorVariable scale(const TaylorVariable& x, const Interval& ivl);
Vector<TaylorVariable> scale(const Vector<TaylorVariable>& x, const Vector<Interval>& bx);

// Evaluate an array of Taylor variables on a vector.
Interval evaluate(const TaylorVariable& x, const Vector<Interval>& sy);
Vector<Interval> evaluate(const Vector<TaylorVariable>& x, const Vector<Interval>& sy);


// Compose an array of Taylor variables with another, after scaling by the interval vectors
TaylorVariable compose(const TaylorVariable& x, const Vector<TaylorVariable>& y);
Vector<TaylorVariable> compose(const Vector<TaylorVariable>& x, const Vector<TaylorVariable>& y);

Vector<TaylorVariable> combine(const Vector<TaylorVariable>& x1, const TaylorVariable& x2);
Vector<TaylorVariable> combine(const Vector<TaylorVariable>& x1, const Vector<TaylorVariable>& x2);

// Antidifferentiation operator
TaylorVariable antiderivative(const TaylorVariable& x, uint k);
Vector<TaylorVariable> antiderivative(const Vector<TaylorVariable>& x, uint k);

// Embed the variable in a space of higher dimension
//TaylorVariable embed(const TaylorVariable& tv, uint as, uint b);
//Vector<TaylorVariable> embed(const Vector<TaylorVariable>& tvs, uint as, uint b);

// Test if a variable refines another
bool refines(const TaylorVariable& tv1, const TaylorVariable& tv2);
bool refines(const Vector<TaylorVariable>& tv1, const Vector<TaylorVariable>& tv2);


/*! \brief A class representing a quantity depending on other quantities.
 *  Based on a power series Expansion, scaled to the unit box.
 *
 * See also TaylorFunction, TaylorSet.
 */
class TaylorVariable
{
    typedef Vector<Interval> DomainType;
    typedef Expansion<Float> ExpansionType;
    typedef Float ErrorType;
    static const Float _zero;
    DomainType _domain;
    ExpansionType _expansion;
    Float _error;
    double _sweep_threshold;
    uint _maximum_degree;
    MultiIndexBound _maximum_index;

  private:
    static double _default_sweep_threshold;
    static uint _default_maximum_degree;
  public:
    static const double em;
    static const double ec;
  public:
    //! \brief The type used for the coefficients.
    typedef Float ScalarType;
    //! \brief The type used to index the coefficients.
    typedef MultiIndex IndexType;
    //! \brief The type used for the coefficients.
    typedef Float ValueType;

    //! \brief An iterator through the (index,coefficient) pairs of the expansion expansion.
    typedef ExpansionType::iterator iterator;
    //! \brief A constant iterator through the (index,coefficient) pairs of the expansion expansion.
    typedef ExpansionType::const_iterator const_iterator;

    //@{
    /*! \name Constructors and destructors. */
    //! \brief Default constructor.
    explicit TaylorVariable();
    //! \brief Construct a TaylorVariable over the domain \a d.
    explicit TaylorVariable(const DomainType& d);
    //! \brief Construct a TaylorVariable over the domain \a d, with scaled power series expansion \a f and error \a e.
    explicit TaylorVariable(const DomainType& d, const ExpansionType& f, const ErrorType& e=0);
    //! \brief Construct a TaylorVariable over the domain \a d from the polynomial \a p.
    template<class X> explicit TaylorVariable(const DomainType& d, const Polynomial<X>& p);
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorVariable& operator=(const Float& c);
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorVariable& operator=(const Interval& c);
    //! \brief Test if the quantity is a better approximation than \a t throughout the domain.
    bool refines(const TaylorVariable& t);
    //@}

    //@{
    /*! \name Data access */
    //! \brief The domain of the quantity.
    const Vector<Interval>& domain() const { return this->_domain; }
    //! \brief The scaled expansion over a unit box.
    const ExpansionType& expansion() const { return this->_expansion; }
    //! \brief The error of the expansion over the domain.
    const ErrorType& error() const { return this->_error; }
    //! \brief A reference to the expansion.
    ExpansionType& expansion() { return this->_expansion; }
    //! \brief A reference to the error of the expansion over the domain.
    ErrorType& error() { return this->_error; }
    //! \brief A reference to the constant term in the expansion.
    Float& value() { return (*this)[MultiIndex::zero(this->argument_size())]; }
    //! \brief The constant term in the expansion.
    const Float& value() const { return (*this)[MultiIndex::zero(this->argument_size())]; }

    //! \brief Set the error of the expansion.
    void set_error(const Float& ne) {
        ARIADNE_ASSERT(ne>=0); this->_error=ne; }
    //! \brief Set the constant term in the expansion.
    void set_value(const Float& c) {
        this->_expansion[MultiIndex::zero(this->argument_size())]=c; }
    //! \brief Set the gradient of \a j term in the expansion.
    void set_gradient(unsigned int j, const Float& c) {
        this->_expansion[MultiIndex::unit(this->argument_size(),j)]=c; }

    //! \brief The coefficient of the term in $x^a$.
    const Float& operator[](const MultiIndex& a) const { return this->_expansion[a]; }
    //! \brief A read/write reference to the coefficient of the term in $x^a$.
    Float& operator[](const MultiIndex& a) { return this->_expansion[a]; }

    //! \brief An iterator to the first term in the expansion expansion.
    iterator begin() { return this->_expansion.begin(); }
    //! \brief A constant iterator to the first term in the expansion expansion.
    const_iterator begin() const { return this->_expansion.begin(); }
    //! \brief An iterator to the end of the expansion expansion.
    iterator end() { return this->_expansion.end(); }
    //! \brief A constant iterator to the end of the expansion expansion.
    const_iterator end() const { return this->_expansion.end(); }
    //! \brief An iterator to the term with index \a.
    iterator find(const MultiIndex& a) { return this->_expansion.find(a); }
    //! \brief A constant iterator to the term with index \a.
    const_iterator find(const MultiIndex& a) const { return this->_expansion.find(a); }

    //! \brief The number of variables in the argument of the quantity.
    uint argument_size() const { return this->_expansion.argument_size(); }
    //! \brief The maximum degree of terms in the expansion expansion.
    uint degree() const { return (--this->_expansion.end())->first.degree(); }
    //! \brief The number of nonzero terms in the expansion expansion.
    uint number_of_nonzeros() const { return this->_expansion.number_of_nonzeros(); }
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct a constant with value c \a c over the interval \a d.
    static TaylorVariable constant(const Interval& d, const Float& c);
    //! \brief Construct the quantity \f$x\f$ over the scalar domain \a d.
    static TaylorVariable scaling(const Interval& d);

    //! \brief Construct the zero quantity in \a as independent variables.
    static TaylorVariable zero(const DomainType& d);
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorVariable constant(const DomainType& d, const Float& c);
    //! \brief Construct the quantity \f$x_j\f$ over the domain \a d.
    static TaylorVariable variable(const DomainType& d, unsigned int j);
    //! \brief Construct the quantity \f$x_j\f$ over the domain \a d.
    static TaylorVariable scaling(const DomainType& d, unsigned int j);

    //! \brief Construct the quantity \f$c+\sum g_jx_j\f$ over the domain \a d.
    static TaylorVariable affine(const DomainType& d, const Float& c, const Vector<Float>& g);
    //! \brief Construct the quantity \f$c+\sum g_jx_j \pm e\f$.
    static TaylorVariable affine(const DomainType& d, const Float& x, const Vector<Float>& g, const Float& e) ;

    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector<TaylorVariable> constants(const DomainType& d, const Vector<Float>& c);
    //! \brief Return the vector of variables with values \a x in \a as arguments.
    static Vector<TaylorVariable> variables(const DomainType& d);
    //! \brief Return the vector of variables with values \a x in \a as arguments.
    static Vector<TaylorVariable> scaling(const DomainType& d);
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    bool operator==(const TaylorVariable& tv) const;
    //! \brief Inequality operator.
    bool operator!=(const TaylorVariable& tv) const { return !(*this==tv); }
    //@}

    //@{
    /*! \name Function operations. */
    //! \brief An over-approximation to the range of the quantity.
    Interval range() const;
    //! \brief Evaluate the quantity at the point \a x.
    Interval evaluate(const Vector<Float>& x) const;
    //! \brief Evaluate the quantity over the interval of points \a x.
    Interval evaluate(const Vector<Interval>& x) const;
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Truncate to the default maximum degree of the quantity.
    TaylorVariable& truncate();
    //! \brief Truncate to degree \a deg.
    TaylorVariable& truncate(uint deg);
    //! \brief Truncate all terms with any coefficient higher than \a a.
    TaylorVariable& truncate(const MultiIndex& a);
    //! \brief Truncate all terms with any coefficient higher than those given by \a a.
    TaylorVariable& truncate(const MultiIndexBound& a);
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorVariable& sweep();
    //! \brief Remove all terms whose coefficient has magnitude less than \a eps.
    TaylorVariable& sweep(double eps);
    //! \brief Remove all terms whose degree is higher than \a deg or
    //! whose coefficient has magnitude less than \a eps.
    TaylorVariable& clean(uint deg, double eps);
    //! \brief Remove all terms which have high degree or small magnitude.
    TaylorVariable& clean();
    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \brief .
    static void set_default_maximum_degree(uint md) { _default_maximum_degree=md; }
    //! \brief .
    static void set_default_sweep_threshold(double me) { ARIADNE_ASSERT(me>=0.0); _default_sweep_threshold=me; }
    //! \brief .
    static uint default_maximum_degree() { return _default_maximum_degree; }
    //! \brief .
    static double default_sweep_threshold() { return _default_sweep_threshold; }
    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \brief .
    void set_maximum_index(MultiIndexBound md) { this->_maximum_index=md; }
    //! \brief .
    void set_maximum_degree(uint md) { this->_maximum_degree=md; }
    //! \brief .
    void set_sweep_threshold(double me) { ARIADNE_ASSERT(me>=0.0); this->_sweep_threshold=me; }
    //! \brief .
    MultiIndexBound maximum_index() const { return this->_maximum_index; }
    //! \brief .
    uint maximum_degree() const { return this->_maximum_degree; }
    //! \brief .
    double sweep_threshold() const { return this->_sweep_threshold; }
    //@}

    //@{
    /*! \name Arithmetic operations. */
    //! \brief Inplace addition of another variable.
    friend TaylorVariable& operator+=(TaylorVariable& x, const TaylorVariable& y);
    //! \brief Inplace subtraction of another variable.
    friend TaylorVariable& operator-=(TaylorVariable& x, const TaylorVariable& y);
    //! \brief Inplace addition of a product of two variables.
    friend TaylorVariable& operator+=(TaylorVariable& x, const Product<TaylorVariable,TaylorVariable>& y);
    //! \brief Inplace addition of an exact floating-point constant.
    friend TaylorVariable& operator+=(TaylorVariable& x, const Float& c);
    //! \brief Inplace addition of an interval constant.
    friend TaylorVariable& operator+=(TaylorVariable& x, const Interval& c);
    //! \brief Inplace subtraction of an exact floating-point constant.
    friend TaylorVariable& operator-=(TaylorVariable& x, const Float& c);
    //! \brief Inplace subtraction of an interval constant.
    friend TaylorVariable& operator-=(TaylorVariable& x, const Interval& c);
    //! \brief Inplace multiplication by an exact scalar.
    friend TaylorVariable& operator*=(TaylorVariable& x, const Float& c);
    //! \brief Inplace multiplication by an approximate scalar.
    friend TaylorVariable& operator*=(TaylorVariable& x, const Interval& c);
    //! \brief Inplace division by an exact scalar.
    friend TaylorVariable& operator/=(TaylorVariable& x, const Float& c);
    //! \brief Inplace division by an approximate scalar.
    friend TaylorVariable& operator/=(TaylorVariable& x, const Interval& c);

    //! \brief Unary plus.
    friend TaylorVariable operator+(const TaylorVariable& x);
    //! \brief Unary minus.
    friend TaylorVariable operator-(const TaylorVariable& x);
    //! \brief Addition.
    friend TaylorVariable operator+(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Subtraction.
    friend TaylorVariable operator-(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Multiplication.
    friend TaylorVariable operator*(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Division.
    friend TaylorVariable operator/(const TaylorVariable& x, const TaylorVariable& y);

    //! \brief Addition of a scakar.
    friend TaylorVariable operator+(const TaylorVariable& x, const Float& c);
    //! \brief Subtraction of a scakar.
    friend TaylorVariable operator-(const TaylorVariable& x, const Float& c);
    //! \brief Multiplication by a scakar.
    friend TaylorVariable operator*(const TaylorVariable& x, const Float& c);
    //! \brief Division by a scakar.
    friend TaylorVariable operator/(const TaylorVariable& x, const Float& c);
    //! \brief Addition of a scakar.
    friend TaylorVariable operator+(const TaylorVariable& x, const Interval& c);
    //! \brief Subtraction of a scakar.
    friend TaylorVariable operator-(const TaylorVariable& x, const Interval& c);
    //! \brief Multiplication by a scakar.
    friend TaylorVariable operator*(const TaylorVariable& x, const Interval& c);
    //! \brief Division by a scakar.
    friend TaylorVariable operator/(const TaylorVariable& x, const Interval& c);
    //! \brief Addition of a scakar.
    friend TaylorVariable operator+(const Float& c, const TaylorVariable& x);
    //! \brief Subtraction from a scakar.
    friend TaylorVariable operator-(const Float& c, const TaylorVariable& x);
    //! \brief Multiplication by a scakar.
    friend TaylorVariable operator*(const Float& c, const TaylorVariable& x);
    //! \brief Division through a scalar.
    friend TaylorVariable operator/(const Float& c, const TaylorVariable& x);
    //! \brief Addition of a scakar.
    friend TaylorVariable operator+(const Interval& c, const TaylorVariable& x);
    //! \brief Subtraction from a scakar.
    friend TaylorVariable operator-(const Interval& c, const TaylorVariable& x);
    //! \brief Multiplication by a scakar.
    friend TaylorVariable operator*(const Interval& c, const TaylorVariable& x);
    //! \brief Division through a scalar.
    friend TaylorVariable operator/(const Interval& c, const TaylorVariable& x);

    //! \brief Maximum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend TaylorVariable max(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Minimum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend TaylorVariable min(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Addition.
    friend TaylorVariable add(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Multiplication.
    friend TaylorVariable mul(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Absolute value. Throws an error if one variable is neither greater than
    //! zero over the entire domain, nor less than zero over the entire domain.
    friend TaylorVariable abs(const TaylorVariable& x);
    //! \brief Negation.
    friend TaylorVariable neg(const TaylorVariable& x);
    //! \brief Reciprocal.
    friend TaylorVariable rec(const TaylorVariable& x);
    //! \brief Square.
    friend TaylorVariable sqr(const TaylorVariable& x);
    //! \brief Power.
    friend TaylorVariable pow(const TaylorVariable& x, int n);
    //! \brief Square root.
    friend TaylorVariable sqrt(const TaylorVariable& x);
    //! \brief Natural exponent.
    friend TaylorVariable exp(const TaylorVariable& x);
    //! \brief Natural logarithm.
    friend TaylorVariable log(const TaylorVariable& x);
    //! \brief Sine.
    friend TaylorVariable sin(const TaylorVariable& x);
    //! \brief Cosine.
    friend TaylorVariable cos(const TaylorVariable& x);
    //! \brief Tangent.
    friend TaylorVariable tan(const TaylorVariable& x);
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const TaylorVariable& x);
    //@}

  public:
    std::string str() const;

    TaylorVariable& clobber();
    TaylorVariable& clobber(uint o);
    TaylorVariable& clobber(uint so, uint to);
};

TaylorVariable max(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable min(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable abs(const TaylorVariable& x);
TaylorVariable neg(const TaylorVariable& x);
TaylorVariable rec(const TaylorVariable& x);
TaylorVariable exp(const TaylorVariable& x);
TaylorVariable log(const TaylorVariable& x);
TaylorVariable sin(const TaylorVariable& x);
TaylorVariable cos(const TaylorVariable& x);
TaylorVariable tan(const TaylorVariable& x);


template<>
class Vector<TaylorVariable>
    : public ublas::vector<TaylorVariable>
{
  public:
    typedef Vector<Interval> DomainType;

    Vector() : ublas::vector<TaylorVariable>() { }
    Vector(uint rs) : ublas::vector<TaylorVariable>(rs) { for(uint i=0; i!=rs; ++i) { (*this)[i]=TaylorVariable(); } }
    Vector(uint rs, const Vector<Interval>& d) : ublas::vector<TaylorVariable>(rs) { for(uint i=0; i!=rs; ++i) { (*this)[i]=TaylorVariable(d); } }
    Vector(uint rs, const TaylorVariable& x) : ublas::vector<TaylorVariable>(rs) { for(uint i=0; i!=rs; ++i) { (*this)[i]=x; } }

    Vector(const Vector<Interval>& d, const Vector< Expansion<Float> >& f)
        : ublas::vector<TaylorVariable>(f.size())
        { for(uint i=0; i!=f.size(); ++i) { (*this)[i]=TaylorVariable(d,f[i],0.0); } }
    Vector(const Vector<Interval>& d, const Vector< Expansion<Float> >& f, const Vector<Float>& e)
        : ublas::vector<TaylorVariable>(f.size())
        { ARIADNE_ASSERT(f.size()==e.size()); for(uint i=0; i!=f.size(); ++i) { (*this)[i]=TaylorVariable(d,f[i],e[i]); } }

    template<class E> Vector(const ublas::vector_expression<E> &ve) : ublas::vector<TaylorVariable>(ve) { }

    uint result_size() const { return this->size(); }
    uint argument_size() const { ARIADNE_ASSERT(this->size()>0); return (*this)[0].argument_size(); }

    DomainType domain() const { return (*this)[0].domain(); }
    Vector< Expansion<Float> > expansion() const;
    Vector<Float> error() const;

    Vector<Float> value() const;
    Vector<Interval> evaluate(const Vector<Float>&) const;
    Vector<Interval> evaluate(const Vector<Interval>&) const;


    void check() const;
};



} // namespace Ariadne

#endif // ARIADNE_TAYLOR_VARIABLE_H
