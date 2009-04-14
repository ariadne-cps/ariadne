/***************************************************************************
 *            taylor_function.h
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

/*! \file taylor_function.h
 *  \brief Over-approximations of functions based on Taylor expansions.
 */

#ifndef ARIADNE_TAYLOR_FUNCTION_H
#define ARIADNE_TAYLOR_FUNCTION_H

#include <iosfwd>
#include "numeric.h"
#include "vector.h"
#include "taylor_model.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Polynomial;

class FunctionInterface;
class MultiIndex;
class TaylorModel;
class TaylorExpression;
class TaylorFunction;


TaylorFunction restrict(const TaylorFunction&, const Vector<Interval>& bx);
TaylorExpression compose(const TaylorExpression&, const TaylorFunction&);
TaylorFunction compose(const TaylorFunction&, const TaylorFunction&);
TaylorFunction compose(const FunctionInterface&, const TaylorFunction&);
TaylorFunction antiderivative(const TaylorFunction&, uint);
TaylorFunction implicit(const TaylorFunction&);
TaylorFunction flow(const TaylorFunction& vf, const Vector<Interval>& d, const Interval& t, uint o);

TaylorExpression unchecked_compose(const TaylorExpression&, const TaylorFunction&);
TaylorFunction unchecked_compose(const TaylorFunction&, const TaylorFunction&);
TaylorFunction unchecked_implicit(const TaylorFunction&);
TaylorFunction unchecked_flow(const TaylorFunction& vf, const Vector<Interval>& d, const Interval& t, uint o);




/*! \brief A taylor_model with multivalued output using the TaylorModel class.
 *
 *  See also TaylorModel, TaylorExpression, TaylorFunction.
 */
class TaylorFunction {
    typedef Float R;
    typedef Interval I;
  public:
    /*! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. */
    TaylorFunction();

    /*! \brief Construct from a domain and the expansion. */
    TaylorFunction(const Vector<Interval>& domain,
                   const Vector< Expansion<Float> >& expansion);

    /*! \brief Construct from a domain, and expansion and errors. */
    TaylorFunction(const Vector<Interval>& domain,
                   const Vector< Expansion<Float> >& expansion,
                   const Vector<Float>& error);

    /*! \brief Construct from a domain and the models. */
    explicit TaylorFunction(const Vector<Interval>& domain, const Vector<TaylorModel>& variables);

    /*! \brief Construct from a domain and a function. */
    TaylorFunction(const Vector<Interval>& domain,
                   const FunctionInterface& function);

    /*! \brief Construct from a domain, a function, and accuracy paramters. */
    TaylorFunction(const Vector<Interval>& domain,
                   const FunctionInterface& function,
                   shared_ptr<TaylorModel::Accuracy> accuracy_ptr);

    /*! \brief Construct from a domain and a polynomial. */
    TaylorFunction(const Vector<Interval>& domain,
                   const Vector< Polynomial<Float> >& polynomial);

    /*! \brief Construct from a domain and a n interval polynomial. */
    TaylorFunction(const Vector<Interval>& domain,
                   const Vector< Polynomial<Interval> >& polynomial);

    /*! \brief Construct from a vector of Taylor variables. */
    TaylorFunction(const Vector<TaylorExpression>& variables);


    /*! \brief Equality operator. */
    bool operator==(const TaylorFunction& p) const;
    /*! \brief Inequality operator. */
    bool operator!=(const TaylorFunction& p) const;

    // Data access
    /*! \brief The accuracy parameter used to control approximation of the Taylor function. */
    shared_ptr<TaylorModel::Accuracy> accuracy_ptr() const;
    /*! \brief Set the accuracy parameter used to control approximation of the Taylor function. */
    void set_accuracy(shared_ptr<TaylorModel::Accuracy> acc);
    /*! \brief The data used to define the domain of the Taylor model. */
    const Vector<Interval>& domain() const;
    /*! \brief The centre of the Taylor model. */
    const Vector<Float> centre() const;
    /*! \brief The range of the Taylor model. */
    const Vector<Interval> range() const;
    /*! \brief The data used to define the centre of the Taylor model. */
    const Vector<TaylorModel>& models() const;

    /*! \brief The size of the argument. */
    uint argument_size() const;
    /*! \brief The size of the result. */
    uint result_size() const;

    /*! \brief The \a ith Taylor variable */
    TaylorExpression operator[](uint i) const;
    /*! \brief Evaluate the Taylor model at the point \a x. */
    Vector<Interval> evaluate(const Vector<Interval>& x) const;
    /*! \brief Evaluate the Taylor model at the point \a x. */
    Vector<Interval> evaluate(const Vector<Float>& x) const;
    /*! \brief Compute an approximation to Jacobian derivative of the Taylor model sat the point \a x. */
    Matrix<Interval> jacobian(const Vector<Interval>& x) const;

    /*! \brief Truncate to a model of lower order and/or smoothness, possibly on a different domain. */
    TaylorFunction truncate(ushort degree) const;

    /*! \brief The constant Taylor model with range \a r and argument domain \a d. */
    static TaylorFunction constant(const Vector<Interval>& d, const Vector<Interval>& r);
    /*! \brief The constant Taylor model with result \a c and argument domain \a d. */
    static TaylorFunction constant(const Vector<Interval>& d, const Vector<Float>& c);
    /*! \brief The identity Taylor model on domain \a d. */
    static TaylorFunction identity(const Vector<Interval>& d);

    /*! \brief Convert to an interval polynomial. */
    Vector< Polynomial<Interval> > polynomial() const;

    /*! \brief Truncate terms higher than \a bd. */
    TaylorFunction& truncate(const MultiIndexBound& bd);

    /*! \brief Write to an output stream. */
    std::ostream& write(std::ostream& os) const;

    /*! \brief Inplace addition. */
    friend TaylorFunction& operator+=(TaylorFunction& f, const TaylorFunction& g);
    /*! \brief Inplace subtraction. */
    friend TaylorFunction& operator-=(TaylorFunction& f, const TaylorFunction& g);
    /*! \brief Inplace addition. */
    friend TaylorFunction& operator+=(TaylorFunction& f, const Vector<Interval>& c);
    /*! \brief Inplace subtraction. */
    friend TaylorFunction& operator-=(TaylorFunction& f, const Vector<Interval>& c);
    /*! \brief Inplace scalar multiplication. */
    friend TaylorFunction& operator*=(TaylorFunction& f, const Float& c);
    /*! \brief Inplace scalar division. */
    friend TaylorFunction& operator/=(TaylorFunction& f, const Float& c);

    /*! \brief Negation. */
    friend TaylorFunction operator-(const TaylorFunction& f);
    /*! \brief Addition. */
    friend TaylorFunction operator+(const TaylorFunction& f1, const TaylorFunction& f2);
    /*! \brief Subtraction. */
    friend TaylorFunction operator-(const TaylorFunction& f1, const TaylorFunction& f2);

    /*! \brief Addition of a constant. */
    friend TaylorFunction operator+(const TaylorFunction& f, const Vector<Float>& c);
    /*! \brief Subtraction of a constant. */
    friend TaylorFunction operator-(const TaylorFunction& f, const Vector<Float>& c);
    /*! \brief Multiplication by a scalar. */
    friend TaylorFunction operator*(const Float& c, const TaylorFunction& f);
    /*! \brief Multiplication by a scalar. */
    friend TaylorFunction operator*(const TaylorFunction& f, const Float& c);
    /*! \brief Division by a scalar. */
    friend TaylorFunction operator/(const TaylorFunction& f, const Float& c);
    /*! \brief Addition of a constant. */
    friend TaylorFunction operator+(const TaylorFunction& f, const Vector<Interval>& c);
    /*! \brief Subtraction of a constant. */
    friend TaylorFunction operator-(const TaylorFunction& f, const Vector<Interval>& c);
    /*! \brief Multiplication by a scalar. */
    friend TaylorFunction operator*(const Interval& c, const TaylorFunction& f);
    /*! \brief Multiplication by a scalar. */
    friend TaylorFunction operator*(const TaylorFunction& f, const Interval& c);
    /*! \brief Division by a scalar. */
    friend TaylorFunction operator/(const TaylorFunction& f, const Interval& c);
    /*! \brief Multiplication by a matrix. */
    friend TaylorFunction operator*(const Matrix<Interval>& A, const TaylorFunction& f);

    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    friend TaylorFunction compose(const FunctionInterface& f, const TaylorFunction& g);
    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    friend TaylorExpression compose(const TaylorExpression& f, const TaylorFunction& g);
    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    friend TaylorFunction compose(const TaylorFunction& f, const TaylorFunction& g);
    //! \brief Antiderivative of \a f with respect to variable \a k.
    friend TaylorFunction antiderivative(const TaylorFunction& f, uint k);
    //! \brief The flow of the vector field \a vf defined over a space domain \a d over a time interval \a t.
    friend TaylorFunction flow(const TaylorFunction& vf, const Vector<Interval>& d, const Interval& t, uint o);
    //! \brief Compute the implicit function of \a f satisfying \f$f(c,h(c))=0\f$,
    //! where \f$c\f$ is the centre of the domain of \f$f\f$.
    friend TaylorFunction implicit(const TaylorFunction& f);
    //! \brief Compute the inverse function of \a f based at the centre of the domain. */
    friend TaylorFunction inverse(const TaylorFunction& f);
    //! \brief Compute the inverse function of \a f based at \f$f(c)\f$. */
    friend TaylorFunction inverse(const TaylorFunction& f, const Vector<Float>& c);
    //! \brief Compute the function \f$(f,g)(x)=(f(x),g(x))\f$.
    friend TaylorFunction join(const TaylorFunction& f, const TaylorFunction& g);
    friend TaylorFunction join(const TaylorFunction& f, const TaylorExpression& g);
    //! \brief Compute the function \f$(f\oplus g)(x,y)=(f(x),g(y))\f$.
    friend TaylorFunction combine(const TaylorFunction& f, const TaylorFunction& g);
    friend TaylorFunction combine(const TaylorFunction& f, const TaylorExpression& g);
    //! \brief Restrict the function \a f to a subdomain \a d.
    friend TaylorFunction restrict(const TaylorFunction& f, const Vector<Interval>& d);
    //! \brief Tests if a function \a f refines another function \a g.
    //! To be a refinement, the domain of \f f must contain the domain of \a g.
    friend bool refines(const TaylorFunction& f, const TaylorFunction& g);

    // For compatibility wit Vector.
    uint size() const { return this->result_size(); }
  private:
    array< array<Interval> > _powers(const Vector<Interval>&) const;
    void _compute_jacobian() const;
    void _set_argument_size(uint n);
    uint _compute_maximum_component_size() const;
    void _resize(uint rs, uint as, ushort d, ushort s);

  private:
    /* Domain of definition. */
    Vector<Interval> _domain;
    Vector<TaylorModel> _models;
};

TaylorFunction join(const TaylorFunction& f, const TaylorFunction& g);
TaylorFunction join(const TaylorFunction& f, const TaylorExpression& g);
TaylorFunction combine(const TaylorFunction& f, const TaylorFunction& g);
TaylorFunction combine(const TaylorFunction& f, const TaylorExpression& g);

std::ostream& operator<<(std::ostream&, const TaylorFunction&);

} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_H
