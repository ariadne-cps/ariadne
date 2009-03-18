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

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;

class FunctionInterface;
class MultiIndex;
class TaylorFunction;

TaylorFunction recentre(const TaylorFunction&, const Vector<Interval>& bx, const Vector<Float>& pt);
TaylorFunction truncate(const TaylorFunction&, const Vector<Interval>&, uint, uint);

TaylorFunction compose(const TaylorFunction&, const TaylorFunction&);
TaylorFunction inverse(const TaylorFunction&);
TaylorFunction implicit(const TaylorFunction&);
TaylorFunction antiderivative(const TaylorFunction&, uint);



/*! \brief A taylor_model with multivalued output using the TaylorVariable class.
 */
class TaylorFunction {
    typedef Float R;
    typedef Interval I;
  public:
    /*! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. */
    TaylorFunction();
    /*! \brief The zero Taylor model in \a as variables with size \a rs image, defined on the unit box. */
    TaylorFunction(uint rs, uint as);
    /*! \brief The zero Taylor model in \a as variables with size \a rs image, order \a o and smoothness \a s, defined on the unit box. */
    TaylorFunction(uint rs, uint as, ushort o, ushort s);
  
    /*! \brief Construct from a domain and the expansionn. */
    TaylorFunction(const Vector<Interval>& domain,
                   const Vector<TaylorVariable>& expansion);
  
    /*! \brief Construct from a domain and a function. */
    TaylorFunction(const Vector<Interval>& domain,
                   const FunctionInterface& function);
  
    /*! \brief Construct from a domain, centre, an order and a function. */
    TaylorFunction(const Vector<Interval>& domain, 
                   ushort order, ushort smoothness,
                   const FunctionInterface& function);
  
  
    /*! \brief Equality operator. */
    bool operator==(const TaylorFunction& p) const;
    /*! \brief Inequality operator. */
    bool operator!=(const TaylorFunction& p) const;
  
    // Data access
    /*! \brief The data used to define the domain of the Taylor model. */
    const Vector<Interval>& domain() const;
    /*! \brief The centre of the Taylor model. */
    const Vector<Float> centre() const;
    /*! \brief The range of the Taylor model. */
    const Vector<Interval> range() const;
    /*! \brief The data used to define the centre of the Taylor model. */
    const Vector<TaylorVariable>& variables() const;
    const Vector<TaylorVariable>& expansion() const;
  
    /*! \brief The size of the argument. */
    uint argument_size() const;
    /*! \brief The size of the result. */
    uint result_size() const;
    
  
    /*! \brief Evaluate the Taylor model at the point \a x. */
    Vector<Interval> evaluate(const Vector<Interval>& x) const;
    /*! \brief Evaluate the Taylor model at the point \a x. */
    Vector<Interval> evaluate(const Vector<Float>& x) const;
  
    /*! \brief Compute an approximation to Jacobian derivative of the Taylor model at the point \a x. */
    Matrix<Float> jacobian(const Vector<Float>& x) const;
  
    /*! \brief Truncate to a model of lower order and/or smoothness, possibly on a different domain. */
    TaylorFunction truncate(ushort degree) const;
  
    /*! \brief The constant Taylor model with result \a c and argument domain \a d. */
    static TaylorFunction constant(const Vector<Interval>& d, const Vector<Float>& c);
    /*! \brief The identity Taylor model on domain \a d. */
    static TaylorFunction identity(const Vector<Interval>& d);
 
    /*! \brief Write to an output stream. */
    std::ostream& write(std::ostream& os) const;
  
    /*! \brief Addition. */
    friend TaylorFunction operator+(const TaylorFunction& f, const TaylorFunction& c);
    /*! \brief Subtraction. */
    friend TaylorFunction operator-(const TaylorFunction& f, const TaylorFunction& c);
  
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
  
    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    friend TaylorFunction compose(const FunctionInterface& f, const TaylorFunction& tg);
    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    friend TaylorFunction compose(const TaylorFunction& f, const TaylorFunction& g);
    //! \brief Antiderivative of \a f with respect to variable \a k.
    friend TaylorFunction antiderivative(const TaylorFunction& f, uint k);
    //! \brief The flow of the vector field \a vf defined over a space domain \a d over a time interval \a t.
    friend TaylorFunction flow(const TaylorFunction& vf, const Vector<Interval>& d, const Interval& t);
    //! \brief Compute the implicit function of \a f satisfying \f$f(c,h(c))=0\f$,
    //! where \f$c\f$ is the centre of the domain of \f$f\f$. 
    friend TaylorFunction implicit(const TaylorFunction& f);
    //! \brief Compute the inverse function of \a f based at the centre of the domain. */
    friend TaylorFunction inverse(const TaylorFunction& f);
    //! \brief Compute the inverse function of \a f based at \f$f(c)\f$. */
    friend TaylorFunction inverse(const TaylorFunction& f, const Vector<Float>& c);
    //! \brief Compute the function \f$(f,g)(x)=(f(x),g(x))\f$.
    friend TaylorFunction join(const TaylorFunction& f, const TaylorFunction& g);
    //! \brief Compute the function \f$(f\oplus g)(x,y)=(f(x),g(y))\f$.
    friend TaylorFunction combine(const TaylorFunction& f, const TaylorFunction& g);
    //! \brief Restrict the function \a f to a subdomain \a d.
    friend TaylorFunction restrict(const TaylorFunction& f, const Vector<Interval>& d);

  private:
    array< array<Interval> > _powers(const Vector<Interval>&) const;
    void _compute_jacobian() const;
    void _set_argument_size(uint n);
    uint _compute_maximum_component_size() const;
    void _resize(uint rs, uint as, ushort d, ushort s);

  private:
    /* Domain of definition. */
    Vector<Interval> _domain;
    Vector<TaylorVariable> _expansion;
};


std::ostream& operator<<(std::ostream&, const TaylorFunction&);

} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_H
