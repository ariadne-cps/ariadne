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
 *  \brief Over-approximations of functions based on Taylor expansions.
 */

#ifndef ARIADNE_TAYLOR_MODEL_H
#define ARIADNE_TAYLOR_MODEL_H

#include <iosfwd>
#include "numeric.h"
#include "vector.h"
#include "sparse_differential.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;

class FunctionInterface;

class MultiIndex;
template<class X> class SparseDifferential;

class TaylorModel;

TaylorModel recentre(const TaylorModel&, const Vector<Interval>& bx, const Vector<Float>& pt);
TaylorModel truncate(const TaylorModel&, const Vector<Interval>&, uint, uint);

TaylorModel compose(const TaylorModel&, const TaylorModel&);
TaylorModel inverse(const TaylorModel&, const Vector<Float>&);
TaylorModel implicit(const TaylorModel&, const Vector<Float>&);
TaylorModel derivative(const TaylorModel&, uint);




/* \brief A taylor_model with multivalued output using the TaylorVariable class.
 */
class TaylorModel {
    typedef Float R;
    typedef Interval I;
  public:
    /*! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. */
    TaylorModel();
    /*! \brief The zero Taylor model in \a as variables with size \a rs image, order \a o and smoothness \a s, defined on the whole space with centre at the origin. */
    TaylorModel(uint rs, uint as, ushort o, ushort s);
  
    /*! \brief Construct from a domain and the expansionn. */
    TaylorModel(const Vector<Interval>& domain,
                const Vector<TaylorVariable>& expansion);
  
    /*! \brief Construct from a domain and the expansionn. */
    TaylorModel(const Vector<Interval>& domain,
                const Vector< SparseDifferential<Float> >& expansion);
  
    /*! \brief Construct from a domain, centre, an order and a function. */
    TaylorModel(const Vector<Interval>& domain, 
                ushort order, ushort smoothness,
                const FunctionInterface& function);
  
  
    /*! \brief Equality operator. */
    bool operator==(const TaylorModel& p) const;
    /*! \brief Inequality operator. */
    bool operator!=(const TaylorModel& p) const;
  
    // Data access
    /*! \brief The data used to define the domain of the Taylor model. */
    const Vector<Interval>& domain() const;
    /*! \brief The centre of the Taylor model. */
    const Vector<Float> centre() const;
    /*! \brief The range of the Taylor model. */
    const Vector<Interval> range() const;
    /*! \brief The data used to define the centre of the Taylor model. */
    const Vector<TaylorVariable>& expansion() const;
  
    /*! \brief The size of the argument. */
    uint argument_size() const;
    /*! \brief The size of the result. */
    uint result_size() const;

  
    /// Resizing
    void resize(uint rs, uint as, ushort d, ushort s);

    /*! \brief Evaluate the Taylor model at the point \a x. */
    Vector<Interval> evaluate(const Vector<Interval>& x) const;
    Vector<Interval> evaluate(const Vector<Float>& x) const;
  
    /*! \brief Compute an approximation to Jacobian derivative of the Taylor model at the point \a x. */
    Matrix<Float> jacobian(const Vector<Float>& x) const;
  
    /*! \brief Truncate to a model of lower order and/or smoothness, possibly on a different domain. */
    TaylorModel truncate(ushort degree) const;
  
    /*! \brief The constant Taylor model with result size 1 and argument size \a as. */
    static TaylorModel constants(uint as, const Vector<Float>& c);
     /*! \brief The constant Taylor model with result size 1 and argument size \a as. */
    static TaylorModel variables(uint as, const Vector<Float>& x);
 
    /*! \brief Write to an output stream. */
    std::ostream& write(std::ostream& os) const;
  
    /*! \brief Addition. */
    friend TaylorModel operator+(const TaylorModel&, const TaylorModel&);
    /*! \brief Subtraction. */
    friend TaylorModel operator-(const TaylorModel&, const TaylorModel&);
  
    /*! \brief Multiplication by a scalar. */
    friend TaylorModel operator*(const Float&, const TaylorModel&);
    /*! \brief Multiplication by a scalar. */
    friend TaylorModel operator*(const TaylorModel&, const Float&);
    /*! \brief Division by a scalar. */
    friend TaylorModel operator/(const TaylorModel&, const Float&);
  
    /*! \brief Composition \f$p\circ q(x)=p(q(x))\f$. */
    friend TaylorModel compose(const TaylorModel&, const TaylorModel&);
    /*! \brief Derivative with respect to variable \a k. */
    friend TaylorModel antiderivative(const TaylorModel&, uint k);
  private:
    array< array<Interval> > _powers(const Vector<Interval>&) const;
    void _compute_jacobian() const;
    void _set_argument_size(uint n);
    uint _compute_maximum_component_size() const;
  private:
    friend TaylorModel recentre(const TaylorModel&, const Vector<Interval>& bx);
    friend TaylorModel inverse(const TaylorModel&, const Vector<Float>&);
    //friend void add(TaylorModel&,const TaylorModel&,const TaylorModel&);
    //friend void sub(TaylorModel&,const TaylorModel&,const TaylorModel&);
    //friend void mul(TaylorModel&,const TaylorModel&,const TaylorModel&);
    //friend void div(TaylorModel&,const TaylorModel&,const TaylorModel&);
    //friend void compose(TaylorModel&,const TaylorModel&,const TaylorModel&);
    //friend void scale(TaylorModel&,const R&);
  private:
    /* Domain of definition. */
    Vector<Interval> _domain;
    Vector<TaylorVariable> _expansion;
};


TaylorModel combine(const TaylorModel&, const TaylorModel&);
TaylorModel join(const TaylorModel&, const TaylorModel&);

std::ostream& operator<<(std::ostream&, const TaylorModel&);

} // namespace Ariadne

#endif // ARIADNE_TAYLOR_MODEL_H
