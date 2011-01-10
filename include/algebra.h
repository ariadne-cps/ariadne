/***************************************************************************
 *            algebra.h
 *
 *  Copyright 2010  Pieter Collins
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

/*! \file algebra.h
 *  \brief Algebras and normed algebras.
 */
#ifndef ARIADNE_ALGEBRA_H
#define ARIADNE_ALGEBRA_H

#include <iosfwd>
#include <iostream>
#include "numeric.h"
#include "pointer.h"
#include "algebra_interface.h"
#include "algebra_mixin.h"

namespace Ariadne {

//! \brief Generic class for elements of unital algebras.
template<class X> class Algebra
    : public AlgebraOperators<Algebra<X>,X>
{
  private:
    boost::shared_ptr< AlgebraInterface<X> > _ptr;
  public:
    typedef X ScalarType;
    typedef typename X::NumericType NumericType;
    Algebra() : _ptr() { }
    Algebra(AlgebraInterface<X>* p) : _ptr(p) { }
    Algebra(const AlgebraInterface<X>& a) : _ptr(a._clone()) { }
    Algebra(const Algebra<X>& a) : _ptr(a._ptr->_clone()) { }
    Algebra<X>& operator=(int c) { *this = this->create(); this->iadd(c); return *this; }
    Algebra<X>& operator=(const X& c) { *this = this->create(); this->iadd(c); return *this; }
    Algebra<X>& operator=(const Algebra<X>& a) { this->_ptr=boost::shared_ptr< AlgebraInterface<X> >(a._ptr->_clone()); return *this; }
    operator const AlgebraInterface<X>& () const { return *_ptr; }
    Algebra<X> create() const { return Algebra<X>(_ptr->_create()); }
    Algebra<X> clone() const { return Algebra<X>(_ptr->_clone()); }
    std::ostream& write(std::ostream& os) const { return _ptr->write(os); }
  public:
    void iadd(const X& c) { _ptr->_iadd(c); }
    void imul(const X& c) { _ptr->_imul(c); }
    void isma(const X& c, const Algebra<X>& x) { _ptr->_isma(c,*x._ptr); }
    void ifma(const Algebra<X>& x1, const Algebra<X>& x2) { _ptr->_ifma(*x1._ptr,*x2._ptr); }
};

template<class X> class NormedAlgebra
    : public NormedAlgebraOperators<NormedAlgebra<X>,X>
{
  private:
    boost::shared_ptr< NormedAlgebraInterface<X> > _ptr;
  public:
    typedef X ScalarType;
    NormedAlgebra(NormedAlgebraInterface<X>* p) : _ptr(p) { }
    NormedAlgebra(const NormedAlgebraInterface<X>& a) : _ptr(a.clone()) { }
    operator const NormedAlgebraInterface<X>& () const { return *_ptr; }
    Float norm() const { return _ptr->norm(); }
    Float average() const { return _ptr->average(); }
    NormedAlgebra<X> create() const { return NormedAlgebra<X>(_ptr->_create()); }
    NormedAlgebra<X> clone() const { return NormedAlgebra<X>(_ptr->_clone()); }
    void clear() { this->imul(0); }
    void iadd(const X& c) { _ptr->_iadd(c); }
    void imul(const X& c) { _ptr->_imul(c); }
    void isma(const X& c, const NormedAlgebra<X>& x) { _ptr->_isma(c,*x._ptr); }
    void ifma(const NormedAlgebra<X>& x1, const NormedAlgebra<X>& x2) { _ptr->_ifma(*x1._ptr,*x2._ptr); }
    std::ostream& write(std::ostream& os) const { return _ptr->write(os); }
};

template<class X> std::ostream& operator<<(std::ostream& os, const Algebra<X>& x) {
    return x.write(os); }
template<class X> std::ostream& operator<<(std::ostream& os, const NormedAlgebra<X>& x) {
    return x.write(os); }

template<class X> Algebra<X> create(const Algebra<X>& a) { return a.create(); }
template<class X> Algebra<X> copy(const Algebra<X>& a) { return a.clone(); }

template<class X> Algebra<X> operator/(const Algebra<X>&,const Algebra<X>&) { ARIADNE_NOT_IMPLEMENTED; }
template<class X> Algebra<X> rec(const Algebra<X>&) { ARIADNE_NOT_IMPLEMENTED; }
template<class X> Algebra<X> pow(const Algebra<X>&, int) { ARIADNE_NOT_IMPLEMENTED; }
template<class X> Algebra<X> sqrt(const Algebra<X>&) { ARIADNE_NOT_IMPLEMENTED; }
template<class X> Algebra<X> exp(const Algebra<X>&) { ARIADNE_NOT_IMPLEMENTED; }
template<class X> Algebra<X> log(const Algebra<X>&) { ARIADNE_NOT_IMPLEMENTED; }
template<class X> Algebra<X> sin(const Algebra<X>&) { ARIADNE_NOT_IMPLEMENTED; }
template<class X> Algebra<X> cos(const Algebra<X>&) { ARIADNE_NOT_IMPLEMENTED; }
template<class X> Algebra<X> tan(const Algebra<X>&) { ARIADNE_NOT_IMPLEMENTED; }



} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_H */