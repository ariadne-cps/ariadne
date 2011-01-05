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
    : public AlgebraMixin<Algebra<X>,X>
{
  private:
    boost::shared_ptr< AlgebraInterface<X> > _ptr;
  public:
    typedef X ScalarType;
    Algebra() : _ptr() { }
    Algebra(AlgebraInterface<X>* p) : _ptr(p) { }
    operator const AlgebraInterface<X>& () const { return *_ptr; }
    Algebra<X> null() const { return Algebra<X>(_ptr->create()); }
    Algebra<X>& operator=(const X& c) { *this = this->create(); this->iadd(c); return *this; }
    Algebra<X>& operator=(int c) { return this->operator=(static_cast<X>(c)); }
    void clear() { this->imul(0); }
    std::ostream& write(std::ostream& os) const { return _ptr->write(os); }
  public:
    void iadd(const X& c) { _ptr->iadd(c); }
    void imul(const X& c) { _ptr->imul(c); }
    void isma(const X& c, const Algebra<X>& x) { _ptr->isma(c,*x._ptr); }
    void ifma(const Algebra<X>& x1, const Algebra<X>& x2) { _ptr->ifma(*x1._ptr,*x2._ptr); }
    Algebra<X> sma(const X& c, const Algebra<X>& x) const { Algebra<X> r(x); r._ptr->isma(c,*x._ptr); return r; }
};

template<class X> class NormedAlgebra
    : public NormedAlgebraMixin<NormedAlgebra<X>,X>
{
  private:
    boost::shared_ptr< NormedAlgebraInterface<X> > _ptr;
  public:
    typedef X ScalarType;
    NormedAlgebra(NormedAlgebraInterface<X>* p) : _ptr(p) { }
    operator const NormedAlgebraInterface<X>& () const { return *_ptr; }
    NormedAlgebra<X> null() const { return NormedAlgebra<X>(_ptr->create()); }
    Float norm() const { return _ptr->norm(); }
    Float average() const { return _ptr->average(); }
    void clear() { this->imul(0); }
    void iadd(const X& c) { _ptr->iadd(c); }
    void imul(const X& c) { _ptr->imul(c); }
    void isma(const X& c, const NormedAlgebra<X>& x) { _ptr->isma(c,*x._ptr); }
    void ifma(const NormedAlgebra<X>& x1, const NormedAlgebra<X>& x2) { _ptr->ifma(*x1._ptr,*x2._ptr); }
    NormedAlgebra<X> sma(const X& c, const NormedAlgebra<X>& x) const { NormedAlgebra<X> r(x); r._ptr->isma(c,*x._ptr); return r; }
    std::ostream& write(std::ostream& os) const { return _ptr->write(os); }
};



template<class X> std::ostream& operator<<(std::ostream& os, const Algebra<X>& x) {
    return x.write(os); }
template<class X> std::ostream& operator<<(std::ostream& os, const NormedAlgebra<X>& x) {
    return x.write(os); }

template<class X> Algebra<X> null(const Algebra<X>& a) { return a.create(); }
template<class X> Algebra<X> copy(const Algebra<X>& a) { return a.clone(); }





} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_H */