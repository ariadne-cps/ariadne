/***************************************************************************
 *            vector_interface.hpp
 *
 *  Copyright 2010-17  Pieter Collins
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
 
/*! \file vector_interface.hpp
 *  \brief Interface for vectors.
 */

#ifndef ARIADNE_VECTOR_INTERFACE_HPP
#define ARIADNE_VECTOR_INTERFACE_HPP

#include "../utility/writable.hpp"

namespace Ariadne {

template<class X, class I=Void> class VectorInterface;

template<class X, class I> class VectorInterface {
    typedef I ScalarInterface;
    typedef X ScalarType;
    typedef typename ScalarType::NumericType NumericType;
  public:
    virtual ScalarInterface* _get(SizeType i) const = 0;
    virtual Void _set(SizeType i, ScalarType const& x) const = 0;
    virtual VectorInterface* _add(VectorInterface const*) const = 0;
    virtual VectorInterface* _sub(VectorInterface const*) const = 0;
    virtual VectorInterface* _mul(ScalarType const&) const = 0;
    virtual VectorInterface* _mul(NumericType const&) const = 0;
};

} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_INTERFACE_HPP */
