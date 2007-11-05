/***************************************************************************
 *            affine_variable.inline.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {

template<class X, class CV> inline
Function::AffineVariable<X,CV>::~AffineVariable() 
{
  ARIADNE_LOG(3,"AffineVariable<X,CV>::~AffineVariable()\n");
}


template<class X, class CV> inline
Function::AffineVariable<X,CV>::AffineVariable() 
  : _x(0), _dx()
{
}

template<class X, class CV> inline
Function::AffineVariable<X,CV>::AffineVariable(dimension_type ad) 
  :  _x(0), _dx(ad)
{
  ARIADNE_LOG(3,"AffineVariable<X,CV>::AffineVariable(dimension_type d)\n");
}


template<class X, class CV> inline
Function::AffineVariable<X,CV>::AffineVariable(const AffineVariable<X,CV>& av) 
  : _x(av._x), _dx(av._dx)
{
  ARIADNE_LOG(3,"AffineVariable<X,CV>::AffineVariable(AffineVariable<X,CV> av)\n");
}


template<class X, class CV> inline
Function::AffineVariable<X,CV>& 
Function::AffineVariable<X,CV>::operator=(const AffineVariable<X,CV>& av) 
{
  ARIADNE_LOG(3,"AffineVariable<X,CV>::operator=(AffineVariable<X,CV> av)\n");
  if(this!=&av) {
    this->_x=av._x;
    this->_dx=av._dx;
  }
  return *this;
}







template<class X, class CV> template<class XT> inline
void 
Function::AffineVariable<X,CV>::set(const XT& x) 
{
  this->_x=x;
}


template<class X, class CV> template<class XT> inline
void 
Function::AffineVariable<X,CV>::set(dimension_type j, const XT& x) 
{
  assert(j<this->_dx.size());
  this->_dx[j]=x;
}




}

