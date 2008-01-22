/***************************************************************************
 *            affine_variable.code.h
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
 
#include "affine_variable.h"

namespace Ariadne {


template<class X>
std::ostream& 
Function::operator<<(std::ostream& os, const AffineVariable<X>& av)
{
  os << "[" << av.value();
  for(uint i=0; i!=av.argument_size(); ++i) {
    os << (i==0?";":",") << av.data()[i+1u]; 
  }
  return os << "]";
}


template<class X>
void
Function::AffineVariable<X>::instantiate()
{
  AffineVariable<X>* av=0;
  std::ostream* os=0;

  operator<<(*os,*av);
}

}
