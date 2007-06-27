/***************************************************************************
 *            lputil.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
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

/*! \file lputil.h
 *  \brief Linear programming utility routines.
 */

#ifndef ARIADNE_LPUTIL_H
#define ARIADNE_LPUTIL_H

#include "exceptions.h"

namespace Ariadne {
  namespace LinearProgramming {
    
    // Returns the index of the permuation list pptr of size sz with requested value
    int getindex(uint* pptr, uint sz, uint value);
      

    template<class R>
    LinearAlgebra::Matrix<R> 
    to_tableau(const LinearAlgebra::Matrix<R>& A, 
               const LinearAlgebra::Vector<R>& b, 
               const LinearAlgebra::Vector<R>& c);


    // Modify the tableau
    // FIXME: Should enter and leave be uint?
    template<class AP>
    void 
    pivot_tableau(uint m, uint n,
                  AP* Aptr, uint Arinc, uint Acinc,
                  AP* bptr, uint binc,
                  AP* cptr, uint cinc,
                  AP* dptr,
                  uint* pptr,
                  int enter, int leave);


    
  } // namespace LinearProgramming
} // namespace Ariadne

#include "lputil.inline.h"
#include "lputil.template.h"

#endif /* ARIADNE_LPUTIL_H */
