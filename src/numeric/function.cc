/***************************************************************************
 *            function.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include <mpfr.h>

#include "numeric/function.h"
#include "numeric/function.tpl"

#include "numeric/numerical_types.h"
#include "numeric/float64.h"
#include "numeric/mpfloat.h"

#include "real_typedef.h"

namespace Ariadne {
  namespace Numeric {

    int 
    mpfr_hypot(mpfr_t y, 
               const __mpfr_struct* x1, 
               const __mpfr_struct* x2, 
               mpfr_rnd_t r)
    {
      mpfr_t z;
      mpfr_init_set(z,y,GMP_RNDN);
      mpfr_mul(y,x1,x1,r);
      mpfr_mul(z,x2,x2,r);
      mpfr_add(y,y,z,r);
      mpfr_clear(z);
      return mpfr_sqrt(y,y,r);
    } 

  } 
}
