/***************************************************************************
 *            <vector>.cc
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
 
#include "real_typedef.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/vector.tpl"

namespace Ariadne {
  namespace LinearAlgebra {
    
    template class Vector<Real>;
    template std::ostream& operator<<(std::ostream&, const Vector<Real>&);
    template std::istream& operator>>(std::istream&, Vector<Real>&);

#ifndef REAL_IS_A_FIELD
    template class Vector<Field>;
    template std::ostream& operator<<(std::ostream&, const Vector<Field>&);
    template std::istream& operator>>(std::istream&, Vector<Field>&);
#endif
  }
}
