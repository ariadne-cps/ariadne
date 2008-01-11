/***************************************************************************
 *            taylor_model.inline.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#include "sorted_index.h"
#include "multi_index.h"

namespace Ariadne {


template<class R> template<class RR> inline
Function::TaylorModel<R>::TaylorModel(const TaylorModel<RR>& p)
  : _domain(p.domain()),
    _centre(p.centre()),
    _smoothness(p.smoothness()),
    _derivatives(p.derivatives())
{
}
        




} // namespace Ariadne
