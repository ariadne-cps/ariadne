/***************************************************************************
 *            grid_paving.cc
 *
 *  Copyright  2008  Ivan S. Zapreev, Pieter Collins
 *  ivan.zapreev@gmail.com, Pieter.Collins@cwi.nl
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "geometry/grid_paving.h"
#include "geometry/grid_paving.code.h"

namespace Ariadne {

#ifdef ENABLE_FLOAT64
	template class GridPavingCell<Float64>;
	template class GridPavingCursor<Float64>;
	template class GridSubPaving<Float64>;
	template class GridPaving<Float64>;
#endif

#ifdef ENABLE_FLOATMP
	template class GridPavingCell<FloatMP>;
	template class GridPavingCursor<FloatMP>;
	template class GridSubPaving<FloatMP>;
	template class GridPaving<FloatMP>;
#endif

} // namespace Ariadne
