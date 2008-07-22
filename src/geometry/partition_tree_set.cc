/***************************************************************************
 *            partition_tree_set.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include "geometry/box.h"
#include "geometry/zonotope.h"
#include "geometry/partition_tree_set.h"
#include "geometry/partition_tree_set.code.h"

#include "numeric/float.h"

namespace Ariadne {
  

    
    
#ifdef ENABLE_FLOAT64
    template class PartitionScheme<Float64>;
    template class PartitionTree<Float64>;
    template class PartitionTreeCell<Float64>;
    template class PartitionTreeSet<Float64>;
#endif

#ifdef ENABLE_FLOATMP
    template class PartitionScheme<FloatMP>;
    template class PartitionTree<FloatMP>;
    template class PartitionTreeCell<FloatMP>;
    template class PartitionTreeSet<FloatMP>;
#endif

  
} // namespace Ariadne
