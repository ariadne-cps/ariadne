/***************************************************************************
 *            partition_tree_set.cc
 *
 *  1 July 2006
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

#include "geometry/parallelotope.h"
#include "geometry/partition_tree_set.h"
#include "geometry/partition_tree_set.tpl"

#include "real_typedef.h"

namespace Ariadne {
  namespace Geometry {

    template class PartitionScheme<Real>;
    template class PartitionTree<Real>;
    template class PartitionTreeCell<Real>;
    template class PartitionTreeSet<Real>;

    template PartitionTreeSet<Real> outer_approximation(const Parallelotope<Real>&, const PartitionScheme<Real>&, const uint);
    template PartitionTreeSet<Real> inner_approximation(const Parallelotope<Real>&, const PartitionScheme<Real>&, const uint);
    template PartitionTreeSet<Real> over_approximation(const Parallelotope<Real>&, const PartitionScheme<Real>&, const uint);
    template PartitionTreeSet<Real> under_approximation(const Parallelotope<Real>&, const PartitionScheme<Real>&, const uint);
    
    template PartitionTreeSet<Real> outer_approximation(const GridMaskSet<Real>&, const PartitionScheme<Real>&, const uint);
    template PartitionTreeSet<Real> inner_approximation(const GridMaskSet<Real>&, const PartitionScheme<Real>&, const uint);
    template PartitionTreeSet<Real> over_approximation(const GridMaskSet<Real>&, const PartitionScheme<Real>&, const uint);
    template PartitionTreeSet<Real> under_approximation(const GridMaskSet<Real>&, const PartitionScheme<Real>&, const uint);
    
    template std::ostream& operator<<(std::ostream&, const PartitionScheme<Real>&);
    template std::ostream& operator<<(std::ostream&, const PartitionTree<Real>&);
    template std::ostream& operator<<(std::ostream&, const PartitionTreeCell<Real>&);
    template std::ostream& operator<<(std::ostream&, const PartitionTreeSet<Real>&);
  }
}
