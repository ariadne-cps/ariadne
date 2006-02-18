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

#include "geometry/partition_tree_set.h"
#include "geometry/partition_tree_set.tpl"

namespace Ariadne {
  namespace Geometry {

    template class PartitionScheme<Dyadic>;
    template class PartitionTree<Dyadic>;
    template class PartitionTreeCell<Dyadic>;
    template class PartitionTreeSet<Dyadic>;

    template class PartitionTreeSetIterator<Dyadic>;

    template std::ostream& operator<<(std::ostream&, const PartitionScheme<Dyadic>&);
    template std::ostream& operator<<(std::ostream&, const PartitionTree<Dyadic>&);
    template std::ostream& operator<<(std::ostream&, const PartitionTreeCell<Dyadic>&);
    template std::ostream& operator<<(std::ostream&, const PartitionTreeSet<Dyadic>&);
  }
}
