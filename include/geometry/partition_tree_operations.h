/***************************************************************************
 *            partition_tree_operations.h
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

/*! \file partition_tree_operations.h
 *  \brief General operations on partition trees.
 */

#ifndef _ARIADNE_PARTITION_TREE_OPERATIONS_H
#define _ARIADNE_PARTITION_TREE_OPERATIONS_H

#include "base/basic_type.h"

#include "geometry/geometry_declarations.h"

namespace Ariadne {
  namespace Geometry {
    SubdivisionSequence default_subdivision_coordinates(dimension_type);
    dimension_type compute_dimension(const SubdivisionSequence&);
  }
}

#endif /* _ARIADNE_PARTITION_TREE_OPERATIONS_H */
