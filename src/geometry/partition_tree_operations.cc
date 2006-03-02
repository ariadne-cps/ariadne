/***************************************************************************
 *            partition_tree_operations.cc
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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

#include "geometry/partition_tree_operations.h"

Ariadne::sequence<Ariadne::dimension_type>
Ariadne::Geometry::default_subdivision_coordinates(dimension_type n) {
  dimension_type coords[n];
  for(dimension_type i=0; i!=n; ++i) {
    coords[i]=i;
  }
  return sequence<dimension_type>(coords,coords,coords+n);
}

Ariadne::dimension_type 
Ariadne::Geometry::compute_dimension(const sequence<dimension_type>& ss) 
{
  Ariadne::dimension_type result=0;
  for(size_type i=0; i!=ss.body_size()+ss.tail_size(); ++i) {
    if(ss[i]>=result) {
      result=ss[i]+1;
    }
  }
  return result;
}
