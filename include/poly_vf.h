/***************************************************************************
 *            poly_vf.h
 *
 *  Mon Feb  7 11:42:22 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
#ifndef _POLY_VF_H
#define _POLY_VF_H

#include <approx_type.h>

namespace Ariadne {	
namespace VectorField{
namespace Affine{

template <typename S>
class PolyAffineIntegrator {	
	
	public:
		typedef S State;
		typedef typename State::Real Real;
		typedef typename Ariadne::Geometry::Polyhedron< State > BasicSet;
		typedef typename Ariadne::Geometry::Polyhedron< State > Polyhedron;
		typedef typename Ariadne::Geometry::DenotableSet< Polyhedron > DenotableSet;
	
		typedef typename Ariadne::Map::Affine::PolyAffineMap< State > Map;
		typedef typename Ariadne::Map::Affine::AffineMap< Map > SolutionMap;
	
		PolyAffineIntegrator() {}
		
		inline BasicSet get_flow_tube_from_to(const BasicSet &A, 
				const BasicSet &B, const SolutionMap &sol_map, 
				const Ariadne::Geometry::ApproxKind &atype) {
					
				return Ariadne::Geometry::convex_hull(A,B);		
		}
};
	
}
}
}

#endif /* _POLY_VF_H */
