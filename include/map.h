/***************************************************************************
 *            map.h
 *
 *  Thu Apr 29 14:33:39 2004
 *  Copyright  2004  Alberto Casagrande
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
 
#ifndef _MAP_H
#define _MAP_H

#include "linear_algebra.h"
#include "approx_type.h"
#include "classes.h"

namespace Ariadne {

/*! \class Map
 *  \brief Represents a (multivalued/nondeterministic) map of the state space..
 * 
 * At the very least, the map f should be able to compute an outer 
 * approximation to a basic closed set \a A with an error of at most \a e.
 */
class Map {
	public:
		virtual ~Map() = 0;
		
		virtual BasicSet* operator() (const BasicSet& A, 
				ApproxType atype) const = 0;
		
		virtual ASet* operator() (const ASet& A) const = 0;
		
		virtual Set* operator() (const Set& A) const = 0;
		
		virtual HSet* operator() (const HSet& A) const = 0;

		virtual bool is_linear() const { return false; }
};

class LinearMap: public Map {
	
		Matrix *A;
		Vector *b;
			
	public:
		LinearMap(Matrix *A, Vector* b); 

		~LinearMap(); 

		BasicSet* operator() (const BasicSet& A, 
				ApproxType atype) const;

		ASet* operator() (const ASet& A) const;

		Set* operator() (const Set& A) const;
		
		HSet* operator() (const HSet& A) const;
		
		const Matrix* get_A() const;
		
		const Vector* get_b() const;
		
		bool is_linear() const { return true; }
};

}

#endif /* _MAP_H */

