/***************************************************************************
 *            map.cpp
 *
 *  Fri Aug  6 10:50:39 2004
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
 

#include "map.h"
#include "basic_set.h"
#include "set.h"

using namespace Ariadne;
	
LinearMap::LinearMap(Matrix *A, Vector* b) {
	this->b=b;
	this->A=A;
}

LinearMap::~LinearMap() {
	delete this->b;
	delete this->A;

	// useless, but safer
	this->b=NULL;
	this->A=NULL;
}

BasicSet* LinearMap::operator() (const BasicSet& A, ApproxType atype) const {
	return A.apply(*this,atype);
}

ASet* LinearMap::operator() (const ASet& A) const {
	return A.apply(*this);
}

Set* LinearMap::operator() (const Set& A) const {
	return A.apply(*this);
}
		
HSet* LinearMap::operator() (const HSet& A) const {
	return A.apply(*this);
}

const Matrix* LinearMap::get_A() const {
	return A;	
}
		
const Vector* LinearMap::get_b() const {
	return b;
}
