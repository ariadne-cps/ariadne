/***************************************************************************
 *            approx_type.cpp
 *
 *  Thu Aug 24 18:08:01 2004
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
 
#include "approx_type.h"
	
using namespace Ariadne;
	

ApproxType invert_approx(ApproxType atype) {
	switch(atype) {
		case OVER:
		default:
			return UNDER;
			break;
		case UNDER:
			return OVER;
			break;
		case NONE:
			return NONE;
			break;
	}

}

