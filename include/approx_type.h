/***************************************************************************
 *            approx_type.h
 *
 *  Wed Apr 21 12:37:01 2004
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
 
#ifndef _APPROX_TYPE_H
#define _APPROX_TYPE_H

namespace Ariadne {
	
/*! \brief Defines the approximation type */ 
typedef enum ApproxType{
	OVER, /*!< Over-approximation */
	UNDER, /*!< Under-approximation */
	NONE /*!< No approximation */
};

ApproxType invert_approx(ApproxType atype); 

}

#endif /* _APPROX_TYPE_H */
