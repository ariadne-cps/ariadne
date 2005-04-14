/***************************************************************************
 *            numerival_type.h
 *
 *  Thu Oct 02 16:27:05 2004
 *  Copyright  2004  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#ifndef _NUMERICAL_TYPE_H
#define _NUMERICAL_TYPE_H

#include <gmpxx.h>

namespace Ariadne {
	
typedef mpz_class AriadneIntegerType;
typedef mpq_class AriadneRationalType;
typedef mpf_class AriadneDyadeticType;

	
inline AriadneRationalType transform_into_rational(
			const AriadneRationalType &num){ return num; }

inline AriadneIntegerType numerator(const AriadneRationalType 
		&num){ return num.get_num(); }

inline AriadneIntegerType denumerator(const AriadneRationalType 
		&num){ return num.get_den();}


		
inline AriadneRationalType epsilon(const AriadneRationalType &a) {
	
	AriadneRationalType e(1,1000);
	
	return e*e*e;
}

		
template <typename NumType>
NumType AriadneGCD(const NumType &a, const NumType &b){
	
	NumType c=a%b;
	
	if (c==0) return b;
		
	return (AriadneGCD(b, c));
}

template <typename NumType>
inline NumType AriadneLCM(const NumType &a, const NumType &b) {
	return ((a*b)/AriadneGCD(a,b));
}

}

#endif
