/***************************************************************************
 *            linear_algebra.cpp
 *
 *  Thu Aug  5 12:11:15 2004
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
 

#include "linear_algebra.h"
	
using namespace Ariadne;

Matrix *exp_Ah(const Matrix &A, const double h, const unsigned int n) {

#ifdef USE_OCTAVE
	DiagMatrix I(A.dim1(),A.dim2(),1.0);
#else 
	Identity_Matrix I(A.size1());
#endif
	
	Matrix tmp(I),e_Ah(I);
	
	/* tmp = \frac{h^{0}}{0!}*A^{0} = I
	 * and e_Ah = \Sum_{j=0}^{0}\frac{h^j}{j!}*A^{j} = I */
	
	for (unsigned int i=1; i< n; i++) {
		/* tmp = \frac{h^{i-1}}{(i-1)!}*A^{i-1}
		 * and e_Ah = \Sum_{j=0}^{i-1}\frac{h^j}{j!}*A^{j} */
#ifdef USE_OCTAVE
		tmp= (h/i)*(tmp * A);
#else 
		tmp *= (h/i);
		tmp = prec_prod(tmp,A);
#endif
		/* tmp =  (h^i/i!)*A^i */
		e_Ah += tmp;
		/*  e_Ah = \Sum_{j=0}^{i}\frac{h^j}{j!}*A^{j} */
	}

	Matrix *out= new Matrix(e_Ah);
	return out;
}

Vector *exp_b(const Matrix &A, const Vector &b, const double h, 
			const unsigned int n) {
	
#ifdef USE_OCTAVE
	DiagMatrix I(A.dim1(),A.dim2(),1.0);
#else 
	Identity_Matrix I(A.size1());
#endif
	Matrix tmp(I),e_b(I);
	
	/* tmp = \frac{h^{0}}{1!}*A^{0}
	 * and e_b = \Sum_{j=0}^{0}\frac{h^j}{(j+1)!}*A^{j} */
	
	for (unsigned int i=1; i< n-1; i++) {
		/* tmp = \frac{h^{i-1}}{i!}*A^{i-1}
		 * and e_b = \Sum_{j=0}^{i-1}\frac{h^j}{(j+1)!}*A^{j} */
#ifdef USE_OCTAVE
		tmp= (h/(i+1))*(tmp * A);
#else 
		tmp *= (h/(i+1));
		tmp = prec_prod(tmp,A);
#endif
		/* tmp =  (h^i/(i+1)!)*A^i */
		e_b += tmp;
		/*  e_b = \Sum_{j=0}^{i}\frac{h^j}{(j+1)!}*A^{j} */
	}

#ifdef USE_OCTAVE
	Vector *out = new Vector(h * (e_b * b));
#else 
	Vector *out = new Vector(h * prec_prod(e_b,b));
#endif	
	/* *out = h ( \Sum_{j=0}^{n-1}\frac{h^j}{(j+1)!}*A^{j} ) b */

	return(out);
}

