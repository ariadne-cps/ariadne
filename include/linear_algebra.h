/***************************************************************************
 *            linear_algebra.h
 *
 *  Mon May  3 12:31:15 2004
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
 
#ifndef _LINEAR_ALGEBRA_H
#define _LINEAR_ALGEBRA_H

#ifdef USE_OCTAVE /* if you have liboctave...*/
#include <octave/config.h>

#include <octave/Matrix.h>
	
#else  /* ...else if you have not liboctave 
	  you should have boost's classes */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
	
using namespace boost::numeric::ublas;

#endif


namespace Ariadne {

#ifdef USE_OCTAVE /* if you have liboctave...*/

/*! \brief The vector type. */
typedef ColumnVector Vector;

#else  /* ...else if you have not liboctave 
	  you should have boost's classes */
	
/*! \brief The vector type. */
typedef vector<double> Vector;

/*! \brief The matrix type. */
typedef matrix<double> Matrix; 

/*! \brief The identity matrix type. */
typedef identity_matrix<double> Identity_Matrix; 

#endif

Matrix *exp_Ah(const Matrix &A, const double h, const unsigned int n); 

Vector *exp_b(const Matrix &A, const Vector &b, const double h, 
			const unsigned int n);

}

#endif /* _LINEAR_ALGEBRA_H */
