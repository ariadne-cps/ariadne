/***************************************************************************
 *            vectorfield.cpp
 *
 *  Mon	Aug 16 15:51:45 2004
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
 
#include "vectorfield.h"
#include "linear_algebra.h"
#include "map.h"

using namespace Ariadne;

VectorField::VectorField() {
	this->M=NULL;
	
	/* TODO: implement a test and check if x'=M(x)
	 * is analiticaly solvable. */
	this->solvable=false;
}


VectorField::VectorField(Map *M) {
	this->M=M;
	
	/* TODO: implement a test and check if x'=M(x)
	 * is analiticaly solvable. */
	this->solvable=false;
}

VectorField::~VectorField() {
	delete this->M;

	this->M=NULL;
}
		
Map *LinearVectorField::runge_kutta(const double h,
		double* error, unsigned int n) const {

	/*TOCHECK: the following IS NOT the Runge-Kutta
	 * method of n-th order. Check if the Runge-Kutta 
	 * method is better for some numerical reason.
	 */

	const LinearMap *m=(LinearMap *)(this->M);
	const Matrix *A=m->get_A();
	const Vector *b=m->get_b();
	
	/* Evalutate the approximation (until the n-th term)
	 * of the exponential matrix that solve the system
	 * y'=Ay
	 */
	Matrix *e_Ah=exp_Ah(*A,h,n);
			
	/* Evalutate the approximation (until the n-th term)
	 * of the vector e_b such that exp_Ah(A,h,n) y + e_b is
	 * a numerical solution of the vector field.
	 */
	Vector *e_b=exp_b(*A,*b,h,n);

	LinearMap *out=new LinearMap(e_Ah,e_b);

	/* Estimate the error done */

	*error=(ETA/h)+pow(h,n);

	return(out);
}
	
Map *LinearVectorField::adaptive_runge_kutta(const ASet &S, double *h,
				double approx_value,
				unsigned int n) const {
	/*TODO: reimplement this method */

	double error;
	
	return this->runge_kutta(*h,&error,n);
	
}

LinearVectorField::LinearVectorField(LinearMap *M) {

	this->M=M;

	/* TODO: implement a test and check if x'=Ax+b
	 * is analiticaly solvable. */
	this->solvable=false;
}


LinearVectorField::~LinearVectorField() {
	delete this->M;

	this->M=NULL;
}	
		
Map *LinearVectorField::get_numerical_solution(const ASet &S, 
		const double h, double* error) const {

	/* The LinearVectorClass evaluates a
	 * numerical solution of the differential equation
	 * inducted by the vector field using the 4-th
	 * order Runge-Kutta method. */
	return(this->runge_kutta(h, error, 4));
}
	
Map *LinearVectorField::get_adaptive_numerical_solution(const ASet &S,
		double *h, const double approx_value) const {
	
	/* The LinearVectorClass evaluates an adaptive
	 * numerical solution of the differential equation
	 * inducted by the vector field using the 4-th
	 * order Runge-Kutta adaptive method. */
	return(this->adaptive_runge_kutta(S, h, approx_value, 4));
	
}
											
Map *LinearVectorField::get_solution(const double t_step) const{
	/* TODO: Reimplement this method */
	if (this->solvable) {
		/* TODO: solve some how the system */
	}	
	return NULL;	
}
		
bool LinearVectorField::exactly_solvable() const {
	return this->solvable;
}

