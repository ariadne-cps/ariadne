/***************************************************************************
 *            vectorfield.h
 *
 *  Wed Apr 28 12:02:45 2004
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
 
#ifndef _VECTORFIELD_H
#define _VECTORFIELD_H

#include "linear_algebra.h"
#include "classes.h"

/* The ETA constant is used to estimate the error done by 
 * the differential equation integrator if it uses the double 
 * precision. (see http://farside.ph.utexas.edu/teaching/329/lectures/node63.html) 
 */
#define ETA 2.22e-16
	
namespace Ariadne {
	
/*! \class VectorField
 *  \brief A class representing a differential equation or inclusion.
 */
class VectorField 
{
	protected:
		/*! \brief Indicates if is analiticaly solvable or not. */
		bool solvable;

		/*! \brief The map M of the differential system x'= M(x) */
		Map *M;
		
	public:
		/*! \brief This is a \a VectorField class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a VectorField.
		 */
		VectorField(); 

		/*! \brief This is a \a VectorField class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a VectorField.
		 */
		VectorField(Map *map); 
			
		/*! \brief This is the destructor of the class 
		 * \a VectorField.
		 *
		 * This destructor deletes in a safe way an object of the class 
		 * \a VectorField.
		 */
		virtual ~VectorField(); 
		
		/*! \brief  Returns a map corresponding to an approximation
		 * of vector field flow.
		 *
		 * This method evaluates a numerical solution of the flow 
		 * using the integration step \a i_step. 
		 * At the end of the computation the parameter \a error 
		 * will be set to the error of the numerical tecnique.
		 * \param S is the flow's initial phase space state set.
		 * \param i_step is the integration step of the numerical 
		 * method.
		 * \param error is an output parameter and will be set to the 
		 * value of the error.
		 * \return A map corresponding to a numerical evalutation 
		 * of the flow.
 		 */
		virtual Map *get_numerical_solution(const ASet &S, 
				const double i_step, 
				double* error) const = 0;
	
		/*! \brief Returns a map corresponding to an adaptive 
		 * approximation of vector field flow.
		 *
		 * This method evaluates an adaptive numerical solution of the 
		 * vector field using as initial integration step \a i_step. 
		 * At the end of the computation the solution error will be 
		 * less or equal to \a approx_value an \a i_step will be set 
		 * to used integration step.
		 * \param S is the flow's initial phase space state set.
		 * \param i_step is the integration step of the numerical 
		 * method.
		 * \param approx_value is the maximal error done by the method.
		 * \return A map corresponding to a numerical evalutation of 
		 * the flow from the set \a S.
 		 */
		virtual Map *get_adaptive_numerical_solution(const ASet &S,
				double *i_step, 
				const double approx_value) const = 0;
											
		/*! \brief Returns a map corresponding to the vector field 
		 * flow.
		 *
		 * This method returns a trasformation corrisponding to a 
		 * \a t_step-timed flow. If the vector field is not 
		 * exactly solvable, a \a NULL pointer should be returned.
		 * \param t_step is the integration time.
		 * \return A trasformation corrisponding to a 
		 * \a t_step-timed flow.
 		 */
		virtual Map *get_solution(const double t_step) const = 0;
		
		/*! \brief Checks if the vector field is exactly solvable.
		 *
		 * This method checks if the vector field of the vector field 
		 * is exactly solvable.
		 * \return \a TRUE if the vector field is exactly solvable, 
		 * \a FALSE otherwise.
 		 */
		virtual bool exactly_solvable() const = 0;
		
};

/*! \brief The linear vector field class
 *
 * This class describes vector fields of the type x'=Ax+b. 
 */
class LinearVectorField: public VectorField
{
		/*! \brief  Returns a map corresponding to an approximation
		 * of vector field flow.
		 *
		 * This method evaluates the numerical evalutation of
		 * the flow from the state \a s using the {\a n}-th order 
		 * Runge-Kutta method and the integration step \a i_step. 
		 * At the end of the computation the parameter \a error will 
		 * be set to the error of the numerical tecnique.
		 * \param i_step is the integration step of the numerical 
		 * method.
		 * \param error is an output parameter and will be set to the 
		 * value of the error.
		 * \param n is the order of the Runge-Kutta method
		 * \return The numerical evalutation of the flow from \a s.
 		 */
		Map *runge_kutta(const double i_step, double* error, 
				 unsigned int n) const;

		/*! \brief  Returns a map corresponding to an adaptive 
		 * approximation of vector field flow.
		 *
		 * This method evaluates an adamptive numerical evaluation of 
		 * the flow from \a s using the n-th order adamptive 
		 * Runge-Kutta method and as initial integration step 
		 * \a i_step. At the end of the computation the solution 
		 * error will be less or equal to \a approx_value an 
		 * \a i_step will be set to used integration step.
		 * \param S is the flow's initial phase space state set.
		 * \param i_step is the integration step of the numerical
		 * method.
		 * \param approx_value is the maximal error done by the method.
		 * \param n is the order of the Runge-Kutta method
		 * \return The numerical solution of the vector field.
 		 */
		Map *adaptive_runge_kutta(const ASet &S, double *i_step,
				double approx_value, 
				unsigned int n) const;
	public:
		/*! \brief This is a \a LinearVectorField class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a LinearVectorField.
		 */
		LinearVectorField(LinearMap *M); 

		/*! \brief This is the destructor of the class 
		 * \a VectorField.
		 *
		 * This destructor deletes in a safe way an object of the class 
		 * \a VectorField.
		 */
		~LinearVectorField(); 
		
		/*! \brief  Returns a map corresponding to an approximation
		 * of vector field flow.
		 *
		 * This method evaluates a numerical solution of the flow 
		 * using the integration step \a i_step. 
		 * At the end of the computation the parameter \a error 
		 * will be set to the error of the numerical tecnique.
		 * \param S is the flow's initial phase space state set.
		 * \param i_step is the integration step of the numerical 
		 * method.
		 * \param error is an output parameter and will be set to the 
		 * value of the error.
		 * \return A map corresponding to a numerical evalutation 
		 * of the flow.
 		 */
		Map *get_numerical_solution(const ASet &S, const double i_step, 
				double* error) const;
	
		/*! \brief Returns a map corresponding to an adaptive 
		 * approximation of vector field flow.
		 *
		 * This method evaluates an adaptive numerical solution of the 
		 * vector field using as initial integration step \a i_step. 
		 * At the end of the computation the solution error will be 
		 * less or equal to \a approx_value an \a i_step will be set 
		 * to used integration step.
		 * \param S is the flow's initial phase space state set.
		 * \param i_step is the integration step of the numerical 
		 * method.
		 * \param approx_value is the maximal error done by the method.
		 * \return A map corresponding to a numerical evalutation of 
		 * the flow from the set \a S.
 		 */
		Map *get_adaptive_numerical_solution(const ASet &S, 
				double *i_step, 
				const double approx_value) const;
		
		/*! \brief Returns a map corresponding to the vector field 
		 * flow.
		 *
		 * This method returns a trasformation corrisponding to a 
		 * \a t_step-timed flow. If the vector field is not 
		 * exactly solvable, a \a NULL pointer should be returned.
		 * \param t_step is the integration time.
		 * \return A trasformation corrisponding to a 
		 * \a t_step-timed flow.
 		 */
		Map *get_solution(const double t_step) const;
		
		/*! \brief Checks if the vector field is exactly solvable.
		 *
		 * This method checks if the vector field of the vector field 
		 * is exactly solvable.
		 * \return \a TRUE if the vector field is exactly solvable, 
		 * \a FALSE otherwise.
 		 */
		bool exactly_solvable() const;
			
};

}

#endif /* _VECTORFIELD_H */
