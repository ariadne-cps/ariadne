/***************************************************************************
 *            solver.h
 *
 *  Mon May  3 13:07:24 2004
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
 
#ifndef _SOLVER_H
#define _SOLVER_H

#include <list>

#include "classes.h"
	
namespace Ariadne {
	
class Solver{
	public:								
		/*! \brief This method evaluates if the 
		 * region \a R_d is reachable from \a R_s.
		 *
		 * This method evaluates if the region \a R_d is reachable from 
		 * the region \a R_s by the hybrid automaton \a A computing an 
		 * approximation of the regions 
		 * reachable in time \a t_max with at most \a n resets. The 
		 * type of the approximation is specificated by the 
		 * parameter \a atype.
		 * If \a R_d is not reachable in time \a t_max with at most
		 * \a n resets, it returns 
		 * \a FALSE, otherwise returns \a TRUE.
		 * After the computation, if \a TRUE is returned, 
		 * the parameter \a bt contains the path from the inital region 
		 * to \a R_d. 
		 * \param A is the hybrid automaton.
		 * \param R_s is the region from which the reachability is 
		 * tested.
		 * \param R_d is the region of which the rechability is tested.
		 * \param t_max is the maximum flow time (\a t_max=0 for 
		 * infinite flow time).
		 * \param n is the maximum number of resets.
		 * \param e is precision.
		 * \param bt is the eventual trace from \a R_s to \a R_d.
		 * \param atype is the type of the approximation.
		 * \return \a TRUE if \a R_d is reachable from \a R_s 
		 * in time \a t_max with at most \n resets, \a FALSE otherwise.
	   	 */
		virtual bool is_reachable(const Automaton &A, HSet &R_s, 
				HSet &R_d, double t_max, unsigned int n, 
				double e,  Trace &bt, ApproxType atype) = 0;
		
		/*! \brief This method evaluates if the 
		 * region \a R_d is chain-reachable from \a R_s.
		 *
		 * This method evaluates if the region \a R_d is chain-reachable 
		 * from the region \a R_s by the hybrid automaton \a A computing 
		 * an approximation of the regions chain-reachable. The 
		 * type of the approximation is specificated by the 
		 * parameter \a atype.
		 * If \a R_d is not chain-reachable it returns 
		 * \a FALSE, otherwise returns \a TRUE.
		 * After the computation, if \a TRUE is returned, 
		 * the parameter \a bt contains the path from the inital region 
		 * to \a R_d. 
		 * \param A is the hybrid automaton.
		 * \param R_s is the region from which the chain-reachability 
		 * is tested.
		 * \param R_d is the region of which the chain-rechability is 
		 * tested.
		 * \param e is precision.
		 * \param bt is the eventual trace from \a R_s to \a R_d.
		 * \param atype is the type of the approximation.
		 * \return \a TRUE if \a R_d is chain-reachable from \a R_s, 
		 * \a FALSE otherwise.
	   	 */
		virtual bool is_chainreachable(const Automaton &A, HSet &R_s, 
				HSet &R_d, double e,  Trace &bt, 
				ApproxType atype) = 0;

		/*! Evolves the set for a time \a t with at most \a n discrete 
		 * transistions and error determined by \a atype (inner- or 
		 * outer- approximation).
		 */
		virtual HSet *evolve_for_time(const Automaton &A, const HSet& R, 
				double t, unsigned int n, ApproxType atype) = 0; 
		
		/*! Evolves a set until immediately after the \a n-th reset, as 
		 * long at this occurs withing time \a t.
		 */
		virtual HSet *evolve_to_event(const Automaton &A, const HSet& R, 
				unsigned int n, double t, ApproxType atype) = 0;
		
		/*! \brief Compute an approximation to the reachable set.
		 *  
		 * This function computes, for a given precision \a e, an 
		 * approximation to the set \f$\mathrm{Reach}(f,a,t,n)\f$, 
		 * which is the set of points reachable starting at the set 
		 * \a a, and continuing for at most time \a t and \a n events.
		 * As long as the evolution of \a f is computable, then this 
		 * problem can always be solved.
		 * Further, by taking \a t and \a n large enough, the 
		 * computation yields an outer approximation to 
		 * \f$\mathrm{Reach}(f,a)\f$.
		 * Unfortunately, it is not possible to determine how 
		 * large \a t and \a n need to be!
		 *
		 * If the system is non-Zeno, then supplying the time 
		 * argument \a t suffices to bound the computation time.
		 * (Alberto: I'm not sure about this statement. Thi should be 
		 * checked)
		 * \param A is the hybrid automaton.
		 * \param a is the initial region.
		 * \param t_max is the maximum flow time (\a t_max=0 for 
		 *infinite flow time)
		 * \param n is the maximum number of reset accepted.
		 * \param e is the precision.
		 * \return The reachability region of the automaton.
		 */
		virtual RSet *reach(const Automaton &A, const HSet &a, 
				double t_max, unsigned int n, double e) = 0;
			
		/*! \brief Compute an outer approximation to the reachable set.
		 *
		 * This function computes a set which is an outer approximation
		 * to the reachable set \f$\mathrm{Reach}(f,a)\f$ by first 
		 * discretising the system to a given tolerance \a e, and then 
		 * computing the reachable set for the approximate system. As 
		 * \f$ e\rightarrow 0^+ \f$, the set computed converges towards
		 * the chain-reachable set, \f$\mathrm{ChainReach}(f,a)\f$. 
		 * Unfortunately, it is not possible to determine the rate of 
		 * convergence, nor whether 
		 * \f$\mathrm{ChainReach}(f,a)=\mathrm{Reach}(f,a)\f$.
		 * \param A is the hybrid automaton.
		 * \param a is the initial region.
		 * \param e is the precision.
		 * \return The chain-reachability region of the automaton.
		 */
		virtual RSet *chainreach(const Automaton &A, const HSet &a, 
				double e) = 0;
			
		/*! Computes a map-based description of the reachability 
		 * relation determined by the function reach.
		 */
		virtual Map *reachability(const Automaton &A, double t, 
				int n, double e) = 0;
			
		/*! Computes a map-based description of the chain-reachability 
		 * relation determined by the function chainreach.
		 */
		virtual Map *chainreachability(const Automaton &A, double t, 
				int n, double e) = 0;
			
};

}

#endif /* _SOLVER_H */
