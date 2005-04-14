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

#include <numerical_type.h>
#include <discrete_location.h>
#include <automaton.h>

#define TIME_PER_UNIT 15

namespace Ariadne {
namespace Evaluation {

template <typename L> class AriadneLocationTrace;
template <typename L> class AriadneTrace;

namespace IO_Operators{
	
template <typename L>
std::ostream& IO_Operators::operator<<(std::ostream &os, 
	const AriadneLocationTrace<L> &t) {
	
	os << "Discrete Location = "<< (t._location).name() << 
		" Time Spent = " << (t._time);
		
	return os;
}

template <typename LT>
std::ostream& IO_Operators::operator<<(std::ostream &os, 
	const  AriadneTrace<LT> &t) {
	
	typename std::list<LT>::iterator i=t.begin();
	size_t f=0;
	
	os << *i;
		
	for (i=t.begin()+1; i!=t.end(); i++) {
		os << std::endl << *i;
	}
	
	return os;
}


}

template <typename LOC>
class AriadneLocationTrace{
	
		typedef LOC DiscreteLocation;
		typedef typename DiscreteLocation::VectorField VectorField;
		typedef typename VectorField::DenotableSet DenotableSet;
		typedef typename DenotableSet::BasicSet BasicSet;
		typedef typename BasicSet::State State;
		typedef typename State::Real Real;
		
	private:
		
		DiscreteLocation &_location;
	
		Real _time;
	
	public:
		
		AriadneLocationTrace(const DiscreteLocation &loc, Real time): 
				_location(loc), _time(time) {}
					
		AriadneLocationTrace(const AriadneLocationTrace<Real> &A): 
				_location(A._location), _time(A._time) {}
		
		inline const AriadneLocationTrace<Real> &operator=(
			const AriadneLocationTrace<Real> &A) { 
				
			this->_time=A._time;
			this->_location=A._location;
				
			return *this;
		}
		
		template <typename L>
		friend std::ostream& IO_Operators::operator<<(std::ostream &os, 
				AriadneTrace<L> &t);
};

template <typename LT>
class AriadneTrace{
	
		typedef LT LocationTrace;
		typedef typename LocationTrace::Real Real;
		typedef typename LocationTrace::DiscreteLocation DiscreteLocation;
	
	private:
		
		std::list<LocationTrace> _trace;
	
	public:
		
		AriadneTrace() {}
			
		inline const std::list<LocationTrace> &trace() const{
			return this->_trace;	
		}
		
		inline const AriadneTrace<Real> &operator=(
			const AriadneTrace<Real> &A) { 
				
			this->_trace=A._trace;
				
			return *this;
		}
		
		inline void add_location_trace(const LocationTrace &A) {
			(this->_trace).push_back(A);
		}
		
		inline void add_location_trace(
					const DiscreteLocation &loc, Real time) {
						
			LocationTrace A(loc,time);
						
			(this->_trace).push_back(A);
		}
		
		inline void del_last_location_trace() {
			(this->_trace).pop_back();
		}
		
		template <typename LocTrace>
		friend std::ostream& IO_Operators::operator<<(std::ostream &os, 
				AriadneTrace<LocTrace> &t);
};

	
template <typename AUTO , typename HDS ,typename TRACE, typename MAINTAIN, 
typename INT , typename EXPORTER>
class AriadneSolver{
	
	public:
		typedef TRACE Trace;
	
		typedef MAINTAIN MaintainDenotableSet;
	
		typedef INT Integrator;
	
		typedef AUTO Automaton;
		typedef typename Automaton::LeavingTrans LeavingTrans;
		typedef typename LeavingTrans::DiscreteLocation DiscreteLocation;
		typedef typename LeavingTrans::ResetMap ResetMap;
		typedef typename DiscreteLocation::VectorField VectorField;
		typedef typename ResetMap::DenotableSet DenotableSet;
		typedef typename DenotableSet::BasicSet BasicSet;
		typedef typename BasicSet::State State;
		typedef typename State::Real Real;
	
		typedef HDS HybridDenotableSet;
		typedef typename HybridDenotableSet::LocationDenotableSet LocationDenotableSet;		
	
		typedef Ariadne::HybridDefinitions::HybridAutomatonID HybridAutomatonID;
		typedef Ariadne::HybridDefinitions::DiscreteLocationID DiscreteLocationID;
		typedef Ariadne::Geometry::ApproxKind ApproxKind;
	
		typedef EXPORTER Exporter;
	
	private:
		bool _used;
	
		HybridAutomatonID _automaton_id;
	
		std::vector<Integrator> _integrators;	
	
		size_t _cache_size;
	
		inline void _empty_cache(Exporter &exporter,
					DenotableSet &maintain_cache,size_t &ds_in_cache) {
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif	
								
			/* if cache is not empty */
			if (ds_in_cache > 0 ) {
				
				#ifdef VERBATIM	
					std::cout << "Flushing cache....";
				#endif				
				
				/* add flow to the reached region of the current location */
				exporter.export_denotableset(maintain_cache);
					
				ds_in_cache=0;
						
				DenotableSet empty_ds;
				maintain_cache=empty_ds;
				
				#ifdef VERBATIM		
					std::cout << "done" <<std::endl;
				#endif
						
			}

			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif			
		
		}
	
		inline void _empty_cache(MaintainDenotableSet &reach,
					DenotableSet &maintain_cache,size_t &ds_in_cache) {
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif	
								
			/* if cache is not empty */
			if (ds_in_cache > 0 ) {
				
				#ifdef VERBATIM	
					std::cout << "Flushing cache....";
				#endif				
				
				/* add flow to the reached region of the current location */
				reach.inplace_union_bbox_intersection(maintain_cache);
					
				ds_in_cache=0;
						
				DenotableSet empty_ds;
				maintain_cache=empty_ds;
				
				#ifdef VERBATIM		
					std::cout << "done" <<std::endl;
				#endif
						
			}

			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif			
		
		}
		
		inline bool _try_to_reset_and_eval_reachability(const Automaton &A, 
				const DiscreteLocation &loc, DenotableSet &flow, 
				Real t_max, Real t, unsigned int n, Real e,
				Ariadne::Geometry::ApproxKind atype,
				const std::vector<LeavingTrans> &leaving_arcs,
				DenotableSet &maintain_cache, size_t &ds_in_cache,
				Exporter &exporter) {
					
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			DenotableSet flow_act;
			bool fix_point=true, thread_fix_point=false;
				
			/* for each leaving arc e_i */
			for (size_t i=0; i < leaving_arcs.size(); i++) {
						
				/* if the flow tube intersect the e_i's activation region */
				if (intersects_interior(flow,leaving_arcs[i].activation())) {
				
					this->_empty_cache(exporter,
								maintain_cache,ds_in_cache);
					
					#ifdef VERBATIM
						std::cout << "Reseting...";
					#endif
					
					/* evaluate the intersection */
					flow_act=closure_of_intersection_of_interior(
							flow,leaving_arcs[i].activation());
					
					/* get e_i's informations */
					const ResetMap &reset=leaving_arcs[i].reset();
					const DiscreteLocationID &dest_id=(leaving_arcs[i].destination()).id();
					
					/* evaluates reachability region on the new location 
					 * from the reset flow  */
					thread_fix_point=this->_reach(A, dest_id, reset(flow_act),
								t_max, t , n-1 , e , atype, exporter, true);
					
					fix_point= fix_point && thread_fix_point;
					
					#ifdef VERBATIM
						std::cout << "done"<<std::endl;
						
						std::cout << "Location "<< loc.name() 
							<< " Flow time left "<< t_max-t <<std::endl;
					#endif

				}
			}
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			return fix_point&&thread_fix_point;
		}
	
		bool _reach(const Automaton &A, 
				const DiscreteLocationID &id, DenotableSet ds0, 
				Real t_max, Real t, unsigned int n, Real e,
				Ariadne::Geometry::ApproxKind atype, 
				Exporter &exporter, const bool from_reset = false) {
					
			Integrator &integrator=this->_integrators[id];
			const DenotableSet &inv=(A._locations[id]).invariant();
			const std::vector<LeavingTrans> &leaving_arcs=A._automaton[id];
			const DiscreteLocation &loc=A._locations[id];
					
			DenotableSet ds1, flow, maintain_cache, empty_ds;
					
					
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
								
			#ifdef VERBATIM	
			uint unit=0;
				
			if (from_reset) {			
				std::cout << "done"<<std::endl;
			}
			
			std::cout << "Location "<< (A._locations[id]).name() 
				<< " Flow time left "<< t_max-t <<std::endl;
			#endif
					
			size_t ds_in_cache=1;
			
			ds0=closure_of_intersection_of_interior(ds0,inv);
			
			if (!ds0.empty()) {
			
				/* not implemented yet
				if (subset_of_interior(ds0,reach_regions[id])) {
					return true;	
				}*/
			
				maintain_cache.inplace_union(ds0);
			}
				
			while (((t_max==0)||(t<t_max))&&(!ds0.empty())) {
				
				#ifdef VERBATIM	
					if (unit>TIME_PER_UNIT) {				
						std::cout << "Time=" << t <<std::endl;
						unit=0;
					}
					
					unit++;
				#endif
	
				/* WARNING!!! this is used to get a faster computation
				 * but if we have not finite math approximation is useless */  
				ds0.set_precision_to_upperapproximating(epsilon(e));
				
				/* evaluates flow and flow-tube */
				flow=integrator.get_flow_tube_from_for_to(ds0,e,ds1,atype);
					
				/* TO DO: we can get more precision considering than the flow can leave 
				 * the invariant. */
				
				/* intersect flow with invariant */
				flow=closure_of_intersection_of_interior(flow,inv);

				maintain_cache.inplace_union(flow);
				ds_in_cache++;
				
				/* if I can reset again */
				if (n>0) {
					
					/* try to reset */
					if (this->_try_to_reset_and_eval_reachability(
							A, loc, flow , t_max, t, n , e, atype,
							leaving_arcs, maintain_cache, ds_in_cache,
							exporter)) {

						return true;
					}
				}
				
				/* intersect the reached reagion with the location's invariant */
				ds0=closure_of_intersection_of_interior(ds1,inv);
				
				/* update the flow time */
				t+=e;
				
				/* if the cache is full */
				if (ds_in_cache >=this->_cache_size) {
					
					/* empty it */
					this->_empty_cache(exporter,
								maintain_cache,ds_in_cache);
					
				}
				
			}
			
			/* if there is something into cache empty it*/
			this->_empty_cache(exporter,maintain_cache,ds_in_cache);
			
			#ifdef VERBATIM
			if (from_reset) {			
				std::cout << "Backtracking...";
			}
			#endif

			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			return false;
		}
		
		inline bool _try_to_reset_and_eval_reachability(const Automaton &A, 
				const DiscreteLocation &loc, DenotableSet &flow, 
				Real t_max, Real t, unsigned int n, Real e,
				Ariadne::Geometry::ApproxKind atype,
				const std::vector<LeavingTrans> &leaving_arcs,
				DenotableSet &maintain_cache, size_t &ds_in_cache,
				std::vector<MaintainDenotableSet> &reach_regions) {
					
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			DenotableSet flow_act;
			bool fix_point=true, thread_fix_point;
					
			/* for each leaving arc e_i */
			for (size_t i=0; i < leaving_arcs.size(); i++) {
						
				/* if the flow tube intersect the e_i's activation region */
				if (intersects_interior(flow,leaving_arcs[i].activation())) {
					
					#ifdef VERBATIM
						std::cout << "Reseting"<<std::endl;
					#endif
					
					this->_empty_cache(reach_regions[loc.id()],
								maintain_cache,ds_in_cache);
					
					/* evaluate the intersection */
					flow_act=closure_of_intersection_of_interior(
							flow,leaving_arcs[i].activation());
							
					/* get e_i's informations */
					const ResetMap &reset=leaving_arcs[i].reset();
					const DiscreteLocationID &dest_id=(leaving_arcs[i].destination()).id();
						
					/* evaluates reachability region on the new location 
					 * from the reset flow  */
					thread_fix_point=this->_reach(A, dest_id, reset(flow_act),
								t_max, t , n-1 , e , atype, reach_regions);
					
					fix_point= fix_point && thread_fix_point;
					
					#ifdef VERBATIM
						std::cout << "Backtracking"<<std::endl;
					#endif
				}
			}
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			return fix_point;
		}
	
		bool _reach(const Automaton &A, 
				const DiscreteLocationID &id, DenotableSet ds0, 
				Real t_max, Real t, unsigned int n, Real e,
				Ariadne::Geometry::ApproxKind atype, 
				std::vector<MaintainDenotableSet> &reach_regions) {
					
			Integrator &integrator=this->_integrators[id];
			const DenotableSet &inv=(A._locations[id]).invariant();
			const std::vector<LeavingTrans> &leaving_arcs=A._automaton[id];
			const DiscreteLocation &loc=A._locations[id];
					
			DenotableSet ds1, flow, maintain_cache, empty_ds;

			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
					
			size_t ds_in_cache=1;
					
			ds0=closure_of_intersection_of_interior(ds0,inv);
					
			if (ds0.empty()) {
				
				#ifdef DEBUG
					std::cout << __FILE__ << ":" << __LINE__ << std::endl;
				#endif
				
				return true;
			}
			
			/* not implemented yet
			if (subset_of_interior(ds0,reach_regions[id])) {
				return true;	
			}*/
			
			maintain_cache.inplace_union(ds0);
			
			while (((t_max==0)||(t<t_max))&&(!ds0.empty())) {
				
				
				/* WARNING!!! this is used to get a faster computation
				 * but if we have not finite math approximation is useless */  
				ds0.set_precision_to_upperapproximating(epsilon(e));
				
				/* evaluates flow and flow-tube */
				flow=integrator.get_flow_tube_from_for_to(ds0,e,ds1,atype);
				
				/* TO DO: we can get more precision considering than the flow can leave 
				 * the invariant. */
				
				/* intersect flow with invariant */
				flow=closure_of_intersection_of_interior(flow,inv);
				
				maintain_cache.inplace_union(flow);
				ds_in_cache++;
				
				/* if I can reset again */
				if (n>0) {
					
					/*try to reset */
					this->_try_to_reset_and_eval_reachability(
								A, loc, flow , t_max, t, n , e, atype,
								leaving_arcs, maintain_cache, ds_in_cache, 
								reach_regions);
					
				}
				
				/* intersect the reached reagion with the location's invariant */
				ds0=closure_of_intersection_of_interior(ds1,inv);
				
				/* update the flow time */
				t+=e;
				
				/* if the cache is full */
				if (ds_in_cache >=this->_cache_size) {
					
					/* empty it */
					this->_empty_cache(reach_regions[id],
								maintain_cache,ds_in_cache);
					
				}
				
			}
			
			/* if there is something into cache empty it*/
			this->_empty_cache(reach_regions[id],
								maintain_cache,ds_in_cache);
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			return false;
			
		}
	
		inline bool _try_to_reset_and_check_reachability(const Automaton &A, 
				const DiscreteLocation &loc, DenotableSet &flow, 
				const HybridDenotableSet &R_d, Real t_max, Real t, Real l_time,
				unsigned int n, Real e, Trace &bt, 
				Ariadne::Geometry::ApproxKind atype,
				const std::vector<LeavingTrans> &leaving_arcs) {
					
			DiscreteLocationID &dest_id;
			DenotableSet flow_act;
			ResetMap &reset;
					
			/* for each leaving arc e_i */
			for (i=0; i < leaving_arcs.size(); i++) {
						
				/* if the flow tube intersect the e_i's activation region */
				if (intersects_interior(flow,leaving_arcs[i].activation())) {
							
					/* evaluate the intersection */
					flow_act=closure_of_intersection_of_interior(
							flow,leaving_arcs[i].activation());
							
					/* get e_i's informations */
					reset=leaving_arcs[i].reset();
					dest_id=(leaving_arcs[i].destination()).id();
						
					/* update the trace */
					bt->add_location_trace(loc, l_time);
						
					/* check reachability on the new location from 
					 * the reset flow  */
					if (this->_is_reachable(A, dest_id, reset(flow_act),
							R_d, t_max, t , n-1 , e , atype)) {
								
						/* if I get reachability return true */
						return true;
					}
				
					/* if I don't get reachability remove last trace entry */
					bt->del_last_location_trace();	
				}
			}						
			
			return false;
		}
	
		bool _is_reachable(const Automaton &A, 
				DiscreteLocationID &id, DenotableSet &ds0, 
				const HybridDenotableSet &R_d, Real t_max, Real t,
				unsigned int n, Real e, Trace &bt, 
				Ariadne::Geometry::ApproxKind atype) {
					
			Integrator &integrator=this->_integrators[id];
			DenotableSet &inv=(A._locations[id]).invariant();
			DenotableSet &dest=R_d[id];
			std::vector<LeavingTrans> &leaving_arcs=A._automaton[id];
			DiscreteLocation &loc=A._locations[id];
					
			DenotableSet ds1, flow;
			Real l_time=0;
			size_t i;
					
			ds0=closure_of_intersection_of_interior(ds0,inv);
					
			if (ds0.empty())
				return false;
			
			if (intersects_interior(ds0, dest)) {
				bs.add_location_trace(loc, t);
				return true;
			}
			
			while ((t_max!=0)||(t<t_max)) {
				
				/* evaluates flow and flow-tube */
				flow=integrator.get_flow_tube_from_for_to(ds0,e,ds1,atype);
				
				/* TO DO: we can get more precision considering than the flow can leave 
				 * the invariant. */

				/* intersect flow with invariant */
				flow=closure_of_intersection_of_interior(flow,inv);
				
				/* if I can reset again */
				if (n>0) {
					
					/*try to reset */
					if (this->_try_to_reset_and_check_reachability(
								A, loc, flow ,R_d, t_max, t, l_time, n,e,atype,
								leaving_arcs)) {
						
						/* if I get reachability return true */
						return true;				
					}
					
				}
				
				/* intersect the reached reagion with the location's invariant */
				ds0=closure_of_intersection_of_interior(ds1,inv);
				
				/* update the flow time */
				l_time+=e;
				t+=e;
					
			}
			
			return false;
			
		}
	
		
		
		inline void _clean_solver() {
			(this->_integrators).clear();
			this->_used=false;
		}
		
		inline void _set_used_by(const Automaton &A) {
			
			this->_clean_solver();
			
			this->_automaton_id=A._id;
			this->_used=true;
				
			for (size_t i=0; i<(A._locations).size(); i++) {
				Integrator i_integrator((A._locations[i]).vector_field());
				
				(this->_integrators).push_back(i_integrator);
			}
		}
		
		inline bool _used_by(const Automaton &A) const {
			
			return (this->_automaton_id==A._id);
		}
		
	public:
		AriadneSolver(const size_t cache = 200): _used(false), _cache_size(cache){}	
	
		~AriadneSolver() {
			this->_clean_solver();
		}
		/*! \brief This method evaluates if the 
		 * region \f$R_d\f$ is reachable from \f$R_s\f$.
		 *
		 * This method evaluates if the region \f$R_d\f$ is reachable from 
		 * the region \f$R_s\f$ by the hybrid automaton \f$A\f$ computing an 
		 * approximation of the regions 
		 * reachable in time \f$t_max\f$ with at most \f$n\f$ resets. The 
		 * type of the approximation is specificated by the 
		 * parameter \a atype.
		 * If \f$R_d\f$ is not reachable in time \f$t_max\f$ with at most
		 * \f$n\f$ resets, it returns 
		 * \a false, otherwise returns \a true.
		 * After the computation, if \a true is returned, 
		 * the parameter \f$bt\f$ contains the path from the inital region 
		 * to \f$R_d\f$. 
		 * \param A is the hybrid automaton.
		 * \param R_s is the region from which the reachability is 
		 * tested.
		 * \param R_d is the region of which the rechability is tested.
		 * \param t_max is the maximum flow time (\f$t_max=0\f$ for 
		 * infinite flow time).
		 * \param n is the maximum number of resets.
		 * \param e is precision.
		 * \param bt is the eventual trace from \f$R_s\f$ to \f$R_d\f$.
		 * \param atype is the type of the approximation.
		 * \return \a true if \f$R_d\f$ is reachable from \f$R_s\f$
		 * in time \f$t_max\f$ with at most \f$n\f$ resets, \a false otherwise.
	   	 */
		inline bool is_reachable(const Automaton &A, 
				const LocationDenotableSet &R_s,
				HybridDenotableSet &R_d, Real t_max, 
				unsigned int n, Real e, Trace &bt, 
				Ariadne::Geometry::ApproxKind atype = Ariadne::Geometry::OVER) {
					
			if 	(atype!= Ariadne::Geometry::OVER) {
				throw std::invalid_argument("Not yet implemented.");	
			}
			
			if ((!this->_used)||(!this->_used_by(A))) {
				this->_set_used_by(A);
			}
			
			if (!R_d.ordered()) {
				R_d.order_using(A._locations);	
			}
			
			DenotableSet &ds0=R_s.set;
			DiscreteLocationID id=(R_s.location).id();
			
			return 	this->_is_reachable(A, id, ds0,R_d,t_max,0,n,e,atype);
					
		}

		/*! \brief This method evaluates if the 
		 * region \f$R_d\f$ is reachable from \f$R_s\f$.
		 *
		 * This method evaluates if the region \f$R_d\f$ is reachable from 
		 * the region \f$R_s\f$ by the hybrid automaton \f$A\f$ computing an 
		 * approximation of the regions 
		 * reachable in time \f$t_max\f$ with at most \f$n\f$ resets. The 
		 * type of the approximation is specificated by the 
		 * parameter \a atype.
		 * If \f$R_d\f$ is not reachable in time \f$t_max\f$ with at most
		 * \f$n\f$ resets, it returns 
		 * \a false, otherwise returns \a true.
		 * After the computation, if \a true is returned, 
		 * the parameter \f$bt\f$ contains the path from the inital region 
		 * to \f$R_d\f$. 
		 * \param A is the hybrid automaton.
		 * \param R_s is the region from which the reachability is 
		 * tested.
		 * \param R_d is the region of which the rechability is tested.
		 * \param t_max is the maximum flow time (\f$t_max=0\f$ for 
		 * infinite flow time).
		 * \param n is the maximum number of resets.
		 * \param e is precision.
		 * \param bt is the eventual trace from \f$R_s\f$ to \f$R_d\f$.
		 * \param atype is the type of the approximation.
		 * \return \a true if \f$R_d\f$ is reachable from \f$R_s\f$
		 * in time \f$t_max\f$ with at most \f$n\f$ resets, \a false otherwise.
	   	 */
		inline bool is_reachable(const Automaton &A, 
				HybridDenotableSet &R_s, 
				HybridDenotableSet &R_d, Real t_max, 
				unsigned int n, Real e, Trace &bt, 
				Ariadne::Geometry::ApproxKind atype = Ariadne::Geometry::OVER) {
					
			if 	(atype!= Ariadne::Geometry::OVER) {
				throw std::invalid_argument("Not yet implemented.");	
			}
			
			if ((!this->_used)||(!this->_used_by(A))) {
				this->_set_used_by(A);
			}

			if (!R_s.ordered()) {
				R_s.order_using(A._locations);	
			}
			
			for (size_t i=0; i< R_s.size(); i++) {
				if (this->is_reachable(A, R_s[i].location, 
					R_s[i].set, R_d, t_max, n, e, bt, atype)) {
				
					return true;
				}
			}				
					
			return false;
		}
		
		/*! \brief This method evaluates if the 
		 * region \a R_d is chain-reachable from \a R_s.
		 *
		 * This method evaluates if the region \a R_d is 
		 * chain-reachable from the region \a R_s by the hybrid 
		 * automaton \a A computing an approximation of the regions 
		 * chain-reachable. The type of the approximation is 
		 * specificated by the parameter \a atype.
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
		 * \return \a TRUE if \a R_d is chain-reachable from \a R_s, 
		 * \a FALSE otherwise.
	   	 */
		bool is_chainreachable(const Automaton &A, 
				const HybridDenotableSet &R_s, 
				HybridDenotableSet &R_d, Real e, Trace &bt);

		/*! Evolves the set for a time \a t with at most \a n discrete 
		 * transistions and error determined by \a atype (inner- or 
		 * outer- approximation).
		 */
		inline HybridDenotableSet evolve_for_time(const Automaton &A, 
				const LocationDenotableSet &R_s,
				HybridDenotableSet &R_d, Real t, 
				unsigned int n, ApproxKind atype);
		
		/*! Evolves the set for a time \a t with at most \a n discrete 
		 * transistions and error determined by \a atype (inner- or 
		 * outer- approximation).
		 */
		inline HybridDenotableSet evolve_for_time(const Automaton &A, 
				const HybridDenotableSet &R, Real t, 
				unsigned int n, ApproxKind atype);				
		
		/*! Evolves a set until immediately after the \a n-th reset, as 
		 * long at this occurs withing time \a t.
		 */
		inline HybridDenotableSet evolve_to_event(const Automaton &A, 
				const HybridDenotableSet &R, unsigned int n, 
				Real t, ApproxKind atype);
		
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
		 * \param A is the hybrid automaton.
		 * \param init is the initial region.
		 * \param t_max is the maximum flow time (\a t_max=0 for 
		 *infinite flow time)
		 * \param n is the maximum number of reset accepted.
		 * \param e is the precision.
		 * \param reach_set is the maintainer for the reach set.
		 * \param atype is the approximation type.
		 * \return \a true if a fixpoint is reached, \a false otherwise.
		 */
		inline bool reach(const Automaton &A, 
				LocationDenotableSet &init, Real t_max, 
				unsigned int n, Real e, Exporter &exporter,
				Ariadne::Geometry::ApproxKind atype = Ariadne::Geometry::OVER) {
					
			if 	(atype!= Ariadne::Geometry::OVER) {
				throw std::invalid_argument("Not yet implemented.");	
			}
			
			if ((!this->_used)||(!this->_used_by(A))) {
				this->_set_used_by(A);
			}
			
			DenotableSet &ds0=init.set;
			DiscreteLocationID id=(init.location).id();
			
			bool fix_point=this->_reach(A, id, ds0,t_max,0,n,e,atype, 
								exporter);
			
			return fix_point;
		}
		
		/*! \brief Compute an approximation to the reachable set.
		 *  
		 * This function computes, for a given precision \a e, an 
		 * approximation to the set \f$\mathrm{Reach}(f,a,t,n)\f$, 
		 * which is the set of points reachable starting at the set 
		 * \f$init\f$, and continuing for at most time \f$t\f$ and 
		 * \f$n\f$ events.
		 * As long as the evolution of \f$f\f$ is computable, then this 
		 * problem can always be solved.
		 * Further, by taking \f$t\f$ and \f$n\f$ large enough, the 
		 * computation yields an outer approximation to 
		 * \f$\mathrm{Reach}(f,a)\f$.
		 * Unfortunately, it is not possible to determine how 
		 * large \f$t\f$ and \f$n\f$ need to be.
		 * \param A is the hybrid automaton.
		 * \param init is the initial region.
		 * \param t_max is the maximum flow time (\a t_max=0 for 
		 *infinite flow time)
		 * \param n is the maximum number of reset accepted.
		 * \param e is the precision.
		 * \param reach_set is the maintainer for the reach set.
		 * \param atype is the approximation type.
		 * \return \a true if a fixpoint is reached, \a false otherwise.
		 */
		inline bool reach(const Automaton &A, 
				HybridDenotableSet &init, Real t_max, 
				unsigned int n, Real e, Exporter &exporter,
				Ariadne::Geometry::ApproxKind atype = Ariadne::Geometry::OVER) {
					
			if 	(atype!= Ariadne::Geometry::OVER) {
				throw std::invalid_argument("Not yet implemented.");	
			}
			
			if ((!this->_used)||(!this->_used_by(A))) {
				this->_set_used_by(A);
			}
			
			if (!init.ordered()) {
				init.order_using(A._locations);	
			}
			bool fix_point=true, local_fix_point;
			
			Real zero=0;
			
			for (size_t i=0; i< init.size(); i++) {
				
				local_fix_point=this->_reach(A, 
					((init[i]).location).id(), 
					(init[i]).set, t_max, zero, n, e, atype, exporter);
				
				fix_point=fix_point && local_fix_point;
			}				
					
			return fix_point;	
		}
		
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
		 * \param A is the hybrid automaton.
		 * \param init is the initial region.
		 * \param t_max is the maximum flow time (\a t_max=0 for 
		 *infinite flow time)
		 * \param n is the maximum number of reset accepted.
		 * \param e is the precision.
		 * \param reach_set is the maintainer for the reach set.
		 * \param atype is the approximation type.
		 * \return \a true if a fixpoint is reached, \a false otherwise.
		 */
		inline bool reach(const Automaton &A, 
				const LocationDenotableSet &init, Real t_max, 
				unsigned int n, Real e, MaintainDenotableSet &reach_set,
				Ariadne::Geometry::ApproxKind atype = Ariadne::Geometry::OVER) {
					
			if 	(atype!= Ariadne::Geometry::OVER) {
				throw std::invalid_argument("Not yet implemented.");	
			}
			
			if ((!this->_used)||(!this->_used_by(A))) {
				this->_set_used_by(A);
			}
			
			std::vector<MaintainDenotableSet> reached_regions;
			
			for (size_t i=0 ; i< A._locations.size() ; i++) {
				reached_regions.push_back(reach_set);
			}
			
			const DenotableSet &ds0=init.set;
			const DiscreteLocationID id=(init.location).id();
			
			bool fix_point=this->_reach(A, id, ds0, t_max,0,n,e,atype, 
								reached_regions);
			
			for (size_t i=0 ; i< A._locations.size() ; i++) {
				reach_set.inplace_union(reached_regions[i]);
			}
			
			return fix_point;
		}
		
		/*! \brief Compute an approximation to the reachable set.
		 *  
		 * This function computes, for a given precision \a e, an 
		 * approximation to the set \f$\mathrm{Reach}(f,a,t,n)\f$, 
		 * which is the set of points reachable starting at the set 
		 * \f$init\f$, and continuing for at most time \f$t\f$ and 
		 * \f$n\f$ events.
		 * As long as the evolution of \f$f\f$ is computable, then this 
		 * problem can always be solved.
		 * Further, by taking \f$t\f$ and \f$n\f$ large enough, the 
		 * computation yields an outer approximation to 
		 * \f$\mathrm{Reach}(f,a)\f$.
		 * Unfortunately, it is not possible to determine how 
		 * large \f$t\f$ and \f$n\f$ need to be.
		 * \param A is the hybrid automaton.
		 * \param init is the initial region.
		 * \param t_max is the maximum flow time (\a t_max=0 for 
		 *infinite flow time)
		 * \param n is the maximum number of reset accepted.
		 * \param e is the precision.
		 * \param reach_set is the maintainer for the reach set.
		 * \param atype is the approximation type.
		 * \return \a true if a fixpoint is reached, \a false otherwise.
		 */
		inline bool reach(const Automaton &A, 
				HybridDenotableSet &init, Real t_max, 
				unsigned int n, Real e, MaintainDenotableSet &reach_set,
				Ariadne::Geometry::ApproxKind atype = Ariadne::Geometry::OVER) {
					
			if 	(atype!= Ariadne::Geometry::OVER) {
				throw std::invalid_argument("Not yet implemented.");	
			}
			
			if ((!this->_used)||(!this->_used_by(A))) {
				this->_set_used_by(A);
			}
			
			if (!init.ordered()) {
				init.order_using(A._locations);	
			}
			bool fix_point=true, local_fix_point;
			
			Real zero=0;
			
			MaintainDenotableSet empty_reach(reach_set);
			std::vector<MaintainDenotableSet> reach;
			
			reach.resize(init.size());
			
			for (size_t i=0; i< init.size(); i++) {
				
				for (size_t j=0; j< init.size(); j++) {
					reach[i]=empty_reach;
				}
				
				local_fix_point=this->_reach(A, 
					((init[i]).location).id(), 
					(init[i]).set, t_max, zero, n, e, atype, reach);
				
				fix_point=fix_point && local_fix_point;
				
				for (size_t j=0; j< init.size(); j++) {
					reach_set.inplace_union(reach[i]);
				}
			}				
					
			return fix_point;	
		}
			
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
		MaintainDenotableSet chainreach(const Automaton &A, 
				const HybridDenotableSet &a, 
				Real e);
			
		/*! Computes a map-based description of the reachability 
		 * relation determined by the function reach.
		 */
/*		Map reachability(const Automaton &A, 
				AriadneRationalType t, 
				int n, AriadneRationalType e) = 0;
			*/
		/*! Computes a map-based description of the chain-reachability 
		 * relation determined by the function chainreach.
		 */
/*		Map *chainreachability(const Automaton &A, 
				AriadneRationalType t, 
				int n, AriadneRationalType e) = 0;
			*/
};

}

}
#endif /* _SOLVER_H */
