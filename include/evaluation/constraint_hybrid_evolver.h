/***************************************************************************
 *            constraint_hybrid_evolver.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
#ifndef ARIADNE_CONSTRAINT_HYBRID_EVOLVER_H
#define ARIADNE_CONSTRAINT_HYBRID_EVOLVER_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "../base/types.h"
#include "../base/tribool.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"
#include "../evaluation/declarations.h"
#include "../evaluation/integrator.h"

namespace Ariadne {  
  namespace Evaluation {
  
    template<class R> class LohnerIntegrator;

    /*! \brief */
    class HybridTime {
     public:
      HybridTime(Numeric::Rational t) : _time(t), _steps(0) { }
      HybridTime(Numeric::Rational t, Numeric::Integer s) : _time(t), _steps(s) { }
      Numeric::Rational& time() { return this->_time; }
      Numeric::Integer& steps() { return this->_steps; }
      const Numeric::Rational& time() const { return this->_time; }
      const Numeric::Integer& steps() const { return this->_steps; }
      bool operator==(const HybridTime& other) const { return this->_time==other._time && this->_steps==other._steps; }
      bool operator!=(const HybridTime& other) const { return !(*this==other); }
      bool operator<(const HybridTime& other) const { return this->_time<other._time; }
      HybridTime operator+(const HybridTime& other) const { return HybridTime(this->_time+other._time,this->_steps+other._steps); }
      HybridTime operator+(const Numeric::Rational& other) const { return HybridTime(this->_time+other,this->_steps); }
     private:
      Numeric::Rational _time;
      Numeric::Integer _steps;
    };

    typedef HybridTime hybrid_time_type;  


    /*! \ingroup Evolve
     *  \brief A class for computing the evolution of a hybrid system.
     */
    template< class R >
    class ConstraintHybridEvolver
    {
      typedef Numeric::Interval<R> I;
     public:
      typedef typename System::ConstraintDiscreteMode<R> mode_type;
      typedef typename System::ConstraintDiscreteTransition<R> transition_type;
      typedef typename Geometry::Zonotope<Numeric::Interval<R> > continuous_basic_set_type;
      typedef Geometry::HybridBasicSet<continuous_basic_set_type> hybrid_basic_set_type;
      typedef Geometry::HybridTimedBasicSet<continuous_basic_set_type> timed_set_type;
      typedef Geometry::HybridListSet<continuous_basic_set_type> hybrid_list_set_type;
      typedef Geometry::ConstraintInterface<R> constraint_type;
      typedef Geometry::Rectangle<R> bounding_box_type;

      typedef boost::shared_ptr< const Geometry::ConstraintInterface<R> > constraint_const_pointer;
     public:

      /*! \brief The type used for real numbers. */
      typedef R real_type;

      /*! \brief Copy constructor. */
      ConstraintHybridEvolver(const ConstraintHybridEvolver<R>& evolver);

      /*! \brief Construct from an applicator and an integrator. */
      ConstraintHybridEvolver(Applicator<R>& applicator, Integrator<R>& integrator);

      /*! \brief Virtual destructor. */
      virtual ~ConstraintHybridEvolver();


      //@{
      //! \name Parameters controlling the evolution.
      /*! \brief The maximum step size for integration. */
      time_type maximum_step_size() const;
      /*! \brief The maximum step size for integration. */
      real_type maximum_basic_set_radius() const;
      /*! \brief The time before the sets are locked to the grid. */
      time_type lock_to_grid_time() const;

      //@}

      //! \name Evolution steps.
      /*! \brief Compute the possible states reached by an unforced jump within time \f$[-h,h]\f$ remaining within \a source_bounding_box and \a destination_bounding_box.
       */
      std::vector<timed_set_type>
      unforced_jump(const transition_type& transition,
                    const timed_set_type& initial_set,
                    const bounding_box_type& source_bounding_box,
                    const bounding_box_type& destination_bounding_box,
                    const time_type& step_size) const;

      /*! \brief Compute the possible states reached by a forced jump within time \f$[-h,h]\f$ remaining within \a source_bounding_box and \a destination_bounding_box.
       */
      std::vector<timed_set_type>
      forced_jump(const transition_type& transition,
                  const timed_set_type& initial_set,
                  const bounding_box_type& source_bounding_box,
                  const bounding_box_type& destination_bounding_box,
                  const time_type& step_size) const;

      /*! \brief Compute the possible states reached by an unforced jump within time \a h.
       *
       * The generators are given by
       * \f[ D\Phi_2 \circ DF \circ D\Phi_1 G; \quad (D\Phi_i\circ DF \circ \dot{\Phi}_1 - \dot{\Phi}_2) (h/2) \f]
       */
      std::vector<timed_set_type>
      forced_jump(const transition_type& transition,
                  const timed_set_type& initial_set,
                  const bounding_box_type& source_bounding_box,
                  const time_type& step_size) const;

      /*! \brief Compute the possible states reached by an unforced jump within time \a h. */
      std::vector<timed_set_type>
      unforced_jump(const transition_type& transition,
                    const timed_set_type& basic_set,
                    const bounding_box_type& source_bounding_box,
                    const time_type& step_size) const;

      timed_set_type
      discrete_step(const transition_type& transition,
                    const timed_set_type& initial_set) const;

      timed_set_type
      integration_step(const mode_type& mode,
                       const timed_set_type& initial_set,
                       const bounding_box_type& bounding_box,
                       const time_type& step_size) const;

      std::vector<timed_set_type>
      lower_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                           const timed_set_type& initial_set,
                           time_type& step_size) const;

      std::vector<timed_set_type>
      upper_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                           const timed_set_type& initial_set,
                           time_type& step_size) const;

      std::vector<timed_set_type>
      lower_reachability_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                              const timed_set_type& initial_set,
                              time_type& step_size) const;

      std::vector<timed_set_type>
      upper_reachability_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                              const timed_set_type& initial_set,
                              time_type& step_size) const;


      /*! \brief Subdivide the set. */
      std::vector<timed_set_type> subdivide(const timed_set_type& ts) const;
      //@}


      //@{
      //! \name Evolution using abstract sets.
      /*! \brief Make a discrete step of the hybrid automaton, starting from initial set. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> discrete_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                           const Geometry::HybridSet<R>& initial_set) const;
     
      /*! \brief Evolve the hybrid automaton within \a bounding_set starting from the \a initial_set respecting invariants and without using discrete transitions. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> continuous_chainreach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                   const Geometry::HybridSet<R>& initial_set,
                                                   const Geometry::HybridSet<R>& bounding_set) const;
     

      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      std::vector<hybrid_basic_set_type> lower_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                  const Geometry::HybridSet<R>&,
                                                  time_type evolution_time,
                                                  size_type maximum_number_of_events) const;
      
      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      std::vector<hybrid_basic_set_type> upper_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                      const Geometry::HybridSet<R>&,
                                                      time_type evolution_time,
                                                      size_type maximum_number_of_events) const;
      
      /*! \brief Compute a lower approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> lower_reach(const System::ConstraintHybridAutomaton<R>&, 
                                         const Geometry::HybridSet<R>&, 
                                         time_type initial_evolution_time, 
                                         time_type final_time, 
                                         size_type maximum_number_of_events) const;
      
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> upper_reach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                         const Geometry::HybridSet<R>& initial_set, 
                                         time_type initial_evolution_time, 
                                         time_type final_time, 
                                         size_type maximum_number_of_events) const;
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. */
      Geometry::HybridSet<R> chainreach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                        const Geometry::HybridSet<R>& initial_set, 
                                        const Geometry::HybridSet<R>& bounding_set) const;

      /*! \brief Compute the viability kernel of \a map within \a bounding_set. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> viable(const System::ConstraintHybridAutomaton<R>& automaton, 
                                         const Geometry::HybridSet<R>& bounding_set) const;
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      tribool verify(const System::ConstraintHybridAutomaton<R>& automaton, 
                     const Geometry::HybridSet<R>& initial_set, 
                     const Geometry::HybridSet<R>& safe_set) const;
      //@}

     public:
      //@{
      //! \name Evolution using concrete sets.
      /*! \brief Make a discrete step of the hybrid automaton, starting from initial set. */
      Geometry::HybridGridMaskSet<R> discrete_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                       const Geometry::HybridGridMaskSet<R>& initial_set) const;
      /*! \brief Evolve the hybrid automaton within \a bounding_set starting from the \a initial_set respecting invariants and without using discrete transitions. */
      Geometry::HybridGridMaskSet<R> continuous_chainreach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                           const Geometry::HybridGridMaskSet<R>& initial_set,
                                                           const Geometry::HybridGridMaskSet<R>& bounding_set) const;


      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using lower semantics. */
      Geometry::HybridListSet<continuous_basic_set_type> lower_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                                      const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                                      time_type evolution_time,
                                                                      size_type maximum_number_of_events) const;
      
      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using upper semantics. */
      Geometry::HybridListSet<continuous_basic_set_type> upper_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                                      const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                                      time_type evolution_time,
                                                                      size_type maximum_number_of_events) const;
      
      /*! \brief Compute the system evolution up to \a time with up to \a maximum_number_of_events using lower semantics. */
      Geometry::HybridListSet<continuous_basic_set_type> lower_reach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                                      const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                                      time_type evolution_time,
                                                                      size_type maximum_number_of_events) const;

      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics.*/
      Geometry::HybridListSet<continuous_basic_set_type> 
      upper_reach(const System::ConstraintHybridAutomaton<R>& automaton, 
                  const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                  time_type evolution_time, 
                  size_type maximum_number_of_events) const;
      

     
     
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
      *  with up to \a maximum_number_of_events using upper semantics. */
      Geometry::HybridGridCellListSet<R> 
      upper_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                   const Geometry::HybridGridCell<R>& initial_set, 
                   const Geometry::HybridGrid<R>& hybrid_grid, 
                   time_type evolution_time, 
                   size_type maximum_number_of_events) const;
      
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridGridCellListSet<R> 
      upper_reach(const System::ConstraintHybridAutomaton<R>& automaton, 
                  const Geometry::HybridGridCell<R>& initial_set, 
                  const Geometry::HybridGrid<R>& hybrid_grid, 
                  time_type evolution_time, 
                  size_type maximum_number_of_events) const;
      
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. */
      Geometry::HybridGridMaskSet<R> upper_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                  const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                  time_type evolution_time, 
                                                  size_type maximum_number_of_events) const;
     
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridGridMaskSet<R> upper_reach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                 const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                 time_type evolution_time, 
                                                 size_type maximum_number_of_events) const;

     
      /*! \brief Compute an over-approximation to the set of points which remain in \a bounding_set under evolution of \a automaton. using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridGridMaskSet<R> viable(const System::ConstraintHybridAutomaton<R>& automaton, 
                                            const Geometry::HybridGridMaskSet<R>& bounding_set) const;
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. */
      Geometry::HybridGridMaskSet<R> chainreach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                const Geometry::HybridGridMaskSet<R>& bounding_set) const;

      /*! \brief Attempt to verify that the reachable set of \a map starting in \a initial_set remains in \a safe_set. */
      tribool verify(const System::ConstraintHybridAutomaton<R>& automaton, 
                     const Geometry::HybridGridMaskSet<R>& initial_set, 
                     const Geometry::HybridGridMaskSet<R>& safe_set) const;
      //@}

      //@{ 
      //! \name Traces and diagnostics
      /*! \brief A trace of the basic sets covered in the last evolution. */
      const std::vector<timed_set_type>& trace() const;
      //@}

     private:
      // Evolve the hybrid automaton within \a domains starting from the initial_set without using discrete transitions (no checking). */
      Geometry::HybridGridMaskSet<R> _discrete_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridGridMaskSet<R>& initial_set,
                                                    const Geometry::HybridGridMaskSet<R>& domain_set) const;
      // Evolve the hybrid automaton within \a domains starting from the initial_set without using discrete transitions (no checking). */
      Geometry::HybridGridMaskSet<R> _continuous_chainreach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                            const Geometry::HybridGridMaskSet<R>& initial_set,
                                                            const Geometry::HybridGridMaskSet<R>& domain_set) const;

     private:

      // Initialisation and finalisation routines
      typedef std::vector<timed_set_type> working_sets_type; 
      working_sets_type _compute_working_sets(const hybrid_list_set_type& set) const;
      hybrid_list_set_type _compute_list_set(const working_sets_type& working_sets, const Geometry::HybridSpace& locations) const;

     private:
      Applicator<R>* _applicator;
      // FIXME: Allow arbitrary integrators
      IntegratorBase< R, System::VectorFieldInterface<R>, Geometry::Zonotope<I,I> >* _integrator;
      
      mutable std::vector<timed_set_type> _trace;
    };


  }
}

#endif /* ARIADNE_HYBRID_EVOLVER_H */
