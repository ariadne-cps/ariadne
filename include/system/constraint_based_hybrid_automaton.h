/***************************************************************************
 *            constraint_based_hybrid_automaton.h
 *
 *  Copyright  2007 Pieter Collins
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
 
#ifndef ARIADNE_HYBRID_AUTOMATON_H
#define ARIADNE_HYBRID_AUTOMATON_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "../base/types.h"
#include "../base/pointer.h"
#include "../base/clonable_container.h"
#include "../base/reference_container.h"
#include "../geometry/constraint_interface.h"
#include "../geometry/set_interface.h"
#include "../geometry/hybrid_space.h"
#include "../system/declarations.h"

namespace Ariadne {  

  namespace System {
  
    template<class R> class MapInterface;
    template<class R> class VectorFieldInterface;
    template<class R> class ConstraintBasedDiscreteMode;
    template<class R> class ConstraintBasedDiscreteTransition;
    template<class R> class ConstraintBasedHybridAutomaton;
  
    using Geometry::ConstraintInterface;
    using Geometry::Constraint;
  
    enum EventKind { invariant_tag, guard_tag, activation_tag };

    /*! \brief A discrete mode of a ConstraintBasedHybridAutomaton. */    
    template<class R>
    class ConstraintBasedDiscreteMode {
      friend class ConstraintBasedHybridAutomaton<R>;
      typedef shared_ptr< const Geometry::ConstraintInterface<R> > constraint_const_pointer;
     public:
      id_type id() const;
      dimension_type dimension() const;
      const VectorFieldInterface<R>& dynamic() const;
      const ConstraintBasedDiscreteTransition<R>& transition(id_type event_id) const;
      const reference_set< const ConstraintBasedDiscreteTransition<R> >& transitions() const;
      std::ostream& write(std::ostream& os) const;
     private:
      ConstraintBasedDiscreteMode(id_type id, const VectorFieldInterface<R>& dynamic);
      id_type _id;
      shared_ptr< const VectorFieldInterface<R> > _dynamic;
      std::map< id_type, shared_ptr< const Geometry::ConstraintInterface<R> > > _invariants;
      std::map< id_type, shared_ptr< const Geometry::ConstraintInterface<R> > > _guards;
      std::map< id_type, shared_ptr< const Geometry::ConstraintInterface<R> > > _activations;
      reference_data_map< id_type, const Geometry::ConstraintInterface<R> > _constraints;
      reference_set< const ConstraintBasedDiscreteTransition<R> > _transitions;
    };
      
    template<class R>
    bool operator<(const ConstraintBasedDiscreteMode<R>& dm1, const ConstraintBasedDiscreteMode<R>& dm2);



    /*! \brief A discrete transition of a ConstraintBasedHybridAutomaton. */
    template<class R>
    class ConstraintBasedDiscreteTransition {
      friend class ConstraintBasedHybridAutomaton<R>;
      typedef shared_ptr< const Geometry::ConstraintInterface<R> > constraint_const_pointer;
     public:
      id_type id() const;
      id_type source_id() const;
      id_type destination_id() const;
      const ConstraintBasedDiscreteMode<R>& source() const;   
      const ConstraintBasedDiscreteMode<R>& destination() const;
      const System::MapInterface<R>& reset() const;
      const Geometry::ConstraintInterface<R>& constraint() const;
      EventKind kind() const;
      bool forced() const;
      std::ostream& write(std::ostream& os) const;
    private:
      ConstraintBasedDiscreteTransition(id_type id, const ConstraintBasedDiscreteMode<R>& source, const Geometry::ConstraintInterface<R>& activation);
      ConstraintBasedDiscreteTransition(id_type id, const ConstraintBasedDiscreteMode<R>& source, const ConstraintBasedDiscreteMode<R>& destination,
                         const System::MapInterface<R>& reset, const Geometry::ConstraintInterface<R>& activation, bool forced);
      id_type _event_id;
      const ConstraintBasedDiscreteMode<R>* _source;   
      const ConstraintBasedDiscreteMode<R>* _destination;
      shared_ptr< const System::MapInterface<R> > _reset;  
      shared_ptr< const Geometry::ConstraintInterface<R> > _constraint;
      EventKind _event_kind;
    };
    
    template<class R>
    bool operator<(const ConstraintBasedDiscreteTransition<R>& dt1, const ConstraintBasedDiscreteTransition<R>& dt2);
 

  
    /*! \ingroup HybridTime
     *  \brief A hybrid automaton, comprising continuous-time behaviour
     *  at each discrete mode (ConstraintBasedDiscreteMode), coupled by instantaneous transitions (ConstraintBasedDiscreteTransition).
     *  The state space is given by a Geometry::HybridSet.  
     *
     * A hybrid automaton is a dynamic system with evolution in both
     * continuous time and discrete time. 
     * The state space is a product \f$X=\bigcup\{q\}\times X_q\f$
     * where \f$q\f$ is the <em>discrete state</em> and \f$X_q\f$
     * is the <em>continuous state space</em> of corresponding to
     * each discrete state.
     *
     * For each %ConstraintBasedDiscreteMode, the dynamics is given by a 
     * %VectorFieldInterface describing the continuous dynamics,
     * and a %Set giving an invariant which must be satisified at
     * all times.
     *
     * The discrete time behaviour is specified by %ConstraintBasedDiscreteTransition
     * objects. 
     * Each discrete transition represents an jump from a \a source
     * mode to a \a destination mode. 
     * There can be at most one discrete transition in an automaton
     * with the same event_id and source_id.
     */
    template< class R >
    class ConstraintBasedHybridAutomaton
    {
     public:
      typedef R real_type;
      typedef ConstraintBasedDiscreteMode<R> mode_type;
      typedef ConstraintBasedDiscreteTransition<R> transition_type;

      typedef typename std::set< mode_type >::iterator discrete_mode_iterator;
      typedef typename std::set< mode_type >::const_iterator discrete_mode_const_iterator;
      typedef typename std::set< transition_type >::const_iterator discrete_transition_const_iterator;
      
      typedef shared_ptr< const Geometry::ConstraintInterface<R> > constraint_const_pointer;
      typedef shared_ptr< const System::MapInterface<R> > map_const_pointer;
      typedef shared_ptr< const System::VectorFieldInterface<R> > vector_field_const_pointer;
      typedef const Geometry::ConstraintInterface<R>& constraint_const_reference;
  
     private:
      /*! \brief The name of the hybrid automaton. */
      std::string _name;

      /*! \brief The list of the hybrid automaton's discrete modes. */
      std::set< mode_type > _modes;
      
      /*! \brief The hybrid automaton's transitions. */
      std::set< transition_type > _transitions;

     public:
      //@{
      //! \name Constructors and destructors
      /*! \brief This is a hybrid automaton class constructor.
       *  
       * This constructor initializes the object of the
       * hybrid automaton class.
       * \param name is the name of the hybrid automaton.
       */
      ConstraintBasedHybridAutomaton(const std::string &name);
      
      /*! \brief  This is the destructor of the class hybrid 
       * automaton.
       *
       * This destructor deletes in a safe way an object of the
       * hybrid automaton class.
       */
      ~ConstraintBasedHybridAutomaton();
      //@}
      
      //@{
      //! \name Methods for adding new objects. (Not currently implemented)
      /*! \brief Adds a new dynamic to the automaton.
       * This method does not introduce a new discrete mode.
       */
     const System::VectorFieldInterface<R>&  new_dynamic(id_type id, const System::VectorFieldInterface<R>& dynamic);
      
      /*! \brief Adds a new reset function to the automaton.
       * This method does not introduce a new discrete transition.
       */
      const System::MapInterface<R>& new_reset(id_type id, const System::MapInterface<R>& reset);
      
      /*! \brief Adds a new reset function to the automaton.
       * This method does not introduce a new discrete transition.
       */
      const Geometry::ConstraintInterface<R>& new_constraint(id_type id, const Geometry::ConstraintInterface<R>& constraint);

      //@}
      
      //@{
      //! \name Methods for adding new modes and transitions based on existing objects. (Not currently implemented)
      /*! \brief Adds a new mode to the automaton with the given dynamic. */
      const mode_type& new_mode(id_type mode_id, id_type dynamic_id);
      /*! \brief Adds a new invariant to the automaton with the given event, mode and constraint. */
      const transition_type& new_invariant(id_type event_id, id_type mode_id, id_type constraint_id);
      /*! \brief Adds a new transition to the automaton with the given event, mode, reset and constraint, which may be forced or unforced. */
      const transition_type& new_transition(id_type event_id, id_type source_id, id_type destination_id, id_type reset_id, id_type constraint_id, bool forced);
      /*! \brief Adds a new force transition to the automaton with the given event, mode, reset and constraint. */
      const transition_type& new_forced_transition(id_type event_id, id_type source_id, id_type destination_id, id_type reset_id, id_type constraint_id);
      /*! \brief Adds a new force transition to the automaton with the given event, mode, reset and constraint. */
      const transition_type& new_unforced_transition(id_type event_id, id_type source_id, id_type destination_id, id_type reset_id, id_type constraint_id);
      //@}

      //@{
      //! \name Methods for building the hybrid automaton based on references.
      
      /*! \brief Adds a discrete mode with a given dynamic but no constraints.
       *
       * This method adds a discrete mode to automaton definition.
       * \param id is the unique key or identifier of the discrete mode.
       * \param dynamic is the discrete mode's vector field.
       */
      const mode_type& new_mode(id_type id,
                                const System::VectorFieldInterface<R>& dynamic);
      
      /*! \brief Adds an invariant to a discrete mode.
       *
       * This method adds an invariant constraint to a discrete mode of a hybrid automaton.
       * \param event_id is the unique key or identifier of the discrete mode.
       * \param mode_id is the unique key or identifier of the discrete mode.
       * \param constraint is the constraint of the invariant.
       */
      const transition_type& new_invariant(id_type event_id,
                                           id_type mode_id, 
                                           const Geometry::ConstraintInterface<R>& constraint);
      
      
      /*! \brief Adds a discrete-time reset with a given event, source mode and destination mode.
       *
       * This method creates a new discrete transition from the mode with  mode to the 
       * destination mode.
       * \param event_id is the identifier of the discrete transition's event.
       * \param source_id is the identifier of the discrete transition's source mode.
       * \param destination_id is the identifier of the discrete transition's destination mode.
       * \param reset is the discrete transition's reset.
       * \param constraint is the discrete transition's activating constraint. The constraint is activated when the activation value is positive.
       * \param forced is \a true if the transition is forced to occur whenever the constraint is violated.
       */
      const transition_type& new_transition(id_type event_id,
                                            id_type source_id, 
                                            id_type destination_id,
                                            const System::MapInterface<R>& reset,
                                            const Geometry::ConstraintInterface<R>& constraint,
                                            bool forced);
      
      /*! \brief Adds a discrete-time reset with a given event, source mode and destination mode..
       *
       * This method creates a new discrete transition from the mode with  mode to the 
       * destination mode.
       * \param event_id is the identifier of the discrete transition's event.
       * \param source_id is the identifier of the discrete transition's source mode.
       * \param destination_id is the identifier of the discrete transition's destination mode.
       * \param reset is the discrete transition's reset.
       * \param activation is the discrete transition's activating constraint. The constraint is activated when the activation value is positive.
       */
      const transition_type& new_unforced_transition(id_type event_id,
                                                     id_type source_id, 
                                                     id_type destination_id,
                                                     const System::MapInterface<R>& reset,
                                                     const Geometry::ConstraintInterface<R>& activation);

      /*! \brief Adds a discrete-time reset with a given event, source mode and destination mode..
       *
       * This method creates a new discrete transition from the mode with  mode to the 
       * destination mode.
       * \param event_id is the identifier of the discrete transition's event.
       * \param source_id is the identifier of the discrete transition's source mode.
       * \param destination_id is the identifier of the discrete transition's destination mode.
       * \param reset is the discrete transition's reset.
       * \param guard is the discrete transition's guard constraint. The transition is activated when the guard value is positive.
       */
      const transition_type& new_forced_transition(id_type event_id,
                                                   id_type source_id, 
                                                   id_type destination_id,
                                                   const System::MapInterface<R>& reset,
                                                   const Geometry::ConstraintInterface<R>& guard);
      
      //@}
      
      //@{
      //! \name Methods for building the hybrid automaton based on shared pointers. (Not currently implemented)
      
      /*! \brief Adds a discrete mode with a given dynamic but no constraints.
       *
       * This method adds a discrete mode to automaton definition.
       * \param id is the unique key or identifier of the discrete mode.
       * \param dynamic is the discrete mode's vector field.
       */
      const mode_type& new_mode(id_type id,
                                shared_ptr< const System::VectorFieldInterface<R> > dynamic);
      
      /*! \brief Adds an invariant to a discrete mode.
       *
       * This method adds an invariant constraint to a discrete mode of a hybrid automaton.
       * \param event_id is the unique key or identifier of the discrete mode.
       * \param mode_id is the unique key or identifier of the discrete mode.
       * \param constraint is the constraint of the invariant.
       */
      const transition_type& new_invariant(id_type event_id,
                                           id_type mode_id, 
                                           shared_ptr< const Geometry::ConstraintInterface<R> > constraint);
      
      
      /*! \brief Adds a discrete-time reset with a given event, source mode and destination mode.
       *
       * This method creates a new discrete transition from the mode with  mode to the 
       * destination mode.
       * \param event_id is the identifier of the discrete transition's event.
       * \param source_id is the identifier of the discrete transition's source mode.
       * \param destination_id is the identifier of the discrete transition's destination mode.
       * \param reset is the discrete transition's reset.
       * \param constraint is the discrete transition's activating constraint. The constraint is activated when the activation value is positive.
       * \param forced is \a true if the transition is forced to occur whenever the constraint is violated.
       */
      const transition_type& new_transition(id_type event_id,
                                            id_type source_id, 
                                            id_type destination_id,
                                            shared_ptr< const System::MapInterface<R> > reset,
                                            shared_ptr< const Geometry::ConstraintInterface<R> > constraint,
                                            bool forced);
      
      /*! \brief Adds a discrete-time reset with a given event, source mode and destination mode..
       *
       * This method creates a new discrete transition from the mode with  mode to the 
       * destination mode.
       * \param event_id is the identifier of the discrete transition's event.
       * \param source_id is the identifier of the discrete transition's source mode.
       * \param destination_id is the identifier of the discrete transition's destination mode.
       * \param reset is the discrete transition's reset.
       * \param activation is the discrete transition's activating constraint. The constraint is activated when the activation value is positive.
       */
      const transition_type& new_unforced_transition(id_type event_id,
                                                     id_type source_id, 
                                                     id_type destination_id,
                                                     shared_ptr< const System::MapInterface<R> > reset,
                                                     shared_ptr< const Geometry::ConstraintInterface<R> > activation);

      /*! \brief Adds a discrete-time reset with a given event, source mode and destination mode..
       *
       * This method creates a new discrete transition from the mode with  mode to the 
       * destination mode.
       * \param event_id is the identifier of the discrete transition's event.
       * \param source_id is the identifier of the discrete transition's source mode.
       * \param destination_id is the identifier of the discrete transition's destination mode.
       * \param reset is the discrete transition's reset.
       * \param guard is the discrete transition's guard constraint. The transition is activated when the guard value is positive.
       */
      const transition_type& new_forced_transition(id_type event_id,
                                                   id_type source_id, 
                                                   id_type destination_id,
                                                   shared_ptr< const System::MapInterface<R> > reset,
                                                   shared_ptr< const Geometry::ConstraintInterface<R> > guard);
      
      //@}
      

      //@{
      //! \name Data access
      /*! \brief Test if the hybrid automaton has a discrete mode with key id. */
      bool has_mode(id_type id) const;
      
      /*! \brief Test if the hybrid automaton has a discrete transition or invariant with \a event_id and \a source_id. */
      bool has_transition(id_type event_id, id_type source_id) const;
      
      /*! \brief A set giving the dimension of the state space for each location identifier. */
      Geometry::HybridSpace locations() const;
      
      /*! \brief The set of discrete modes. */
      const std::set< mode_type >& modes() const;
      
      /*! \brief The discrete mode with given id. */
      const mode_type& mode(id_type id) const;
      
      /*! \brief The set of discrete transitions. */
      const std::set< transition_type >& transitions() const;
      /*! \brief The set of discrete transitions with a given \a mode_id. */
      const reference_set<const transition_type>& transitions(id_type mode_id) const;
      /*! \brief The discrete transition with given \a event_id and \a source id. */
      const transition_type& transition(id_type event_id, id_type source_id) const;
      
      /*! \brief Returns the hybrid automaton's name. */
      const std::string& name() const;
      //@}

      //@{
      //! \name Output methods.
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;  
      //@}
     private:
      /* The discrete mode with given id. */
      mode_type& _mode(id_type id);
    };




    template<class R> inline
    bool operator<(const ConstraintBasedDiscreteMode<R>& dm1, const ConstraintBasedDiscreteMode<R>& dm2) {
      return dm1.id() < dm2.id();
    }


    template<class R> inline
    bool operator<(const ConstraintBasedDiscreteTransition<R>& dt1, const ConstraintBasedDiscreteTransition<R>& dt2) {
      return dt1.id() < dt2.id() or (dt1.id()==dt2.id() and dt1.source_id()<dt2.source_id());
    }


    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const ConstraintBasedDiscreteMode<R>& dm) {
      return dm.write(os);
    }

    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const ConstraintBasedDiscreteTransition<R>& dt) {
      return dt.write(os);
    }

    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const ConstraintBasedHybridAutomaton<R>& ha) {
      return ha.write(os);
    }

  }
}

#endif /* ARIADNE_HYBRID_AUTOMATON_H */
