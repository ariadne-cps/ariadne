/***************************************************************************
 *            constraint_hybrid_automaton.h
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
 
#ifndef ARIADNE_CONSTRAINT_HYBRID_AUTOMATON_H
#define ARIADNE_CONSTRAINT_HYBRID_AUTOMATON_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "../base/stlio.h"
#include "../geometry/constraint.h"
#include "../geometry/hybrid_space.h"
#include "../system/declarations.h"



namespace Ariadne {  

  namespace System {
  
    template<class R> class MapInterface;
    template<class R> class VectorFieldInterface;
    template<class R> class ConstraintHybridAutomaton;
  
    using Geometry::ConstraintInterface;
    using Geometry::Constraint;
  
    template<class R>
    class ConstraintDiscreteMode {
      friend class ConstraintHybridAutomaton<R>;
      typedef typename boost::shared_ptr< const ConstraintInterface<R> > constraint_const_pointer;
     public:
      id_type id() const;
      dimension_type dimension() const;
      const VectorFieldInterface<R>& dynamic() const;
      const ConstraintInterface<R>& invariant(size_type k) const;
      const ConstraintInterface<R>& activation(id_type event_id) const;
      const ConstraintInterface<R>& guard(id_type event_id) const;
      std::ostream& write(std::ostream& os) const;
     private:
      ConstraintDiscreteMode(id_type id, const VectorFieldInterface<R>& dynamic, 
                             const std::vector< boost::shared_ptr< const ConstraintInterface<R> > >& invariant);
      id_type _id;
      boost::shared_ptr< const VectorFieldInterface<R> > _dynamic;
      std::vector< boost::shared_ptr< const ConstraintInterface<R> > > _invariant;
      std::map< id_type, boost::shared_ptr< const ConstraintInterface<R> > > _guards;
      std::map< id_type, boost::shared_ptr< const ConstraintInterface<R> > > _activations;
    };
      
    template<class R>
    bool operator<(const ConstraintDiscreteMode<R>& dm1, const ConstraintDiscreteMode<R>& dm2);


    template<class R>
    class ConstraintDiscreteTransition {
      friend class ConstraintHybridAutomaton<R>;
      typedef typename boost::shared_ptr< const ConstraintInterface<R> > constraint_const_pointer;
     public:
      id_type id() const;
      id_type source_id() const;
      const ConstraintDiscreteMode<R>& source() const;   
      const ConstraintDiscreteMode<R>& destination() const;
      const MapInterface<R>& reset() const;
      const ConstraintInterface<R>& activation() const;
      bool forced() const;
      std::ostream& write(std::ostream& os) const;
    private:
      ConstraintDiscreteTransition(id_type id, const ConstraintDiscreteMode<R>& source, const ConstraintDiscreteMode<R>& destination,
                                   const MapInterface<R>& reset, const ConstraintInterface<R>& activation, bool forced);
      id_type _event_id;
      const ConstraintDiscreteMode<R>* _source;   
      const ConstraintDiscreteMode<R>* _destination;
      boost::shared_ptr< const MapInterface<R> > _reset;  
      boost::shared_ptr< const ConstraintInterface<R> > _activation;
      bool _forced;
    };
    
    template<class R>
    bool operator<(const ConstraintDiscreteTransition<R>& dt1, const ConstraintDiscreteTransition<R>& dt2);


  
    /*! \ingroup HybridTime
     *  \brief A hybrid automaton, comprising continuous-time behaviour
     *  at each DiscreteMode, coupled by instantaneous DiscreteTransition events.
     *  The state space is given by a Geometry::HybridSet.  
     *
     * A hybrid automaton is a dynamic system with evolution in both
     * continuous time and discrete time. 
     * The state space is a product \f$X=\bigcup\{q\}\times X_q\f$
     * where \f$q\f$ is the <em>discrete state</em> and \f$X_q\f$
     * is the <em>continuous state space</em> of corresponding to
     * each discrete state.
     *
     * For each %DiscreteMode, the dynamics is given by a 
     * %VectorFieldInterface describing the continuous dynamics,
     * and a %Set giving an invariant which must be satisified at
     * all times.
     *
     * The discrete time behaviour is specified by %DiscreteTransition
     * objects. 
     * Each discrete transition represents an jump from a \a source
     * mode to a \a destination mode. 
     * There can be at most one discrete transition in an automaton
     * with the same event_id and source_id.
     */
    template< class R >
    class ConstraintHybridAutomaton
    {
      
     public:
      typedef R real_type;
      typedef ConstraintDiscreteMode<R> mode_type;
      typedef ConstraintDiscreteTransition<R> transition_type;

      typedef typename std::set< mode_type >::iterator discrete_mode_iterator;
      typedef typename std::set< mode_type >::const_iterator discrete_mode_const_iterator;
      typedef typename std::set< transition_type >::const_iterator discrete_transition_const_iterator;
      
      typedef typename boost::shared_ptr< const ConstraintInterface<R> > constraint_const_pointer;
      typedef typename boost::shared_ptr< const MapInterface<R> > map_const_pointer;
      typedef typename boost::shared_ptr< const VectorFieldInterface<R> > vector_field_const_pointer;
      typedef const ConstraintInterface<R>& constraint_const_reference;
  
     private:
      /*! \brief The name of the hybrid automaton. */
      std::string _name;

      /*! \brief The list of the hybrid automaton's discrete modes. */
      std::set< mode_type > _modes;
      
      /*! \brief The hybrid automaton's transitions. */
      std::set< transition_type > _transitions;

     public:
      
      /*! \brief This is a hybrid automaton class constructor.
       *  
       * This constructor initializes the object of the
       * hybrid automaton class.
       * \param name is the name of the hybrid automaton.
       */
      ConstraintHybridAutomaton(const std::string &name);
      
      /*! \brief  This is the destructor of the class hybrid 
       * automaton.
       *
       * This destructor deletes in a safe way an object of the
       * hybrid automaton class.
       */
      ~ConstraintHybridAutomaton();
      
      /*! \brief Adds a discrete mode with a given dynamic but no constraints.
       *
       * This method adds a discrete mode to automaton definition.
       * \param id is the unique key or identifier of the discrete mode.
       * \param dynamic is the discrete mode's vector field.
       */
      const mode_type& new_mode(id_type id,
                                const VectorFieldInterface<R>& dynamic);
      
      /*! \brief Adds a discrete mode with a given dynamic but no constraints.
       *
       * This method adds a discrete mode to automaton definition.
       * \param id is the unique key or identifier of the discrete mode.
       * \param dynamic is the discrete mode's vector field.
       * \param invariant is the discrete mode's invariant. 
       */
      const mode_type& new_mode(id_type id,
                                const VectorFieldInterface<R>& dynamic,
                                const ConstraintInterface<R>& invariant);
      
      /*! \brief Adds a discrete mode with a given dynamic but no constraints.
       *
       * This method adds a discrete mode to automaton definition.
       * \param id is the unique key or identifier of the discrete mode.
       * \param dynamic is the discrete mode's vector field.
       * \param invariant is the discrete mode's invariant.
       */
      const ConstraintInterface<R>& new_invariant(id_type mode_id, 
                                                  const ConstraintInterface<R>& constraint);
      
      
      /*! \brief Adds a discrete-time reset with a given event, source mode and destination mode..
       *
       * This method creates a new discrete transition from the mode with  mode to the 
       * destination mode.
       * \param event_id is the identifier of the discrete transition's event.
       * \param source_id is the identifier of the discrete transition's source mode.
       * \param destination_id is the identifier of the discrete transition's destination mode.
       * \param reset is the discrete transition's reset.
       */
      const transition_type& new_transition(id_type event_id,
                                            id_type source_id, 
                                            id_type destination_id,
                                            const MapInterface<R>& reset,
                                            const ConstraintInterface<R>& activation);
      
      /*! \brief Adds a discrete-time reset with a given event, source mode and destination mode..
       *
       * This method creates a new discrete transition from the mode with  mode to the 
       * destination mode.
       * \param event_id is the identifier of the discrete transition's event.
       * \param source_id is the identifier of the discrete transition's source mode.
       * \param destination_id is the identifier of the discrete transition's destination mode.
       * \param reset is the discrete transition's reset.
       */
      const transition_type& new_forced_transition(id_type event_id,
                                                   id_type source_id, 
                                                   id_type destination_id,
                                                   const MapInterface<R>& reset,
                                                   const ConstraintInterface<R>& activation);
      
      /*! \brief Test if the hybrid automaton has a discrete mode with key id. */
      bool has_mode(id_type id) const;
      
      /*! \brief Test if the hybrid automaton has a discrete transition with \a event_id and \a source_id. */
      bool has_transition(id_type event_id, id_type source_id) const;
      
      /*! \brief A set giving the dimension of the state space for each location identifier. */
      Geometry::HybridSpace locations() const;
      
      /*! \brief The set of discrete modes. */
      const std::set< mode_type >& modes() const;
      
      /*! \brief The discrete mode with given id. */
      const mode_type& mode(id_type id) const;
      
      /*! \brief The set of discrete transitions. */
      const std::set< transition_type >& transitions() const;
      
      /*! \brief The discrete transition with given \a event_id and \a source id. */
      const transition_type& transition(id_type event_id, id_type source_id) const;
      
      const System::VectorFieldInterface<R>& reset(id_type mode_id) const;
      const System::MapInterface<R>& reset(id_type event_id, id_type source_id) const;

      const std::vector< constraint_const_pointer >& invariants(id_type mode_id) const;
      const std::map< id_type, constraint_const_pointer >& activations(id_type source_id) const;
      const std::map< id_type, constraint_const_pointer >& guards(id_type source_id) const;
      constraint_const_reference activation(id_type event_id, id_type source_id) const;
      constraint_const_reference guard(id_type event_id, id_type source_id) const;

      
      /*! \brief Returns the hybrid automaton's name. */
      const std::string &name() const;
      
      std::ostream& write(std::ostream& os) const;  
     private:
      /* The discrete mode with given id. */
      mode_type& _mode(id_type id);
    };




    template<class R> inline
    bool operator<(const ConstraintDiscreteMode<R>& dm1, const ConstraintDiscreteMode<R>& dm2) {
      return dm1.id() < dm2.id();
    }


    template<class R> inline
    bool operator<(const ConstraintDiscreteTransition<R>& dt1, const ConstraintDiscreteTransition<R>& dt2) {
      return dt1.id() < dt2.id() or (dt1.id()==dt2.id() and dt1.source_id()==dt2.source_id());
    }


    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const ConstraintDiscreteMode<R>& dm) {
      return dm.write(os);
    }

    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const ConstraintDiscreteTransition<R>& dt) {
      return dt.write(os);
    }

    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const ConstraintHybridAutomaton<R>& ha) {
      return ha.write(os);
    }

  }
}

#endif /* ARIADNE_CONSTRAINT_HYBRID_AUTOMATON_H */
