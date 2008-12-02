/***************************************************************************
 *            evolver_interface.h
 *
 *  Copyright  2008  Pieter Collins
 * 
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
 
/*! \file evolver_interface.h
 *  \brief Interface for computing a single time step of the evolution of a system.
 */

#ifndef ARIADNE_EVOLVER_INTERFACE_H
#define ARIADNE_EVOLVER_INTERFACE_H


namespace Ariadne {

using std::pair;


template<class ES> class ListSet;
template<class ES> class Orbit;

enum Semantics { LOWER_SEMANTICS, UPPER_SEMANTICS }; 

  
/*! \brief Interface for evolving a dynamic system.
 *
 * \sa \link Ariadne::ToolboxInterface \c CalculusInterface<S,M,F> \endlink
 *   , \link Ariadne::ReachabilityAnalyserInterface \c ReachabilityAnalyserInterface<SYS> \endlink
 */
template<class SYS, class ES> 
class EvolverInterface 
{
  public:
    typedef SYS SystemType;
    typedef ES EnclosureType;
    typedef typename SystemType::TimeType TimeType;
    typedef ListSet<EnclosureType> EnclosureListType;

    //! \brief Virtual destructor. 
    virtual ~EvolverInterface() {};

    //! \brief Cloning operator.
    virtual EvolverInterface<SYS,ES>* clone() const = 0;

    //! \brief Write to an output stream. 
    virtual std::ostream& write(std::ostream& os) const = 0;

  public:
    //! \brief Compute an approximation to the evolved set under the given semantics. 
    virtual 
    Orbit<EnclosureType>
    orbit(const SystemType& system, 
          const EnclosureType& initial_set, 
          const TimeType& time, 
          Semantics semantics) const = 0;

    //! \brief Compute an approximation to the evolved set under the given semantics. 
    virtual 
    EnclosureListType 
    evolve(const SystemType& system, 
           const EnclosureType& initial_set, 
           const TimeType& time, 
           Semantics semantics) const = 0;

    //! \brief Compute an approximation to the reachable set under the given semantics. 
    virtual 
    EnclosureListType 
    reach(const SystemType& system, 
          const EnclosureType& initial_set, 
          const TimeType& time, 
          Semantics semantics) const = 0;

    //! \brief Compute an approximation to the evolved and reachable sets under the given semantics. 
    virtual 
    pair<EnclosureListType,EnclosureListType> 
    reach_evolve(const SystemType& system, 
                 const EnclosureType& initial_set, 
                 const TimeType& time, 
                 Semantics semantics) const = 0;
  
  
    //! \brief Compute an approximation to the evolved set under the given semantics. 
    virtual 
    void 
    evolution(EnclosureListType& final, 
              const SystemType& system, 
              const EnclosureType& initial, 
              const TimeType& time, 
              Semantics semantics) const = 0;

    //! \brief Compute an approximation to the evolved and reachable sets 
    //! under the given semantics. 
    virtual void evolution(EnclosureListType& final, 
                           EnclosureListType& intermediate, 
                           const SystemType& system, 
                           const EnclosureType& initial, 
                           const TimeType& time, 
                           Semantics semantics) const = 0;
  

    //! \brief Compute an approximation to the evolved set under the given semantics, 
    //! starting from a list of enclosure sets. 
    virtual 
    void 
    evolution(EnclosureListType& final, 
              const SystemType& system, 
              const EnclosureListType& initial, 
              const TimeType& time, 
              Semantics semantics) const = 0;

    //! \brief Compute an approximation to the evolved and reachable sets 
    //! under the given semantics starting from a list of enclosure sets. 
    virtual 
    void 
    evolution(EnclosureListType& final, 
              EnclosureListType& intermediate, 
              const SystemType& system, 
              const EnclosureListType& initial, 
              const TimeType& time, 
              Semantics semantics) const = 0;
  

};


template<class SYS, class ES> inline
std::ostream& 
operator<<(std::ostream& os, const EvolverInterface<SYS,ES>& e) {
    return e.write(os); 
}


/*! \brief Interface for evolving a system and discretising the result on a grid.
 */
template<class SYS, class BS>
class DiscreteEvolverInterface
{
  public:
    typedef SYS SystemType;
    typedef typename SYS::TimeType TimeType;
    typedef BS BasicSetType;
 
    /*! \brief Destructor. */
    virtual ~DiscreteEvolverInterface() { };
  
    /*! \brief Make a dynamically-allocated copy. */
    virtual DiscreteEvolverInterface<SystemType,BasicSetType>* clone() const = 0;
  
    /*! \brief Compute a lower-approximation to the the reachable and evolved sets under the system evolution. */
    virtual Orbit<BasicSetType> 
    lower_evolution(const SystemType& f, const BasicSetType& s, const TimeType& t, const int accuracy) const = 0;
  
    /*! \brief Compute a lower-approximation to the the reachable and evolved sets under the system evolution. */
    virtual Orbit<BasicSetType> 
    upper_evolution(const SystemType& f, const BasicSetType& s, const TimeType& t, const int accuracy) const = 0;
};


} // namespace Ariadne



#endif // ARIADNE_EVOLVER_INTERFACE_H
