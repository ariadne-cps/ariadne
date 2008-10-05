#ifndef ARIADNE_EVOLVER_INTERFACE_H
#define ARIADNE_EVOLVER_INTERFACE_H

#include "semantics.h"

template<class ES> class ListSet;

namespace Ariadne {
  
    /*! \ingroup EvaluatorInterfaces \ingroup Evolvers
     *  \brief Interface for evolving a dynamic system.
     *  
     */
    template<class Sys, class ES> class EvolverInterface 
    {
      typedef typename Sys::time_type T;
      typedef ListSet<ES> ESL;
     public:
      /*! \brief Destructor. */
      virtual ~EvolverInterface<Sys,ES>() {};
      /*! \brief Make a dynamically-allocated copy. */
      virtual EvolverInterface<Sys,ES>* clone() const = 0;
     public:
      /*! \brief Compute an approximation to the evolved set under the given semantics. */
      virtual ESL evolve(const Sys& system, const ES& initial_set, const T& time, Semantics semantics) const = 0;
      /*! \brief Compute an approximation to the reachable set under the given semantics. */
      virtual ESL reach(const Sys& system, const ES& initial_set, const T& time, Semantics semantics) const = 0;
      /*! \brief Compute an approximation to the evolved and reachable sets under the given semantics. */
      virtual std::pair<ESL,ESL> reach_evolve(const Sys& system, const ES& initial_set, const T& time, Semantics semantics) const = 0;


      /*! \brief Compute an approximation to the evolved set under the given semantics. */
      virtual void evolution(ESL& final, const Sys& system, const ES& initial, const T& time, Semantics semantics) const = 0;
      /*! \brief Compute an approximation to the evolved and reachable sets under the given semantics. */
      virtual void evolution(ESL& final, ESL& intermediate, const Sys& system, const ES& initial, const T& time, Semantics semantics) const = 0;

      /*! \brief Compute an approximation to the evolved set under the given semantics, starting from a list of enclosure sets. */
      virtual void evolution(ESL& final, const Sys& system, const ESL& initial, const T& time, Semantics semantics) const = 0;
      /*! \brief Compute an approximation to the evolved and reachable sets under the given semantics starting from a list of enclosure sets. */
      virtual void evolution(ESL& final, ESL& intermediate, const Sys& system, const ESL& initial, const T& time, Semantics semantics) const = 0;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const = 0;
    };

    template<class Sys, class ES> inline
    std::ostream& operator<<(std::ostream& os, const EvolverInterface<Sys,ES>& e) {
      return e.write(os); }

} // namespace Ariadne



#endif /* ARIADNE_EVOLVER_INTERFACE_H */
