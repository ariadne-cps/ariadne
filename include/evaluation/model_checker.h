/***************************************************************************
 *            model_checker.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file model_checker.h
 *  \brief Methods for computing the images of sets under maps.
 */

#ifndef ARIADNE_MODEL_CHECKER_H
#define ARIADNE_MODEL_CHECKER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "system/transition_system_interface.h"
#include "evaluation/evolver_interface.h"

#include "evaluation/evolution_parameters.h"

namespace Ariadne {
  

    class TimedLogicFormula;

    /*! \ingroup Analysers
     *  \brief A class for checking validity of logical formulae for transition systems.
     */
    template<class T, class Aprx>
    class ModelChecker {
      typedef typename Aprx::real_type R;
      typedef typename Aprx::PartitionTreeSet PartitionTreeSet;
      typedef TransitionSystemInterface<T,Aprx> TransitionSystem;
     private:
      boost::shared_ptr< EvolutionParameters<R> > _parameters;
     public:

      //@{
      //! \name Constructors and destructors

      /*! \brief Virtual destructor. */
      virtual ~ModelChecker() { }
      /*! \brief Construct from evolution parameters. */
      ModelChecker(const EvolutionParameters<R>& parameters) : _parameters(parameters.clone()) { }
      /*! \brief Make a dynamically-allocated copy. */
      ModelChecker<T,Aprx>* clone() const { return new ModelChecker<T,Aprx>(*this); }
  
      //@}

      public:
      //@{ 
      //! \name Verification of logical formulae
      /*! \brief Attempt to verify that the reachable set of \a system starting in \a initial_set remains in \a safe_set. */
      virtual tribool verify(const TransitionSystem& system, const PartitionTreeSet& initial_set, const PartitionTreeSet& safe_set) const;
      /*! \brief Attempt to verify that the reachable set of \a system starting in \a initial_set remains in \a safe_set. (Not currently implemented) */
      virtual tribool verify(const TransitionSystem& system, const PartitionTreeSet& initial_set, const TimedLogicFormula& formula) const;
      //@}

     private:
      typedef typename Aprx::Paving Pv;
      typedef typename Aprx::BasicSet Bx;
      typedef typename Aprx::BasicSet BS;
      typedef typename Aprx::CoverListSet CLS;
      typedef typename Aprx::PartitionListSet PLS;
      typedef TransitionSystem TSI; 
      typedef PartitionTreeSet PTS;

      const int verbosity() const { return this->_parameters->verbosity(); }
      const EvolutionParameters<R>& parameters() const { return *this->_parameters; }
   };



  
} // namespace Ariadne

#endif /* ARIADNE_MODEL_CHECKER_H */
