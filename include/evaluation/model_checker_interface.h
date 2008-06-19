/***************************************************************************
 *            model_checker_interface.h
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
 
/*! \file model_checker_interface.h
 *  \brief Interface for model checking class.
 */

#ifndef ARIADNE_MODEL_CHECKER_INTERFACE_H
#define ARIADNE_MODEL_CHECKER_INTERFACE_H

namespace Ariadne {
  

    class TimedLogicFormula;

    /*! \brief Interface for checking validity of logical formulae for transition systems.
     *  \ingroup EvaluatorInterfaces \ingroup Analysers
     */
    template<class T, class Aprx>
    class ModelCheckerInterface {
      typedef typename Aprx::PartitionTreeSet PartitionTreeSet;
      typedef TransitionSystemInterface<T,Aprx> TransitionSystem;
     public:
      /*! \brief Virtual destructor. */
      ~ModelChecker() { }
      /*! \brief Make a dynamically-allocated copy. */
      ModelChecker<T,Aprx>* clone() const = 0;
      /*! \brief Attempt to verify that the reachable set of \a system starting in \a initial_set remains in \a safe_set. (Not currently implemented) */
      virtual tribool verify(const TransitionSystem& system, const PartitionTreeSet& initial_set, const TimedLogicFormula& formula) const;
   };



  }
}

#endif /* ARIADNE_MODEL_CHECKER_INTERFACE_H */
