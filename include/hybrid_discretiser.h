/***************************************************************************
 *            hybrid_discretiser.h
 *
 *  Copyright  2006-11  Alberto Casagrande, Pieter Collins
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

/*! \file hybrid_discretiser.h
 *  \brief Methods for computing the evolution of hybrid systems on grids/pavings.
 */

#ifndef ARIADNE_HYBRID_DISCRETISER_H
#define ARIADNE_HYBRID_DISCRETISER_H

#include <boost/smart_ptr.hpp>

#include "numeric.h"
#include "evolver_interface.h"
#include "discretiser_interface.h"
#include "hybrid_automaton_interface.h"
#include "hybrid_grid.h"
#include "vector_field.h"

#include "logging.h"


namespace Ariadne {

class VectorField;
class Grid;
class GridCell;
class GridTreeSet;

class HybridEvolverInterface;
class HybridAutomatonInterface;
class DiscreteLocation;
template<class BS> class HybridBasicSet;

class HybridScalingInterface;
class HybridGrid;
class HybridGridCell;
class HybridGridTreeSet;



/*!  \brief A class for computing the evolution of a discrete-time autonomous system.
 */
template<class HES>
class HybridDiscretiser
    : public DiscretiserInterface<HybridAutomatonInterface,HybridGridCell>
    , public Loggable
{
    typedef typename HES::ContinuousStateSetType ES;
    typedef int AccuracyType;
    typedef HybridTime TimeType;
    typedef HybridAutomatonInterface SystemType;
    typedef HybridGridTreeSet PavingType;
    typedef HybridGridCell CellType;
    typedef HybridGridCell BasicSetType;
    typedef HybridGridTreeSet DenotableSetType;
    typedef ES ContinuousEnclosureType;
    typedef HES EnclosureType;
  private:
    boost::shared_ptr< HybridEvolverInterface > _evolver_ptr;
    boost::shared_ptr< HybridScalingInterface > _scaling_ptr;
  public:
    //@{
    //! \name Constructors and destructors

    //! \brief Construct from evolution parameters and a method for evolving basic sets,
    //!  and a scheme for approximating sets.
    HybridDiscretiser(const HybridEvolverInterface& evolver);

    //! \brief Make a dynamically-allocated copy.
    HybridDiscretiser<HES>* clone() const { return new HybridDiscretiser<HES>(*this); }

    //! \brief Set the length scales for the variables in each locations.
    void set_scaling(const HybridScalingInterface& hsc);
    //@}

    //@{
    //! \name Evaluation on basic sets.

    //! \brief Compute approximations to the reachable and evolved sets
    //! of \a system starting in \a initial_set over \a time. */
    virtual Orbit<BasicSetType>
    evolution(const SystemType& system,
              const BasicSetType& initial_set,
              const TimeType& time,
              const AccuracyType accuracy,
              const Semantics semantics) const;

    //! \brief Compute approximations to the reachable set
    //! of \a system starting in \a initial_set over \a time. */
    virtual DenotableSetType
    reach(const SystemType& system,
                const BasicSetType& initial_set,
                const TimeType& time,
                const AccuracyType accuracy,
                const Semantics semantics) const;

    //! \brief Compute approximations to the evolved set
    //! of \a system starting in \a initial_set over \a time. */
    virtual DenotableSetType
    evolve(const SystemType& system,
                 const BasicSetType& initial_set,
                 const TimeType& time,
                 const AccuracyType accuracy,
                 const Semantics semantics) const;

    /*! \brief Compute a lower-approximation to the the reachable and evolved sets under the system evolution. */
    virtual Orbit<BasicSetType>
    lower_evolution(const SystemType& system,
                    const BasicSetType& initial_set,
                    const TimeType& time,
                    const AccuracyType accuracy) const;

    /*! \brief Compute a lower-approximation to the the reachable and evolved sets under the system evolution. */
    virtual Orbit<BasicSetType>
    upper_evolution(const SystemType& system,
                    const BasicSetType& initial_set,
                    const TimeType& time,
                    const AccuracyType accuracy) const;


  private:
    EnclosureType _enclosure(const BasicSetType& bs) const;
    Orbit<BasicSetType> _discretise(const Orbit<EnclosureType>& orb,
                                    const BasicSetType& initial_set,
                                    const HybridGrid& system_grid,
                                    const AccuracyType accuracy) const;

    DenotableSetType _discretise(const ListSet<EnclosureType>& ls,
                                 const BasicSetType& initial_set,
                                 const HybridGrid& system_grid,
                                 const AccuracyType accuracy) const;


    HybridGrid _hybrid_grid(const SystemType& system) const;
};

} // namespace Ariadne

#endif /* ARIADNE_DISCRETISER_H */
