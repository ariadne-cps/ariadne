/***************************************************************************
 *            evolution_parameters.h
 *
 *  Copyright  2007-8  Davide Bresolim, Alberto Casagrande, Pieter Collins
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
 
/*! \file evolution_parameters.h
 *  \brief Parameters for controlling the accuracy of evaluation methods.
 */

#ifndef ARIADNE_EVOLUTION_PARAMETERS_H
#define ARIADNE_EVOLUTION_PARAMETERS_H

#include <cstddef>
#include <boost/smart_ptr.hpp>

#include "grid_set.h"
#include "hybrid_set.h"

namespace Ariadne {

//! \brief Parameters for controlling the accuracy of continuous evolution methods.
class ContinuousEvolutionParameters {
  public:
    typedef uint UnsignedIntType;
    typedef double RealType;

    //! \brief Default constructer gives reasonable values. 
    ContinuousEvolutionParameters();

    //! \brief A suggested order for the representation of enclosure sets.
    UnsignedIntType spacial_order;

    //! \brief A suggested order for the time-stepping method used for numerical integration.
    UnsignedIntType temporal_order;

    //! \brief A suggested minimum step size for integration. 
    //! This value may be ignored if an integration step cannot be performed without reducing the step size below this value. 
    RealType minimum_step_size;

    //! \brief The maximum allowable step size for integration. 
    //! Decreasing this value increases the accuracy of the computation. 
    RealType maximum_step_size;

    //! \brief A suggested minimum cell of a basic set after a subdivision (not a strict bound). 
    Vector<RealType> minimum_enclosure_cell;

    //! \brief The maximum allowable cell of a basic set during integration. 
    //! Decreasing the volume of the cell increases the accuracy of the computation of an over-approximation. 
    Vector<RealType> maximum_enclosure_cell;

	//! \brief The interleaving of events before set model reduction.
    //! \details The value represents the number of events between two set model reductions: a value of zero implies that a set model is always reduced.
	//! Used in conjunction with enable_set_model_reduction. 
	UnsignedIntType set_model_events_size_interleaving;

	//! \brief The maximum overall volume of working sets that is allowed inside a location.
	//! \details Used in conjunction with enable_working_sets_pruning.
	std::map<DiscreteState,Float> hybrid_maximum_working_sets_volume;
    
    //! \brief Enable subdivision of basic sets (false by default).
    bool enable_subdivisions;

    //! \brief Terminate evolution if basic sets became too large (true by default).
    bool enable_premature_termination;

	//! \brief Reduces a set model to the equivalent of its bounding box, every set_model_events_size_interleaving events (false by default).
	bool enable_set_model_reduction;

	//! \brief Enables the pruning of the working sets when too large (false by default).
    //! \details The pruning is done probabilistically, as soon as the volume of the working sets is larger than the hybrid_maximum_working_sets_volume.
	//! <br>
    //! This parameter is used only under lower semantics.
	bool enable_working_sets_pruning;
};


//! \brief Parameters for controlling the accuracy of discretised evolution methods and reachability analysis.
class DiscreteEvolutionParameters {
  public:
    //! \brief The integer type.
    typedef int IntType;
    //! \brief The unsigned integer type.
    typedef uint UnsignedIntType;
    //! \brief The real type.
    typedef double RealType;
  
    //! \brief Default constructer gives reasonable values. 
    DiscreteEvolutionParameters();

    //! \brief The time after which infinite-time upper-evolution routines
    //! may approximate computed sets on a grid. 
    //! \details
    //! This value should be set to the time after which the transient
    //! behaviour has mostly died out. If there are no transients (i.e. the system evolves
    //! for a certain time and then stops), then this parameter should be set to a value 
    //! slightly higher than the maximum running time.
    //! <br> 
    //! Setting this value to the time at which transients die out can improve the
    //! speed and accuracy of the computations. 
    //!  <br> 
    //! This parameter is only used by chain_reach routine.
    RealType transient_time;

    //! \brief The number of discrete steps after which infinite-time upper evolution 
    //! routines may approximate computed sets on the grid. 
    //! \details
    //! Note that the transients are assumed to have died out after <em>either</em> 
    //! transient_time or transient_steps has been reached. 
    //! <br> 
    //! For example, if a hybrid system makes at most three steps before settling down
    //! into its limiting behaviour, then this parameter should be set to three. 
    //! <br> 
    //! Setting this value to the number of steps at which transients die out can improve the
    //! speed and accuracy of the computations. 
    //! <br> 
    //! This parameter is only used by chain_reach() routine.
    //! \sa #transient_time
    UnsignedIntType transient_steps;
    // (See the #transient_time parameter.) 

    //! \brief The time after which an upper evolution or reachability analysis routine
    //! may approximate computed sets on a grid, in order to use previously cached 
    //! integration results for the grid. 
    //! \details
    //! Increasing this parameter improves the accuracy of the computations. 
    //! Setting this parameter too low usually results in meaningless computations. 
    //! As a rule of thumb, a typical system trajectory should move at least four 
    //! times the grid size between locking to the grid. <br>
    //! For forced oscillators, this parameter should be set to the forcing time, 
    //! or a multiple or fraction thereof. 
    //! <br>
    //! This parameter is only used for continuous-time computation.
    RealType lock_to_grid_time;

    //! \brief The time after which an evolver may approximate computed sets on a grid,
    //! in order to use previously cached results for the grid. 
    //! \details
    //! Increasing this parameter may improve the accuracy of the computations.  
    //! If there is recurrence in the system, then this parameter should be set to 
    //! the average recurrence time, if known. 
    //!  <br>
    //! This parameter is only used for discrete-time computation.
    UnsignedIntType lock_to_grid_steps;

    //! \brief Set the depth used for approximation on a grid for the initial set in computations using lower semantics.
    //! \details
    //! Setting this value to \a d will on a grid with lengths \a l will result in the use of initial boxes
    //! with sides of length \f$l/2^{d}\f$.
    //! If the initial set is an open set, then this parameter may be unused; instead the initial sets are points,
    //! spaced according to the initial_grid_density.
    //!  <br> 
    //! Increasing this value increases the accuracy of the computation of lower evolution.
    //!  <br> 
    //! This parameter is only used in the lower_evolve() and lower_reach() routines.
    IntType initial_grid_depth;

    //! \brief Set the density of initial values on a grid for the initial set in computations using lower semantics.
    //! \details
    //! Setting this value to \a g will on a grid with lengths \a l will result in the use of one initial box
    //! (one simulation) for each cell of size \f$l/2^g\f$. If the #initial_grid_depth parameter is higher, then
    //! the initial sets will be smaller than the specified cell. 
    //!  <br> 
    //! Increasing this value increases the number of initial sets used in the computation of lower_evolve() and lower_reach().
    //! and decreases their spacing.
    //!  <br> 
    //! This parameter is only used in the lower_evolve() and lower_reach() routines.
    //! \internal Pieter: I don't like the name of this parameter very much, any other suggestions?
    IntType initial_grid_density;

    //! \brief Set the depth used for approximation on a grid for computations using upper semantics.
    //! \details
    //! Increasing this value increases the accuracy of the computation. 
    //!  <br> 
    //! This parameter is only used in upper_evolve(), upper_reach() and chain_reach() routines.
    IntType maximum_grid_depth;

    //! \brief Set the maximum height used for approximation on a grid for chain reachability computations.
    //! \details
    //! Increasing this value increases domain over which computation is performed. 
    //!  <br> 
    //! This parameter is only used in the chain_reach() routine.
    IntType maximum_grid_height;

    //! \brief Set the allowed bounding domain for chain reachability computations.
	//! \details
    //! Please note that the box is ultimately outer approximated on the grid, therefore the HybridBoxes does not represent a strict bounding on the HybridGridTreeSet resulting from computation.
    //! <br>
	//! This parameters is only used in the chain_reach() routine.
    HybridBoxes bounding_domain;

};

//! \brief Parameters for controlling the accuracy of evolution methods and reachability analysis.
class EvolutionParameters
    : public ContinuousEvolutionParameters, public DiscreteEvolutionParameters 
{ };


inline
ContinuousEvolutionParameters::ContinuousEvolutionParameters() 
    : spacial_order(1),
      temporal_order(4),
      minimum_step_size(0.0),
      maximum_step_size(1.0),
      minimum_enclosure_cell(Vector<RealType>(0)),
      maximum_enclosure_cell(Vector<RealType>(0)),
	  set_model_events_size_interleaving(0),
	  hybrid_maximum_working_sets_volume(),
      enable_subdivisions(false),
      enable_premature_termination(true),
	  enable_set_model_reduction(false),
	  enable_working_sets_pruning(false)
{ }

inline
DiscreteEvolutionParameters::DiscreteEvolutionParameters() 
    : transient_time(0.0),
      transient_steps(0),
      lock_to_grid_time(1.0),
      lock_to_grid_steps(1),
      initial_grid_depth(10),
      initial_grid_density(8),
      maximum_grid_depth(6),
      maximum_grid_height(16)
{ }


inline
std::ostream& 
operator<<(std::ostream& os, const ContinuousEvolutionParameters& p) 
{
    os << "ContinuousEvolutionParameters"
       << "(\n  spacial_order=" << p.spacial_order
       << ",\n  temporal_order=" << p.temporal_order
       << ",\n  minimum_step_size=" << p.minimum_step_size
       << ",\n  maximum_step_size=" << p.maximum_step_size
       << ",\n  minimum_enclosure_cell=" << p.minimum_enclosure_cell
       << ",\n  maximum_enclosure_cell=" << p.maximum_enclosure_cell
       << ",\n  set_model_events_size_interleaving=" << p.set_model_events_size_interleaving
       << ",\n  hybrid_maximum_working_sets_volume=" << p.hybrid_maximum_working_sets_volume
       << ",\n  enable_subdivisions=" << p.enable_subdivisions
       << ",\n  enable_premature_termination=" << p.enable_premature_termination
       << ",\n  enable_set_model_reduction=" << p.enable_set_model_reduction
       << ",\n  enable_working_sets_pruning=" << p.enable_working_sets_pruning
       << "\n)\n";
    return os;
}


inline
std::ostream& 
operator<<(std::ostream& os, const DiscreteEvolutionParameters& p) 
{
    os << "DiscreteEvolutionParameters"
       << "(\n  lock_to_grid_steps=" << p.lock_to_grid_steps
       << ",\n  lock_to_grid_time=" << p.lock_to_grid_time

       << ",\n  transient_time=" << p.transient_time
       << ",\n  transient_steps=" << p.transient_steps


       << ",\n  initial_grid_depth=" << p.initial_grid_depth
       << ",\n  initial_grid_density=" << p.initial_grid_density
       << ",\n  maximum_grid_depth=" << p.maximum_grid_depth
       << ",\n  maximum_grid_height=" << p.maximum_grid_height
       << ",\n  bounding_domain=" << p.bounding_domain

       << "\n)\n";
    return os;
}


inline
std::ostream& 
operator<<(std::ostream& os, const EvolutionParameters& p) 
{
    os << "EvolutionParameters"
       << "(\n  spacial_order=" << p.spacial_order
       << ",\n  temporal_order=" << p.temporal_order
       << ",\n  minimum_step_size=" << p.minimum_step_size
       << ",\n  maximum_step_size=" << p.maximum_step_size
       << ",\n  minimum_enclosure_cell=" << p.minimum_enclosure_cell
       << ",\n  maximum_enclosure_cell=" << p.maximum_enclosure_cell
       << ",\n  set_model_events_size_interleaving=" << p.set_model_events_size_interleaving
       << ",\n  hybrid_maximum_working_sets_volume=" << p.hybrid_maximum_working_sets_volume

       << ",\n\n  lock_to_grid_steps=" << p.lock_to_grid_steps
       << ",\n  lock_to_grid_time=" << p.lock_to_grid_time

       << ",\n  transient_time=" << p.transient_time
       << ",\n  transient_steps=" << p.transient_steps

       << ",\n  initial_grid_depth=" << p.initial_grid_depth
       << ",\n  initial_grid_density=" << p.initial_grid_density
       << ",\n  maximum_grid_depth=" << p.maximum_grid_depth
       << ",\n  maximum_grid_height=" << p.maximum_grid_height
       << ",\n  bounding_domain=" << p.bounding_domain

       << "\n)\n";
    return os;
}

} //!namespace Ariadne

#endif //!ARIADNE_EVOLUTION_PARAMETERS_H
