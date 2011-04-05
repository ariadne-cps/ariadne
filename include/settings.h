/***************************************************************************
 *            settings.h
 *
 *  Copyright  2007-8  Davide Bresolin, Alberto Casagrande, Pieter Collins, Luca Geretti
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
 
/*! \file settings.h
 *  \brief Settings for controlling the accuracy of evaluation methods.
 */

#ifndef ARIADNE_SETTINGS_H
#define ARIADNE_SETTINGS_H

#include <cstddef>
#include <boost/smart_ptr.hpp>

#include "grid_set.h"
#include "hybrid_set.h"
#include "variables.h"
#include "hybrid_automaton.h"

typedef std::map<RealConstant,int,ConstantComparator<Real> > RealConstantIntMap;

namespace Ariadne {

enum EvolutionDirection { FORWARD, BACKWARD };

//! \brief Settings for controlling the accuracy of continuous evolution methods.
class ContinuousEvolutionSettings {
  public:
    typedef uint UnsignedIntType;
    typedef double RealType;

    //! \brief Default constructer gives reasonable values. 
    ContinuousEvolutionSettings();

	//! \brief The direction of continuous evolution. */
	EvolutionDirection direction;

    //! \brief A suggested order for the representation of enclosure sets.
    UnsignedIntType spacial_order;

    //! \brief A suggested order for the time-stepping method used for numerical integration.
    UnsignedIntType temporal_order;

    //! \brief A suggested minimum step size for integration. 
    //! This value may be ignored if an integration step cannot be performed without reducing the step size below this value. 
    RealType minimum_step_size;

    //! \brief The maximum allowable step size for integration, different for each location.
    //! Decreasing the values increases the accuracy of the computation.
    std::map<DiscreteState,RealType> hybrid_maximum_step_size;

    //! \brief A suggested minimum cell of a basic set after a subdivision (not a strict bound). 
    Vector<RealType> minimum_enclosure_cell;

    //! \brief The maximum allowable cell of a basic set during integration. 
    //! Decreasing the volume of the cell increases the accuracy of the computation of an over-approximation. 
    Vector<RealType> maximum_enclosure_cell;
    
    //! \brief Enable subdivision of basic sets (false by default).
    bool enable_subdivisions;

    //! \brief Terminate evolution if basic sets became too large (true by default).
	//! \details In the case of upper semantics, if true and no subdivisions are present, the set is put into the final sets. In the case of lower semantics, the set is discarded.
    bool enable_premature_termination;
};


//! \brief Settings for controlling the accuracy of discretised evolution methods and reachability analysis.
class DiscreteEvolutionSettings {
  public:
    //! \brief The integer type.
    typedef int IntType;
    //! \brief The unsigned integer type.
    typedef uint UnsignedIntType;
    //! \brief The real type.
    typedef double RealType;
  
    //! \brief Default constructer gives reasonable values. 
    DiscreteEvolutionSettings();

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

    //! \brief Set the highest allowed value of maximum_grid_depth to be used by an analysis method that iteratively refines the grid depth.
    //! \details
    //! Increasing this value (greater or equal to zero) reduces the computation time when very high accuracy is required for verification,
	//! but can increase the computation time when very low accuracy is sufficient instead.
    //!  <br> 
    //! This parameter is only used in the iterative routines such as verify_iterative() and the *_parametric() routines.
	IntType lowest_maximum_grid_depth;

    //! \brief Set the highest allowed value of maximum_grid_depth to be used by an analysis method that iteratively refines the grid depth.
    //! \details
    //! Increasing this value increases the maximum accuracy of the computation for iterative methods. 
    //!  <br> 
    //! This parameter is only used in the verify_iterative() routine.
	IntType highest_maximum_grid_depth;

    //! \brief Set the maximum height used for approximation on a grid for chain reachability computations.
    //! \details
    //! Increasing this value increases domain over which computation is performed. 
    //!  <br> 
    //! This parameter is only used in the chain_reach() routine.
    IntType maximum_grid_height;

    //! \brief Set the allowed bounding domain for chain reachability computations.
    //! <br>
	//! This parameters is only used in the chain_reach() routines.
    HybridBoxes domain_bounds;

    //! \brief The reached region for constraining a chain outer reach.
    HybridGridTreeSet outer_approx_constraint;

    //! \brief The constants that must not be automatically split inside a system.
    //! \details The actual intervals values of the constants are irrelevant.
    RealConstantSet locked_constants;

    //! \brief The target ratio of derivatives width to obtain when splitting constants.
    RealType splitting_constants_target_ratio;

	//! \brief Enable the pruning of the trajectories when too large (false by default).
    //! \details The pruning is done probabilistically.
	//! <br>
    //! This parameter is used only under lower semantics.
	bool enable_lower_pruning;

};

//! \brief Settings for controlling the accuracy of evolution methods and reachability analysis.
class EvolutionSettings
    : public ContinuousEvolutionSettings, public DiscreteEvolutionSettings 
{ };


//! \brief Settings for controlling the accuracy of continuous evolution methods.
class VerificationSettings {
  public:
	typedef uint UnsignedIntType;
	typedef double RealType;

	//! \brief Default constructor gives reasonable values.
	VerificationSettings();

    /*! \brief Whether the analysis results must be plotted. */
	bool plot_results;
	/*! \brief The maximum depth of parameter range splitting.
	 * \details A value of zero means that the parameter space is not splitted at all. */
	uint maximum_parameter_depth;
	/*! \brief Whether we allow to skip further outer reach calculation in safety as soon as proving has failed. */
	bool allow_quick_safety_proving;
	/*! \brief Whether we allow to skip further disproving in safety as soon as a counterexample is found. */
	bool allow_quick_safety_disproving;
	/*! \brief Whether we allow to skip further outer reach calculation in dominance as soon as proving has failed. */
	bool allow_quick_dominance_proving;
	/*! \brief Whether we allow to skip further outer reach calculation in dominance as soon as disproving has failed. */
	bool allow_quick_dominance_disproving;
	/*! \brief Whether to substitute midpoints of parameter boxes when proving.
	 * \details Defaults to false. A value of true would not yield a formal result for the parameter box
	 * but would be useful for quick pre-analysis. */
	bool use_param_midpoints_for_proving;
	/*! \brief Whether to substitute midpoints of parameter boxes when disproving.
	 * \details Defaults to true. Indeed, if we use a value of false and successfully disprove, we gain no additional insight.
	 * Choosing false has the benefit of exploring the whole parameter box, but the drawback of possibly be unable to successfully disprove at all
	 * due to error radii. */
	bool use_param_midpoints_for_disproving;
};

inline
ContinuousEvolutionSettings::ContinuousEvolutionSettings() 
    : direction(FORWARD),
      spacial_order(1),
      temporal_order(4),
      minimum_step_size(0.0),
      minimum_enclosure_cell(Vector<RealType>(0)),
      maximum_enclosure_cell(Vector<RealType>(0)),
      enable_subdivisions(false),
      enable_premature_termination(true)
{ }

inline
DiscreteEvolutionSettings::DiscreteEvolutionSettings() 
    : transient_time(0.0),
      transient_steps(0),
      lock_to_grid_time(1.0),
      lock_to_grid_steps(1),
      initial_grid_depth(10),
      initial_grid_density(8),
      maximum_grid_depth(6),
	  lowest_maximum_grid_depth(0),
	  highest_maximum_grid_depth(9),
      maximum_grid_height(16),
      splitting_constants_target_ratio(0.1),
	  enable_lower_pruning(true)
{ }

inline
VerificationSettings::VerificationSettings() :
		plot_results(false),
    	maximum_parameter_depth(3),
    	allow_quick_safety_proving(true),
    	allow_quick_safety_disproving(true),
    	allow_quick_dominance_proving(true),
    	allow_quick_dominance_disproving(true),
    	use_param_midpoints_for_proving(false),
    	use_param_midpoints_for_disproving(true)
{ }


inline
std::ostream& 
operator<<(std::ostream& os, const ContinuousEvolutionSettings& p) 
{
    os << "ContinuousEvolutionSettings"
       << "(\n  direction=" << p.direction
       << ",\n  spacial_order=" << p.spacial_order
       << ",\n  temporal_order=" << p.temporal_order
       << ",\n  minimum_step_size=" << p.minimum_step_size
       << ",\n  hybrid_maximum_step_size=" << p.hybrid_maximum_step_size
       << ",\n  minimum_enclosure_cell=" << p.minimum_enclosure_cell
       << ",\n  maximum_enclosure_cell=" << p.maximum_enclosure_cell
       << ",\n  enable_subdivisions=" << p.enable_subdivisions
       << ",\n  enable_premature_termination=" << p.enable_premature_termination
       << "\n)\n";
    return os;
}


inline
std::ostream& 
operator<<(std::ostream& os, const DiscreteEvolutionSettings& p) 
{
    os << "DiscreteEvolutionSettings"
       << "(\n  lock_to_grid_steps=" << p.lock_to_grid_steps
       << ",\n  lock_to_grid_time=" << p.lock_to_grid_time

       << ",\n  transient_time=" << p.transient_time
       << ",\n  transient_steps=" << p.transient_steps

       << ",\n  initial_grid_depth=" << p.initial_grid_depth
       << ",\n  initial_grid_density=" << p.initial_grid_density
       << ",\n  maximum_grid_depth=" << p.maximum_grid_depth
       << ",\n  lowest_maximum_grid_depth=" << p.lowest_maximum_grid_depth
       << ",\n  highest_maximum_grid_depth=" << p.highest_maximum_grid_depth
       << ",\n  maximum_grid_height=" << p.maximum_grid_height
       << ",\n  bounding_domain=" << p.domain_bounds
       << ",\n  constraint_reach=" << p.outer_approx_constraint
       << ",\n  splitting_constants_target_ratio=" << p.splitting_constants_target_ratio
       << ",\n  enable_lower_pruning=" << p.enable_lower_pruning
       << "\n)\n";
    return os;
}


inline
std::ostream& 
operator<<(std::ostream& os, const EvolutionSettings& p) 
{
    os << "EvolutionSettings"
       << "(\n  direction=" << p.direction
       << ",\n  spacial_order=" << p.spacial_order
       << ",\n  temporal_order=" << p.temporal_order
       << ",\n  minimum_step_size=" << p.minimum_step_size
       << ",\n  hybrid_maximum_step_size=" << p.hybrid_maximum_step_size
       << ",\n  minimum_enclosure_cell=" << p.minimum_enclosure_cell
       << ",\n  maximum_enclosure_cell=" << p.maximum_enclosure_cell

       << ",\n\n  lock_to_grid_steps=" << p.lock_to_grid_steps
       << ",\n  lock_to_grid_time=" << p.lock_to_grid_time

       << ",\n  transient_time=" << p.transient_time
       << ",\n  transient_steps=" << p.transient_steps

       << ",\n  initial_grid_depth=" << p.initial_grid_depth
       << ",\n  initial_grid_density=" << p.initial_grid_density
       << ",\n  maximum_grid_depth=" << p.maximum_grid_depth
       << ",\n  lowest_maximum_grid_depth=" << p.lowest_maximum_grid_depth
       << ",\n  highest_maximum_grid_depth=" << p.highest_maximum_grid_depth
       << ",\n  maximum_grid_height=" << p.maximum_grid_height
       << ",\n  bounding_domain=" << p.domain_bounds
       << ",\n  constraint_reach=" << p.outer_approx_constraint
       << ",\n  splitting_constants_target_ratio=" << p.splitting_constants_target_ratio
       << ",\n  enable_lower_pruning=" << p.enable_lower_pruning
       << "\n)\n";
    return os;
}


inline
std::ostream&
operator<<(std::ostream& os, const VerificationSettings& p)
{
    os << "VerificationSettings"
       << "(\n  plot_results=" << p.plot_results
       << ",\n  maximum_parameter_depth=" << p.maximum_parameter_depth
       << ",\n  allow_quick_safety_proving=" << p.allow_quick_safety_proving
       << ",\n  allow_quick_safety_disproving=" << p.allow_quick_safety_disproving
       << ",\n  allow_quick_dominance_proving=" << p.allow_quick_dominance_proving
       << ",\n  allow_quick_dominance_disproving=" << p.allow_quick_dominance_disproving
       << ",\n  use_param_midpoints_for_proving=" << p.use_param_midpoints_for_proving
       << ",\n  use_param_midpoints_for_disproving=" << p.use_param_midpoints_for_disproving
       << "\n)\n";
    return os;
}

} //!namespace Ariadne

#endif //!ARIADNE_SETTINGS_H
