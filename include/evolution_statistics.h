/***************************************************************************
 *            evolution_statistics.h
 *
 *  Copyright  2010  Luca Geretti
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
 
/*! \file evolution_statistics.h
 *  \brief Statistics on the execution of the evolution methods and the reachability analyses.
 */

#ifndef ARIADNE_EVOLUTION_STATISTICS_H
#define ARIADNE_EVOLUTION_STATISTICS_H

#include <cstddef>
#include <boost/smart_ptr.hpp>

#include "grid_set.h"
#include "hybrid_set.h"

namespace Ariadne {

//! \brief Statistics related to the execution of the evolution methods.
//! \details The statistics are not reset between the execution of evolution methods. This choice 
//! allows to gather statistics between subsequent executions stemming from the analysis methods.
//! <br>
//! The statistics are reset at the beginning of any public HybridReachabilityAnalyser analysis method.
class ContinuousEvolutionStatistics {
  public:
    //! \brief The unsigned integer type.
    typedef uint UnsignedIntType;
    //! \brief The real type.
    typedef double RealType;

    //! \brief Default constructor gives "empty" values.
    ContinuousEvolutionStatistics();

	//! \brief The largest evolution time among the working sets
	RealType largest_evol_time;

	//! \brief The largest evolution steps among the working sets
	UnsignedIntType largest_evol_steps;

	//! \brief The largest enclosure cell among the working sets
	Vector<RealType> largest_enclosure_cell;

	//! \brief The largest total number of working sets
	UnsignedIntType largest_working_sets_total;

	//! \brief Whether the maximum enclosure cell has been reached by the actual enclosure bounding box
	bool has_max_enclosure_been_reached;

	//! \brief Resets the statistics to their initial values
	void reset() {
		largest_evol_time = 0.0;
		largest_evol_steps = 0;
		largest_enclosure_cell = Vector<Float>(0);
		largest_working_sets_total = 0;
		has_max_enclosure_been_reached = false;
	}
};


//! \brief Statistics on the execution of the evolution methods and the reachability analyses.
//! \details The statistics are reset at the beginning of any public HybridReachabilityAnalyser analysis method.
class DiscreteEvolutionStatistics {
  public:
    //! \brief The unsigned integer type.
    typedef uint UnsignedIntType;
    //! \brief The real type.
    typedef double RealType;
  
    //! \brief Default constructor gives "empty" values. 
    DiscreteEvolutionStatistics();

	//! \brief The total number of grid locks.
	UnsignedIntType total_locks;

	//! \brief The total number of recurrent locks (for chain_reach only).
	UnsignedIntType total_recurrent_locks;

	//! \brief Whether a reach region restriction occurred (for chain_reach only).
	bool has_restriction_occurred;	

	//! \brief Resets the statistics to their initial values
	void reset() {
		total_locks = 0;
		total_recurrent_locks = 0;
		has_restriction_occurred = false;
	}
};

//! \brief Statistics for controlling the accuracy of evolution methods and reachability analysis.
class EvolutionStatistics
    : public ContinuousEvolutionStatistics, public DiscreteEvolutionStatistics 
{ };

//! \brief Constructs the default continuous evolution statistics
inline
ContinuousEvolutionStatistics::ContinuousEvolutionStatistics() { 
	reset();
}

//! \brief Constructs the default discrete evolution statistics
inline
DiscreteEvolutionStatistics::DiscreteEvolutionStatistics() { 
	reset();
}

//! \brief Outputs the continuous evolution statistics
inline
std::ostream& 
operator<<(std::ostream& os, const ContinuousEvolutionStatistics& p) 
{
    os << "ContinuousEvolutionStatistics"
       << "(\n  largest_evol_time=" << p.largest_evol_time
       << ",\n  largest_evol_steps=" << p.largest_evol_steps
       << ",\n  largest_enclosure_cell=" << p.largest_enclosure_cell
       << ",\n  largest_working_sets_total=" << p.largest_working_sets_total
       << ",\n  has_max_enclosure_been_reached=" << p.has_max_enclosure_been_reached
       << "\n)\n";
    return os;
}

//! \brief Outputs the discrete evolution statistics
inline
std::ostream& 
operator<<(std::ostream& os, const DiscreteEvolutionStatistics& p) 
{
    os << "DiscreteEvolutionStatistics"
       << "(\n  total_locks=" << p.total_locks
       << ",\n  total_recurrent_locks=" << p.total_recurrent_locks
       << ",\n  has_restriction_occurred=" << p.has_restriction_occurred
       << "\n)\n";
    return os;
}

//! \brief Outputs the evolution statistics
inline
std::ostream& 
operator<<(std::ostream& os, const EvolutionStatistics& p) 
{
    os << "EvolutionStatistics"       
	   << "(\n  largest_evol_time=" << p.largest_evol_time
       << ",\n  largest_evol_steps=" << p.largest_evol_steps
       << ",\n  largest_enclosure_cell=" << p.largest_enclosure_cell
       << ",\n  largest_working_sets_total=" << p.largest_working_sets_total
       << ",\n  has_max_enclosure_been_reached=" << p.has_max_enclosure_been_reached
       << ",\n  total_locks=" << p.total_locks
       << ",\n  total_recurrent_locks=" << p.total_recurrent_locks
       << ",\n  has_restriction_occurred=" << p.has_restriction_occurred
       << "\n)\n";
    return os;
}

} //!namespace Ariadne

#endif //!ARIADNE_EVOLUTION_STATISTICS_H
