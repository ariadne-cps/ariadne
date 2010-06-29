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

//! \brief Statistics on the execution of the evolution methods and the reachability analyses using lower semantics.
//! \details The statistics are reset at the beginning of any public HybridReachabilityAnalyser analysis method using lower semantics.
class LowerDiscreteEvolutionStatistics {
  public:
  
    //! \brief Default constructor gives "empty" values. 
    LowerDiscreteEvolutionStatistics();
	//! \brief Constructs the lower discrete evolution statistics from existing statistics.
	LowerDiscreteEvolutionStatistics(const LowerDiscreteEvolutionStatistics& statistics);

	//! \brief The reached region.
	HybridGridTreeSet reach;

};

//! \brief Statistics on the execution of the evolution methods and the reachability analyses using upper semantics.
//! \details The statistics are reset at the beginning of any public HybridReachabilityAnalyser analysis method using upper semantics.
class UpperDiscreteEvolutionStatistics {
  public:
  
    //! \brief Default constructor gives "empty" values. 
    UpperDiscreteEvolutionStatistics();
	//! \brief Constructs the upper discrete evolution statistics from existing statistics.
	UpperDiscreteEvolutionStatistics(const UpperDiscreteEvolutionStatistics& statistics);

	//! \brief The reached region.
	HybridGridTreeSet reach;

};

//! \brief Statistics on the execution of the evolution methods and the reachability analyses, separated into lower or upper semantics.
//! \details The statistics are reset at the beginning of any public HybridReachabilityAnalyser analysis method, depending on the semantics.
class DiscreteEvolutionStatistics {
  public:
	//! \brief The type of the lower statistics
	typedef LowerDiscreteEvolutionStatistics LowerDiscreteEvolutionStatisticsType;
	//! \brief The type of the upper statistics
	typedef UpperDiscreteEvolutionStatistics UpperDiscreteEvolutionStatisticsType;
	
	//! \brief Default constructor gives "empty" values.
	DiscreteEvolutionStatistics();
	//! \brief Constructor from existing statistics.
	DiscreteEvolutionStatistics(const DiscreteEvolutionStatistics& statistics);

    //! \brief A reference to the statistics for lower semantics.
    LowerDiscreteEvolutionStatisticsType& lower() { return *this->_lower; }
	//! \brief A constant reference to the statistics for lower semantics.
    const LowerDiscreteEvolutionStatisticsType& lower() const { return *this->_lower; }

    //! \brief A reference to the statistics for upper semantics.
    UpperDiscreteEvolutionStatisticsType& upper() { return *this->_upper; }
	//! \brief A constant reference to the statistics for upper semantics.
    const UpperDiscreteEvolutionStatisticsType& upper() const { return *this->_upper; }

  protected:
	boost::shared_ptr< LowerDiscreteEvolutionStatisticsType > _lower;
	boost::shared_ptr< UpperDiscreteEvolutionStatisticsType > _upper;
};

//! \brief Constructs the default lower discrete evolution statistics.
inline
LowerDiscreteEvolutionStatistics::LowerDiscreteEvolutionStatistics() { 
	reach = HybridGridTreeSet();
}

//! \brief Constructs the lower discrete evolution statistics from existing statistics.
inline
LowerDiscreteEvolutionStatistics::LowerDiscreteEvolutionStatistics(const LowerDiscreteEvolutionStatistics& statistics) : 
	reach(statistics.reach)
{ }

//! \brief Constructs the default upper discrete evolution statistics.
inline
UpperDiscreteEvolutionStatistics::UpperDiscreteEvolutionStatistics() { 
	reach = HybridGridTreeSet();
}

//! \brief Constructs the upper discrete evolution statistics from existing statistics.
inline
UpperDiscreteEvolutionStatistics::UpperDiscreteEvolutionStatistics(const UpperDiscreteEvolutionStatistics& statistics) : 
	reach(statistics.reach)
{ }

//! \brief Constructs the default discrete evolution statistics.
inline
DiscreteEvolutionStatistics::DiscreteEvolutionStatistics(): _lower(new LowerDiscreteEvolutionStatisticsType()), 
															_upper(new UpperDiscreteEvolutionStatisticsType())
{ }

//! \brief Constructor from existing statistics.
inline 
DiscreteEvolutionStatistics::DiscreteEvolutionStatistics(const DiscreteEvolutionStatistics& statistics): _lower(new LowerDiscreteEvolutionStatisticsType(statistics.lower())),
																								   		 _upper(new UpperDiscreteEvolutionStatisticsType(statistics.upper()))
{ }

} //!namespace Ariadne

#endif //!ARIADNE_EVOLUTION_STATISTICS_H
