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

//! \brief Statistics related to the execution of the evolution methods, common to lower or upper semantics.
//! \details The statistics are not reset between the executions of evolution methods. This choice 
//! allows to gather statistics between subsequent executions stemming from the analysis methods.
//! <br>
//! The statistics are reset at the beginning of any public HybridReachabilityAnalyser analysis method.
class CommonContinuousEvolutionStatistics {
  public:
    //! \brief The real type.
    typedef double RealType;

    //! \brief Default constructor gives "empty" values.
    CommonContinuousEvolutionStatistics();
	//! \brief Constructs the common continuous evolution statistics from existing statistics.
	CommonContinuousEvolutionStatistics(const CommonContinuousEvolutionStatistics& statistics);

	//! \brief The largest enclosure cell among the working sets.
	Vector<RealType> largest_enclosure_cell;

	//! \brief Whether the maximum enclosure cell has been reached by the actual enclosure bounding box.
	bool has_max_enclosure_been_reached;

	//! \brief Resets the statistics to their initial values.
	void reset() {
		largest_enclosure_cell = Vector<RealType>(0);
		has_max_enclosure_been_reached = false;
	}
};

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

	//! \brief Resets the statistics to their initial values,
	void reset() {
		reach = HybridGridTreeSet();
	}
};

//! \brief Statistics on the execution of the evolution methods and the reachability analyses using upper semantics.
//! \details The statistics are reset at the beginning of any public HybridReachabilityAnalyser analysis method using upper semantics.
class UpperDiscreteEvolutionStatistics {
  public:
  
    //! \brief Default constructor gives "empty" values. 
    UpperDiscreteEvolutionStatistics();
	//! \brief Constructs the upper discrete evolution statistics from existing statistics.
	UpperDiscreteEvolutionStatistics(const UpperDiscreteEvolutionStatistics& statistics);

	//! \brief Whether a reach region restriction occurred (for any chain_reach only).
	bool has_restriction_occurred;

	//! \brief The reached region.
	HybridGridTreeSet reach;

	//! \brief Resets the statistics to their initial values,
	void reset() {
		has_restriction_occurred = false;
		reach = HybridGridTreeSet();
	}
};

//! \brief Statistics related to the execution of the evolution methods, separated into lower or upper semantics.
//! \details The statistics are not reset between the executions of evolution methods. This choice 
//! allows to gather statistics between subsequent executions stemming from the analysis methods.
//! <br>
//! The statistics are selectively reset at the beginning of any public HybridReachabilityAnalyser analysis method, depending on the semantics.
class ContinuousEvolutionStatistics {
  public:	
	//! \brief The type of the lower statistics
	typedef CommonContinuousEvolutionStatistics LowerContinuousEvolutionStatisticsType;
	//! \brief The type of the upper statistics
	typedef CommonContinuousEvolutionStatistics UpperContinuousEvolutionStatisticsType;

	//! \brief Default constructor gives "empty" values.
	ContinuousEvolutionStatistics();
	//! \brief Constructor from existing statistics.
	ContinuousEvolutionStatistics(const ContinuousEvolutionStatistics& statistics);

    //! \brief A reference to the statistics for lower semantics.
    LowerContinuousEvolutionStatisticsType& lower() { return *this->_lower; }
	//! \brief A constant reference to the statistics for lower semantics.
    const LowerContinuousEvolutionStatisticsType& lower() const { return *this->_lower; }

    //! \brief A reference to the statistics for upper semantics.
    UpperContinuousEvolutionStatisticsType& upper() { return *this->_upper; }
	//! \brief A constant reference to the statistics for upper semantics.
    const UpperContinuousEvolutionStatisticsType& upper() const { return *this->_upper; }

  protected:
	boost::shared_ptr< LowerContinuousEvolutionStatisticsType > _lower;
	boost::shared_ptr< UpperContinuousEvolutionStatisticsType > _upper;
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


//! \brief Aggregate class for continuous and discrete statistics.
//! \details Should be used exclusively to store the whole statistics against overwriting.
class EvolutionStatistics
{ 
  public:	
	//! \brief Default constructor with "empty" statistics.
	EvolutionStatistics();
	//! \brief Default constructor from existing statistics.
	EvolutionStatistics(const ContinuousEvolutionStatistics& continuous, const DiscreteEvolutionStatistics& discrete);

	//! \brief A constant reference to the continuous statistics.
    const ContinuousEvolutionStatistics& continuous() const { return *this->_continuous; }
	//! \brief A constant reference to the discrete statistics.
    const DiscreteEvolutionStatistics& discrete() const { return *this->_discrete; }

  protected:
	boost::shared_ptr< ContinuousEvolutionStatistics > _continuous;
	boost::shared_ptr< DiscreteEvolutionStatistics > _discrete;
};

//! \brief Constructs the default common continuous evolution statistics.
inline
CommonContinuousEvolutionStatistics::CommonContinuousEvolutionStatistics() { 
	reset();
}

//! \brief Constructs the common continuous evolution statistics from existing statistics.
inline
CommonContinuousEvolutionStatistics::CommonContinuousEvolutionStatistics(const CommonContinuousEvolutionStatistics& statistics) : 
	largest_enclosure_cell(statistics.largest_enclosure_cell),
	has_max_enclosure_been_reached(statistics.has_max_enclosure_been_reached)
{ }

//! \brief Constructs the default lower discrete evolution statistics.
inline
LowerDiscreteEvolutionStatistics::LowerDiscreteEvolutionStatistics() { 
	reset();
}

//! \brief Constructs the lower discrete evolution statistics from existing statistics.
inline
LowerDiscreteEvolutionStatistics::LowerDiscreteEvolutionStatistics(const LowerDiscreteEvolutionStatistics& statistics) : 
	reach(statistics.reach)
{ }

//! \brief Constructs the default upper discrete evolution statistics.
inline
UpperDiscreteEvolutionStatistics::UpperDiscreteEvolutionStatistics() { 
	reset();
}

//! \brief Constructs the upper discrete evolution statistics from existing statistics.
inline
UpperDiscreteEvolutionStatistics::UpperDiscreteEvolutionStatistics(const UpperDiscreteEvolutionStatistics& statistics) : 
	has_restriction_occurred(statistics.has_restriction_occurred),
	reach(statistics.reach)
{ }

//! \brief Constructs the default continuous evolution statistics.
inline
ContinuousEvolutionStatistics::ContinuousEvolutionStatistics(): _lower(new LowerContinuousEvolutionStatisticsType()), 
															    _upper(new UpperContinuousEvolutionStatisticsType())
{ }

//! \brief Constructor from existing statistics.
inline 
ContinuousEvolutionStatistics::ContinuousEvolutionStatistics(const ContinuousEvolutionStatistics& statistics): _lower(new LowerContinuousEvolutionStatisticsType(statistics.lower())),
 																								   		 	   _upper(new UpperContinuousEvolutionStatisticsType(statistics.upper()))
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

//! \brief Default constructor with "empty" statistics.
inline
EvolutionStatistics::EvolutionStatistics(): _continuous(new ContinuousEvolutionStatistics()),
											_discrete(new DiscreteEvolutionStatistics())
{ }

//! \brief Default constructor from existing statistics.
inline
EvolutionStatistics::EvolutionStatistics(const ContinuousEvolutionStatistics& continuous, 
										 const DiscreteEvolutionStatistics& discrete): _continuous(new ContinuousEvolutionStatistics(continuous)), 
																					   _discrete(new DiscreteEvolutionStatistics(discrete))
{ }


//! \brief Outputs the continuous evolution statistics.
inline
std::ostream& 
operator<<(std::ostream& os, const CommonContinuousEvolutionStatistics& p) 
{
    os << "CommonContinuousEvolutionStatistics"
       << ",\n  largest_enclosure_cell=" << p.largest_enclosure_cell
       << ",\n  has_max_enclosure_been_reached=" << p.has_max_enclosure_been_reached
       << "\n)\n";
    return os;
}

//! \brief Outputs the lower discrete evolution statistics.
inline
std::ostream& 
operator<<(std::ostream& os, const LowerDiscreteEvolutionStatistics& p) 
{
    os << "LowerDiscreteEvolutionStatistics"
	   << ",\n  reach=" << p.reach
       << "\n)\n";
    return os;
}

//! \brief Outputs the upper discrete evolution statistics.
inline
std::ostream& 
operator<<(std::ostream& os, const UpperDiscreteEvolutionStatistics& p) 
{
    os << "UpperDiscreteEvolutionStatistics"
       << ",\n  has_restriction_occurred=" << p.has_restriction_occurred
	   << ",\n  reach=" << p.reach
       << "\n)\n";
    return os;
}

} //!namespace Ariadne

#endif //!ARIADNE_EVOLUTION_STATISTICS_H
