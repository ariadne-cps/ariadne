/***************************************************************************
 *            cluster_list.h
 *
 *  Fri Apr 16 11:56:12 2004
 *  Copyright  2004  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
#ifndef _CLUSTER_LIST_H
#define _CLUSTER_LIST_H

#include <list>
#include <vector>
	
namespace Ariadne {
	
/*! \brief Describes an interval in the real space.
 *
 *  This class describes a closed interval in the real space 
 *  domain.
 */
class Interval
{
		/*! \brief The upper bound of the interval */
		double upper;
	
		/*! \brief The upper bound of the interval */
		double lower;
	public:
		/*! \brief This is a interval object constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a Inteval.
		 * \param lower is the lower bound of the Interval.
		 * \param upper is the upper bound of the Interval.
		 */
		Interval(double lower, double upper);

		/*! \brief This is the destructor of the class 
		 * \a Interval.
		 *
		 * This destructor deletes in a safe way an object of the class 
		 * \a Interval.
		 */
		~Interval();

		/*! \brief Returns the interval's lower bound.
		 *
		 * \return The interval's lower bound.
		 */
		double get_lower();
		
		/*! \brief Returns the interval's upper bound.
		 *
		 * \return The interval's upper bound.
		 */
		double get_upper();

		/*! \brief Sets the interval's lower bound.
		 *
		 * \param lbound is the new lower bound.
		 */
		void set_lower(double lbound);
		
		/*! \brief Sets the interval's upper bound.
		 *
		 * \param ubound is the new upper bound.
		 */
		void set_upper(double ubound);
};

/*! \typedef SpaceBox
 *  \brief Represents boxes in the real \a n-dimentional space.
 *
 *  The objects of this class represent space boxes in the 
 *  \a n-dimentional real space.
 */
typedef std::vector<Interval> SpaceBox;

/*! \typedef ClusterList
 *  \brief Represents a list of boxes in the real \a n-dimentional space.
 *
 * The objects of this class represent a list of boxes in the real 
 * \a n-dimentional space.
 */
typedef std::list<SpaceBox *> ClusterList;

}

#endif /* _CLUSTER_LIST_H */
