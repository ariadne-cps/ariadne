/***************************************************************************
 *            maintain.h
 *
 *  Wed Apr 14 10:04:12 2004
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
 
#ifndef _MAINTAIN_H
#define _MAINTAIN_H

#include <list>
	
#include "classes.h"

namespace Ariadne {

typedef struct {
	_Maintain_Changes *over;
	_Maintain_Changes *under;
} _Maintain_Approx_Changes;

/*! \brief This is the class of the phase space maintain system.
 *
 * Its methods are used to memorize phase space regions that 
 * are already reached in a particular location. This information 
 * is used both to avoid useless recomputations of the flow and 
 * to remember reached regions for reachability analysis output.	
 */
class Maintainer
{
		/*! \brief The over-approximating maintain structure. 
		 *
		 * This member of the class is an over-approximatin 
		 * basic maintainer.
		 */
		BasicMaintainer *over; 
	
		/*! \brief The under-approximating maintain structure. 
		 *
		 * This member of the class is an under-approximatin 
		 * basic maintainer.
		 */
		BasicMaintainer *under;
	
	
	public:
		/*! \brief This is a \a Maintainer class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a Maintainer. 
		 * \param over is the \a BasicMaintainer that will maintain 
		 * the over-approximation of the reached set.
		 * \param under is the \a BasicMaintainer that will maintain 
		 * the under-approximation of the reached set.
		 */
		Maintainer(BasicMaintainer* over, BasicMaintainer *under);

		/*! \brief This is a \a Maintainer class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a Maintainer using as template the object \a templ.
		 * \param templ is the template object for the constructor. 
		 */
		Maintainer(const Maintainer &templ); 
		
		/*! \brief This is the destructor of the class 
		 * \a Maintainer.
		 *
		 * This destructor deletes in a safe way an object of the class 
		 * \a Maintainer.
		 */
		~Maintainer(); 

		/*! \brief Memoizes a phase space region \a R.
		 *
		 * This method memoizes a phase space region \a R into the 
		 * maintain system's structure.
		 * \param R is the phase space region which should be memorized.
		 */
		void memoize(const Set &R);

		/*! \brief Checks if all the phase space states of \a R are 
		 * already memorized.
		 *
		 * This method checks if all the phase space states of \a R are 
		 * already memorized or not.
		 * \param R is the phase space region of which the complete 
		 * memorization should be tested.
		 * \return \a TRUE if all the states of \a R are already 
		 * memorized and \a FALSE otherwise.
		 */
		bool fully_memoized(const Set &R);
		
		/*! \brief Returns the region of phase space states of \a R 
		 * that are not previously memorized and memorizes \a R.
		 *
		 * This method returns the region of the phase space states 
		 * included in \a R that are not previously memorized and 
		 * memoizes \a R.
		 * \param R is the phase space state region which should 
		 * memorized.
		 * \return The region of states of R that are not previously 
		 * memorized.
		 */
		Set& get_unmemoized(const Set &R);
		
		/*! \brief Returns the reached regions.
		 *
		 * This method returns the phase space states memoized. 
		 * \return The phase space region memoized.
		 */
		RSet& get_memoized();
		
		/*! \brief Returns the changes of both maintain containers.
		 *
		 * This method returns the changes of both
		 * maintain containers since the last 
		 * \a clean_changes.
		 * \return The list of changes of both
		 * maintain container since the last 
		 * \a clean_changes.
		 */
		_Maintain_Approx_Changes *get_changes();
		
		/*! \brief Applies changes to both maintain 
		 * containers.
		 *
		 * This method applies the changes stored into
		 * \a changes to both maintain containers. 
		 */
		void apply_changes(_Maintain_Approx_Changes *changes);
		
		/*! \brief Cleans the maintain system's changes.*/
		void clean_alteration_lists();
};

/*! \brief This is the class of the phase space maintain system plus location.
 *
 * It manages a mantain system and the identificator of the location 
 * associated to it. Its methods are used to memorize phase space regions 
 * that are already reached in a particular location. This information is 
 * used both to avoid useless recomputations of the flow and to remember
 * reached regions for reachability analysis output.
 * @see MaintainSystem	
 */
class HMaintainer : public Maintainer
{
		/*! \brief The location to which the maintain 
		 * system is refered. */
		Location *l;
	
	public:
		/*! \brief This is a \a HMaintainer class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a HMaintainer.
		 * \param over is the \a BasicMaintainer that will maintain 
		 * the over-approximation of the reached set.
		 * \param under is the \a BasicMaintainer that will maintain 
		 * the under-approximation of the reached set.
		 */
		HMaintainer(BasicMaintainer *over, BasicMaintainer *under);
		
		/*! \brief This is a \a HMaintainer class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a HMaintainer and set the locations \a l. 
		 * \param over is the \a BasicMaintainer that will maintain 
		 * the over-approximation of the reached set.
		 * \param under is the \a BasicMaintainer that will maintain 
		 * the under-approximation of the reached set.
		 * \param l is the location associated to this maintain system.
		 */
		HMaintainer(BasicMaintainer *over, BasicMaintainer *under, 
				const Location *l);
	
		/*! \brief This is a \a HMaintainer class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a HMaintainer using as template the object 
		 * \a templ.
		 * \param templ is the template object for the constructor. 
		 */
		HMaintainer(const HMaintainer &templ);
		
		/*! \brief This is the destructor of the class 
		 * \a HMaintainer.
		 *
		 * This destructor deletes in a safe way an object of the class 
		 * \a HMaintainer.
		 */
		~HMaintainer();
	
		/*! \brief Set the location.
		 *
		 * This method set the location to which the object refers to.
		 * \param l is the new location. 
		 */
		void set_location(const Location *l);
		
		/*! \brief Get the location.
		 *
		 * This method return a pointer to the location to which 
		 * the object refers to.
		 * \return The location to which the object refers to. 
		 */
		Location* get_location() const;
		
};

/*! \brief This is the class of the state space maintain system
 *
 * Its methods are used to memorize regions that are already 
 * reached in a particular location. This information is used both 
 * to avoid useless recomputations of the flow and to remember 
 * reached regions for reachability analysis output. 
 * @see MaintainSystem
 */
class MaintainSystem
{
		/*! \brief Contains a maintain system for each 
		 * location. 
		 */
		std::list<HMaintainer *> maintain_list;
		//HMaintainer* maintain_list;

		/*! \brief Number of maintained locations. 
		 */
		unsigned int loc_num;
	public:
		/*! \brief This is a \a MaintainSystem class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a MaintainSystem. Objects initializated with this method 
		 * can maintain reachability region of one location only.
		 * \param over is the \a BasicMaintainer that will maintain 
		 * the over-approximation of the reached set.
		 * \param under is the \a BasicMaintainer that will maintain 
		 * the under-approximation of the reached set.
		 */
		MaintainSystem(BasicMaintainer *over, BasicMaintainer *under);
	
		/*! \brief This is a \a MaintainSystem class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a MaintainSystem using as template the object \a templ.
		 * \param templ is the template object for the constructor. 
		 */
		MaintainSystem(const MaintainSystem &templ);
	
		/*! \brief This is a \a MaintainSystem class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a MaintainSystem using as template the structures 
		 * used for the location \a l by the object \a templ.
		 * \param templ is the template object for the constructor.
		 * \param l is the location of which the memorization 
		 * structures should be copied.			 
		 */
		MaintainSystem(const MaintainSystem &templ, const Location *l);
	
		/*! \brief This is a \a MaintainSystem class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a MaintainSystem. The parameter indicates the number 
		 * of locations of which reachability region should be 
		 * maintained.
		 * \param loc_num is the number of locations of which 
		 * reachability region should be maintained.
		 */
		MaintainSystem(const unsigned int loc_num);
		
		/*! \brief This is a \a MaintainSystem class constructor.
		 *
		 * This constructor initializes the object of the class
		 * \a MaintainSystem for the location graph \a loc_graph.
		 * \param automaton is the automaton of which reachability
		 * region should be maintained.
		 */
		MaintainSystem(const Automaton &automaton);
		
		/*! \brief This is the destructor of the class 
		 * \a MaintainSystem.
		 *
		 * This destructor deletes in a safe way an object of the class 
		 * \a MaintainSystem.
		 */
		~MaintainSystem();
		
		/*! \brief Memoizes a region \a R.
		 * 
		 * This method memoizes a region \a R into the maintain system
		 * including the continuous part of the region in the reached 
		 * states of the region's location.
		 * \param R is the region which should be memorized.
		 */
		void memoize(const HSet &R);
	
		/*! \brief Memoizes a phase space region \a R.
		 *
		 * This method memoizes a phase space region \a R into the 
		 * maintain system structure corrisponding to the location \a l.
		 * \param R is the phase space region which should be memorized.
		 * \param l is the location for which \a R is memorized.
		 */
		void memoize(const Set &R, const Location *l);
		
		/*! \brief Checks if all the states of \a R are already 
		 * memorized.
		 *
		 * This method checks if all the states of \a R are already 
		 * memorized or not in the object.
		 * \param R is the region of which the complete memorization 
		 * should be tested.
		 * \return \a TRUE if all the states of \a R are already 
		 * memorized and \a FALSE otherwise.
		 */
		bool fully_memoized(const HSet &R);
		
		/*! \brief Checks if all the phase space states of \a R are 
		 * already memorized.
		 *
		 * This method checks if all the phase space states of \a R are 
		 * already memorized or not in location \a l.
		 * \param R is the phase space region of which the complete 
		 * memorization should be tested.
		 * \param l is the location for which memorization should be 
		 * tested.
		 * \return \a TRUE if all the states of \a R are already 
		 * memorized for \a l and \a FALSE otherwise.
		 */
		bool fully_memoized(const Set &R, const Location *l);
		
		/*! \brief Returns the region of states of \a R that are not 
		 * previously memorized and memoizes \a R.
		 *
		 * This method returns the region of the states included in 
		 * \a R that are not previously memorized and memoizes \a R.
		 * \param R is the region which should memoized.
		 * \return The region of states of \a R that are not previously 
		 * reached.
		 */
		HSet* get_unmemoized(const HSet &R);
		
		/*! \brief Returns the region of phase space states of \a R 
		 * that are not previously reached in \a l and includes \a R 
		 * into the reached regions.
		 *
		 * This method returns the region of the phase space states 
		 * included in \a R that are not previously reached in the 
		 * location \a l and includes \a R into the reached regions.
		 * \param R is the phase space state region which should 
		 * included into the reached region.
		 * \param l is the location to which R is refered.
		 * \return The region of states of \a R that are not previously 
		 * reached in \a l.
		 */
		Set* get_unmemoized(const Set &R, const Location *l);
		
		/*! \brief Returns the memorizated regions.
		 *
		 * This method returns the region memorized for location \a l. 
		 * \param l is the location for which reached region should be 
		 * returned.
		 * \return The region of states that are memorized for 
		 * location \a l.
		 */
		RSet* get_memoized(const Location* l);
		
		/*! \brief Returns the memorizated regions.
		 *
		 * This method returns the phase space states memoized. 
		 * \return The phase space region memoized.
		 */
		RSet* get_memoized();
		
		/*! \brief Returns the changes of
		 * maintain containers.
		 *
		 * This method returns the changes of both 
		 * maintain containers of 
		 * the location \a l since the last 
		 * \a clean_changes.
		 * \return The changes of 
		 * both maintain containers of 
		 * the location \a l since the last 
		 * \a clean_changes.
		 */
		_Maintain_Approx_Changes* get_changes(const Location *l);

		/*! \brief Returns the changes of
		 * maintain containers.
		 *
		 * This method returns the list of changes of 
		 * maintain containers since the last 
		 * \a clean_changes.
		 * \return The list of changes of  
		 * maintain containers since the last \a clean_changes.
		 */
		std::list<_Maintain_Approx_Changes *> *get_changes();
		
		/*! \brief Applies the changes of
		 * maintain containers.
		 *
		 * This method applies the changes stored into 
		 * \a changes to both basic maintainers of location \a l.
		 * \param l is the location of the maintainer to which 
		 * changes should be applyed.
		 * \param changes are the changes that should be applyed.
		 */
		void apply_changes(const Location *l,
					_Maintain_Approx_Changes *changes);

		/*! \brief Applies the changes of basic maintainers.
		 *
		 * This method applies the changes stored into 
		 * \a changes to both basic maintainer associated to each 
		 * the automaton's location.
		 * \param changes are the changes that should be applyed
		 */
		void apply_changes(std::list<_Maintain_Approx_Changes *> 
				*changes);
		
		/*! \brief Cleans the changes.
		 *
		 * This method cleans the maintain system's changes
		 * of the locatation \a l.
		 * \param l is the location for which alteretion lists should
		 * be cleaned.
		 */
		void clean_changes(const Location *l);
		 
		/*! \brief Cleans the changes.
		 *
		 * This method cleans the maintain system's changes
		 * for all locations.
		 */
		void clean_changes();
};

}

#endif /* _MAINTAIN_H */
