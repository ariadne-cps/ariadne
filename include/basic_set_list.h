/***************************************************************************
 *            basic_set_list.h
 *
 *  Thu Aug 19 17:31:59 2004
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
 
#ifndef _BASIC_SET_LIST_H
#define _BASIC_SET_LIST_H

#include <list>

#include "approx_type.h"
#include "state.h"
#include "classes.h"

namespace Ariadne {

/*! \class BasicSetList
 * \brief Rappresents lists of basic set. 
 */
class BasicSetList {

	protected:
		/*! \brief The basic set list */
		std::list<BasicSet *> l;
				
		/*! \brief Makes the union of two basic set 
		 * list.
		 *
		 * This method attach a copy of a 
		 * basic set \a l1 to the current one.
		 * \param l1 is the basic set list which should
		 * be concatenated to the current list.
		 * \return A reference to the new list.
		 */
		BasicSetList& attach_copy(const BasicSetList& l1);

#ifdef FASTER_BUT_DIRTIER
		/*! \brief Attaches two list.
		 *
		 * This method attaches the basic set \a l1 to 
		 * the current one.
		 * \param l1 is the basic set list which should
		 * be attached to the current list.
		 * \return A reference to the new list.
		 */
		BasicSetList& attach(BasicSetList& l1);
#endif

		/*! \brief Removes the sets such that they are 
		 * subsets of other objects into the list.
		 *
		 * This method removes from the list and deletes
		 * the basic sets such that they are included by 
		 * other sets contained into the list.
		 * \return A reference to the new list.
		 */
		BasicSetList& remove_doubly_included_set();
		
	public:
		
		/*! \brief Clears the list 
		 *
		 * This method clears the list deleting
		 * all the objects contined in it.
		 */
		void clear();
		
		/*! \brief This is a \a BasicSetList class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * {\a BasicSetList}.
		 */
		BasicSetList();

		/*! \brief This is a \a BasicSetList class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * {\a BasicSetList}.
		 * \param tl is a template for the new object.
		 */
		BasicSetList(const BasicSetList& tl);
		
		/*! \brief This is the destructor of the class 
		 * {\a BasicSetList}.
		 *
		 * This destructor deletes in a safe way an object 
		 * of the class {\a BasicSetList}.
		 */
		~BasicSetList();

		/*! \brief Makes the union of two basic set lists.
		 *
		 * This method concatenates the current object's list 
		 * with the one of \a l1 and returns the resulting list 
		 * on a new object.
		 * The current object does not change.
		 * \param l1 is the basic set list which 
		 * should be concatenated to the current object's 
		 * list.
		 * \return A \a BasicSetList containing the union
		 * of thecurrent object and {\a l1}.
		 */
		BasicSetList join(const BasicSetList& l1) const;
		
#ifdef FASTER_BUT_DIRTIER
			
		/*! \brief Makes the union of a basic set 
		 * list with the current one.
		 *
		 * This method concatenates a basic set \a l1
		 * to the current one.
		 * The object is set to the union.
		 * \param l1 is the basic set list which should
		 * be concatenated to the current list.
		 * \return A reference to the new list.
		 */
		BasicSetList& join_on_obj(const BasicSetList& l1);
#endif
		/*! \brief Makes the intersection of a basic set list 
		 * with a basic set.
		 *
		 * This method intersects the current object's list 
		 * and the basic set \a bset returning the resulting 
		 * list on a new object.
		 * \param bset is the basic set which 
		 * should be intersected with the current object.
		 * \param atype is the type of approximation 
		 * requested. If \a atype=NONE and the 
		 * intersection is not exacly representable,
		 * an over-approximation is done.
		 * \return A \a BasicSetList containing the
		 * intersection of the current object
		 * with {\a bset}.
		 */
		BasicSetList intersect(const BasicSet& bset, 
				ApproxType atype) const;

		
		/*! \brief Makes the intersection of two basic 
		 * set lists.
		 *
		 * This method intersects the sets represented
		 * by the current object's list and \a l2 
		 * returning the resulting list on a new object.
		 * \param l1 is the basic set list which 
		 * should be intersected with the current object.
		 * \param atype is the type of approximation 
		 * requested. If \a atype=NONE and the 
		 * intersection is not exacly representable,
		 * an over-approximation is done.
		 * \return A \a BasicSetList containing the
		 * intersection of the current object's list 
		 * with {\a l1}.
		 */
		BasicSetList intersect(const BasicSetList& l1, 
				ApproxType atype) const;

#ifdef FASTER_BUT_DIRTIER
		/*! \brief Makes the intersection of a basic set list 
		 * with a basic set.
		 *
		 * This method intersects the current object's list 
		 * and the basic set \a bset.
		 * The object is set to the intersection.
		 * \param bset is the basic set which 
		 * should be intersected with the current object.
		 * \param atype is the type of approximation 
		 * requested. If \a atype=NONE and the 
		 * intersection is not exacly representable,
		 * an over-approximation is done.
		 * \return A reference to the new list.
		 */
		BasicSetList& intersect_on_obj(const BasicSet& bset,
				const ApproxType atype) const;
#endif

		/*! \brief Copies a basic set list on the current 
		 * object.
		 *
		 * This method copies a basic set list on the current 
		 * object. The basic sets initially contained by the 
		 * current object are removed from the list and 
		 * deleted.
		 * \param tl is the basic set list which should be
		 * copied on the current object.
		 * \return A reference to the new list.
		 */
		BasicSetList& operator=(const BasicSetList& tl);

		/*! \brief Checks whenever the current basic set list
		 * is contained by another list.
		 *
		 * This method checks whenever the current basic set 
		 * list is contained by another list {\a l1}.
		 * \param l1 is the list which should include the
		 * current object.
		 * \return {\a true}, if the current list represents
		 * a subset of {\a l1}, \a false otherwise.
		 */
		bool is_subset_of(const BasicSetList& l1) const;

		/*! \brief Checks whenever the current basic set list
		 * intersects another list.
		 *
		 * This method checks whenever the current basic set 
		 * list intersects another list {\a l1}.
		 * \param l1 is the list which should intersect the
		 * current object.
		 * \return {\a true}, if the current list intersects 
		 * {\a l1}, \a false otherwise.
		 */
		bool does_intersect(const BasicSetList& l1) const;
		
		/*! \brief Checks whenever the current basic set list
		 * intersects a basic set.
		 *
		 * This method checks whenever the current basic set 
		 * list intersects the basic set {\a bset}.
		 * \param bset is the basic set which should intersect 
		 * the current object.
		 * \return {\a true}, if the current list intersects 
		 * {\a bset}, \a false otherwise.
		 */
		bool does_intersect(const BasicSet& bset) const;
	
		/*! \brief Includes into the current list a basic set.
		 *
		 * This method includes into the current list the 
		 * basic set {\a bset}.
		 * \param bset is the basic set which should be 
		 * included into the current list.
		 * \return A reference to the new list.
		 */
		BasicSetList& include(const BasicSet &bset);
		
		/*! \brief Includes into the current list a basic set.
		 *
		 * This method includes into the current list the 
		 * basic set {\a bset}.
		 * \param bset is the basic set which should be 
		 * included into the current list.
		 * \return A reference to the new list.
		 */
		BasicSetList& include(BasicSet *bset);

		/*! \brief Checks whenever the current basic set list
		 * include a state.
		 *
		 * This method checks whenever the current basic set list
		 * include the state {\a s}.
		 * \param s is the state of which inclusion 
		 * into the current list should be tested.
		 * \return  {\a true}, if {\a s} is contained into the 
		 * current list, \a false otherwise.
		 */
		bool contains(const State &s) const;

		/*! \brief Checks whenever the current basic set list
		 * include a basic set.
		 *
		 * This method checks whenever the current basic set list
		 * include the basic set {\a bset}.
		 * \param bset is the basic set of which inclusion 
		 * into the current list should be tested.
		 * \return  {\a true}, if {\a bset} is contained into the 
		 * current list, \a false otherwise.
		 */
		bool contains(const BasicSet &bset) const;

		/*! \brief Compacts the current list.
		 *
		 * Thi method approximate the current list with a list 
		 * containing less then \a max_elem+1 basic sets.
		 * The parameter \a atype indicated the type of 
		 * approximation wanted. If \a atype=NONE and the 
		 * compact list can not represent exacly the 
		 * current list, an over-approxiamtion is done.
		 * \param max_elem is the compact list's maximum 
		 * number of basic set.
		 * \param atype is the type of approximation done.
		 * \return A reference to the new list.
		 */
		BasicSetList& compact_list(const unsigned int max_elem, 
				ApproxType atype);
		
		/*! \brief Compacts the current list.
		 *
		 * Thi method approximate the current list with a list 
		 * containing less then \a max_elem+1 basic sets.
		 * If the compact list can not represent exacly the 
		 * current list, an over-approxiamtion is done.
		 * \param max_elem is the compact list's maximum 
		 * number of basic set.
		 * \return A reference to the new list.
		 */
		BasicSetList& compact_list(const unsigned int max_elem);

		/*! \brief Checks whenever the current list is empty
		 * or contains only empty sets.
		 *
		 * \return \a true if either the current list is empty
		 * or contains only empty sets, \a false otherwise.
		 */
		bool empty();

		/*! \brief Applies the map \a M to the list's basic 
		 * sets.
		 *
		 * This method applies the trasformation \a M the list's basic 
		 * sets.
		 * The object does not change.
		 * \param M is the trasformation to be applied.
		 * \param atype is the type of the approximation.
		 * \return The trasformed list.
 		 */
		BasicSetList apply(const Map &M, 
				const ApproxType atype) const;
		
		/*! \brief Applies the map \a M to the list's basic 
		 * sets.
		 *
		 * This method applies the trasformation \a M the list's basic 
		 * sets.
		 * The object does not change.
		 * \param M is the trasformation to be applied.
		 * \param atype is the type of the approximation.
		 * \return A pointer to the trasformed list.
 		 */
		BasicSetList* apply_p(const Map &M, 
				const ApproxType atype) const;
		
#ifdef FASTER_BUT_DIRTIER
		/*! \brief Applies the map \a M to the list's basic 
		 * sets.
		 *
		 * This method applies the trasformation \a M the list's basic 
		 * sets.
		 * The object is set to the trasformed list.
		 * \param M is the trasformation to be applied.
		 * \param atype is the type of the approximation.
		 * \return A reference to the transformed list.
 		 */
		BasicSetList& apply_on_obj(const Map &M, 
				const ApproxType atype);
#endif

		/*! \brief Evaluates the flow slice from the list's basic
		 * sets. 
		 *
		 * This method evaluates the list of the \f$\delta\f$-timed 
		 * flow slices of each basic sets in the current object's list.
		 * The current object does not change.
		 * \param field is the vector field.
		 * \param delta is the time of the slice.
		 * \param atype is the approximation type.
		 * \return A pointer to the list of the flow slices.
 		 */
		BasicSetList* flow_slice_p(const VectorField &field,
				const double delta,
				const ApproxType atype) const;

#ifdef FASTER_BUT_DIRTIER
		/*! \brief Evaluates the flow slice from the list's basic
		 * sets. 
		 *
		 * This method evaluates the list of the \f$\delta\f$-timed 
		 * flow slices of each basic sets in the current object's list.
		 * The current object is set to the list of slices.
		 * \param field is the vector field.
		 * \param delta is the time of the slice.
		 * \param atype is the approximation type.
		 * \return A pointer to the list of the flow slices.
 		 */
		BasicSetList& flow_slice_on_obj(const VectorField &field,
				const double delta,
				const ApproxType atype);

#endif
		/* TODO: comment the following methods*/
		BasicSetList expand_of(const double delta) const;
		BasicSetList* expand_of_p(const double delta) const;
#ifdef FASTER_BUT_DIRTIER
		BasicSetList& expand_of_on_obj(const double delta);
#endif
		BasicSetList shrink_of(const double delta) const;
		BasicSetList* shrink_of_p(const double delta) const;
#ifdef FASTER_BUT_DIRTIER
		BasicSetList& shrink_of_on_obj(const double delta);
#endif

};

class Ex_BSL;

/*! \brief Represents sets listing of basic sets included into such sets. */
class In_BSL : public BasicSetList {
	public:
		In_BSL& operator=(const BasicSetList &tl);
		
		In_BSL& compact_list(Ex_BSL &ex, const unsigned int max_elem,
				const ApproxType atype);

		friend class Ex_BSL;
};

/*! \brief Represents sets listing of basic sets that does not intersect
 * such sets. */
class Ex_BSL : public BasicSetList {
	public:
		Ex_BSL& operator=(const BasicSetList &tl);
		
		Ex_BSL& compact_list(In_BSL &in, const unsigned int max_elem,
				const ApproxType atype);
		
		friend class In_BSL;
};

}

#endif /* _BASIC_SET_LIST_H */
