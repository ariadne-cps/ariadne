/***************************************************************************
 *            set.h
 *
 *  Wed Apr  7 11:57:59 2004
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
 
#ifndef _SET_H
#define _SET_H

#include <list>
	
#include "classes.h"
#include "approx_type.h"
#include "basic_set_list.h"

#define DEFAULT_MAX_OBJECTS 20

namespace Ariadne {
	
/*! \brief Rappresents an approximation of a phase space region. 
 *
 * Each object of this class maintains two lists of basic set, 
 * \a in_list and \a ex_list of at most \a max_basic_sets elements.
 * The set \a S represented by an ASet object is 
 * \f$ S = \left( \cup_{i \in in\_list } i \right) \setminus \left( \cup_{j \in ex\_list } j \right)\f$. 
 * The type of approximation is maintained by the member 
 * {\a atype}.  
 */
class ASet
{
	protected:
		/*! \brief The maximum number of basic sets representing the 
		 * region.
	   	 *
		 * This member indicates the maximum number of basic sets in
		 * each of the class's lists.
		 */
		unsigned int max_basic_sets;
	
		/*! \brief Type of approximating region.
		 *
		 * This member specifies if the approximating set if
		 * over-approximating or a under-approximating set.
		 */
		ApproxType atype;
	
		/*! \brief The list of basic sets of the region.
		 *
		 * This member is a list of basic sets representing a part
		 * of the region.
		 */	
		In_BSL in_list;
	
		/*! \brief The list of basic sets not included in the region.
		 *
		 * This member is a list of basic sets representing a set of 
		 * states not included in the region.
		 */	
		Ex_BSL ex_list;

	protected:

		/*! \brief Expands or shrinks the set to condider an error.
		 *
		 * This method expands the set of \a e if the set
		 * is an OVER approximation, otherwise it shrink it of 
		 * \a e .
		 * \param e is the error done durring the computation
		 * that should be taken in account by the set.
		 * \return A reference to the expanded/shrinked set.
		 */
		ASet& take_in_account_error(const double e);

	public:
		/*! \brief This is a ASet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * ASet.
		 */
		ASet(); 
		
		/*! \brief This is a ASet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * ASet using as template the parameter \a t.
		 * \param t is the template object.
		 */
		ASet(const ASet &t); 
		
		/*! \brief This is a ASet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * ASet.
		 * \param max_basic_sets is the maximum number of basic sets 
		 * usable to represents the region.
		 */
		ASet(unsigned int max_basic_sets); 
		
		/*! \brief This is a ASet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * ASet setting the type to \a type.
		 * \param type is the type of approximation made by the 
		 * object.
		 */
		ASet(ApproxType type); 
		
		/*! \brief This is a ASet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * ASet setting the type to \a type.
		 * \param max_basic_sets is the maximum number of basic sets 
		 * usable to represents the region.
		 * \param type is the type of approximation made by the 
		 * object.
		 */
		ASet(unsigned int max_basic_sets, ApproxType type); 

		/*! \brief This is the destructor of the class ASet.
		 *
		 * This destructor deletes in a safe way an object of the class 
		 * ASet.
		 */
		~ASet(); 
	
		/*! \brief Sets the number of basic sets used in the region
		 * reppresentation.
		 *
		 * This method sets the maximum number of basic sets used in 
		 * the representation of the region.
		 * \param max_basic_sets is the number of basic sets used in 
		 * the region's reppresentation.
		 * \return A reference to the new object.
 		 */
		ASet& set_max_basic_sets(unsigned int max_basic_sets); 
	
		/*! \brief Gets the number of basic sets used in the region's
		 * reppresentation.
		 *
		 * This method gets the maximum number of basic sets used in 
		 * the representation of the region.
		 * \return The number of basic sets used in the region's 
		 * reppresentation.
 		 */
		unsigned int get_max_basic_sets(); 
	
		/*! \brief Includes a basic set in the region.
		 *
		 * This method includes a basic set in the 
		 * region.		 
		 * \param bset is the basic set that should be 
		 * included.
		 * \return A reference to the new object.
 		 */
		ASet& include(const BasicSet &bset); 

		/*! \brief Includes a basic set in the region.
		 *
		 * This method includes a basic set in the 
		 * region.		 
		 * \param bset is the basic set that should be 
		 * included.
		 * \return A reference to the new object.
 		 */
		ASet& include(BasicSet *bset); 
		
		/*! \brief Excludes a basic set from the region.
		 *
		 * This method excludes a set of phase space states from the 
		 * set adding it to the \a ex_list list.
		 * \param set is the basic set that should be
		 * excluded.
 		 * \return A reference to the new object.
		 */
		ASet& exclude(const BasicSet &set); 
		
		/*! \brief Excludes a basic set from the region.
		 *
		 * This method excludes a set of phase space states from the 
		 * set adding it to the \a ex_list list.
		 * \param set is the basic set that should be
		 * excluded.
 		 * \return A reference to the new object.
		 */
		ASet& exclude(BasicSet *set); 
		
		/*! \brief Compacts the lists of basic sets.
		 *
		 * This method tries to compact the lists of basic sets.
		 * \return A reference to the new object.
		 */
		ASet& compact_lists();

		/*! \brief Checks if the intersection of two set is not null.
		 *
		 * This method checks if the intersection between the current 
		 * set and the set R is not null. If it is not null, the 
		 * method returns TRUE and FALSE otherwise.		 
		 * \param R is the set of which the intersection is to be 
		 * tested.
		 * \return TRUE if the intersection between the set and R 
		 * is not null, FALSE otherwise.
 		 */
		bool does_intersect(const ASet &R) const;

		/*! \brief Tests whether the set is empty. */
		bool empty();

		/*! \brief Tests whether the set is a subset of other.*/
		bool is_subset_of(const ASet& other) const;
		
		/*! \brief Tests whether the set intersects other. */
		bool is_disjoint_from(const ASet& other) const;
		
		/*! \brief Tests whether the set contains the specified state.*/
		bool contains(const State &s) const;

		/*! \brief Copies a \a ASet on the current object.
		 *
		 * This method copies \a ts on the current object. 
		 * The set initially contained by the current 
		 * object is deleted.
		 * \param ts is the set which should be
		 * copied on the current object.
		 * \return A reference to the new object.
		 */
		ASet& operator=(const ASet& ts);
		
		/*! \brief Intersects the set with a set of phase space 
		 * states.
		 *
		 * This method intersects the set with a set of phase space 
		 * states.
		 * The current object does not change.
		 * \param bset is the set of phase space states to be 
		 * intersected.
 		 * \return A new \a ASet representing the intersection.
		 */
		ASet operator&&(const BasicSet &bset) const;
		
		/*! \brief Intersects the set with a set of phase space 
		 * states.
		 *
		 * This method intersects the set with a set of phase space 
		 * states.
		 * The current object does not change.
		 * \param bset is the set of phase space states to be 
		 * intersected.
 		 * \return A pointer to a new \a ASet representing the 
		 * intersection.
		 */
		ASet* intersect(const BasicSet &bset) const; 
		
		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the current 
		 * set and the set \a R.		 
		 * \param R is the set of which the intersection is to be 
		 * computed.
		 * \return The intersection between the set and \a R.
 		 */
		ASet operator&&(const ASet &R) const;
		
		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the current 
		 * set and the set \a R.		 
		 * \param R is the set of which the intersection is to be 
		 * computed.
		 * \return A pointer to the intersection between the set 
		 * and \a R.
 		 */
		ASet* intersect(const ASet &R) const;

		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the current 
		 * set and the set \a R.		 
		 * \param R is the set of which the intersection is to be 
		 * computed.
		 * \param type is the wanted type of approximation.
		 * \return A pointer to the intersection between the set 
		 * and \a R.
 		 */
		ASet* intersect(const ASet &R, const ApproxType type) const;

		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return The union of the current set and \a R.
		 */
		ASet operator||(const ASet &R) const;

		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return A pointer to the union of the current set and 
		 * \a R.
		 */
		ASet* join(const ASet &R) const;

		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \param type is the wanted type of approximation.
		 * \return A pointer to the union of the current set and 
		 * \a R.
		 */
		ASet* join(const ASet &R, const ApproxType type) const;
		
		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the current 
		 * set and the set \a R.		 
		 * \param R is the set of which the intersection is to be 
		 * computed.
		 * \return The intersection between the set and \a R.
 		 */
		Set operator&&(const Set &R) const;
		
		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the current 
		 * set and the set \a R.		 
		 * \param R is the set of which the intersection is to be 
		 * computed.
		 * \return A pointer to the intersection between the set 
		 * and \a R.
 		 */
		Set* intersect(const Set &R) const;

		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return The union of the current set and \a R.
		 */
		Set operator||(const Set &R) const;

		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return A pointer to the union of the current set and 
		 * \a R.
		 */
		Set* join(const Set &R) const;
	
		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the current 
		 * set and the set \a R.		 
		 * \param R is the set of which the intersection is to be 
		 * computed.
		 * \return The intersection between the set and \a R.
 		 */
		HSet operator&&(const HSet &R) const;
		
		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the current 
		 * set and the set \a R.		 
		 * \param R is the set of which the intersection is to be 
		 * computed.
		 * \return A pointer to the intersection between the set 
		 * and \a R.
 		 */
		HSet* intersect(const HSet &R) const;

		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return The union of the current set and \a R.
		 */
		HSet operator||(const HSet &R) const;

		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return A pointer to the union of the current set and 
		 * \a R.
		 */
		HSet* join(const HSet &R) const;

		/*! \brief Applies the map \a M to the current set.
		 *
		 * This method applies the trasformation \a M to the set.
		 * The object does not change.
		 * \param M is the trasformation to be applied.
		 * \return The trasformed set.
 		 */
		ASet operator<<(const Map &M) const;
		
		/*! \brief Applies the map \a M to the current set.
		 *
		 * This method applies the trasformation \a M to the set.
		 * The object does not change.
		 * \param M is the trasformation to be applied.
		 * \return A pointer to the trasformed set.
 		 */
		ASet* apply(const Map &M) const;

#ifdef FASTER_BUT_DIRTIER
		
		/*! \brief Applies the map \a M to the current set.
		 *
		 * This method applies the trasformation \a M to the set.
		 * The object is set to the trasformed object.
		 * \param M is the trasformation to be applied.
		 * \return A reference to the transformed object.
 		 */
		ASet& apply_on_obj(const Map &M);
#endif
		
		/*! \brief Evaluates the flow slice from the current set. 
		 *
		 * This method evaluates the \f$\delta\f$-timed flow slice 
		 * from the current set dued to the vector field \a field .
		 * \param field is the vector field.
		 * \param delta is the time of the slice.
		 * \return The flow slice.
 		 */
		ASet* flow_slice_p(const VectorField &field,
				const double delta) const;

		/*! \brief Evaluates the flow. 
		 *
		 * This method evaluates the set reachable in exactely 
		 * \f$delta\f$-timed from the current set with a flow 
		 * dued to the vector field \a field .
		 * \param field is the vector field.
		 * \param delta is the time of the flow.
		 * \return The set reached after a \f$delta\f$-timed flow 
		 * from the current set.
 		 */
		ASet* flow_p(const VectorField &field,
				const double delta) const;
		
		/*! \brief Evaluates the flow slice from the current set. 
		 *
		 * This method evaluates the \f$\delta\f$-timed flow slice 
		 * from the current set dued to the vector field \a field .
		 * \param field is the vector field.
		 * \param delta is the time of the slice.
		 * \return The flow slice.
 		 */
		ASet flow_slice(const VectorField &field,
				const double delta) const;

		/*! \brief Evaluates the flow. 
		 *
		 * This method evaluates the set reachable in exactely 
		 * \f$delta\f$-timed from the current set with a flow 
		 * dued to the vector field \a field .
		 * \param field is the vector field.
		 * \param delta is the time of the flow.
		 * \return The set reached after a \f$delta\f$-timed flow 
		 * from the current set.
 		 */
		ASet flow(const VectorField &field,
				const double delta) const;

#ifdef FASTER_BUT_DIRTIER
		/*! \brief Evaluates the flow slice from the current set. 
		 *
		 * This method evaluates the \f$\delta\f$-timed flow slice 
		 * from the current set dued to the vector field \a field .
		 * \param field is the vector field.
		 * \param delta is the time of the slice.
		 * \return The flow slice.
 		 */
		ASet& flow_slice_on_obj(const VectorField &field,
				const double delta);

		/*! \brief Evaluates the flow. 
		 *
		 * This method evaluates the set reachable in exactely 
		 * \f$delta\f$-timed from the current set with a flow 
		 * dued to the vector field \a field .
		 * \param field is the vector field.
		 * \param delta is the time of the flow.
		 * \return The set reached after a \f$delta\f$-timed flow 
		 * from the current set.
 		 */
		ASet& flow_on_obj(const VectorField &field,
				const double delta);
#endif		
};


/*! \brief Rappresents an over-approximation of a phase 
 * space region. 
 */
class OverSet: public ASet {
	public:
		/*! \brief This is a OverSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * OverSet.
		 * \param max_basic_sets is the maximum number of basic 
		 * sets usable to represents the region.
		 */
		OverSet(unsigned int max_basic_sets);
		
		/*! \brief This is a OverSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * OverSet.
		 * \param t is the template for the new OverSet.
		 */
		OverSet(const ASet &t);
			
		/*! \brief This is a OverSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * OverSet.
		 */
		OverSet();

};

/*! \brief Rappresents an under-approximation of a phase 
 * space region. 
 */
class UnderSet: public ASet {
	public:
		/*! \brief This is a UnderSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * UnderSet.
		 * \param max_basic_sets is the maximum number of basic 
		 * sets usable to represents the region.
		 */
		UnderSet(unsigned int max_basic_sets); 

		/*! \brief This is a UnderSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * UnderSet.
		 * \param t is the template for the new UnderSet.
		 */
		UnderSet(const ASet &t);
		
		/*! \brief This is a UnderSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * UnderSet.
		 */
		UnderSet(); 

};

/*! \brief This is the class of set of phase space states.
 *
 * It maintains both over and under approximations.
 */
class Set 
{
	protected:
		/*! \brief The over-approximation of the set.
		 *
		 * This member maintains the over-approximation of the 
		 * current set. 
		 */ 
		OverSet *over_appr;
	
		/*! \brief The under-approximation of the set.
		 *
		 * This member maintains the under-approximation of the 
		 * current set. 
		 */ 
		UnderSet *under_appr;
		
	public:
		/*! \brief This is a Set class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * Set.
		 */
		Set();
		
		/*! \brief This is a \a Set class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a Set using as template the parameter \a t.
		 * \param t is the template object.
		 */
		Set(const Set &t);
		
		/*! \brief This is a \a Set class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a Set using the parameters \a over and \a under 
		 * as templates for set's over and under
		 * approximation respectevly.
		 * \param over is the template for the over 
		 * approximation of the set.
		 * \param under is the template for the under 
		 * approximation of the set.
		 */
		Set(const ASet &over,const ASet &under);

		/*! \brief This is a \a Set class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a Set using the parameters \a over and \a under 
		 * as set's over and under approximation respectevly.
		 * \param over is the over approximation of the set.
		 * \param under is the under approximation of the set.
		 */
		Set(ASet *over, ASet *under);
		
		/*! \brief This is a Set class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * Set.
		 * \param max_basic_sets is the maximum number of basic sets 
		 * usable to represents the region.
		 */
		Set(unsigned int max_basic_sets); 
			
		/*! \brief This is the destructor of the class Set.
		 *
		 * This destructor deletes in a safe way an object of the class 
		 * Set.
		 */
		~Set();
			
		/*! \brief Sets the number of s used in the region's
		 * reppresentation.
		 *
		 * This method sets the maximum number of basic sets used in 
		 * the representation of the set.
		 * \param max_basic_sets is the number of basic sets used in 
		 * the region's reppresentation.
 		 */
		void set_max_basic_sets(unsigned int max_basic_sets); 

		/*! \brief Gets the number of basic sets used in the region's
		 * reppresentation.
		 *
		 * This method gets the maximum number of basic sets used in 
		 * the representation of the region.
		 * \return is the number of basic sets used in the region's
		 * reppresentation.
 		 */
		unsigned int get_max_basic_sets() const; 

		/*! \brief Copies a \a Set on the current object.
		 *
		 * This method copies \a ts on the current object. 
		 * The set initially contained by the current 
		 * object is deleted.
		 * \param ts is the set which should be
		 * copied on the current object.
		 * \return A reference to the new object.
		 */
		Set& operator=(const Set& ts);
		
		/*! \brief Includes a basic set into the region.
		 *
		 * This method includes a basic set in both 
		 * over and under approximations of the region.		 
		 * \param set is the basic set that should be 
		 * included.
 		 */
		Set& include(BasicSet *set);
	  	
		/*! \brief Excludes a basic set from the region.
		 *
		 * This method excludes a basic set from both 
		 * over and under approximations of the region
		 * adding it to its \a ex_list list.
		 * \param set is the basic set that should be 
		 * excluded.
 		 */
		Set& exclude(BasicSet *set); 
		
		/*! \brief Checks if the intersection of two set is not null.
		 *
		 * This method checks if the intersection between the object
		 * and the set \a set is not null. If it is not null, the 
		 * method returns \a TRUE and \a FALSE otherwise.
		 * \param set is the set of which the intersection should be 
		 * tested.
		 * \param this_type is the current object's approximation 
		 * that should be used to test the interesection.
		 * \param set_type is the \a set 's approximation 
		 * that should be used to test the interesection.
		 * \return \a TRUE if the intersection between the 
		 * current set's approximation of type \a this_type 
		 * and the \a set 's approximation of type \a set_type 
		 * is not null, \a FALSE otherwise.
 		 */
		bool does_intersect(const Set &set, ApproxType this_type,
				ApproxType set_type) const; 
	
		/*! \brief Checks if the intersection of two set is not null.
		 *
		 * This method checks if the intersection between the 
		 * over-approximations of the current object and the set 
		 * \a set is not null. If it is not null, the 
		 * method returns \a TRUE and \a FALSE otherwise.
		 * \param set is the set of which the intersection should be 
		 * tested.
		 * \return \a TRUE if the intersection between the 
		 * over approximations of set and of 
		 * \a set is not null, \a FALSE otherwise.
 		 */
		bool does_intersect(const Set &set) const;

		/*! \brief Checks if the intersection of two set is not null.
		 *
		 * This method checks if the intersection between the object
		 * and the set \a set is not null. If it is not null, the 
		 * method returns \a TRUE and \a FALSE otherwise.
		 * \param set is the set of which the intersection should be 
		 * tested.
		 * \param atype is the current object's approximation 
		 * that should be used to test the interesection.
		 * \return \a TRUE if the intersection between the 
		 * current set's approximation of type \a atype 
		 * and \a set is not null, \a FALSE otherwise.
 		 */
		bool does_intersect(const ASet &set, ApproxType atype) const; 
		
		/*! \brief Checks if the intersection of two set is not null.
		 *
		 * This method checks if the intersection between the 
		 * over-approximations of the current object and the set 
		 * \a set is not null. If it is not null, the 
		 * method returns \a TRUE and \a FALSE otherwise.
		 * \param set is the set of which the intersection should be 
		 * tested.
		 * \return \a TRUE if the intersection between the 
		 * over approximations of set and of 
		 * \a set is not null, \a FALSE otherwise.
 		 */
		bool does_intersect(const ASet &set) const; 
	
		/*! \brief Tests whether the current set is a subset 
		 * of another.
		 *
		 * \param set is the set which should be a superset of 
		 * the current object.
		 * \param this_type is the current object's approximation 
		 * that should be used to test the inclusion.
		 * \param set_type is the \a set 's approximation 
		 * that should be used to test the inclusion.
		 * \return \a TRUE if the 
		 * current set's approximation of type \a this_type 
		 * is a subset of the \a set 's approximation of 
		 * type \a set_type , \a FALSE otherwise.
		 */
		bool is_subset_of(const Set& set, ApproxType this_type,
				ApproxType set_type) const;
		
		/*! \brief Tests whether the set is a subset of other.
		 *
		 * \param set is the set which should be a superset of 
		 * the current object.
		 * \return \a TRUE if the 
		 * current set's over-approximation is a subset of 
		 * the \a set 's over-approximation, \a FALSE otherwise.
		 */
		bool is_subset_of(const Set& set) const;
	
		/*! \brief Tests whether the current set is a subset 
		 * of another.
		 *
		 * \param set is the set which should be a superset of 
		 * the current object.
		 * \param atype is the current object's approximation 
		 * that should be used to test the inclusion.
		 * \return \a TRUE if the current set's approximation 
		 * of type \a this_type is a subset of \a set , 
		 * \a FALSE otherwise.
		 */
		bool is_subset_of(const ASet& set, ApproxType atype) const;
		
		/*! \brief Tests whether the set is a subset of other.
		 *
		 * \param set is the set which should be a superset of 
		 * the current object.
		 * \return \a TRUE if the 
		 * current set's over-approximation is a subset of 
		 * \a set , \a FALSE otherwise.
		 */
		bool is_subset_of(const ASet& set) const;
		
		/*! \brief Tests whether the set intersects other. 
		 *
		 * \param set is the set of which the disjointness from
		 * the current object should be tested.
		 * \param this_type is the current object's approximation 
		 * that should be used to test the disjointness.
		 * \param set_type is the \a set 's approximation 
		 * that should be used to test the disjointness.
 		 * \return \a TRUE if the current set's approximation 
		 * of type \a this_type is disjoint from the \a set 's 
		 * approximation of type \a set_type , \a FALSE otherwise.
		 */
		bool is_disjoint_from(const Set& set, ApproxType this_type,
				ApproxType set_type) const;

		/*! \brief Tests whether the set intersects other. 
		 *
		 * \param set is the set of which the disjointness from
		 * the current object should be tested.
		 * \return \a TRUE if the current set's over-approximation 
		 * is disjoint from the \a set 's over-approximation, 
		 * \a FALSE otherwise.
		 */
		bool is_disjoint_from(const Set& set) const;
	
		/*! \brief Tests whether the set contains the specified state.
		 *
		 * \param s is the state which should be cointaned into 
		 * the current object.
		 * \param type is the current object's approximation 
		 * that should be used to test the inclusion of \a s .
		 * \return \a TRUE if the current set's approximation of
		 * type \a type contains \a s , \a FALSE otherwise*/
		bool contains(const State &s, ApproxType type) const;

		/*! \brief Tests whether the set contains the specified state.
		 *
		 * \param s is the state which should be cointaned into 
		 * the current object.
		 * \return \a TRUE if the current set's over-approximation 
		 * contains \a s , \a FALSE otherwise*/
		bool contains(const State &s) const;
	
		/*! \brief Intersects the set with a basic set.
		 *
		 * This method intersects both over and under approximations 
		 * of the region with a basic set.
		 * \param set is the basic set that intersects the object.
		 * \return A pointer to the intersection of the current 
		 * object and the basic set \a set .
 		 */
		Set* intersect(const BasicSet &set);

		/*! \brief Intersects the set with a basic set.
		 *
		 * This method intersects both over and under approximations 
		 * of the region with a basic set.
		 * \param set is the basic set that intersects the object.
		 * \return The intersection of the current object and the
		 * basic set \a set .
 		 */
		Set operator&&(const BasicSet &set);
	
		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the object
		 * and the set \a set.		 
		 * \param set is the set of which the intersection should be
		 * computed.
		 * \return A pointer to the intersection between the object 
		 * and \a set.
 		 */
		Set* intersect(const ASet &set) const;

		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the object
		 * and the set \a set.		 
		 * \param set is the set of which the intersection should be
		 * computed.
		 * \return The object intersection between the object and 
		 * \a set.
 		 */
		Set operator&&(const ASet &set) const;
	
		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return A pointer to the union of the current set and 
		 * \a R.
		 */
		Set* join(const ASet &R) const;
		
		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return The union of the current set and \a R.
		 */
		Set operator||(const ASet &R) const;
		
		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the object
		 * and the set \a set.		 
		 * \param set is the set of which the intersection should be
		 * computed.
		 * \return A pointer to the intersection between the object 
		 * and \a set.
 		 */
		Set* intersect(const Set &set) const;

		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the object
		 * and the set \a set.		 
		 * \param set is the set of which the intersection should be
		 * computed.
		 * \return The object intersection between the object and 
		 * \a set.
 		 */
		Set operator&&(const Set &set) const;
	
		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return A pointer to the union of the current set and 
		 * \a R.
		 */
		Set* join(const Set &R) const;
		
		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return The union of the current set and \a R.
		 */
		Set operator||(const Set &R) const;

		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the object
		 * and the set \a set.		 
		 * \param set is the set of which the intersection should be
		 * computed.
		 * \return A pointer to the intersection between the object 
		 * and \a set.
 		 */
		HSet* intersect(const HSet &set) const;

		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the object
		 * and the set \a set.		 
		 * \param set is the set of which the intersection should be
		 * computed.
		 * \return The object intersection between the object and 
		 * \a set.
 		 */
		HSet operator&&(const HSet &set) const;
	
		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return A pointer to the union of the current set and 
		 * \a R.
		 */
		HSet* join(const HSet &R) const;
		
		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return The union of the current set and \a R.
		 */
		HSet operator||(const HSet &R) const;
		
		/*! \brief Applies the map \a M to the current set.
		 *
		 * This method applies the trasformation \a M to the set.
		 * The original object does not change.
		 * \param M is the trasformation to be applied.
		 * \return The trasformed set.
 		 */
		Set operator<<(const Map &M) const;

		/*! \brief Applies the map \a M to the current set.
		 *
		 * This method applies the trasformation \a M to the set.
		 * The object does not change.
		 * \param M is the trasformation to be applied.
		 * \return A pointer to a trasformed set.
 		 */
		Set* apply(const Map &M) const;

#ifdef FASTER_BUT_DIRTIER
		
		/*! \brief Applies the map \a M to the current set.
		 *
		 * This method applies the trasformation \a M to the set.
		 * The current object is set to the trasformed set.
		 * \param M is the trasformation to be applied.
		 * \return A reference to the trasformed set.
 		 */
		Set& apply_on_obj(const Map &M);
#endif
		
		/*! \brief Evaluates the flow slice from the current set. 
		 *
		 * This method evaluates the \a delta-timed flow slice from the 
		 * current set dued to the vector field \a field.
		 * \param field is the vector field.
		 * \param delta is the time of the slice.
		 * \return The flow slice.
 		 */
		Set flow_slice(const VectorField &field, 
				const double delta) const;
		 
		/*! \brief Evaluates the flow from the current set. 
		 *
		 * This method evaluates the set reachable with 
		 * a flow of exactely \a delta-timed from the 
		 * current set dued to the vector field \a field.
		 * \param field is the vector field.
		 * \param delta is the time of the flow.
		 * \return The set reached after a \a delta-timed flow from 
		 * the current set.
 		 */
		Set flow(const VectorField &field,
				const double delta) const;

		/*! \brief Evaluates the flow slice from the current set. 
		 *
		 * This method evaluates the \a delta-timed flow slice from the 
		 * current set dued to the vector field \a field.
		 * \param field is the vector field.
		 * \param delta is the time of the slice.
		 * \return The flow slice.
 		 */
		Set* flow_slice_p(const VectorField &field, 
				const double delta) const;
		 
		/*! \brief Evaluates the flow from the current set. 
		 *
		 * This method evaluates the set reachable with 
		 * a flow of exactely \a delta-timed from the 
		 * current set dued to the vector field \a field.
		 * \param field is the vector field.
		 * \param delta is the time of the flow.
		 * \return The set reached after a \a delta-timed flow from 
		 * the current set.
 		 */
		Set* flow_p(const VectorField &field,
				const double delta) const;

#ifdef FASTER_BUT_DIRTIER
		/*! \brief Evaluates the flow slice from the current set. 
		 *
		 * This method evaluates the \a delta-timed flow slice from the 
		 * current set dued to the vector field \a field.
		 * \param field is the vector field.
		 * \param delta is the time of the slice.
		 * \return The flow slice.
 		 */
		Set& flow_slice_on_obj(const VectorField &field, 
				const double delta);
		 
		/*! \brief Evaluates the flow from the current set. 
		 *
		 * This method evaluates the set reachable with 
		 * a flow of exactely \a delta-timed from the 
		 * current set dued to the vector field \a field.
		 * \param field is the vector field.
		 * \param delta is the time of the flow.
		 * \return The set reached after a \a delta-timed flow from 
		 * the current set.
 		 */
		Set& flow_on_obj(const VectorField &field,
				const double delta);
#endif

		/*! \brief Returns the over-approximation of the current 
		 * set. 
		 *
		 * This method returns the over-approximation of the current 
		 * set.
		 * \return The over-approximation of the set.
 		 */
		ASet* over_approximation() const;
		 
		/*! \brief Returns the under-approximation of the current 
		 * set. 
		 *
		 * This method returns the under-approximation of the current 
		 * set.
		 * \return The under-approximation of the set.
 		 */
		ASet* under_approximation() const;

		friend class HSet;
};

/*! \brief This is the class of set of hybrid space states.
 *
 * It associates a location to a set of phase space states.
 */ 
class HSet: public Set 
{
		/*! \brief The set's location. */
		Location *location;

	public:
		/*! \brief This is a HSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * HSet.
		 * \param ts is the template set for the new object.
		 */
		HSet(const HSet &ts);
		
		/*! \brief This is a HSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * HSet.
		 * \param l is the location of the set.
		 */
		HSet(Location *l);
		
		/*! \brief This is a HSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * HSet.
		 * \param set is the phase space set.
		 * \param l is the location of the set.
		 */
		HSet(const Set &set, Location *l);

		/*! \brief This is a HSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * HSet.
		 * \param over is the set's over-approximation.
		 * \param under is the set's under-approximation.
		 * \param l is the set's location.
		 */
		HSet(ASet *over, ASet *under, Location *l);
			
		/*! \brief This is a HSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * HSet.
		 * \param l is the location of the set.
		 * \param max_basic_sets is the maximum number of basic_sets 
		 * usable to represents the continuous part of the set.
		 */
		HSet(int max_basic_sets, Location *l);

		/*! \brief This is the destructor of the class HSet.
		 *
		 * This destructor deletes in a safe way an object of the class 
		 * HSet.
		 */
		~HSet();
	
		/*! \brief Copies a \a HSet on the current object.
		 *
		 * This method copies \a ts on the current object. 
		 * The set initially contained by the current 
		 * object is deleted.
		 * \param ts is the set which should be
		 * copied on the current object.
		 * \return A reference to the new object.
		 */
		HSet& operator=(const HSet& ts);
	
		/*! \brief Checks if two sets have the same location.
		 *
		 * This method checks if the object set and \a hset have the 
		 * same location.
		 * \param hset is the set of which should be checked if is has
		 * the same location of the object.
		 * \return \a TRUE if the object set and \a hset have the same 
		 * location, \a FALSE otherwise.
		 */
		bool same_location(const HSet &hset) const; 
	
		/*! \brief Intersects the set with a basic set.
		 *
		 * This method intersects both over and under approximations 
		 * of the region with a basic set.
		 * \param set is the basic set that intersects the object.
		 * \return A pointer to the intersection of the current 
		 * object and the basic set \a set .
 		 */
		HSet* intersect(const BasicSet &set);

		/*! \brief Intersects the set with a basic set.
		 *
		 * This method intersects both over and under approximations 
		 * of the region with a basic set.
		 * \param set is the basic set that intersects the object.
		 * \return The intersection of the current object and the
		 * basic set \a set .
 		 */
		HSet operator&&(const BasicSet &set);	
	
		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the object
		 * and the set \a set.		 
		 * \param set is the set of which the intersection should be
		 * computed.
		 * \return A pointer to the intersection between the object 
		 * and \a set.
 		 */
		HSet* intersect(const ASet &set) const;

		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the object
		 * and the set \a set.		 
		 * \param set is the set of which the intersection should be
		 * computed.
		 * \return The object intersection between the object and 
		 * \a set.
 		 */
		HSet operator&&(const ASet &set) const;
	
		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return A pointer to the union of the current set and 
		 * \a R.
		 */
		HSet* join(const ASet &R) const;
		
		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param R is the set of which the union should to be 
		 * computed.
		 * \return The union of the current set and \a R.
		 */
		HSet operator||(const ASet &R) const;
		
		/*! \brief Evaluates the intersection of two sets.
		 *
		 * This method computes the intersection between the current 
		 * set and the set \a hset.		 
		 * \param set is the set of which the intersection should be
		 * computed.
		 * \return The intersection between the set and \a hset.
 		 */
		HSet* intersect(const Set &set) const;
		
		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the object
		 * and the set \a set.		 
		 * \param set is the set of which the intersection should be
		 * computed.
		 * \return The object intersection between the object and 
		 * \a set.
 		 */
		HSet operator&&(const Set &set) const;
	
		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param set is the set of which the union with the
		 * current object should to be computed.
		 * \return A pointer to the union of the current set and 
		 * \a R.
		 */
		HSet* join(const Set &set) const;
		
		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param set is the set of which the union with the
		 * current object should to be computed.
		 * \return The union of the current set and \a R.
		 */
		HSet operator||(const Set &set) const;

		/*! \brief Evaluates the intersection of two sets.
		 *
		 * This method computes the intersection between the current 
		 * set and the set \a hset.		 
		 * \param set is the set of which the intersection should be
		 * computed.
		 * \return The intersection between the set and \a hset.
 		 */
		HSet* intersect(const HSet &set) const;
		
		/*! \brief Evaluates the intersection of two set.
		 *
		 * This method computes the intersection between the object
		 * and the set \a set.		 
		 * \param set is the set of which the intersection should be
		 * computed.
		 * \return The object intersection between the object and 
		 * \a set.
 		 */
		HSet operator&&(const HSet &set) const;
	
		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param set is the set of which the union with the
		 * current object should to be computed.
		 * \return A pointer to the union of the current set and 
		 * \a R.
		 */
		HSet* join(const HSet &set) const;
		
		/*! \brief Evaluates the union of two sets. 
		 *
		 * This method computes the union of the current 
		 * set and \a R.		 
		 * \param set is the set of which the union with the
		 * current object should to be computed.
		 * \return The union of the current set and \a R.
		 */
		HSet operator||(const HSet &set) const;
		
		/*! \brief Evaluates the flow slice from the set. 
		 *
		 * This method evaluates the \a delta-timed flow slice
		 * from the object dued to the vector field of the 
		 * current location.
		 * \param delta is the time of the slice.
		 * \return The flow slice.
 		 */
		HSet* flow_slice_p(const double delta) const;
		 
		/*! \brief Evaluates the flow from the object. 
		 *
		 * This method evaluates the set reachable with 
		 * a flow of exactely \a delta-timed from the 
		 * object dued to the vector field of the current 
		 * location.
		 * \param delta is the time of the flow.
		 * \return The set reached after a \a delta-timed
		 * flow from the object.
 		 */
		HSet* flow_p(const double delta) const;

		/*! \brief Evaluates the flow slice from the set. 
		 *
		 * This method evaluates the \a delta-timed flow 
		 * slice from the object dued to the vector field 
		 * of the current location.
		 * \param delta is the time of the slice.
		 * \return The flow slice.
 		 */
		HSet flow_slice(const double delta) const;
		 
		/*! \brief Evaluates the flow from the object. 
		 *
		 * This method evaluates the set reachable with 
		 * a flow of exactely \a delta-timed from the 
		 * object dued to the vector field of the current
		 * location.
		 * \param delta is the time of the flow.
		 * \return The set reached after a \a delta-timed
		 * flow from the object.
 		 */
		HSet flow(const double delta) const;
		
#ifdef FASTER_BUT_DIRTIER
		/*! \brief Evaluates the flow slice from the set. 
		 *
		 * This method evaluates the \a delta-timed flow 
		 * slice from the object dued to the vector field 
		 * of the current location.
		 * \param delta is the time of the slice.
		 * \return The flow slice.
 		 */
		HSet& flow_slice_on_obj(const double delta);
		 
		/*! \brief Evaluates the flow from the object. 
		 *
		 * This method evaluates the set reachable with 
		 * a flow of exactely \a delta-timed from the 
		 * object dued to the vector field of the current 
		 * location.
		 * \param delta is the time of the flow.
		 * \return The set reached after a \a delta-timed 
		 * flow from the object.
 		 */
		HSet& flow_on_obj(const double delta);
#endif
		
		/*! \brief Applies the map \a M to the current set.
		 *
		 * This method applies the trasformation \a M to the set.
		 * The original object does not change.
		 * \param M is the trasformation to be applied.
		 * \return The trasformed set.
 		 */
		HSet operator<<(const Map &M) const;
		
		/*! \brief Applies the map \a M to the current set.
		 *
		 * This method applies the trasformation \a M to the set.
		 * The object does not change.
		 * \param M is the trasformation to be applied.
		 * \return A pointer to a trasformed set.
 		 */
		HSet* apply(const Map &M) const;

#ifdef FASTER_BUT_DIRTIER
		
		/*! \brief Applies the map \a M to the current set.
		 *
		 * This method applies the trasformation \a M to the set.
		 * The current object is set to the trasformed set.
		 * \param M is the trasformation to be applied.
		 * \return A reference to the trasformed set.
 		 */
		HSet& apply_on_obj(const Map &M);
#endif
		
		/*! \brief Returns a pointer to the set's location.
		 *
		 * This method returns a pointer to the set's location.
		 * \return The set's location.
		 */
		const Location* get_location() const; 
		 
		/*! \brief Returns the corresponding phase space set.
		 *
		 * This method returns the projection the object 
		 * in the phase space.
		 * \return The projection the object 
		 * in the phase space.
		 */
		Set* get_Set() const;
		 
		/*! \brief Sets the set's location.
		 *
		 * This method sets the set's location.
		 * \param l is the new location of the set.
		 * \return A reference to the new object.
		 */
		HSet& set_location(Location *l); 

		/*! \brief Checks whenever a set activates
		 * an automaton's arc.
		 *
		 * This methods checks if the set, which 
		 * if represented by the current object,
		 * activates the leaving arc \a e or not.
		 * \param e is the leaving arc which may
		 * be activated by the current set.
		 * \return \a TRUE , if the current set 
		 * activates the leaving arc \a e , 
		 * \a FALSE otherwise.
		 */
		bool does_activate(const LeavingArc &e) const;

		/*! \brief Checks whenever a set activates
		 * an automaton's arc.
		 *
		 * This methods checks if the set, which 
		 * if represented by the current object,
		 * activates the arc \a e or not.
		 * \param e is the arc which may
		 * be activated by the current set.
		 * \return \a TRUE , if the current set 
		 * activates the arc \a e , \a FALSE 
		 * otherwise.
		 */
		bool does_activate(const Arc &e) const;

		/*! \brief Makes a set jump over a 
		 * leaving arc.
		 *
		 * This methods intersects the current set
		 * with the activation region and apply 
		 * the reset of the leaving arc \a e.
		 * The current object does not change.
		 * \param e is the leaving arc that 
		 * should be crossed.
		 * \return A pointer to the region 
		 * reached crossing the leaving arc \e ,
		 * if \a e is activated, a pointer a copy 
		 * of the current set otherwise.
		 */
		HSet* jump_p(const LeavingArc &e) const;
		
		/*! \brief Makes a set jump over an arc.
		 *
		 * This methods intersects the current 
		 * set with the activation region and 
		 * apply the reset of the arc \a e.
		 * The current object does not change.
		 * \param e is the arc that should be 
		 * crossed.
		 * \return A pointer to the region 
		 * reached crossing the arc \e , if \a e 
		 * is activated, a pointer a copy 
		 * of the current set otherwise.
		 */
		HSet* jump_p(const Arc &e) const;

		/*! \brief Makes a set jump over a 
		 * leaving arc.
		 *
		 * This methods intersects the current set
		 * with the activation region and apply 
		 * the reset of the leaving arc \a e.
		 * The current object does not change.
		 * \param e is the leaving arc that 
		 * should be crossed.
		 * \return The region reached crossing 
		 * the leaving arc \e , if \a e is 
		 * activated, a pointer a copy 
		 * of the current set otherwise.
		 */
		HSet jump(const LeavingArc &e) const;
		
		/*! \brief Makes a set jump over an arc.
		 *
		 * This methods intersects the current 
		 * set with the activation region and 
		 * apply the reset of the arc \a e.
		 * The current object does not change.
		 * \param e is the arc that should be 
		 * crossed.
		 * \return The region reached crossing 
		 * the arc \e , if \a e is activated, 
		 * a pointer a copy of the current 
		 * set otherwise.
		 */
		HSet jump(const Arc &e) const;
		
};

class RSet
{
	
};

}

#endif /* _SET_H */
