/***************************************************************************
 *            basic_set.h
 *
 *  Thu Apr 15 15:18:51 2004
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
 
#ifndef _BASIC_SET_H
#define _BASIC_SET_H

#include <string>

#include "approx_type.h"
#include "io_error.h"
#include "state.h"
#include "classes.h"
	
namespace Ariadne {	

typedef unsigned int BasicSetID;
	
/*! \brief The simplest structure representing a set of phase space states.
 *
 * It's a virtual class for convex closed geometrical shapes. A shape \a S 
 * represents the set of all that phase space states \a x 
 * such that \a x stays in \a S. The class methods supply basic operations
 * on geometrical shapes such intersection and union and, moreover, should
 * provide methods to parse a shape definition from a buffer, to 
 * insert a basic set into a basic maintainer and to check if
 * a basic set is fully included into a basic maintainer.
 */
class BasicSet
{
	protected:
		/*! \brief The basic set id */
		BasicSetID id;

	public:
		/*! \brief This is a \a BasicSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a BasicSet to an empty basic set.
		 */	
		BasicSet() { this->id=0; }

		/*! \brief This is a \a BasicSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a BasicSet to an empty basic set.
		 * \param dim is the space's dimension of the new basic
		 * set.
		 */	
		//BasicSet(unsigned int dim);
		
		/*! \brief This is a \a BasicSet class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a BasicSet.
		 * \param t is a template for the new basic set.
		 */	
		//BasicSet(const BasicSet &t);

		/*! \brief This is the destructor of the class 
		 * \a BasicSet.
		 *
		 * This destructor deletes in a safe way an object of the class 
		 * \a BasicSet.
		 */
		//virtual ~BasicSet() = 0;
		
#ifdef FAST_BUT_DIRTIER		
		/*! \brief Evalutes the geometric sum of two basic sets.
		 *
		 * This method makes the geometric sum of two sets. If 
		 * the sum is not exactely representable by a 
		 * set, an over approximation is evaluated insteed. 
		 * The object is set to the geometric sum.
		 * \param set is the set that should be summed.
		 * \return A reference to the geometric sum of the sets.
 		 */
		virtual BasicSet& operator+= (const BasicSet &set) = 0; 
#endif

		/*! \brief Evalutes the geometric sum of two basic sets.
		 *
		 * This method makes the geometric sum of two sets. If 
		 * the sum is not exactely representable by a 
		 * set, an over approximation is evaluated insteed. 
		 * The object does not change.
		 * \param set is the set that should be summed.
		 * \return The geometric sum of the two sets.
 		 */
		virtual BasicSet* operator+ (const BasicSet &set) const = 0;

		/* TODO: comment
 		 */
		virtual BasicSet* expand_of(const double e) const = 0; 
		virtual BasicSet* shrink_of(const double e) const = 0;

		/*! \brief Evaluates the set difference of two sets.
		 *
		 * This method evaluates the set difference between 
		 * the current basic set and the basic set {\a set}.
		 * The type of approximation wanted is indicated by 
		 * {\a atype}. If \a atype=NONE and the difference
		 * is not exactely representable by a list of basic 
		 * sets, an over approximation is evaluated insteed. 
		 * The current object does not change.
		 * \param set is the basic set which should be subtracted
		 * from the current object.
		 * \param atype is the type of approximation wanted.
		 * \return A list of basic sets such that their union
		 * represent the difference of the two sets.
		 */
		virtual BasicSetList* subtract(const BasicSet &set,
				const ApproxType atype) const = 0;
		
		/*! \brief Evaluates the set difference of two sets.
		 *
		 * This method evaluates the set difference between 
		 * the current basic set and the basic set {\a set}.
		 * The type of approximation wanted is indicated by 
		 * {\a atype}. If \a atype=NONE and the difference
		 * is not exactely representable by a list of basic 
		 * sets, an over approximation is evaluated insteed. 
		 * The current object does not change.
		 * \param set is the basic set which should be subtracted
		 * from the current object.
		 * \param atype is the wanted type of approximation. 
		 * \param max_elem is the maximum number of basic sets
		 * in the output list.
		 * \return A list of basic sets such that their union
		 * represent the difference of the two sets.
		 */
		virtual BasicSetList* subtract(const BasicSet &set,
				const ApproxType atype,
				const unsigned int max_elem) const = 0;
		
		/*! \brief Makes the union of two sets.
		 *
		 * This method makes the union of two sets. If the union 
		 * is not exactely representable by a basic set, an 
		 * over approximation is evaluated insteed. 
		 * The object does not change.
		 * \param set is the basic set that should be united.
		 * \return The union of the two sets.
 		 */
		virtual BasicSet* operator||(const BasicSet &set) const = 0; 

		/*! \brief Makes the union of two sets.
		 *
		 * This method makes the union of two sets. If the union 
		 * is not exactely representable by a set, an 
		 * approximation is evaluated insteed. The type of the 
		 * incidental approximation is indicated by the second 
		 * parameter.
		 * The object does not change.
		 * \param set is the set that should be united to 
		 * the object.
		 * \param atype is the type of the approximation.
		 * \return The union of the two sets.
 		 */
		virtual BasicSet* join(const BasicSet &set, 
					const ApproxType atype) const= 0;
	
#ifdef FASTER_BUT_DIRTIER
		/*! \brief Makes the union of two sets.
		 *
		 * This method makes the union of two sets. If the union 
		 * is not exactely representable by a basic set, an 
		 * over approximation is evaluated insteed. 
		 * The object is set to the union.
		 * \param set is the basic set that should be united.
		 * \return A reference to the union of the sets.
 		 */
		virtual BasicSet& operator|=(const BasicSet &set) = 0; 

		/*! \brief Makes the union of two sets.
		 *
		 * This method makes the union of two sets. If the union 
		 * is not exactely representable by a set, an 
		 * approximation is evaluated insteed. The type of the 
		 * incidental approximation is indicated by the second 
		 * parameter.
		 * The object is set to the union of the sets.
		 * \param set is the set that should be united to 
		 * the object.
		 * \param atype is the type of the approximation.
 		 */
		virtual BasicSet& join_on_obj(const BasicSet &set, 
					const ApproxType atype) = 0; 
#endif		
		/*! \brief Create an object using as template the 
		 * current one.
		 */
		virtual BasicSet* copy() const = 0;

		/*! \brief Copy a basic set on the object
		 * 
		 * \param bs is the template for the current object.
		 * \return A reference to the copied object.
		 */
		virtual BasicSet& operator=(const BasicSet &bs) = 0;
				
		/*! \brief Evalutes the intersections of two sets.
		 *
		 * This method makes the intersections of two sets. If 
		 * the intersection is not exactely representable by a 
		 * set, an over approximation is evaluated insteed.
		 * The object does not change.
		 * \param set is the set that should be intersected.
		 * \return The intersection of the two sets.
 		 */
		virtual BasicSet* operator&& (const BasicSet &set) const = 0; 

		/*! \brief Evalutes the intersections of two sets.
		 *
		 * This method makes the intersections of two sets. If 
		 * the intersection is not exactely representable by a 
		 * set, an approximation is evaluated insteed. The type 
		 * of the incidental approximation is indicated by the second 
		 * parameter.
		 * The object does not change.
		 * \param set is the set that should be intersected 
		 * to the object.
		 * \param atype is the type of the approximation.
		 * \return The intersection of the two sets.
 		 */
		virtual BasicSet* intersect(const BasicSet &set, 
					const ApproxType atype) const = 0; 

#ifdef FASTER_BUT_DIRTIER
		/*! \brief Evalutes the intersections of two sets.
		 *
		 * This method makes the intersections of two sets. If 
		 * the intersection is not exactely representable by a 
		 * set, an over approximation is evaluated insteed. 
		 * The object is set to the intersection.
		 * \param set is the set that should be intersected.
		 * \return A reference to the intersection of the sets.
 		 */
		virtual BasicSet& operator&= (const BasicSet &set) = 0; 

		/*! \brief Evalutes the intersections of two sets.
		 *
		 * This method makes the intersections of two sets. If 
		 * the intersection is not exactely representable by a 
		 * set, an approximation is evaluated insteed. The type 
		 * of the incidental approximation is indicated by the second 
		 * parameter.
		 * The object is set to the intersection of the sets.
		 * \param set is the set that should be intersected 
		 * to the object.
		 * \param atype is the type of the approximation.
 		 */
		virtual BasicSet& intersect_on_obj(const BasicSet &set, 
					const ApproxType atype) = 0; 
#endif
	
		/*! \brief Checks if the intersection of two set is not 
		 * null.
		 *
		 * This method checks if the intersection between the current 
		 * set and the set S is not null. If it is not 
		 * null, the method returns true and false otherwise.
		 * \param set is the set of which the intersection 
		 * should be tested.
		 * \return true if the intersection between the set and 
		 * S is not null, false otherwise.
 		 */
		virtual bool does_intersect(const BasicSet &set) const = 0;
		
		/*! \brief Checks if the type of \a set is known.
		 *
		 * This method checks if the type of \a set is known. 
		 * The method returns true, if it is so, and false otherwise.
		 * This method will be used to test if intersection or union
		 * between the current object and the object \a set can
		 * be done.
		 * \param set is the set of which the knowlge of 
		 * type should be checked.
		 * \return \a true if the type of \a set is known, 
		 * \a false otherwise.
 		 */
		virtual bool known_basicset(const BasicSet &set) const = 0;
		 
		/*! \brief Checks if the set is an empty set.
		 *
		 * This method checks if the object represents an empty set.
		 * \return true if the set is an empty set, false 
		 * otherwise.
 		 */
		virtual bool empty() const = 0;
	
		/* \brief Tests whether the set is a subset of {\a other}.*/
		virtual bool is_subset_of(const BasicSet& other) const = 0;
		
		/* \brief Tests whether the set is disjoint from {\a other}. */
		virtual bool is_disjoint_from(const BasicSet& other) const = 0;
		
		/* \brief Tests whether the set contains the specified state.*/
		virtual bool contains(const State &s) const = 0;
		
		/*! \brief Evaluates the distance between two sets.
		 *
		 * This method computes the distance between two sets.
		 * \param set is the set from which the distance 
		 * should be evaluated.
		 * \return The distance from the object set and the 
		 * set \a cont.
 		 */
		virtual double distance(const BasicSet &set) const = 0;
		
		/*! \brief Applies the map \a M on the set.
		 *
		 * This method applies the map \a M on the set.
		 * If the trasformed region is not exactely representable 
		 * by a set, an over approximation is evaluated insteed. 
		 * The object does not change.
		 * \param M is the trasformation to be applied.
		 * \return The transformed basic set.
 		 */
		virtual BasicSet* apply(const Map &M) const = 0;
		
		/*! \brief Applies the map \a M on the set.
		 *
		 * This method applies the map \a M on the set.
		 * If the trasformed region is not exactely representable 
		 * by a set, an approximation is evaluated insteed. The 
		 * type of the incidental approximation is indicated by a 
		 * parameter.
		 * The object does not change.
		 * \param M is the trasformation to be applied.
		 * \param atype is the type of the approximation.
		 * \return The transformed basic set.
 		 */
		virtual BasicSet* apply(const Map &M, 
					const ApproxType atype) const = 0;
	
#ifdef FASTER_BUT_DIRTIER
		/*! \brief Applies the map \a M on the set.
		 *
		 * This method applies the map \a M on the set.
		 * If the trasformed region is not exactely representable 
		 * by a set, an over approximation is evaluated insteed. 
		 * The object is set to the transformed basic set.
		 * \param M is the trasformation to be applied.
		 * \return A reference to the transformed basic set.
 		 */
		virtual BasicSet& apply_on_obj(const Map &M) = 0;
		
		/*! \brief Applies the map \a M on the set.
		 *
		 * This method applies the map \a M on the set.
		 * If the trasformed region is not exactely representable 
		 * by a set, an approximation is evaluated insteed. The 
		 * type of the incidental approximation is indicated by a 
		 * parameter.
		 * The object is set to the transformed basic set.
		 * \param M is the trasformation to be applied.
		 * \param atype is the type of the approximation.
		 * \return A reference to the transformed basic set.
 		 */
		virtual BasicSet& apply_on_obj(const Map &M, 
					const ApproxType atype) = 0;
#endif
		
		/*! \brief Evaluates the flow slice from the current set. 
		 *
		 * This method evaluates the \f$\delta\f$-timed flow slice 
		 * from the current region dued to the vector field \a field. 
		 * If the flow slice is not exactely representable 
		 * by a set, an over approximation is evaluated insteed. 
		 * The object does not change.
		 * \param field is the vector field.
		 * \param delta is the time of the slice.
		 * \return The \f$\delta\f$-timed flow slice from the 
		 * current set.
 		 */
		virtual BasicSet* flow_slice(const VectorField &field, 
					const double delta) const = 0;
		
		/*! \brief Evaluates the flow slice from the current set. 
		 *
		 * This method evaluates the \f$\delta\f$-timed flow slice
		 * from the current region dued to the vector field \a field. 
		 * If the flow slice is not exactely representable 
		 * by a set, an over approximation is evaluated insteed. 
		 * The object does not change.
		 * \param field is the vector field.
		 * \param delta is the time of the slice.
		 * \param atype is the approximation type.
		 * \return The \f$\delta\f$-timed flow slice from the current 
		 * set.
 		 */
		virtual BasicSet* flow_slice(const VectorField &field, 
					const double delta, 
					const ApproxType atype) const = 0;
		
		/*! \brief Includes the object into a basic maintainer.
		 *
		 * This method inserts the object into a basic maintainer.
		 * \param maintainer is the basic maintainer into which the
		 * object is stored.
		 * \param atype is the approximation type that can be used
		 * in the insertion.
		 * \return The type of approximation used in the insertion.
		 */
		virtual ApproxType maintain(BasicMaintainer &maintainer,
				const ApproxType atype) const = 0;
		
		/*! \brief Checks if the object is fully included into the 
		 * maintain system.
		 *
		 * This method checks if the object is fully included into the 
		 * maintain system or not.
		 * \param maintainer is the basic maintainer in which the 
		 * inclusion should be tested.
		 * \return \a true if the object is fully included into the 
		 * maintain system, \a false otherwise.
		 */
		virtual bool fully_maintained(BasicMaintainer &maintainer) 
			const = 0;
	
		/*! \brief Gets an approximation of the unmaintained part
		 * of the object.
		 *
		 * This method gets an approximation of the unmaintained part
		 * of the object. The resulting basic set is stored 
		 * into the current object.
		 * \param maintainer is the basic maintainer into which the
		 * object could be stored.
		 * \param atype is the approximation type.
		 * \return The type of approximation used.
		 */
		virtual ApproxType get_unmaintained(
				BasicMaintainer &maintainer,
		 		ApproxType atype) = 0;

		/*! \brief Reads a set definition and stores it into the 
		 * object.
		 *
		 * This method parses a set definition maintained into 
		 * a buffer and stores the corresponding set into the 
		 * current object.
		 * The language used to describe classes of sets does not 
		 * depend on the automaton specification language and can be 
		 * freely choosen by set's developers.
		 * \see write_set
		 * \param buf is the buffer which stores the set 
		 * definition.
		 * \return An error code. 		 
		 */
		virtual Read_Error read_set(Buffer &buf) = 0;
		
		/*! \brief Writes the set definition of the 
		 * object.
		 *
		 * This method writes the set definition corresponding 
		 * to the current object into a buffer.
		 * The language used to describe classes of sets does not 
		 * depend on the automaton specification language and can be 
		 * freely choosen by set's developers.
		 * \see read_set
		 * \param buf is the buffer in which the set 
		 * definition should be stored.
		 * \return An error code. 		 
		 */
		virtual Write_Error write_set(Buffer &buf) const = 0;
		
		/*! \brief Returns the name of the basic set type.
		 *
		 * This methods returns the name of the basic set type and 
		 * it will be used by the tool parser to identify the 
		 * specification of an object of the current class.
		 * \return The name of the basic set type.
		 */
		virtual std::string get_name() = 0;

		/*! \brief Returns the id of the basic set type.
		 *
		 * This methods returns the id of the basic set type and 
		 * it will be used by the tool parser to identify the 
		 * specification of an object of the current class.
		 * \return The name of the basic set type.
		 */
		BasicSetID get_id() { return this->id; }
};

}


#endif /* _BASIC_SET_H */
