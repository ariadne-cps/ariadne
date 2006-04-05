/***************************************************************************
 *            approximate_set.h
 *
 *  Thu Oct  2 10:32:15 2004
 *  Copyright  2004  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#ifndef _APPROXIMATE_SET_H
#define _APPROXIMATE_SET_H

#include "../geometry/approx_type.h"
#include "../geometry/denotable_set.h"
#include "../geometry/point.h"
#include "../geometry/classes.h"

namespace Ariadne {

namespace Geometry {

/*! \class ApproximateSet
 * \brief Represents an approximating denotable set. 
 */
class ApproximateSet {
  /* FOR PIETER: This class should be inherited from some 
   * kind of denotable set. 
   * But, in that case, we should 
   * define as many approximate set classes as denotable 
   * set (e.g. ApprosimateSetList<T>, ApproximateSetGrid 
   * and so on). What do you prefer? */

  /*! \brief The denotable set used for the approximation. */
  AbstractDenotableSet *_set;

  /*! \brief The approximation type. */
  ApproximationType _type;
 public:

  /*! \brief This is a \a ApproximateSet class constructor.
   *
   * This constructor initializes the object of the class 
   * \a ApproximateSet.
   * \param set is the denotable set used for the approximation.
   * \param type is the approximation type.
   */
  ApproximateSet(AbstractDenotableSet *set, ApproximationType type);

  /*! \brief This is a \a ApproximateSet class constructor.
   *
   * This constructor initializes the object of the class 
   * \a ApproximateSet.
   * \param set is the set which is used as model for the
   * new object.
   */
  ApproximateSet(const ApproximateSet& set);
  
  /*! \brief This is the destructor of the class 
   * \a ApproximateSet.
   *
   * This destructor deletes in a safe way an object 
   * of the class \a ApproximateSet.
   */
  ~ApproximateSet();

  /*! \brief Makes the union of two approximate sets.
   *
   * This method tries to make the union of two approximate sets.
   * If the set's types do not agree, this methods throws
   * a \a NotCompatibleApproximateSets exception.
   * \param set the set .
   * \return A new object representing the
   * union of the two sets.
   */
  ApproximateSet* join(const ApproximateSet& set) const;

  /*! \brief Makes the intersection of two denotable sets.
   *
   * This method intersects the sets represented
   * by the current object and \a set returning 
   * a new object.
   * If the set's types do not agree, this methods throws
   * a \a NotCompatibleApproximateSets exception.
   * \param set is the approximate set which should be 
   * intersected with the current object.
   * \return A \a ApproximateSet containing the
   * intersection of the current set with \a set.
   */
  ApproximateSet* intersect(const ApproximateSet& set) const;

  /*! \brief Evaluates the set difference of two sets.
   *
   * This method evaluates the set difference between 
   * the current set and the set \a set.
   * \param set is the set which should be subtracted
   * from the current object.
   * \return A pointer to a \a AbstractDenotableSet 
   * which represents the difference of the two sets.
   */
  AbstractDenotableSet* subtract(const AbstractDenotableSet &set) const;

  /*! \brief Copies an object on the current one.
   *
   * This method copies an object on the current 
   * one. The current object is deleted.
   * \param set is the set which should be
   * copied on the current object.
   * \return A reference to the new set.
   */
  ApproximateSet& operator=(const ApproximateSet& set);

  /*! \brief Checks if a denotable set includes a state.
   *
   * This method checks whenever the current denotable set
   * includes the state \a s.
   * \param s is the state of which inclusion 
   * into the current denotable set should be tested.
   * \return  \a true, if \a s is contained into the 
   * current set, \a false otherwise.
   */
  bool contains(const State &s) const;

  /*! \brief Checks if the intersection of two set is not null.
   *
   * This method checks if the intersection between the current 
   * set and the set \a set is not null. If it is not null, the 
   * method returns \a true and \a false otherwise.
   * \param set is the set of which the intersection is to be 
   * tested.
   * \return \a true if the intersection between the set and 
   * \a set is not null, \a false otherwise.
   */
  bool does_intersect(const ApproximateSet &set) const;

  /*! \brief Tests whether a set is a subset of an other.
   *
   * \return \a true if the current set is s subset of \a set,
   * \a false otherwise.*/
  bool is_subset_of(const ApproximateSet &set) const;
  
  /*! \brief Checks whenever a set is empty.
   *
   * \return \a true if the current set is empty,
   * \a false otherwise.
   */
  bool empty() const;

  /*! \brief Create an object using as original the 
   * current one.
   * \return A pointer to a copy of the current set.
   */
  ApproximateSet* copy() const;

};

}
}

#endif /* _APPROXIMATE_SET_H */
