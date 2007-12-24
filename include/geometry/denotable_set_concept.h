/***************************************************************************
 *            denotable_set_concept.h
 *
 *  Copyright 2007  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file denotable_set_concept.h
 *  \brief Documentation for the BasicSet concept.
 */


namespace Ariadne {
  namespace Geometry {


    /*! \ingroup DenotableSet
     *  \brief The specification of the DenotableSet concept.
     *
     *  A denotable set is any set which can be exactly expressed as a finite
     *  union of sets of some type satisfying the BasicSetConcept.
     *
     *  See also the SetInterface interface.
     */
    class DenotableSetConcept 
    {
     public:
      /*! \brief The type of real number used to describle the set. */
      typedef RealConcept real_type;
      /*! \brief The type of point that the set contains. */
      typedef PointConcept state_type;
      /*! \brief The type of basic set that makes up the denotable set. */
      typedef BasicSetConcept basic_set_type;
      /*! \brief A modifying iterator through the basic set elements of the denoatable set. (Optional) */
      typedef ForwardIterator iterator;
      /*! \brief A modifying iterator through the basic set elements of the denoatable set. (Optional) */
      typedef ConstForwardIterator const_iterator;
     public:
      //@{ 
      //! \name Constructors and assignment operators
      /*! \brief Copy constructor. */
      DenotableSetConcept(const DenotableSetConcept&);
    
      /*! \brief Copy assignment operator. */
      DenotableSetConcept& operator=(const DenotableSetConcept&);

      //@}
      
      
      //@{
      //! \name Required geometric operations
      /*! \brief The dimension of the Euclidean space the set lies in. */
      size_type dimension() const;
      
      /*! \brief Tests if a point is an element of the set. */
      tribool contains(const state_type& pt) const;
      
      /*! \brief A rectangle containing the denotable set. */
      Box bounding_box() const;
      //@}

      //@{
      //! \name List operations
      /*! \brief The number of basic sets comprising the set. */
      size_type size() const;
      
      /*! \brief A constant forward iterator to the first basic set in the list. */
      const_iterator begin() const;
      
      /*! \brief A constant forward iterator to the last basic set in the list. */
      const_iterator end() const;
      
      /*! \brief A constant reference to the \a ith basic set in the list. */
      const basic_set_type& operator[](size_type i) const;
      
      /*! \brief Adjoin a basic set to the list. */
      void adjoin(const basic_set_type& bs);
      /*! \brief Adjoin all elements of a denotable set to the list. */
      void adjoin(const DenotableSetConcept& ds);
      //@}

      //@{
      //! \name Optional geometric operations
      /*! \brief Tests if the set is empty. (Optional) */
      tribool empty() const;
      
      /*! \brief Tests if the set is bounded. (Optional) */
      tribool bounded() const;
      
      /*! \brief An approximation to the volume. (Optional) */
      real_type volume() const;

      //@}
     
#ifdef DOXYGEN
      //@{ 
      //! \name Required geometric predicates and operations
      /*! \brief Tests if a denotable set contains a point. */
      friend tribool contains(const DenotableSetConcept& ds, const state_type>& pt);
      /*! \brief Tests if a denotable set is disjoint from a rectangle. */
      friend tribool disjoint(const DenotableSetConcept& ds, const Box<real_type>& r);
      /*! \brief Tests if a denotable set is a subset of a rectangle. */
      friend tribool subset(const DenotableSetConcept& ds, const Box<real_type>& r);
      /*! \brief The union of a denotable set and a basic set. (Note that union is a reserved word in C++) */
      friend DenotableSetConcept join(const DenotableSetConcept& ds, const basic_set_type& r);
      /*! \brief The union of two denotable sets. (Note that union is a reserved word in C++) */
      friend DenotableSetConcept join(const DenotableSetConcept& ds, const DenotableSetConcept& r);
      //@}

      //@{ 
      //! \name Optional binary geometric predicates and operations
      /*! \brief Tests if a rectangle is a subset of the basic set. (Optional) */
      friend tribool subset(const Box<real_type>& s, const DenotableSetConcept& ds);
      /*! \brief The intersection of \a ds1 and \a ds2. (Optional) */
      friend DenotableSetConcept intersection(const DenotableSetConcept& ds1, const DenotableSetConcept& ds2); 
      /*! \brief The closure of the intersection of the interiors of \a ds1 and \a ds2. (Optional) */
      friend DenotableSetConcept regular_intersection(const DenotableSetConcept& ds1, const DenotableSetConcept& ds2); 

      //@}
#endif
      
      //@{ 
      //! \name Input/output operations
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. (Optional) */
      std::istream& read(std::istream& is);
      //@}
    };
  }
}
