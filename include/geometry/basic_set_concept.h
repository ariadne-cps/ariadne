/***************************************************************************
 *            basic_set_concept.h
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
 
/*! \file basic_set_concept.h
 *  \brief Documentation for the BasicSet concept
 */


namespace Ariadne {
  namespace Geometry {


    /*! \ingroup BasicSet
     *  \brief The specification of the BasicSet concept.
     *
     *  Geometric predicates return tribool values, with the \a indeterminate
     *  value indicating that the result could not be determined by computing at
     *  the given precision, or is not robust with respect to perturbation of
     *  the parameters.
     *
     *  See also the DenotableSetConcept and the SetInterface interface.
     */
    class BasicSetConcept 
    {
     public:
      /*! \brief The type of real number used to describle the set. */
      typedef RealConcept real_type;
      /*! \brief The type of point that the set contains. */
      typedef PointConcept state_type;
     public:
      //@{
      //! \name Constructors and assignment operators
    
      /*! \brief Copy constructor. */
      BasicSetConcept(const BasicSetConcept&);
    
      /*! \brief Copy assignment operator. */
      BasicSetConcept& operator=(const BasicSetConcept&);

      //@}
      
      
      //@{
      //! \name Required geometric operations
      /*! \brief The dimension of the Euclidean space the set lies in. */
      size_type dimension() const;
      
      /*! \brief Tests if a point is an element of the set. */
      tribool contains(const state_type& pt) const;
      
      /*! \brief A rectangle containing the given rectangle; returns a copy. */
      Box bounding_box() const;
      //@}

      //@{
      //! \name Optional geometric operations
      /*! \brief Tests if the set is empty. May return \a indeterminate if emptiness is not sufficiently robust. (Optional) */
      tribool empty() const;
      
      /*! \brief Tests if the set is bounded. (Optional) */
      tribool bounded() const;
      
      /*! \brief A point in the set, typically the centre. (Optional) */
      Point<real_type> centre() const;
      
      /*! \brief The radius in the supremum norm. (Optional) */
      real_type radius() const;
      
      /*! \brief An approximation to the volume. (Optional) */
      real_type volume() const;

      //@{
      //! \name Polyhedral operations
      /*! \brief The number of vertices. (Polyhedral sets only; Optional) */
      size_type number_of_vertices() const;
      /*! \brief The \a i th vertex. (Polyhedral sets only; Optional) */
      Point<real_type> vertex(size_type i) const;
      /*! \brief A constant iterator to the first vertex of the set. (Polyhedral sets only; Optional) */
      vertices_const_iterator vertices_begin() const;
      /*! \brief A constant iterator to the end vertex of the rectangle. (Polyhedral sets only; Optional) */
      vertices_const_iterator vertices_end() const;
      //@}
      //@}
     
#ifdef DOXYGEN
      //@{ 
      //! \name Required geometric predicates and operations
      /*! \brief Tests if the set contains a point. */
      friend tribool contains(const BasicSetConcept& bs, const state_type>& pt);
      /*! \brief Tests disjointness with a rectangle. */
      friend tribool disjoint(const BasicSetConcept& bs, const Box<real_type>& r);
      /*! \brief Tests if the basic set is a subset of a rectangle. */
      friend tribool subset(const BasicSetConcept& bs, const Box<real_type>& r);
      /*! \brief Tests if a rectangle is a subset of the basic set. */
      friend tribool subset(const Box<real_type>& s, const BasicSetConcept& bs);
      //@}

      //@{ 
      //! \name Optional binary geometric predicates and operations
      /*! \brief %Set equality operator. (Optional) */
      friend tribool equal(const BasicSetConcept& bs1, const BasicSetConcept& bs2);

      /*! \brief The intersection of \a bs1 and \a bs2. (Optional) */
      friend BasicSetConcept intersection(const BasicSetConcept& bs1, const BasicSetConcept& bs2); 
      /*! \brief The closure of the intersection of the interiors of \a bs1 and \a bs2. (Optional) */
      friend BasicSetConcept regular_intersection(const BasicSetConcept& bs1, const BasicSetConcept& bs2); 

      /*! \brief The smallest set containing \a bs1 and \a bs2. (Optional) */
      friend BasicSetConcept hull(const BasicSetConcept& bs2, const BasicSetConcept& bs2); 

      /*! \brief The componentwise sum of basic sets \a bs1 and \a bs2. (Optional) */
      friend BasicSetConcept minkowski_sum(const BasicSetConcept& bs2, const BasicSetConcept& bs2); 
      /*! \brief The componentwise difference of basic sets \a bs1 and \a bs2. (Optional) */
      friend BasicSetConcept minkowski_difference(const BasicSetConcept& bs2, const BasicSetConcept& bs2); 
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
