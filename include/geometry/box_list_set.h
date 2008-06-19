/***************************************************************************
 *            box_list_set.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file list_set.h
 *  \brief Unions of boxes implemented as lists.
 */

#ifndef ARIADNE_BOX_LIST_SET_H
#define ARIADNE_BOX_LIST_SET_H

#include <iosfwd>
#include <iostream>
#include <exception>
#include <stdexcept>

#include <vector>

#include "linear_algebra/declarations.h"
#include "base/utility.h"
#include "base/tribool.h"
#include "geometry/box.h"
#include "geometry/set_interface.h"

namespace Ariadne {
  

    template<class R> class Box;
    
    class basic_set_tag;
    class denotable_set_tag;
  
    /*!\ingroup DenotableSet
     * \ingroup List
     * \brief A finite union of boxed, represented as a sequence.
     */
    template<class R>
    class BoxListSet 
      : public SetInterface< Box<R> >
    {
     private:
      /* List of basic sets. Note that std::vector provides a
       * reserve(size_type) method to increase the capacity.
       */
      std::vector< Box<R> > _vector;
     public:
      /*! \brief A tag describing the type of set. */
      typedef denotable_set_tag set_category;
      /*!\brief The type of denotable real number used to represent points in the space. */
      typedef R real_type;
      /*!\brief The type of point contained by the set. */
      typedef Point<R> state_type;
      /*!\brief The type of basic set in the list of sets. */
      typedef Box<R> basic_set_type;
      /*!\brief The type of basic set in the list of sets. */
      typedef Box<R> value_type;

      typedef typename std::vector< Box<R> >::const_iterator const_iterator;
      typedef typename std::vector< Box<R> >::iterator iterator;

     public:
      /*! \brief An empty list set which can hold sets of an unspecified dimension. */
      BoxListSet();

      /*! \brief An empty list set which can hold sets of a dimension \a d. */
      BoxListSet(dimension_type d);

      /*! \brief Returns the number of basic sets forming this object. */
      size_type size() const;

      /*! \brief Adjoins a basic set to the back of the list. */
      void push_back(const Box<R>& A);

      /*! \brief Removes the box at the back of the list. */
      void pop_back();

      /*! \brief Accesses the i-th box in the list. */
      const Box<R>& get(size_type index) const;

      /*! \brief Assigns to the i-th box. */
      void set(size_type index, const Box<R>& set);

      /*! \brief Accesses the i-th box. */
      const Box<R>& operator[](size_type index) const;

      /*! \brief Returns the denotable set's space dimension. */
      dimension_type dimension() const;



 
      //@{ 
      //! \name Standard geometry operators
      /*! \brief Returns the denotable set's space dimension. */
      virtual BoxListSet<R>* clone() const;

      /*! \brief Returns the denotable set's space. */
      virtual EuclideanSpace space() const;

      /*!\brief Checks if a denotable set includes a point. */
      virtual tribool contains(const Point<R>& p) const;

      /*! \brief Tests for superset of a box. 
       *
       * Currently always returns \a indeterminate, since the 
       * test is difficult for general list sets. 
       */
      virtual tribool superset(const Box<R>& r) const;

      /*! \brief Tests for intersection with a box. */
      virtual tribool intersects(const Box<R>& r) const;

      /*! \brief Tests for disjointness with a box. */
      virtual tribool disjoint(const Box<R>& r) const;

      /*! \brief Tests for subset of a box. */
      virtual tribool subset(const Box<R>& r) const;

      /*! \brief Return a rectangle containing the set. */
      virtual Box<R> bounding_box() const;
        
      //@}

      /*! \brief Checks for emptyness. Returns \a true if the denotable set is empty, \a false otherwise. */
      tribool empty() const;

      /*! \brief Checks for boundedness. */
      tribool bounded() const;

      /*! \brief Make the set empty. */
      void clear();
      
      /*! \brief An iterator to the beginning of the list of basic sets. */
      iterator begin();

      /*! \brief An iterator to the end of the list of basic sets. */
      iterator end();

      /*! \brief A constant iterator to the beginning of the list of basic sets. */
      const_iterator begin() const;

      /*! \brief A constant iterator to the end of the list of basic sets. */
      const_iterator end() const;

      /*! \brief Adjoins (makes union with) a basic set. */
      void adjoin(const Box<R>& bx);

      /*! \brief Adjoins (makes union with) another list set. */
      void adjoin(const BoxListSet<R>& bxls);
      
      //@{
      //! \name Input/output operators
      /*! \brief Write to an output stream. */
      std::string summary() const;
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}

#ifdef DOXYGEN
      //@{ 
      //! \name Geometric binary predicates
      /*! \brief Tests disjointness.
       */
      friend tribool disjoint(const BoxListSet<BS>& A,
                              const BoxListSet<BS>& B);

      /*! \brief Tests inclusion of \a A in \a B.
       */
      friend tribool subset(const BoxListSet<BS>& A,
                            const BoxListSet<BS>& B);
      //@}
      
      
      //@{
      //! \name Geometric binary operations
      /*! \brief The union of \a A and \a B.
       * <br><br>
       * Note that 'union' is a reserved word in C++.
       */
      //@}
#endif      
     private:
      static void _instantiate();

    };

  
    template<class R>
    std::ostream& operator<<(std::ostream& os, const BoxListSet<R>& A);


    template<class R>
    std::istream&
    operator>>(std::istream& is, BoxListSet<R>& A);




    template<class R>
    BoxListSet<R>
    open_intersection(const BoxListSet<R>& A,
                      const BoxListSet<R>& B);
  
    template<class R>
    tribool
    disjoint(const BoxListSet<R>& A,
             const BoxListSet<R>& B);

    template<class R>
    tribool
    subset(const BoxListSet<R>& A,
           const BoxListSet<R>& B);
    
    
    template<class R>
    tribool
    subset(const BoxListSet<R>& A,
           const Box<R>& B);
    
    
    
    template<class R>
    BoxListSet<R>
    join(const BoxListSet<R>& A,
         const BoxListSet<R>& B);
    

    template<class R>
    BoxListSet<R>
    open_intersection(const BoxListSet<R>& A,
                      const BoxListSet<R>& B);
  
    template<class R>
    BoxListSet<R>
    inner_intersection(const BoxListSet<R>& A,
                       const SetInterface<typename R::real_type>& B);
  
    template<class R>
    BoxListSet<R>
    lower_intersection(const BoxListSet<R>& A,
                       const SetInterface<typename R::real_type>& B);
  
    template<class R>
    BoxListSet<R>
    outer_intersection(const BoxListSet<R>& A,
                       const SetInterface<typename R::real_type>& B);
  
  
} // namespace Ariadne

#include "box_list_set.inline.h"

#endif /* ARIADNE_BOX_LIST_SET_H */
