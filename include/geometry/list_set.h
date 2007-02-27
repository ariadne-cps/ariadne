/***************************************************************************
 *            list_set.h
 *
 *  23 June 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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
 *  \brief Denotable sets implemented as lists.
 */

#ifndef _ARIADNE_LIST_SET_H
#define _ARIADNE_LIST_SET_H

#include <iosfwd>
#include <iostream>
#include <exception>
#include <stdexcept>

#include <vector>

#include "../declarations.h"
#include "../base/utility.h"
#include "../base/tribool.h"
#include "../geometry/set.h"

namespace Ariadne {
  namespace Geometry {

    template<class R> class Rectangle;
      
    /*!\ingroup DenotableSet
     * \ingroup List
     * \brief A finite union of basic sets, represented as a sequence.
     *
     * A list set is the simplest type of denotable set class. It can hold
     * arbitrary lists of basic sets of the same type. Hence, the %ListSet class
     * also takes a template parameter which is the type of basic set contained
     * in the list.
     *
     * A list set is ordered by the order of insertion. Hence, as well as the
     * standard adjoin() method for denotable sets, a %ListSet also provides
     * the STL methods push_back() and pop_back().
     */
    template<class R, template<class> class BS>
    class ListSet 
      : public Set<R>
    {
     private:
      /* List of basic sets. Note that std::vector provides a
       * reserve(size_type) method to increase the capacity.
       */
      dimension_type _dimension;
      std::vector< BS<R> > _vector;

     public:
      /*!\brief The type of denotable real number used to represent the sets in the list. */
      typedef R real_type;
      /*!\brief The type of point contained by the set. */
      typedef typename BS<R>::state_type state_type;
      /*!\brief The type of basic set making up the denotable set. */
      typedef BS<R> basic_set_type;
      /*!\brief The type of basic set in the list of sets. */
      typedef BS<R> value_type;

      typedef typename std::vector<basic_set_type>::const_iterator const_iterator;
      typedef typename std::vector<basic_set_type>::iterator iterator;

     public:
      /*! \brief An empty list set which can hold sets of an unspecified dimension. */
      ListSet();

      /*! \brief An empty ListSet which can only hold sets of dimension \a n. */
      ListSet(size_type n);

      /*! \brief A denotable set constructor. */
      ListSet(const BS<R>& A);

      /*! \brief The copy constructor. */
      ListSet(const ListSet<R,BS>& A);

      /*! \brief The destructor. */
      ~ListSet();

      /*! \brief Returns the number of basic sets forming this object. */
      size_type size() const;

      /*! \brief Adjoins a basic set to the back of the list. */
      void push_back(const BS<R>& A);

      /*! \brief Removes the basic set at the back of the list. */
      void pop_back();

      /*! \brief Accesses the i-th basic set in the list. */
      const BS<R>& get(size_type index) const;

      /*! \brief Assigns to the i-th BasicSet. */
      void set(size_type index, const BS<R>& set);

      /*! \brief Accesses the i-th BasicSet. */
      const BS<R>& operator[](size_type index) const;


      /*! \brief Copy assignment. */
      const ListSet<R,BS>& operator=(const ListSet<R,BS>& A);

      /*!\brief Convert to a list set BSt<Rl>. */
      template<class Rl, template<class> class BSt> 
      operator ListSet<Rl,BSt> () const;
      
      //@{
      //! \name Set methods
      /*! \brief Tests for disjointness with a Rectangle. */
      virtual ListSet<R,BS>* clone() const;

      /*! \brief Returns the denotable set's space dimension. */
      virtual dimension_type dimension() const;

      /*!\brief Checks if a denotable set includes a point. */
      virtual tribool contains(const Point<R>& p) const;

      /*! \brief Tests for disjointness with a Rectangle. */
      virtual tribool disjoint(const Rectangle<R>& r) const;

      /*! \brief Tests for superset of a Rectangle. 
       *
       * Currently always returns \a indeterminate, since the 
       * test is difficult for general list sets. 
       */
      virtual tribool superset(const Rectangle<R>& r) const;

      /*! \brief Tests for subset of a Rectangle. */
      virtual tribool subset(const Rectangle<R>& r) const;

      /*! \brief Return a rectangle containing the set. */
      virtual Rectangle<R> bounding_box() const;
        
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

      /*! \brief Adjoins (makes union with) another denotable set. */
      void adjoin(const ListSet<R,BS>& A);
      
      /*! \brief Adjoins (makes union with) another denotable set. */
      void inplace_union(const ListSet<R,BS>& A);

      /*! \brief Adjoins (makes union with) a basic set. */
      void adjoin(const BS<R>& A);

      /*! \brief Adjoins (makes union with) a basic set. */
      void inplace_union(const BS<R>& A);
      
      //@{
      //! \name Input/output operators
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
      friend tribool disjoint(const ListSet<R,BS>& A,
                           const ListSet<R,BS>& B);

      /*! \brief Tests inclusion of \a A in \a B.
       */
      friend tribool subset(const ListSet<R,BS>& A,
                         const ListSet<R,BS>& B);
      //@}
      
      
      //@{
      //! \name Geometric binary operations
      /*! \brief The union of \a A and \a B.
       * <br><br>
       * Note that 'union' is a reserved word in C++.
       */
      //FIXME: Compiler doesn't like this
      //friend ListSet<R,BS> join<> (const ListSet<R,BS>& A,
      //                             const ListSet<R,BS>& B);
      friend ListSet<R,BS> join(const ListSet<R,BS>& A,
                                const ListSet<R,BS>& B);

      /*! \brief The closure of intersection of the interior of \a A with the interior of \a B.
       */
      friend ListSet<R,BS> regular_intersection(const ListSet<R,BS>& A,
                                                const ListSet<R,BS>& B);
      //@}
#endif      
     private:
      static void _instantiate_geometry_operators();

    };

  
    template<class R, template<class> class BS>
    std::ostream& operator<<(std::ostream& os, const ListSet<R,BS>& A);


    template<class R, template<class> class BS>
    std::istream&
    operator>>(std::istream& is, ListSet<R,BS>& A);




    template<class R, template<class> class BS>
    ListSet<R,BS>
    open_intersection(const ListSet<R,BS>& A,
                      const ListSet<R,BS>& B);
  


    template<class R, template<class> class BS>
    tribool
    disjoint(const ListSet<R,BS>& A,
             const ListSet<R,BS>& B);


    template<class R>
    tribool
    subset(const ListSet<R,Rectangle>& A,
           const ListSet<R,Rectangle>& B);
    
    
    template<class R, template<class> class BS>
    tribool
    disjoint(const ListSet<R,BS>& A,
             const Rectangle<R>& B);
    
    
    template<class R, template<class> class BS>
    tribool
    subset(const ListSet<R,BS>& A,
           const Rectangle<R>& B);
    
    
    template<class R, template<class> class BS>
    ListSet<R,BS>
    join(const ListSet<R,BS>& A,
         const ListSet<R,BS>& B);
    

    template<class R, template<class> class BS>
    ListSet<R,BS>
    open_intersection(const ListSet<R,BS>& A,
                         const ListSet<R,BS>& B);
  
  }
}

#include "list_set.inline.h"

#endif /* _LIST_SET_H */
