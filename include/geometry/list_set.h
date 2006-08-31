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
#include <exception>
#include <stdexcept>

#include <vector>

#include "../declarations.h"
#include "../base/utility.h"

namespace Ariadne {
  namespace Geometry {

    template<typename R> class Rectangle;
    template<typename R, template<typename> class BS> class ListSet;
  
    template<typename R, template<typename> class BS>
    std::ostream& operator<<(std::ostream&, const ListSet<R,BS>&);
    template<typename R, template<typename> class BS>
    std::istream& operator>>(std::istream&, ListSet<R,BS>&);
    
    /*! \brief A finite union of basic sets, represented as a sequence.
     *  \ingroup DenotableSet
     *  \ingroup List
     */
    template<typename R, template<typename> class BS>
    class ListSet {
     private:
      /* List of basic sets. Note that std::vector provides a
       * reserve(size_type) method to increase the capacity.
       */
      dimension_type _dimension;
      std::vector< BS<R> > _vector;

     public:
      typedef R real_type;
      typedef typename BS<R>::state_type state_type;
      typedef BS<R> basic_set_type;
      
      typedef BS<R> value_type;
      typedef typename std::vector<basic_set_type>::const_iterator const_iterator;
      typedef typename std::vector<basic_set_type>::iterator iterator;

      /*! \brief An empty list set which can hold sets of an unspecified dimension. */
      ListSet() : _dimension(0), _vector() { }

      /*! \brief An empty ListSet which can only hold sets of dimension \a n. */
      ListSet(size_type n) : _dimension(n), _vector() { }

      /*! \brief A denotable set constructor. */
      ListSet(const BS<R>& A) : _dimension(A.dimension()), _vector() {
        if (A.empty()) {
            return;
        }
        _vector.push_back(A);
      }

      /*! \brief The copy constructor. */
      ListSet(const ListSet<R,BS>& A) : _dimension(A.dimension()), _vector(A._vector) { }

      /*! \brief The destructor. */
      ~ListSet() {
        this->_vector.clear();
      }

      /*! \brief Return the number of basic sets forming this object.
      *
      * \return The number of basic sets forming this object.
      */
      const size_type size() const {
          return this->_vector.size();
      }

      /*! \brief Adjoins a basic set to the back of the list. */
      void push_back(const BS<R>& A) {
        if (this->dimension()==0) { this->_dimension=A.dimension(); }

        if (A.dimension()!=this->dimension()) {
          throw std::invalid_argument("The denotable set has a different space dimension to the list.");
        }

        _vector.push_back(A);
      }

      /*! \brief Return the denotable set's space dimension. 
      *
      * \return The space dimension of the ListSet.
      */
      const size_type dimension() const {
          return this->_dimension;
      }

      /*! \brief Accesses the i-th BasicSet.
      *
      * \param index is the index of the returned basic set.
      * \return The i-th basic set maitained by the ListSet.
      */
      const BS<R>& get(size_type index) const {
        if (this->size()<=index)
          throw std::invalid_argument("Invalid index in ListSet::get(size_type).");

        return this->_vector[index];
      }

      /*! \brief Assigns to the i-th BasicSet.
      *
      * \param index is the index of the returned basic set.
      * \param set is the new set.
      */
      void set(size_type index, const BS<R>& set) {
        if (this->size()<=index) 
          throw std::invalid_argument("Invalid index in ListSet::set(size_type, const BasicSet& ).");
        
        this->_vector[index]=set;
      }

      /*! \brief Accesses the i-th BasicSet.
      *
      * \param index is the index of the returned basic set.
      * \return The i-th basic set maitained by the ListSet.
      */
      const BS<R>& operator[](size_type index) const {
        if (this->size()<=index)
          throw std::invalid_argument("Invalid index in ListSet::operator[](size_type).");

        return this->_vector[index];
      }

      /*! \brief What does this do? */
      ListSet<R,BS> operator+(const ListSet<R,BS>& A) const{

        #ifdef DEBUG
          std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif

        ListSet<R,BS> sum(A.dimension());

        for (size_type i=0; i< this->size(); i++) {
          for (size_type j=0; j< A.size(); j++) {
              // sum.inplace_union((this->_vector[i])+A[j]);
          }
        }

        #ifdef DEBUG
          std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif

        return sum;
      }


      /*! \brief Copy assignment. */
      const ListSet<R,BS>& 
      operator=(const ListSet<R,BS>& A) {
        #ifdef DEBUG
          std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif

        if(this !=& A) {
            this->_dimension = A._dimension;
            this->_vector = A._vector;
        }

        #ifdef DEBUG
          std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif

        return *this;
      }

      template <template<typename> class BS2>
      operator ListSet<R,BS2> () const {
         ListSet<R,BS2> result(this->dimension());

         for (size_type i=0; i<this->size(); i++) 
           result.push_back((BS2<R>)((*this)[i]));
      
         return result;
      }
      
      /*! \brief Checks if a denotable set includes a point.
      *
      * This method checks whenever the current denotable set
      * includes the point \a p.
      * \param p is the point of which inclusion
      * into the current denotable set should be tested.
      * \return  \a true, if \a s is contained into the
      * current set, \a false otherwise.
      */
      bool contains(const Point<R>& p) const {

        for (size_type i=0; i<this->size(); i++) {
          if ((this->_vector[i]).contains(p))
            return true;
        }

        return false;
      }

      /*! \brief Checks if the interior of the denotable set includes a point.
       * FIXME! Incorrect since interior of union may be larger than union of interiors.
       *
       * This method checks whenever the interior of the current denotable set
       * includes the point \a s.
       * \param p is the point of which inclusion
       * into the current denotable set should be tested.
       * \return  \a true, if \a s is contained into the
       * current set, \a false otherwise.
       */
      bool interior_contains(const Point<R>& p) const {
        throw(std::domain_error("Not implemented"));
        for (size_type i=0; i<this->size(); i++) {
          if ((this->_vector[i]).interior_contains(p))
            return true;
        }

        return false;
      }

      /*! \brief Checks for emptyness.
       *
       * \return \a true if the denotable set is empty,
       * \a false otherwise.
       */
      bool empty() const {
        for (size_type i=0; i<this->size(); ++i) {
          if (!(this->_vector[i]).empty())
            return false;
        }

        return true;
      }

      /*! \brief Return a rectangle containing the set. */
      Rectangle<R> bounding_box() const {
        assert(!this->empty());
        Rectangle<R> result=(*this)[0].bounding_box();
        for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
          Rectangle<R> bb=iter->bounding_box();
          for(size_type i=0; i!=result.dimension(); ++i) {
            if(bb.lower_bound(i) < result.lower_bound(i)) {
              result.set_lower_bound(i,bb.lower_bound(i));
            }
            if(bb.upper_bound(i) > result.upper_bound(i)) {
              result.set_upper_bound(i,bb.upper_bound(i));
            }
          }
        }
        return result;
      }
        
      /*! \brief Make the set empty. */
      void clear() { 
        this->_vector.clear();
      }
      
      /*! \brief An iterator to the beginning of the list of basic sets.
       *
       * \return The begin of the maintained basic set vector.
       */
      iterator begin() {
          return _vector.begin();
      }

      /*! \brief An iterator to the end of the list of basic sets.
       *
       * \return The end of the maintained basic set vector.
       */
      iterator end() {
        return _vector.end();
      }

      /*! \brief A constant iterator to the beginning of the list of basic sets.
       *
       * \return The begin of the maintained basic set vector.
       */
      const_iterator begin() const {
          return _vector.begin();
      }

      /*! \brief A constant iterator to the end of the list of basic sets.
       *
       * \return The end of the maintained basic set vector.
       */
      const_iterator end() const {
        return _vector.end();
      }

      /*! \brief Adjoins (makes union with) another denotable set.
       *
       * Makes the union of the current denotable set with another of the same type.
       * The result is stored into the current object.
       * \param A is a ListSet.
       */
      void adjoin(const ListSet<R,BS>& A) {
        if(this->dimension()==0) { 
          this->_dimension=A.dimension(); 
        }

        if(A.dimension()!=this->dimension()) {
          throw std::invalid_argument("The two denotable sets have different space dimensions.");
        }

        this->_vector.reserve(A.size());
        for(typename ListSet<R,BS>::const_iterator iter=A.begin(); iter!=A.end(); ++iter) {
          this->_vector.push_back(*iter);
        }
      }
      
      /*! \brief Adjoins (makes union with) another denotable set.
       *
       * Makes the union of the current denotable set with another of the same type.
       * The result is stored into the current object.
       * \param A is a ListSet.
       */
      void inplace_union(const ListSet<R,BS>& A) {
        this->adjoin(A);
      }

      /*! \brief Adjoins (makes union with) a basic set.
       *
       * Makes the union of the current denotable set with a basic set.
       * The result is stored into the current object.
       * \param A is a basic set.
       */
      void adjoin(const BS<R>& A) {
        if(this->dimension()==0) { 
          this->_dimension=A.dimension(); 
        }
        if (A.dimension()!=this->dimension()) {
          throw std::invalid_argument("The denotable set and the basic set have different space dimensions.");
        }
        if(!A.empty()) {
          this->_vector.push_back(A);
        }
      }

      /*! \brief Adjoins (makes union with) a basic set. */
      void inplace_union(const BS<R>& A) {
        this->adjoin(A);
      }

#ifdef DOXYGEN
      //@{ 
      //! \name Geometric binary predicates
      /*! \brief Tests disjointness.
       */
      friend bool disjoint(const ListSet<R,BS>& A,
                           const ListSet<R,BS>& B);

      /*! \brief Tests if the interior of \a A intersects the interior of \a B.
       */
      friend bool interiors_intersect(const ListSet<R,BS>& A,
                                      const ListSet<R,BS>& B);

      /*! \brief Tests inclusion of \a A in the interior of \a B.
       */
      friend bool inner_subset(const ListSet<R,BS>& A,
                               const ListSet<R,BS>& B);

      /*! \brief Tests inclusion of \A in \B.
       */
      friend bool subset(const ListSet<R,BS>& A,
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


      //@{
      //! \name Input/output operators
      /*! \brief Stream insertion operator. */
      friend std::ostream& operator<<(std::ostream& is, const ListSet<R,BS>& S);

      /*! \brief Stream extraction operator. */
      friend std::istream& operator>>(std::istream& is, ListSet<R,BS>& S);
      //@}
#endif      
      
      friend std::ostream& operator<< <> (std::ostream& is, const ListSet<R,BS>& S);

      friend std::istream& operator>> <> (std::istream& is, ListSet<R,BS>& S);
      //@}
      
    };



    template <typename R, template<class> class BS>
    bool
    disjoint (const ListSet<R,BS>& A,
              const ListSet<R,BS>& B);

    template <typename R, template <typename> class BS>
    bool
    interiors_intersect(const ListSet<R,BS>& A,
                        const ListSet<R,BS>& B);


    template <typename R, template <typename> class BS>
    bool
    inner_subset(const ListSet<R,BS>& A,
                 const ListSet<R,BS>& B);




    template <typename R, template <typename> class BS>
    bool
    subset(const ListSet<R,BS>& A,
           const ListSet<R,BS>& B);
    
    
    template <typename R, template <typename> class BS>
    ListSet<R,BS>
    join(const ListSet<R,BS>& A,
         const ListSet<R,BS>& B);
    

    template <typename R, template <typename> class BS>
    ListSet<R,BS>
    regular_intersection(const ListSet<R,BS>& A,
                         const ListSet<R,BS>& B);
  
  }
}

#endif /* _LIST_SET_H */
