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
 
#ifndef _LIST_SET_H
#define _LIST_SET_H

#include <vector>
#include <iosfwd>
#include <iostream>

namespace Ariadne {
  namespace Geometry {

    template<typename R> class Rectangle;
    template<typename R, template<typename> class BS> class ListSet;

    /*! \brief Stream extraction operator. */
    template <typename R, template<typename> class BS>
    std::istream&
    operator>> (std::istream &is, const ListSet<R,BS> &S);

    /*! \brief A finite union of basic sets, represented as a sequence.
     */
    template<typename R, template<typename> class BS>
    class ListSet {
     private:
      /* List of basic sets. Note that std::vector provides a
       * reserve(size_t) method to increase the capacity.
       */
      size_t _dimension;
      std::vector< BS<R> > _vector;

     public:
      typedef BS<R> BasicSet;
      typedef typename BasicSet::State State;
      typedef typename BasicSet::State::Real Real;

      typedef typename std::vector<BasicSet>::const_iterator const_iterator;
      typedef typename std::vector<BasicSet>::iterator iterator;

      /*! \brief A denotable set constructor. */
      ListSet() : _dimension(0), _vector() {
      }

      /*! \brief A denotable set constructor. */
      ListSet(const BasicSet &A) : _dimension(A.dimension()), _vector() {
        if (A.empty()) {
            return;
        }
        _vector.push_back(A);
      }

      /*! \brief The copy constructor. */
      ListSet(const ListSet<R,BS>& A) : _dimension(A.dimension()), _vector(A._vector) { }

      /*! \brief Construct a ListSet<BasicSet> to hold sets of dimension \a n. */
      ListSet(size_t n) : _dimension(n), _vector() { }

      /*! \brief The destructor. */
      ~ListSet() {
        this->_vector.clear();
      }

      /*! \brief Return the number of basic sets forming this object.
      *
      * \return The number of basic sets forming this object.
      */
      inline const size_t size() const {
          return this->_vector.size();
      }

      /*! \brief Adjoins a basic set to the back of the list. */
      inline void push_back(const BS<R>& A) {
        if (this->dimension()==0) { this->_dimension=A.dimension(); }

        if (A.dimension()!=this->dimension()) {
          throw std::invalid_argument("The two denotable set has a different space dimension to the list.");
        }

        _vector.push_back(A);
      }

      /*! \brief Return the denotable set's space dimension. 
      *
      * \return The space dimension of the ListSet.
      */
      inline const size_t dimension() const {
          return this->_dimension;
      }

      /*! \brief Accesses the i-th BasicSet.
      *
      * \param index is the index of the returned basic set.
      * \return The i-th basic set maitained by the ListSet.
      */
      inline const BasicSet& get(size_t index) const {
        if (this->size()<=index)
          throw std::invalid_argument("Index overlaps vector bounds.");

        return this->_vector[index];
      }

      /*! \brief Assigns to the i-th BasicSet.
      *
      * \param index is the index of the returned basic set.
      * \param set is the new set.
      */
      inline void set(size_t index, const BasicSet& set) {
        if (this->size()<=index)
          throw std::invalid_argument("Index overlaps vector bounds.");

        this->_vector[index]=set;
      }

      /*! \brief Accesses the i-th BasicSet.
      *
      * \param index is the index of the returned basic set.
      * \return The i-th basic set maitained by the ListSet.
      */
      inline const BasicSet& operator[](size_t index) const {
        if (this->size()<=index)
          throw std::invalid_argument("Index overlaps vector bounds.");

        return this->_vector[index];
      }

      /*! \brief What does this do? */
      inline ListSet<R,BS> operator+(const ListSet<R,BS>& A) const{

        #ifdef DEBUG
          std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif

        ListSet<R,BS> sum(A.dimension());

        for (size_t i=0; i< this->size(); i++) {
          for (size_t j=0; j< A.size(); j++) {
              // sum.inplace_union((this->_vector[i])+A[j]);
          }
        }

        #ifdef DEBUG
          std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif

        return sum;
      }

      /*! \brief Expands set by \a delta. */
      /*  \internal This is dangerous since it modifies the current set.
      */
      inline void expand_by(const Real &delta) {
        for (size_t i=0; i< this->size(); i++) {
          (this->_vector[i]).expand_by(delta);
        }
      }

      /*! \brief Replaces set be an over-approximation by at most \a delta. */
      /*  \internal This is dangerous since it modifies the current set.
      */
      inline void
      set_precision_to_upperapproximating(const Real &delta)  {
        for (size_t i=0; i< this->size(); i++) {
            //				(this->_vector[i]).set_precision_to_upperapproximating(delta);
        }
      }

      /*! \brief Copy assignment. */
      inline const ListSet<R,BS> &
      operator=(const ListSet<R,BS> &A) {
        #ifdef DEBUG
          std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif

        if(this != &A) {
            this->_dimension = A._dimension;
            this->_vector = A._vector;
        }

        #ifdef DEBUG
          std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif

        return *this;
      }

      /*! \brief Checks if a denotable set includes a state.
      *
      * This method checks whenever the current denotable set
      * includes the state \a s.
      * \param s is the state of which inclusion
      * into the current denotable set should be tested.
      * \return  \a true, if \a s is contained into the
      * current set, \a false otherwise.
      */
      inline bool contains(const State &s) const {

        for (size_t i=0; i<this->size(); i++) {
          if ((this->_vector[i]).contains(s))
            return true;
        }

        return false;
      }

      /*! \brief Checks if the interior of the denotable set includes a state.
       * FIXME! Incorrect since interior of union may be larger than union of interiors.
       *
       * This method checks whenever the interior of the current denotable set
       * includes the state \a s.
       * \param s is the state of which inclusion
       * into the current denotable set should be tested.
       * \return  \a true, if \a s is contained into the
       * current set, \a false otherwise.
       */
      inline bool interior_contains(const State & state) const {
        throw(std::domain_error("Not implemented"));
        for (size_t i=0; i<this->size(); i++) {
          if ((this->_vector[i]).interior_contains(state))
            return true;
        }

        return false;
      }

      /*! \brief Checks for emptyness.
       *
       * \return \a true if the denotable set is empty,
       * \a false otherwise.
       */
      inline bool empty() const {
        for (size_t i=0; i<this->size(); i++) {
          if (!(this->_vector[i]).empty())
            return false;
        }

        return true;
      }

      /*! \brief A constant iterator to the beginning of the list of basic sets.
       *
       * \return The begin of the maintained basic set vector.
       */
      inline const_iterator begin() const {
          return _vector.begin();
      }

      /*! \brief A constant iterator to the end of the list of basic sets.
       *
       * \return The end of the maintained basic set vector.
       */
      inline const_iterator end() const {
        return _vector.end();
      }

      /*! \brief Adjoins (makes union with) another denotable set.
       *
       * Makes the union of the current denotable set with another of the same type.
       * The result is stored into the current object.
       * \param A is a ListSet.
       */
      inline void inplace_union(const ListSet<R,BS>& A) {
        #ifdef DEBUG
          std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif

        if (this->dimension()==0) { this->_dimension=A.dimension(); }

        if (A.dimension()!=this->dimension()) {
          throw std::invalid_argument("The two denotable set have different space dimensions.");
        }

        if (A.empty()) {
          #ifdef DEBUG
            std::cout << __FILE__ << ":" << __LINE__ << std::endl;
          #endif

          return;
        }

        this->_vector.reserve(A.size());
        for (size_t i = 0; i < A.size(); i++) {
          this->_vector.push_back(A[i]);
        }

        #ifdef DEBUG
          std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif
      }

      /*! \brief Adjoins (makes union with) a basic set.
       *
       * Makes the union of the current denotable set with a basic set.
       * The result is stored into the current object.
       * \param A is a basic set.
       */
      inline void inplace_union(const BasicSet &A) {
        #ifdef DEBUG
          std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif

        if (this->dimension()==0) { this->_dimension=A.dimension(); }

        if (A.dimension()!=this->dimension()) {
          throw std::invalid_argument("The denotable set the basic set have different space dimensions.");
        }

        if (!A.empty()) {
            this->_vector.push_back(A);
        }

        #ifdef DEBUG
          std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif
      }

      /*! \brief Stream extraction operator. */
      friend std::istream& operator>> <> (std::istream &is,
                                          ListSet<R,BS> &S);

      /*! \brief Tests disjointness.
       *
       * Tests disjointness of two denotable sets.
       * \param A is a denotable set.
       * \param B is a denotable set.
       * \return \a true if A and B are disjoint, \a false otherwise.
       */
      friend bool disjoint <> (const ListSet<R,BS> &A,
                               const ListSet<R,BS> &B);

      /*! \brief Tests intersection of interiors.
       *
       * Tests intersection of a denotable set with the interior of an other
       * denotable set.
       * \param A is a denotable set.
       * \param B is a denotable set.
       * \return \a true if A intersects the interior of B,
       * \a false otherwise.
       */
      friend bool interiors_intersect<> (const ListSet<R,BS> &A,
                                         const ListSet<R,BS> &B);

      /*! \brief Tests intersection of interiors.
       *
       * Tests intersection of a denotable set with the interior of a rectangle.
       * \param A is a rectangle.
       * \param B is a denotable set.
       * \return \a true if A intersects the interior of B,
       * \a false otherwise.
       */
      friend bool interiors_intersect<> (const Rectangle<R> &rect,
                                          const ListSet<R,BS> &A);

      /*! \brief Tests inclusion of interiors.
       *
       * Tests inclusion of a denotable set into the interior of a denotable set.
       * \param A is a denotable set.
       * \param B is a denotable set.
       * \return \a true if A is a subset of the interior of B,
       * \a false otherwise.
       */
      friend bool inner_subset<> (const ListSet<R,BS> &A,
                                  const ListSet<R,BS> &B);

      /*! \brief Tests inclusion of interiors.
       *
       * Tests inclusion of a rectangle into the interior of a denotable set.
       * \param A is a rectangle.
       * \param B is a denotable set.
       * \return \a true if A is a subset of the interior of B,
       * \a false otherwise.
       */
      friend bool inner_subset<> (const Rectangle<R> &rect,
                                  const ListSet<R,BS> &A);

      /*! \brief Tests inclusion of interiors.
       *
       * Tests inclusion of a denotable set into the interior of a rectangle.
       * \param A is a denotable set.
       * \param B is a rectangle.
       * \return \a true if A is a subset of the interior of B,
       * \a false otherwise.
       */
      friend bool inner_subset<> (const ListSet<R,BS> &A,
                                  const Rectangle<R> &rect);

      /*! \brief Makes union of two interiors.
       *
       * Evalutates the union of two denotable sets.
       * \param A is a denotable set.
       * \param B	is a denotable set.
       * \return The union of A and B.
       */
      friend ListSet<R,BS> join<> (const ListSet<R,BS> &A,
                                   const ListSet<R,BS> &B);

      /*! \brief Makes intersection of two interiors.
       *
       * Evalutates the closure of the intersection of a denotable set
       * with the interiors of an other denotable sets.
       * \param A is a denotable set.
       * \param B is a denotable set.
       * \return The closure of the intersection of A with the interiors of B.
       */
      friend ListSet<R,BS>
      regular_intersection <> (const ListSet<R,BS> &A,
                               const ListSet<R,BS> &B);
    };



    template <typename R, template<class> class BS>
    bool
    disjoint (const ListSet<R,BS> &A,
              const ListSet<R,BS> &B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::invalid_argument("The two denotable set have different space dimensions.");
      }

      size_t i,j;
      for (i=0; i<A.size() ; i++) {
        for (j=0; j<B.size() ; j++) {
          if (!disjoint(A[i],B[j])) { return false; }
        }
      }

      return true;
    }

    template <typename R, template <typename> class BS>
    bool
    interiors_intersect(const ListSet<R,BS> &A,
                        const ListSet<R,BS> &B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::invalid_argument("The two denotable set have different space dimensions.");
      }

      size_t i,j;

      for (i=0; i<A.size() ; i++) {
        for (j=0; j<B.size() ; j++) {
          if (intersects_interior(A[i],B[j])) { return true; }
        }
      }

      return false;
    }

    template <typename R, template <typename> class BS>
    bool
    interiors_intersect(const Rectangle<R> &rect,
                        const ListSet<R,BS> &A)
    {
      if (A.dimension()!=rect.dimension()) {
        throw std::invalid_argument(
          "The denotable set and the rectangle have different space dimensions.");
      }

      for (size_t i=0; i<A.size() ; i++) {
        if (intersects_interior(rect,A[i])) { return true; }
      }

      return false;
    }

    template <typename R, template <typename> class BS>
    bool
    interiors_intersect(const ListSet<R,BS> &A,
                        const Rectangle<R> &rect)
    {
      if (A.dimension()!=rect.dimension()) {
        throw std::invalid_argument(
          "The denotable set and the rectangle have different space dimensions.");
      }

      for (size_t i=0; i<A.size() ; i++) {
        if (intersects_interior(A[i],rect)) { return true; }
      }

      return false;
    }

    template <typename R, template <typename> class BS>
    bool
    inner_subset(const ListSet<R,BS> &A,
                 const ListSet<R,BS> &B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::invalid_argument("The two denotable set have different space dimensions.");
      }

      for (size_t i=0; i<A.size() ; i++) {
        if (!subset_of_open_cover(A[i], B._vector)) { return false; }
      }

      return true;
    }

    template <typename R, template <typename> class BS>
    bool
    inner_subset(const Rectangle<R> &rect,
                 const ListSet<R,BS> &A)
    {
      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif

      if (A.dimension()!=rect.dimension()) {
        throw std::invalid_argument(
          "The denotable set and the rectangle have different space dimensions.");
      }

      /* FIXME: TO REIMPLEMENT */
      for (size_t i=0; i<A.size() ; i++) {
        if (subset_of_interior(rect,A[i])) {
          #ifdef DEBUG
            std::cout << __FILE__ << ":" << __LINE__ << std::endl;
          #endif
          return true;
        }
      }

      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif

      return false;
    }

    template <typename R, template <typename> class BS>
    bool
    inner_subset(const ListSet<R,BS> &A,
                 const Rectangle<R> &rect)
    {
      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif

      if (A.dimension()!=rect.dimension()) {
        throw std::invalid_argument(
          "The denotable set and the rectangle have different space dimensions.");
      }

      for (size_t i=0; i<A.size() ; i++) {
        if (!subset_of_interior(A[i], rect)) {
          #ifdef DEBUG
            std::cout << __FILE__ << ":" << __LINE__ << std::endl;
          #endif
          return false;
        }
      }

      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif

      return true;
    }


    template <typename R, template <typename> class BS>
    bool
    subset(const ListSet<R,BS> &A,
           const ListSet<R,BS> &B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::invalid_argument("The two denotable set have different space dimensions.");
      }

      throw std::runtime_error("Not implemented");
    }

    template <typename R, template <typename> class BS>
    ListSet<R,BS>
    join(const ListSet<R,BS> &A,
         const ListSet<R,BS> &B)
    {
      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif

      if (A.dimension()!=B.dimension()) {
        throw std::invalid_argument("The two denotable set have different space dimensions.");
      }

      ListSet<R,BS> ds_union(A);

      ds_union.inplace_union(B);

      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif

      return ds_union;
    }

    template <typename R, template <typename> class BS>
    ListSet<R,BS>
    regular_intersection(const ListSet<R,BS> &A,
                         const ListSet<R,BS> &B)
    {
      ListSet<R,BS> ds_inter;
      std::vector< BS<R> > &vector=ds_inter._vector;

      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif

      if (A.dimension()!=B.dimension()) {
        throw std::invalid_argument("The two denotable set have different space dimensions.");
      }

      for (size_t i=0; i<A.size(); i++) {
        for (size_t j=0; j<B.size(); j++) {
          if (intersects_interior(A[i],B[j])) {
              vector.push_back(closure_of_intersection_of_interior(A[i],B[j]));
          }
        }
      }

      ds_inter._dimension=A._dimension;

      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif

      return ds_inter;

    }

    template <typename R, template<typename> class BS>
    std::ostream& operator<<(std::ostream &os,
                             const ListSet<R,BS> &A)
    {
      os << "{ class: ListSet<" << name<R>() << ",BS>,\n  basic_set_list: [\n    ";
      if (A.size() >0 ) {
        os << A[0];
      }
      for (size_t i=1; i<A.size(); i++) {
        os << ",\n    " << A[i];
      }
      os << "\n  ]\n}" << std::endl;
      return os;
    }

    template <typename R, template<typename> class BS>
    std::istream& operator>>(std::istream &is,
                             ListSet<R,BS> &A)
    {
      std::vector< BS<R> >& vec(A._vector);
      is >> vec;

      if(vec.size()==0) {
        A._dimension = 0;
      }
      else {
        A._dimension=vec[0].dimension();
      }
      return is;
    }


  }
}

#endif /* _LIST_SET_H */
