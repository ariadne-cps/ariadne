/***************************************************************************
 *            set_operations.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file set_operations.h
 *  \brief Operations for general sets.
 */

#ifndef ARIADNE_SET_OPERATIONS_H
#define ARIADNE_SET_OPERATIONS_H

#include <iosfwd>

#include "base/types.h"
#include "base/tribool.h"


namespace Ariadne {
  

    template<class R> class Point;
    template<class R> class Box;
      
    //! \ingroup SetInterface
    /*! \brief An class representing the intersection of two sets. */
    template<class R>
    class IntersectionSet {
     public:
      typedef R real_type;
      typedef Point<R> state_type;
     
      /*! \brief Construct the intersection of the sets \a s1 and s2. */
      IntersectionSet(const SetInterface<R>& s1, const SetInterface<R>& s2);

      /*! \brief A dynamically-allocated copy of the set. */
      virtual SetInterface<R>* clone() const;
     
      /*! \brief The dimension of the Euclidean space the set lies in. */
      virtual dimension_type dimension() const;
      
      /*! \brief Tests if the set contains a point. */
      virtual tribool contains(const Point<R>&) const;
     
      /*! \brief Tests if the set is a superset of a rectangle. */
      virtual tribool superset(const Box<R>&) const;
      /*! \brief Tests if the set intersects a rectangle. */
      virtual tribool intersects(const Box<R>&) const;
      /*! \brief Tests if the set is disjoint from a rectangle. */
      virtual tribool disjoint(const Box<R>&) const;
      /*! \brief Tests if the set is a subset of a rectangle. */
      virtual tribool subset(const Box<R>&) const;
      
      /*! \brief A rectangle containing the set. Throws UnboundedSet exception if the set is unbounded. */
      virtual Box<R> bounding_box() const;

      /*! \brief Write to an output stream. 
       *  Called by operator<<(std::ostream&, const SetInterface<R>&) to dynamically dispatch stream output. 
       */
      virtual std::ostream& write(std::ostream& os) const;
     private:
      boost::shared_ptr< const SetInterface<R> > _set1;
      boost::shared_ptr< const SetInterface<R> > _set2;
    };




    //! \ingroup SetInterface
    /*! \brief An class representing the union of two sets. */
    template<class R>
    class UnionSet {
     public:
      typedef R real_type;
      typedef Point<R> state_type;
     
       /*! \brief Construct the union of the sets \a s1 and s2. */
      UnionSet(const SetInterface<R>& s1, const SetInterface<R>& s2);

      /*! \brief A dynamically-allocated copy of the set. */
      virtual SetInterface<R>* clone() const;
     
      /*! \brief The dimension of the Euclidean space the set lies in. */
      virtual dimension_type dimension() const;
      
      /*! \brief Tests if the set contains a point. */
      virtual tribool contains(const Point<R>&) const;
     
      /*! \brief Tests if the set is a superset of a rectangle. */
      virtual tribool superset(const Box<R>&) const;
      /*! \brief Tests if the set intersects a rectangle. */
      virtual tribool intersects(const Box<R>&) const;
      /*! \brief Tests if the set is disjoint from a rectangle. */
      virtual tribool disjoint(const Box<R>&) const;
      /*! \brief Tests if the set is a subset of a rectangle. */
      virtual tribool subset(const Box<R>&) const;
      
      /*! \brief A rectangle containing the set. Throws UnboundedSet exception if the set is unbounded. */
      virtual Box<R> bounding_box() const;

      /*! \brief Write to an output stream. 
       *  Called by operator<<(std::ostream&, const SetInterface<R>&) to dynamically dispatch stream output. 
       */
      virtual std::ostream& write(std::ostream& os) const;
     private:
      boost::shared_ptr< const SetInterface<R> > _set1;
      boost::shared_ptr< const SetInterface<R> > _set2;
    };


    //! \ingroup SetInterface
    /*! \brief An class representing the complement of a set. */
    template<class R>
    class ComplementSet {
     public:
      typedef R real_type;
      typedef Point<R> state_type;
     
      /*! \brief Construct the complement of the set \a s. */
      ComplementSet(const SetInterface<R>& s);

      /*! \brief Construct the complement of the set \a s. */
      ComplementSet(boost::shared_ptr< const SetInterface<R> > s);

      /*! \brief A dynamically-allocated copy of the set. */
      virtual SetInterface<R>* clone() const;
     
      /*! \brief The dimension of the Euclidean space the set lies in. */
      virtual dimension_type dimension() const;
      
      /*! \brief Tests if the set contains a point. */
      virtual tribool contains(const Point<R>&) const;
     
      /*! \brief Tests if the set is a superset of a rectangle. */
      virtual tribool superset(const Box<R>&) const;
      /*! \brief Tests if the set intersects a rectangle. */
      virtual tribool intersects(const Box<R>&) const;
      /*! \brief Tests if the set is disjoint from a rectangle. */
      virtual tribool disjoint(const Box<R>&) const;
      /*! \brief Tests if the set is a subset of a rectangle. */
      virtual tribool subset(const Box<R>&) const;
      
      /*! \brief A rectangle containing the set. Throws UnboundedSet exception if the set is unbounded. */
      virtual Box<R> bounding_box() const;

      /*! \brief Write to an output stream. 
       *  Called by operator<<(std::ostream&, const SetInterface<R>&) to dynamically dispatch stream output. 
       */
      virtual std::ostream& write(std::ostream& os) const;
     private:
      boost::shared_ptr< const SetInterface<R> > _set;
    };


    template<class R>
    IntersectionSet<R>
    intersection(const SetInterface<R>&, const SetInterface<R>&);
  
    template<class R>
    UnionSet<R>
    join(const SetInterface<R>&, const SetInterface<R>&);
  
    template<class R>
    ComplementSet<R>
    complement(const SetInterface<R>&);
  
  }
}


namespace Ariadne {
  

template<class R>
IntersectionSet<R>::IntersectionSet(const SetInterface<R>& s1, 
                                              const SetInterface<R>& s2) 
  : _set1(s1.clone()), _set2(s2._clone()) 
{
}

template<class R>
IntersectionSet<R>*
IntersectionSet<R>::clone() const
{
  return new IntersectionSet<R>(*this);
}

template<class R>
tribool
IntersectionSet<R>::dimension() const {
  return this->_set1->dimension();
}

template<class R>
tribool
IntersectionSet<R>::contains(const Point<R>& x) const {
  return this->_set1->contains(x) && this->_set2->contains(x); 
}

template<class R>
tribool
IntersectionSet<R>::superset(const Box<R>& x) const {
  return this->_set1->superset(x) && this->_set2->superset(x); 
}

template<class R>
tribool
IntersectionSet<R>::disjoint(const Box<R>& x) const {
  if(this->_set1->disjoint(x) || this->_set2->disjoint(x)) {
    return true;
  } else {
    return indeterminate;
  }
}

template<class R>
tribool
IntersectionSet<R>::intersects(const Box<R>& x) const {
  return !this->disjoint(x); 
}

template<class R>
Box<R>
IntersectionSet<R>::bounding_box() const {
  return intersection(this->_set1->bounding_box(),this->_set2->bounding_box());
}




template<class R>
UnionSet<R>::UnionSet(const SetInterface<R>& s1, 
                                const SetInterface<R>& s2) 
  : _set1(s1.clone()), _set2(s2._clone()) 
{
}

template<class R>
UnionSet<R>*
UnionSet<R>::clone() const {
  return new UnionSet<R>(*this);
}

template<class R>
dimension_type
UnionSet<R>::dimension() const {
  return this->_set1->dimension();
}

template<class R>
tribool
UnionSet<R>::contains(const Point<R>& x) const {
  return this->_set1->contains(x) || this->_set2->contains(x); 
}

template<class R>
tribool
UnionSet<R>::superset(const Box<R>& x) const {
  if(this->_set1->superset(x) || this->_set2->superset(x)) {
    return true;
  } else {
    return indeterminate;
  }
}

template<class R>
tribool
UnionSet<R>::disjoint(const Box<R>& x) const {
  return this->_set1->disjoint(x) && this->_set2->disjoint(x);
}

template<class R>
tribool
UnionSet<R>::intersects(const Box<R>& x) const {
  return !this->disjoint(x); 
}

template<class R>
Box<R>
UnionSet<R>::bounding_box() const {
  return rectangular_hull(this->_set1->bounding_box(),this->_set2->bounding_box());
}




template<class R>
ComplementSet<R>::ComplementSet(boost::shared_ptr<const SetInterface<R> > s)
  : _set(s)
{
}

template<class R>
ComplementSet<R>::ComplementSet(const SetInterface<R>& s) const
  : _set(s.clone())
{
}

template<class R>
ComplementSet<R>*
ComplementSet<R>::clone() const {
  return new ComplementSet<R>(*this);
}

template<class R>
dimension_type
ComplementSet<R>::dimension() const {
  return this->_set1->dimension();
}

template<class R>
tribool
ComplementSet<R>::contains(const Point<R>& x) const {
  return !this->_set->contains(x); 
}

template<class R>
tribool
ComplementSet<R>::superset(const Box<R>& x) const {
  return this->_set->disjoint(x); 
}

template<class R>
tribool
ComplementSet<R>::disjoint(const Box<R>& x) const {
  return this->_set->superset();
}

template<class R>
tribool
ComplementSet<R>::intersects(const Box<R>& x) const {
  return !this->disjoint(x); 
}

template<class R>
Box<R>
ComplementSet<R>::bounding_box() const {
  Interval<R> r(-inf<R>(),inf<R>());
  return Box(this->dimension(),r);
}




template<class R> inline
IntersectionSet<R>
intersection(const SetInterface<R>& s1, const SetInterface<R>& s2) {
  return IntersectionSet(s1,s2);
}

template<class R>
UnionSet<R>
join(const SetInterface<R>& s1, const SetInterface<R>& s2) {
  return UnionSet(s1,s2);
}

template<class R>
ComplementSet<R>
complement(const SetInterface<R>& s) {
  return ComplementSet(s);
}
  



} // namespace Ariadne



#endif /* ARIADNE_SET_OPERATIONS_H */
