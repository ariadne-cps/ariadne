/***************************************************************************
 *            set_interface.h
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
 
/*! \file set_interface.h
 *  \brief Interface for general sets.
 */

#ifndef ARIADNE_SET_INTERFACE_H
#define ARIADNE_SET_INTERFACE_H

#include <iosfwd>

#include "base/types.h"
#include "base/tribool.h"


namespace Ariadne {
  namespace Geometry {

    template<class R> class Point;
    template<class R> class Box;
      
    //! \ingroup SetInterface 
    /*! %Base class for subsets of Euclidean space. */
    template<class R>
    class SetBaseInterface 
    {
     public:
      /*! \brief The type used to describe real numbers. */
      typedef R real_type;
      /*! \brief Virtual destructor. */
      virtual ~SetBaseInterface() { }
      /*! \brief Make a dynamically-allocated copy. */
      virtual SetBaseInterface<R>* clone() const = 0;
      /*! \brief The dimension of the space the set lies in. */
      virtual dimension_type dimension() const = 0;
      /*! \brief Test is the set contains a point. */
      virtual tribool contains(const Point<R>& pt) const = 0;
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const = 0;
#ifdef DOXYGEN
      /*! \brief Write to an output stream. */
      friend std::ostream& operator<<(std::ostream& os, const SetBaseInterface<R>& s);
#endif
     protected:
      // Only base classes may create the interface
      SetBaseInterface() { }
    };

    //! \ingroup SetInterface 
    /*! \brief Interface for open sets. */
    template<class R>
    class OpenSetInterface 
      : public virtual SetBaseInterface<R>
    {
     public:
      virtual OpenSetInterface<R>* clone() const = 0;
      /*! \brief Tests if an open set is a superset of a closed box. */
      virtual tribool superset(const Box<R>&) const = 0;
#ifdef DOXYGEN
      /*! \brief Tests if the set \a s is a subset of the box \a r. */
      friend tribool superset(const OpenSetInterface<R>& s, const Box<R>& r);
      /*! \brief Tests if the box \a r is a subset of the set \a s. */
      friend tribool subset(const Box<R>& r, const OpenSetInterface<R>& s);
#endif 
    };

    //! \ingroup SetInterface 
    /*! \brief Interface for the closed sets described by the lower Fell topology. */
    template<class R>
    class LowerClosedSetInterface 
      : public virtual SetBaseInterface<R>
    {
     public:
      virtual LowerClosedSetInterface<R>* clone() const = 0;
      /*! \brief Tests if a closed set intersects an open box. */
      virtual tribool intersects(const Box<R>&) const = 0;
#ifdef DOXYGEN
      /*! \brief Tests if the set \a s is intersects the box \a r. */
      friend tribool intersect(const LowerClosedSetInterface<R>& s, const Box<R>& r);
      /*! \brief Tests if the box \a r intersects the set \a s. */
      friend tribool intersect(const Box<R>& r, const LowerClosedSetInterface<R>& s);
#endif
    };

    //! \ingroup SetInterface 
    /*! \brief Interface for the closed sets described by the upper Fell topology. */    
    template<class R>
    class UpperClosedSetInterface 
      : public virtual SetBaseInterface<R>
    {
     public:
      virtual UpperClosedSetInterface<R>* clone() const = 0;
      /*! \brief Tests if a closed set is disjoint from a closed box. */
      virtual tribool disjoint(const Box<R>&) const = 0;
#ifdef DOXYGEN
      /*! \brief Tests if the set \a s is disjoint from the box \a r. */
      friend tribool disjoint(const UpperClosedSetInterface<R>& s, const Box<R>& r);
      /*! \brief Tests if the box \a r is disjoint from the set \a s. */
      friend tribool disjoint(const Box<R>& r, const UpperClosedSetInterface<R>& s);
#endif
    };

    //! \ingroup SetInterface 
    /*! \brief Interface for bounded sets. */    
    template<class R>
    class BoundedSetInterface 
      : public virtual SetBaseInterface<R>
   {
     public:
      virtual BoundedSetInterface<R>* clone() const = 0;
      /*! \brief Tests if a bouded set is a subset of from an open box. */
      virtual tribool subset(const Box<R>&) const = 0;
      /*! \brief Returns a bounding box for the set. */
      virtual Box<R> bounding_box() const = 0;
#ifdef DOXYGEN
      /*! \brief Tests if the set \a s is a subset of the box \a r. */
      friend tribool subset(const BoundedSetInterface<R>& s, const Box<R>& r);
      /*! \brief Tests if the box \a r is a superset of the set \a s. */
      friend tribool superset(const Box<R>& r, const BoundedSetInterface<R>& s);
      /*! \brief Tests if the box \a r is a superset of the set \a s. */
      friend Box<R> bounding_box(const BoundedSetInterface<R>& s);
#endif
    };

    //! \ingroup SetInterface 
    /*! Interface for closed sets for which both lower and upper descriptions are available. */
    template<class R>
    class ClosedSetInterface 
      : public virtual LowerClosedSetInterface<R>
      , public virtual UpperClosedSetInterface<R>
    { 
     public:
      virtual ClosedSetInterface<R>* clone() const = 0;
    };
      
    //! \ingroup SetInterface 
    /*! Interface for compact sets for which an upper description and a bounding box are available. */
    template<class R>
    class UpperCompactSetInterface 
      : public virtual UpperClosedSetInterface<R>
      , public virtual BoundedSetInterface<R>
    { 
     public:
      virtual UpperCompactSetInterface<R>* clone() const = 0;
    };
      
    //! \ingroup SetInterface 
    /*! Interface for closed sets for which both lower and upper descriptions, and a bounding box are available. */
    template<class R>
    class CompactSetInterface 
      : public virtual LowerClosedSetInterface<R>
      , public virtual UpperCompactSetInterface<R>
    { 
     public:
      virtual CompactSetInterface<R>* clone() const = 0;
    };
      
    //! \ingroup SetInterface 
    /*! \brief Interface for regular sets for which both lower and upper descriptions, and a bounding box are available. 
     *  A set is regular if its interior is the interior of its closure, and its closure is the closure of its interior.
     *  
     *  Note that a valid lower description of the closure can be obtained from a lower description of the interior.
     */
    template<class R>
    class RegularSetInterface 
      : public virtual OpenSetInterface<R>
      , public virtual ClosedSetInterface<R>
    {
     public:
      virtual RegularSetInterface<R>* clone() const = 0;
      virtual tribool intersects(const Box<R>& bx) const { return this->superset(bx); }
    };
      
    //! \ingroup SetInterface 
    /*! \brief Interface for regular sets for which both lower and upper descriptions, and a bounding box are available. */
    template<class R>
    class SetInterface 
      : public virtual RegularSetInterface<R>
      , public virtual CompactSetInterface<R>
    {
     public:
      typedef R real_type;
      virtual SetInterface<R>* clone() const = 0;
    };




    template<class R> inline tribool disjoint(const UpperClosedSetInterface<R>& s, const Box<R>& r) {
      return s.disjoint(r);
    }
    
    template<class R> inline tribool disjoint(const Box<R>& r, const UpperClosedSetInterface<R>& s) {
      return s.disjoint(r);
    }
    
    template<class R> inline tribool intersects(const UpperClosedSetInterface<R>& s, const Box<R>& r) {
      return s.intersects(r);
    }
    
    template<class R> inline tribool intersects(const Box<R>& r, const UpperClosedSetInterface<R>& s) {
      return s.intersects(r);
    }
    
    template<class R> inline tribool superset(const OpenSetInterface<R>& s, const Box<R>& r) {
      return s.superset(r);
    }
    
    template<class R> inline tribool subset(const Box<R>& r, const OpenSetInterface<R>& s) {
      return s.superset(r);
    }
    
    template<class R> inline tribool subset(const BoundedSetInterface<R>& s, const Box<R>& r) {
      return s.subset(r);
    }
    
    template<class R> inline tribool superset(const Box<R>& r, const BoundedSetInterface<R>& s) {
      return s.subset(r);
    }
    
    template<class R> inline Box<R> bounding_box(const BoundedSetInterface<R>& s) {
      return s.bounding_box();
    }
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const SetInterface<R>& s) {
      return s.write(os);
    }
    

  }
}


#endif /* ARIADNE_SET_INTERFACE_H */
