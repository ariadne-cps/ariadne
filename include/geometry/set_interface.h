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

//#include <iosfwd>

#include "base/types.h"
#include "base/tribool.h"


namespace Ariadne {
  

    template<class R> class Point;
    template<class R> class Box;
      
    //! \ingroup SetInterface 
    /*! %Base class for subsets of Euclidean space. */
    template<class BS>
    class SetBaseInterface 
    {
      typedef typename BS::real_type R;
      typedef typename BS::state_type Pt;
      typedef typename BS::space_type Spc;
     public:
      /*! \brief The type used to describe the space. */
      typedef Spc space_type;
      /*! \brief The type used to represent real numbers. */
      typedef R real_type;
      /*! \brief The type used to represent points. */
      typedef Pt state_type;
      /*! \brief The type used to represent basic sets. */
      typedef BS basic_set_type;

      /*! \brief Virtual destructor. */
      virtual ~SetBaseInterface() { }
      /*! \brief Make a dynamically-allocated copy. */
      virtual SetBaseInterface<BS>* clone() const = 0;
      /*! \brief The dimension of the space the set lies in. */
      virtual Spc space() const = 0;
      /*! \brief Test is the set contains a point. */
      virtual tribool contains(const Pt& pt) const = 0;
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const = 0;
#ifdef DOXYGEN
      /*! \brief Write to an output stream. */
      friend std::ostream& operator<<(std::ostream& os, const SetBaseInterface<BS>& s);
#endif
     protected:
      // Only base classes may create the interface
      SetBaseInterface() { }
    };

    //! \ingroup SetInterface 
    /*! \brief Interface for open sets. */
    template<class BS>
    class OpenSetInterface 
      : public virtual SetBaseInterface<BS>
    {
     public:
      virtual OpenSetInterface<BS>* clone() const = 0;
      /*! \brief Tests if an open set is a superset of a closed box. */
      virtual tribool superset(const BS&) const = 0;
#ifdef DOXYGEN
      /*! \brief Tests if the set \a s is a subset of the box \a r. */
      friend tribool superset(const OpenSetInterface<BS>& s, const BS& r);
      /*! \brief Tests if the box \a r is a subset of the set \a s. */
      friend tribool subset(const BS& r, const OpenSetInterface<BS>& s);
#endif 
    };

    //! \ingroup SetInterface 
    /*! \brief Interface for the closed sets described by the lower Fell topology. */
    template<class BS>
    class OvertSetInterface 
      : public virtual SetBaseInterface<BS>
    {
     public:
      virtual OvertSetInterface<BS>* clone() const = 0;
      /*! \brief Tests if a closed set intersects an open box. */
      virtual tribool intersects(const BS&) const = 0;
#ifdef DOXYGEN
      /*! \brief Tests if the set \a s is intersects the box \a r. */
      friend tribool intersect(const OvertSetInterface<BS>& s, const BS& r);
      /*! \brief Tests if the box \a r intersects the set \a s. */
      friend tribool intersect(const BS& r, const OvertSetInterface<BS>& s);
#endif
    };

    //! \ingroup SetInterface 
    /*! \brief Interface for the closed sets described by the upper Fell topology. */    
    template<class BS>
    class ClosedSetInterface 
      : public virtual SetBaseInterface<BS>
    {
     public:
      virtual ClosedSetInterface<BS>* clone() const = 0;
      /*! \brief Tests if a closed set is disjoint from a closed box. */
      virtual tribool disjoint(const BS&) const = 0;
#ifdef DOXYGEN
      /*! \brief Tests if the set \a s is disjoint from the box \a r. */
      friend tribool disjoint(const ClosedSetInterface<BS>& s, const BS& r);
      /*! \brief Tests if the box \a r is disjoint from the set \a s. */
      friend tribool disjoint(const BS& r, const ClosedSetInterface<BS>& s);
#endif
    };

    //! \ingroup SetInterface 
    /*! \brief Interface for bounded sets. */    
    template<class BS>
    class BoundedSetInterface 
      : public virtual SetBaseInterface<BS>
   {
     public:
      virtual BoundedSetInterface<BS>* clone() const = 0;
      /*! \brief Tests if a bouded set is a subset of from an open box. */
      virtual tribool subset(const BS&) const = 0;
      /*! \brief Returns a bounding box for the set. */
      virtual BS bounding_box() const = 0;
#ifdef DOXYGEN
      /*! \brief Tests if the set \a s is a subset of the box \a r. */
      friend tribool subset(const BoundedSetInterface<BS>& s, const BS& r);
      /*! \brief Tests if the box \a r is a superset of the set \a s. */
      friend tribool superset(const BS& r, const BoundedSetInterface<BS>& s);
      /*! \brief Tests if the box \a r is a superset of the set \a s. */
      friend BS bounding_box(const BoundedSetInterface<BS>& s);
#endif
    };

    //! \ingroup SetInterface 
    /*! Interface for compact sets for which an upper description and a bounding box are available. */
    template<class BS>
    class CompactSetInterface 
      : public virtual ClosedSetInterface<BS>
      , public virtual BoundedSetInterface<BS>
    { 
     public:
      virtual CompactSetInterface<BS>* clone() const = 0;
    };
      
    //! \ingroup SetInterface 
    /*! \brief Interface for regular sets for which both lower and upper descriptions, and a bounding box are available. 
     *  A set is regular if its interior is the interior of its closure, and its closure is the closure of its interior.
     *  
     *  Note that a valid lower description of the closure can be obtained from a lower description of the interior.
     */
    template<class BS>
    class RegularSetInterface 
      : public virtual OpenSetInterface<BS>
      , public virtual ClosedSetInterface<BS>
    {
     public:
      virtual RegularSetInterface<BS>* clone() const = 0;
      virtual tribool intersects(const BS& bx) const { return this->superset(bx); }
    };
      
    //! \ingroup SetInterface 
    /*! \brief Interface for regular sets for which both lower and upper descriptions, and a bounding box are available. */
    template<class BS>
    class SetInterface 
      : public virtual RegularSetInterface<BS>
      , public virtual CompactSetInterface<BS>
    {
     public:
      typedef typename BS::real_type real_type;
      virtual SetInterface<BS>* clone() const = 0;
    };

    /*! \ingroup SetInterface
     *  \brief Interface for sets in Euclidean space with boxes as basic sets. 
     */
    template< class R >
    class EuclideanSetInterface
      : public SetInterface< Box<R> >
    { };




    template<class BS> inline tribool disjoint(const ClosedSetInterface<BS>& s, const BS& r) {
      return s.disjoint(r);
    }
    
    template<class BS> inline tribool disjoint(const BS& r, const ClosedSetInterface<BS>& s) {
      return s.disjoint(r);
    }
    
    template<class BS> inline tribool intersects(const OvertSetInterface<BS>& s, const BS& r) {
      return s.intersects(r);
    }
    
    template<class BS> inline tribool intersects(const BS& r, const OvertSetInterface<BS>& s) {
      return s.intersects(r);
    }
    
    template<class BS> inline tribool superset(const OpenSetInterface<BS>& s, const BS& r) {
      return s.superset(r);
    }
    
    template<class BS> inline tribool subset(const BS& r, const OpenSetInterface<BS>& s) {
      return s.superset(r);
    }
    
    template<class BS> inline tribool subset(const BoundedSetInterface<BS>& s, const BS& r) {
      return s.subset(r);
    }
    
    template<class BS> inline tribool superset(const BS& r, const BoundedSetInterface<BS>& s) {
      return s.subset(r);
    }
    
    template<class BS> inline BS bounding_box(const BoundedSetInterface<BS>& s) {
      return s.bounding_box();
    }
    
    template<class BS> inline
    std::ostream& operator<<(std::ostream& os, const SetInterface<BS>& s) {
      return s.write(os);
    }
    

} // namespace Ariadne


#endif /* ARIADNE_SET_INTERFACE_H */
