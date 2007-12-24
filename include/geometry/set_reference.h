/***************************************************************************
 *            set_reference.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
/*! \file set_reference.h
 *  \brief Arbitrary sets.
 */

#ifndef ARIADNE_SET_REFERENCE_H
#define ARIADNE_SET_REFERENCE_H

#include <iosfwd>
#include <boost/shared_ptr.hpp>

#include "set_interface.h"
#include "rectangle.h"
#include "rectangular_set.h"
#include "polyhedral_set.h"


namespace Ariadne {
  namespace Geometry {


    //! \ingroup ExactSet
    /*! \brief A class which can store an object of any SetInterface type. 
     *
     * Useful in containers or when returning a SetInterface object from a function. 
     * 
     * TODO: Make thread safe
     */
  /*
    template<class R>
    class SetReference
    {
     public:
      typedef R real_type;
      typedef Point<R> state_type;
     
      //! \brief Construct by making a clone of a given set.
      SetReference(const SetInterface<R>& set)
        : _ref_count(new int(1)), _ptr(set.clone()) { }
      
      //! \brief Assign from another set by making a clone.
      SetReference<R>& operator=(const SetInterface<R>& set) { 
        SetInterface<R>* new_ptr=set.clone();
        --*this->_ref_count; 
        if(*this->_ref_count==0) { delete this->_ptr; *this->_ref_count=1; } else { this->_ref_count=new int(1); }
        this->_ptr=new_ptr; 
        return *this;
      }

      //! \brief Copy constructor. Stores a reference; does not make a clone of the set.
      SetReference(const SetReference<R>& set)
        : _ref_count(set._ref_count), _ptr(set._ptr) { ++this->_ref_count; }

     //! \brief Copy assignment operator.  Stores a reference; does not make a clone of the set.
      SetReference<R>& operator=(const SetReference<R>& set) { 
        if(&set!=this) {  
          --*this->_ref_count; 
          if(*this->_ref_count==0) { delete this->_ref_count; delete this->_ptr; } 
          this->_ref_count=set._ref_count;
          this->_ptr=set._ptr;
          ++*this->_ref_count;
        }
        return *this;
      }
      
      operator const SetInterface<R>& () const { return *this->_ptr; }

      ~SetReference()  { --*this->_ref_count; if(this->_ref_count==0) { delete this->_ref_count; delete this->_ptr; } }
      SetInterface<R>* clone() const { return this->_ptr->clone(); }
      dimension_type dimension() const { return this->_ptr->dimension(); }
      tribool contains(const Point<R>& pt) const { return this->_ptr->contains(pt); }
      tribool disjoint(const Box<R>& r) const { return this->_ptr->disjoint(r); }
      tribool superset(const Box<R>& r) const { return this->_ptr->superset(r); }
      tribool subset(const Box<R>& r) const { return this->_ptr->subset(r); }
      Box<R> bounding_box() const { return this->_ptr->bounding_box(); }
      std::ostream& write(std::ostream& os) const { return this->_ptr->write(os); }
     private:
      int* _ref_count;
      SetInterface<R>* _ptr;
    };
*/

    //! \ingroup ExactSet
    /*! \brief A reference to a set with automatic memory management
     *
     * Useful in containers or when returning a SetInterface object from a function.
     * 
     * Note that this class is derived from SetInterface<R>. This may seem unnecessary
     * since there is a conversion operator to const SetInterface<R>& . However, due
     * to the face that the class is a template, C++ is unable to automatically make
     * the conversion to an ordinary reference when passing the class to a function.
     * Hence the requirement that the class be a template.
     *  
     * \internal Maybe we don't need this class to be derived from SetInterface.
     */
    template<class R>
    class SetReference
      : public SetInterface<R>
    {
     private:
      boost::shared_ptr< const SetInterface<R> > _ptr;
     public:
      typedef R real_type;
      typedef Point<R> state_type;
     
      //! \brief Copy constructor.  Stores a reference; does not make a clone of the set.
      SetReference(const SetReference<R>& set) : _ptr(set._ptr) { }
      //! \brief Copy assignment operator.  Stores a reference; does not make a clone of the set.
      SetReference<R>& operator=(const SetReference<R>& set) { this->_ptr=set._ptr; return *this; }

      //! \brief Construct by making a clone of a given set. 
      //! 
      //! Note that a copy of the reference is constructed; the argument is not deleted when the set goes out of scope. 
      //! This is unlike the Boost shared pointers, which take a dynamically-allocated object.
      SetReference(const SetInterface<R>& set) : _ptr() { (*this)=set; }
      //! \brief Assign from another set by making a clone.
      SetReference<R>& operator=(const SetInterface<R>& set) { 
        const SetReference<R>* set_ref=dynamic_cast<const SetReference<R>*>(&set); 
        if(set_ref) { this->_ptr=set_ref->_ptr; } 
        else { this->_ptr=boost::shared_ptr< const SetInterface<R> >(set.clone()); }
        return *this; 
      }

      //! \brief Construct by making a copy of a Box. 
      SetReference(const Box<R>& r) : _ptr(new RectangularSet<R>(r)) { }

      //! \brief Construct by making a copy of a Box. 
      SetReference(const Polyhedron<R>& p) : _ptr(new PolyhedralSet<R>(p)) { }

      //! \brief Convert to an ordinary reference. 
      operator SetInterface<R>& () { return *this->_ptr; }

      //! \brief Convert to an ordinary const reference. 
      operator const SetInterface<R>& () const { return *this->_ptr; }

      //! \brief Destructor. 
      ~SetReference() { }
      //! \brief Create a dynamically-allocated copy of the set.
      SetInterface<R>* clone() const { return this->_ptr->clone(); }
      dimension_type dimension() const { return this->_ptr->dimension(); }
      tribool contains(const Point<R>& pt) const { return this->_ptr->contains(pt); }
      tribool superset(const Box<R>& r) const { return this->_ptr->superset(r); }
      tribool intersects(const Box<R>& r) const { return this->_ptr->intersects(r); }
      tribool disjoint(const Box<R>& r) const { return this->_ptr->disjoint(r); }
      tribool subset(const Box<R>& r) const { return this->_ptr->subset(r); }
      tribool bounded() const { return this->_ptr->bounded(); }
      Box<R> bounding_box() const { return this->_ptr->bounding_box(); }
      std::ostream& write(std::ostream& os) const { return this->_ptr->write(os); }
    };

    template<class R> 
    std::ostream& operator<<(std::ostream& os, const SetReference<R>& set) {
      return set.write(os);
    }

  }
}


#endif /* ARIADNE_SET_REFERENCE_H */
