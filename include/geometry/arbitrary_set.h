/***************************************************************************
 *            arbitrary_set.h
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
 
/*! \file set.h
 *  \brief General sets.
 */

#ifndef _ARIADNE_ARBITRARY_SET_H
#define _ARIADNE_ARBITRARY_SET_H

#include <iosfwd>

#include "set.h"


namespace Ariadne {
  namespace Geometry {


    //! \ingroup ExactSet
    /*! \brief A class which can store an object of any Set type. 
     *
     * Useful in containers or when returning a Set object from a function. 
     * 
     * TODO: Make thread safe
     */
    template<class R>
    class ArbitrarySet
      : public Set<R>
    {
     public:
      typedef R real_type;
      typedef Point<R> state_type;
     
      /*! \brief Construct by making a clone of a given set. */
      ArbitrarySet(const Set<R>& set)
        : _ref_count(new int(1)), _ptr(set.clone()) { }
      
      /*! \brief Assign from another set by making a clone. */
      ArbitrarySet<R>& operator=(const Set<R>& set) { 
        Set<R>* new_ptr=set.clone();
        --*this->_ref_count; 
        if(*this->_ref_count==0) { delete this->_ptr; *this->_ref_count=1; } else { this->_ref_count=new int(1); }
        this->_ptr=new_ptr; 
        return *this;
      }

      /*! \brief Copy constructor. Stores a reference; does not make a clone of the set. */
      ArbitrarySet(const ArbitrarySet<R>& set)
        : _ref_count(set._ref_count), _ptr(set._ptr) { ++this->_ref_count; }

      /*! \brief Copy assignment operator.  Stores a reference; does not make a clone of the set. */
      ArbitrarySet<R>& operator=(const ArbitrarySet<R>& set) { 
        if(&set!=this) { 
          --*this->_ref_count; 
          if(*this->_ref_count==0) { delete this->_ref_count; delete this->_ptr; } 
          this->_ref_count=set._ref_count;
          this->_ptr=set._ptr;
          ++*this->_ref_count;
        }
        return *this;
      }

      virtual ~ArbitrarySet()  { --*this->_ref_count; if(this->_ref_count==0) { delete this->_ref_count; delete this->_ptr; } }
      virtual Set<R>* clone() const { return this->_ptr->clone(); }
      virtual dimension_type dimension() const { return this->_ptr->dimension(); }
      virtual tribool contains(const Point<R>& pt) const { return this->_ptr->contains(pt); }
      virtual tribool disjoint(const Rectangle<R>& r) const { return this->_ptr->disjoint(r); }
      virtual tribool superset(const Rectangle<R>& r) const { return this->_ptr->superset(r); }
      virtual tribool subset(const Rectangle<R>& r) const { return this->_ptr->subset(r); }
      virtual Rectangle<R> bounding_box() const { return this->_ptr->bounding_box(); }
      virtual std::ostream& write(std::ostream& os) const { return this->_ptr->write(os); }
     private:
      int* _ref_count;
      Set<R>* _ptr;
    };
    
  }
}


#endif /* _ARIADNE_ARBITRARY_SET_H */
