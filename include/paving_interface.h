/***************************************************************************
 *            paving_interface.h
 *
 *  Copyright 2008  Pieter Collins
 * 
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
 
/*! \file paving_interface.h
 *  \brief Storage sets based on pavings and covers
 */

#ifndef ARIADNE_PAVING_INTERFACE_H
#define ARIADNE_PAVING_INTERFACE_H

#include <iosfwd>
#include "set_interface.h"

namespace Ariadne {

class Box;
    

template<class T> class ForwardConstIteratorInterface
{
  private:
    virtual ForwardConstIteratorInterface<T>* clone() const = 0;
    virtual void increment() const = 0;
    virtual T dereference() const = 0;
    virtual bool equals(const ForwardConstIteratorInterface<T>&) const = 0;  
};

template<class T> class ForwardConstIteratorGeneric
{
    ForwardConstIteratorInterface<T>* _ptr;
  public:
    ~ForwardConstIteratorGeneric() { delete _ptr; }
    ForwardConstIteratorGeneric(ForwardConstIteratorInterface<T>* p) : _ptr(p) { }
    ForwardConstIteratorGeneric(const ForwardConstIteratorInterface<T>& p) : _ptr(p.clone()) { }
    ForwardConstIteratorGeneric(const ForwardConstIteratorGeneric<T>& other) : _ptr(other._ptr->clone()) { }
    ForwardConstIteratorGeneric<T>& operator=(const ForwardConstIteratorGeneric<T>& other) { 
        if(_ptr!=other._ptr) { delete _ptr; _ptr=other._ptr->clone(); } return *this; }

    T operator*() const { return _ptr->dereference(); }
    ForwardConstIteratorGeneric<T>& operator++() { _ptr->increment(); return *this; }
    ForwardConstIteratorGeneric<T> operator++(int) { ForwardConstIteratorGeneric<T> tmp(*this); this->_ptr->increment(); return tmp; }
    bool operator==(const ForwardConstIteratorInterface<T>& other) const { return this->_ptr->equals(*other._ptr); }
    bool operator!=(const ForwardConstIteratorInterface<T>& other) const { return !this->_ptr->equals(*other._ptr); }
};


template<class SET> class PavingInterface
    : public LocatedSetInterface
{
  public:
    typedef ForwardConstIteratorGeneric<SET> const_iterator;
  public:
    ForwardConstIteratorGeneric<SET> begin() const { return this->_begin(); }
    ForwardConstIteratorGeneric<SET> end() const { return this->_end(); }
  public:
    virtual PavingInterface<SET>* clone() const = 0;
    virtual Box bounding_box() const = 0;
    virtual tribool disjoint(const Box& bx) const = 0;
    virtual tribool overlaps(const Box& bx) const = 0;
    virtual std::ostream& write(const std::ostream& os) const = 0;
    virtual void adjoin_outer_approximation(const CompactSetInterface& set, uint depth) = 0;
    virtual void adjoin_inner_approximation(const OpenSetInterface& set, uint depth) = 0;
    virtual void mince(uint depth) = 0;
    virtual void recombine() = 0;
  protected:
    virtual ForwardConstIteratorInterface<SET>* _begin() const = 0;
    virtual ForwardConstIteratorInterface<SET>* _end() const = 0;
        
};

template<class SET> class CoverInterface
    : public RegularSetInterface
{
  public:
    typedef ForwardConstIteratorGeneric<SET> const_iterator;
  public:
    ForwardConstIteratorGeneric<SET> begin() const { return this->_begin(); }
    ForwardConstIteratorGeneric<SET> end() const { return this->_end(); }
  public:
    virtual PavingInterface<SET>* clone() const = 0;
    virtual tribool subset(const Box& bx) const = 0;
    virtual tribool superset(const Box& bx) const = 0;
    virtual std::ostream& write(const std::ostream& os) = 0;
    virtual void adjoin_lower_approximation(const OvertSetInterface& set, uint depth) = 0;
    virtual void adjoin_over_approximation(const CompactSetInterface& set, uint depth) = 0;
    virtual void mince(uint depth) = 0;
    virtual void recombine() = 0;
  protected:
    virtual ForwardConstIteratorInterface<SET>* _begin() const = 0;
    virtual ForwardConstIteratorInterface<SET>* _end() const = 0;
        
};


} // namespace Ariadne


#endif // ARIADNE_PAVING_INTERFACE_H
