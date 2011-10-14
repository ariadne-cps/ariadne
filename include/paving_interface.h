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
class Grid;
class GridCell;

#define ARIADNE_ABSTRACT_METHOD { throw std::runtime_error(std::string("ERROR: Unimplemented abstract method ")+__PRETTY_FUNCTION__); }

// A continuous predicate taking values in a three-valued logic.
// NOTE: Corresponds to your SetCheckerInterface; I think that this name is better.
template<class BS> class PredicateInterface {
  public:
    virtual tribool check(const BS&) = 0;
};

template<class T> class ForwardConstantIteratorHandle;

template<class T> class ForwardConstantIteratorInterface
{
    friend class ForwardConstantIteratorHandle<T>;
  private:
    virtual ForwardConstantIteratorInterface<T>* clone() const = 0;
    virtual void increment() const = 0;
    virtual const T dereference() const = 0;
    virtual bool equals(const ForwardConstantIteratorInterface<T>&) const = 0;
};

//! \brief A generic forward iterator through constant data.
//! \details Since the iterator is through constant data, the data may be safely be returned by value rather than by reference,
//! though this is not strictly standards-conforming, and means that operator->() cannot be provided.
template<class T> class ForwardConstantIteratorHandle
{
    ForwardConstantIteratorInterface<T>* _ptr;
  public:
    ~ForwardConstantIteratorHandle() { delete _ptr; }
    ForwardConstantIteratorHandle(ForwardConstantIteratorInterface<T>* p) : _ptr(p) { }
    ForwardConstantIteratorHandle(const ForwardConstantIteratorInterface<T>& p) : _ptr(p.clone()) { }
    ForwardConstantIteratorHandle(const ForwardConstantIteratorHandle<T>& other) : _ptr(other._ptr->clone()) { }
    ForwardConstantIteratorHandle<T>& operator=(const ForwardConstantIteratorHandle<T>& other) {
        if(_ptr!=other._ptr) { delete _ptr; _ptr=other._ptr->clone(); } return *this; }

    const T operator*() const { return _ptr->dereference(); }
    ForwardConstantIteratorHandle<T>& operator++() { _ptr->increment(); return *this; }
    ForwardConstantIteratorHandle<T> operator++(int) { ForwardConstantIteratorHandle<T> tmp(*this); this->_ptr->increment(); return tmp; }
    bool operator==(const ForwardConstantIteratorHandle<T>& other) const { return this->_ptr->equals(*other._ptr); }
    bool operator!=(const ForwardConstantIteratorHandle<T>& other) const { return !this->_ptr->equals(*other._ptr); }
};


//! \brief A prototype interface for sets
//! \details Since the iterator is through constant data, the data may be safely be returned by value rather than by reference,
//! though this is not strictly standards-conforming, and means that operator->() cannot be provided.
template<class BS> class DenotableSetInterface
{
  public:
    inline ForwardConstantIteratorHandle<BS> begin() const;
    inline ForwardConstantIteratorHandle<BS> end() const;
    virtual void adjoin(const BS&) = 0;
    virtual void adjoin(const DenotableSetInterface<BS>&) = 0;
  protected:
    virtual ForwardConstantIteratorInterface<BS>* _begin() const = 0;
    virtual ForwardConstantIteratorInterface<BS>* _end() const = 0;
};

class SubPavingHandle;

//! \brief A prototype interface for subpavings.
//! \details A sub-paving is a paving-type set which is the restriction of another paving to a cell, called the \em root_cell of the subpaving.
//! In general, a sub-paving supports the non-modifying operations of a paving.
//!
//! The only modifying operations on a sub-paving which are logically allowable are those which apply to sub-cells of the root cell.
//! These operations modify the paving itself.
//!
//! \internal
//! NOTE: One difficulty with sub-pavings is that returning a sub-paving requires using either a raw pointer to SubPavingInterface with corresponding memory issues,
//! a shared pointer, which requires using pointer syntax, or a handle, which requires forwarding functions but is cleanest for users.
//! NOTE: Since operations modifying a sub-paving modify the original paving, we need to be careful with constantness, in particular when making copies.
//! This may be an issue in the branch method or when implementing Handle classes.
class SubPavingInterface
    : public virtual LocatedSetInterface
{
  public:
    typedef ForwardConstantIteratorHandle<GridCell> const_iterator;
  public:
    //! \brief A constant iterator through the cells of the paving.
    inline ForwardConstantIteratorHandle<GridCell> begin() const { return this->_begin(); }
    inline ForwardConstantIteratorHandle<GridCell> end() const { return this->_end(); }
    //! \brief A sub-paving obtained by considering one of the immediate subcells of the root cell.
    inline SubPavingHandle branch(bool left_or_right);
  public:
    //! \brief The underlying grid of the paving.
    virtual const Grid& grid() const = 0;

    //! \brief The root cell of the (sub)paving. All enabled cells are subcells of the root cell.
    virtual GridCell root_cell() const = 0;
    //! \brief Tests if the cell \a c is contained in one of the cells of the subpaving.
    virtual bool superset(const GridCell& c) const = 0; // Could also be called "contains"

    //! \brief Tests if the paving is a subset of another paving on the same Grid.
    virtual bool subset(const SubPavingInterface&) const { ARIADNE_ABSTRACT_METHOD; }
    //! \brief Tests if the paving intersects another paving on the same Grid.
    //! Tests intersection of interiors, so returns false if the pavings touch but do not overlap.
    virtual bool intersects(const SubPavingInterface&) const { ARIADNE_ABSTRACT_METHOD; }

    //! \brief Modify the underlying paving so that the root cell of the subpaving is either enabled or disabled.
    virtual void set_root_cell(bool onoff) = 0;

    virtual Box bounding_box() const = 0; // Inherited from CompactSetInterface
    virtual tribool inside(const Box& bx) const = 0; // Inherited from CompactSetInterface
    virtual tribool separated(const Box& bx) const = 0; // Inherited from ClosedSetInterface
    virtual tribool overlaps(const Box& bx) const = 0; // Inherited from OvertSetInterface
    virtual std::ostream& write(std::ostream& os) const = 0;
  protected:
    virtual SubPavingInterface* _branch(bool left_or_right) { ARIADNE_ABSTRACT_METHOD; } // Could also be called "child" or "subpaving"
    virtual ForwardConstantIteratorInterface<GridCell>* _begin() const { ARIADNE_ABSTRACT_METHOD; }
    virtual ForwardConstantIteratorInterface<GridCell>* _end() const { ARIADNE_ABSTRACT_METHOD; }
};

class PavingInterface;

// Stub concrete class implementing SubPavingInterface by forwarding all functions onto top-level paving.
class SubPaving
{
    PavingInterface* _base_paving;
    // GridCell _root_cell;
    // void* _paving_subtree_data;
  public:
    // Public interface here
};


//! \brief A paving, which is a set constructed as a finite union of non-overlapping cells, where the possible cells are obtained by successive subdivision.
class PavingInterface
    : public virtual SubPavingInterface
{
  public:
    virtual PavingInterface* clone() const = 0;
    virtual GridCell smallest_enclosing_primary_cell(const Box& bx) const = 0; // Useful query, but can also be implemented at the Grid level.
    virtual void adjoin_cells(const PredicateInterface<Box>& predicate, const uint depth) { ARIADNE_ABSTRACT_METHOD; }
    virtual void adjoin_outer_approximation(const CompactSetInterface& set, const uint depth) = 0;
    virtual void adjoin_inner_approximation(const OpenSetInterface& set, const uint height, const uint depth) = 0;
    virtual void adjoin(const SubPavingInterface& paving) = 0;
    virtual void restrict(const SubPavingInterface& paving) = 0;
    virtual void remove(const SubPavingInterface& paving) = 0;
    virtual void adjoin(const GridCell& cell) = 0; // Deprecated?
};


} // namespace Ariadne


#endif // ARIADNE_PAVING_INTERFACE_H
