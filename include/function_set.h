/***************************************************************************
 *            function_set.h
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
 
/*! \file function_set.h
 *  \brief Images and preimages of boxes in Euclidean space.
 */

#ifndef ARIADNE_FUNCTION_SET_H
#define ARIADNE_FUNCTION_SET_H

#include <iosfwd>

#include "boost/shared_ptr.hpp"
#include "macros.h"
#include "numeric.h"
#include "vector.h"
#include "set_interface.h"
#include "function.h"
#include "graphics_interface.h"

#include "box.h"

namespace Ariadne {

class Zonotope;
class Polyhedron;

//! \brief A set defined as the image of a box under a continuous function.
class ImageSet
    : public LocatedSetInterface
    , public DrawableInterface
{
    Vector<Interval> _domain;
    VectorFunction _function;
  public:
    //! \brief Default constructor constructs the singleton in \f$\R^0\f$.
    ImageSet();
    //! \brief Construct the image of \a dom under the identity function.
    ImageSet(const Vector<Interval>& dom);
    //! \brief Construct the image of \a dom under the function \a fn.
    ImageSet(const Vector<Interval>& dom, const VectorFunction& fn);
    //! \brief The box used to define the set.
    const Vector<Interval>& domain() const { return this->_domain; }
    //! \brief The function used to define the set.
    const VectorFunction& function() const { return this->_function; }
    //! \brief Equality operator. Compares functions by referential equality.
    bool operator==(const ImageSet& ims) const {
        return this->_domain==ims._domain && this->_function.pointer()==ims._function.pointer(); }

    ImageSet* clone() const;
    uint dimension() const;
    tribool empty() const;
    tribool disjoint(const Box&) const;
    tribool overlaps(const Box&) const;
    tribool inside(const Box&) const;
    Box bounding_box() const;
    void draw(CanvasInterface&) const;
    std::ostream& write(std::ostream&) const;
};

template<class Mdl>
class ModelSet
    : public CompactSetInterface
{
    Mdl _model;
  public:
    ModelSet(const Mdl& mdl) : _model(mdl) { }
    const Vector<Interval>& domain() const { return this->_model.domain(); }
    ModelSet* clone() const { return new ModelSet<Mdl>(*this); }
    uint dimension() const { return this->_model.result_size(); }
    tribool disjoint(const Box& bx) const { 
        return Ariadne::disjoint(this->_model,bx); }
    tribool inside(const Box& bx) const { 
        return Ariadne::inside(this->_model.range(),bx) || indeterminate; }
    Box bounding_box() const { 
        return this->_model.range(); }
    std::ostream& write(std::ostream& os) const {
        return os << "ModelSet( " << this->_model << ")"; }
};
 


//! \brief A set defined as the preimage of a box (the \em codomain) under a continuous function. The set is described as \f$S=f^{-1}(B) = \{ x \mid f(x)\in B\}\f$ where \f$B\f$ is the codomain and \f$f\f$ the function.
class ConstraintSet
    : public RegularSetInterface
{
    Vector<Interval> _codomain;
    VectorFunction _function;
  public:
    //! \brief Construct the preimage of \a codom under \a fn.
    ConstraintSet(const Vector<Interval>& codom, const VectorFunction& fn);
    //! \brief Construct from a polyhedron.
    ConstraintSet(const Polyhedron& p);
    //! \brief The codomain of the set.
    const Vector<Interval>& codomain() const { return this->_codomain; }
    //! \brief The function used to define the set.
    const VectorFunction& function() const { return this->_function; };

    ConstraintSet* clone() const;
    uint dimension() const;
    tribool disjoint(const Box&) const;
    tribool overlaps(const Box&) const;
    tribool covers(const Box&) const;
    std::ostream& write(std::ostream&) const;
};

} // namespace Ariadne 

#endif
