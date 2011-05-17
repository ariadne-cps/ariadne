/***************************************************************************
 *            enclosure.h
 *
 *  Copyright  2011  Pieter Collins
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

/*! \file enclosure.h
 *  \brief Enclosure sets for continuous systems
 */

#ifndef ARIADNE_ENCLOSURE_H
#define ARIADNE_ENCLOSURE_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "logging.h"
#include "tribool.h"
#include "container.h"
#include "graphics_interface.h"

namespace Ariadne {

class Interval;
template<class X> class Vector;
typedef Vector<Interval> IntervalVector;
class Box;
class GridTreeSet;
template<class X> class ScalarFunction;
template<class X> class VectorFunction;
typedef ScalarFunction<Interval> IntervalScalarFunction;
typedef VectorFunction<Interval> IntervalVectorFunction;

class EnclosureInterface
{
    virtual EnclosureInterface* _clone() const = 0;

};

//! \brief A class representing an enclosure for a hybrid evolution.
//! Handles progress, activation and guard constraints internally.
//! The set is represented as the image of a box \f$D\f$ under a function model \f$\hat{f}(s)\f$, under the constraints
//! \f$\hat{c}(s) \leq 0\f$ and \f$\hat{e}(s)=0\f$. Also keeps track of the current time \f$\hat{t}(s)\f$.
//!
//! In other words,
//! \f[ S=\{ \hat{f}(s);\  \hat{t}(s) \mid s\in D \mid \hat{c}(s) \leq 0 \ \wedge \hat{e}(s)=0 \} . \f]
//! In the following documentation, we sometimes write \f$\xi(s)\f$ for \f$\hat{f}(s)\f$ and \f$\tau(s)\f$ for \f$\hat{t}(s)\f$.
//! yielding \f[ S=\{ \xi(s);\  \tau(s) \mid s\in D \mid \rho(s) \in C \} . \f]
//!
//! <b>Rationale:</b><\par>
//! <b>Parameterisation of the time specifier</b>
//! When computing a flow step, the evolution time can be given either as a
//! function \f$\varepsilon(x)\f$ of the state, or as an absolute final time
//! \f$\omega(s)\f$ of the parameters.
//! The rationale for this difference is that we often use the
//! \f$\varepsilon\f$ form when computing the flow up to a crossing with a
//! guard, in which case the time depends purely on the state, and it may
//! be possible to simplify the constraints by expressing them in terms of
//! the current states (in an implementation in which intermediate states are
//! considered explicitly).
//! The rationale for using a parameterised form with independent variable
//! \f$s\f$ is that this is needed to distinguished different states based on
//! crossings with the guard sets. The rationale for using a final time
//! \f$\omega(s)\f$ as a final time is that this allows us to remove
//! errors in the elapsed time \f$\tau(s)\f$.
//!
//! A possible extension would be to allow \f$\omega(s)\f$ to depend on
//! intermediate variables \f$x_i\f$ at previous stages. The reason for not
//! doing this is that it would enforce a policy for tracking intermediate
//! states, which we do not want to insist on during an implementation.
//!
//! <b>Replacing events</b>
//! If \f$g_e\f$ is increasing along trajectories
//! it is safe to remove a constraint \f$g_e(x_{i-1})\leq0\f$ and
//! replace it with the constraint \f$g_e(x_i)\leq 0\f$.
class Enclosure
    : public Loggable
{
  private:
    friend class EnclosureFactory;
    EnclosureInterface* _handle;
  public:
    //! \brief The dimension of the set.
    uint dimension() const;
    void apply_map(IntervalVectorFunction f);
    void new_negative_constraint(IntervalScalarFunction g);
    void new_zero_constraint(IntervalScalarFunction h);
    void new_parameter(Interval d);
    void new_variable(Interval d);

    inline Enclosure() : _handle(0) { }

    //! The domain of the parameters describing the set.
    IntervalVector domain() const;
    //! \brief Returns a bounding box for the set. Computed by a simple interval evaluation of \f$f(D)\f$.
    Box bounding_box() const;
    //! \brief Tests whether the set is empty.
    tribool empty() const;
    //! \brief Tests whether the set is disjoint from the box \a hbx.
    tribool disjoint(const Box& bx) const;
    //! \brief Tests whether the set is a subset of the box \a hbx.
    tribool subset(const Box& hbx) const;
    //! \brief Adjoins an outer approximation of the set to the grid-based set \a paving, with accuracy given by
    //! \a depth subdivisions in each component.
    void adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const;

    //! \brief Simplifies the representation.
    void recondition();

    //! \brief Splits into smaller subsets.
    List<Enclosure> split() const;

    //! \brief Draws onto a canvas.
    void draw(CanvasInterface&) const;
  private:
    inline Enclosure(EnclosureInterface* e) : _handle(e) { }
};

class EnclosureFactory {
    virtual EnclosureInterface* _create(const Box& bx) const = 0;
  public:
    inline Enclosure create(const Box& bx) const { return Enclosure(this->_create(bx)); }
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_ENCLOSURE_H
