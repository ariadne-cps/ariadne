/***************************************************************************
 *            polyhedron.h
 *
 *  Copyright  2005-8  Alberto Casagrande, Pieter Collins
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

/*! \file polyhedron.h
 *  \brief Polyhedra.
 */
 
#ifndef ARIADNE_POLYHEDRON_H
#define ARIADNE_POLYHEDRON_H

#include <iosfwd>
#include <vector>

#include "tribool.h"

#include "vector.h"
#include "matrix.h"

#include "function.h"

#include "set_interface.h"
#include "function_set.h"

namespace Ariadne {  
  

class Box;
class Polytope;
class Polyhedron;
    



/*! \ingroup BasicSet
 *  \brief A polyhedron (not necessarily bounded polyhedral set) described by a system of linear inequalities.
 *
 *  The set is described as
 *  \f$ \{ x\in\mathbb{R}^d \mid Ax \leq b \} \f$
 *  where \f$A\f$ is a \f$n\times d\f$ matrix and \f$b\f$ is a vector of size \f$n\f$.
 */ 
class Polyhedron 
    : public ConstraintSet
{
  public:
    //@{ 
    //! \name Constructors and destructor

    //! \brief Default constructor constructs a polytope in zero dimensions with no constraints.
    explicit Polyhedron();
  
    //! \brief Construct full Euclidean space of dimension \a n.
    explicit Polyhedron(uint n=0u);
  
    //! \brief Construct a polyhedron of dimension \a d with \a nc constraints from the data in the
    //! array beginning at \a data. The jth element of the ith constraint is stored in position i*(d+1)+j, 
    //! and the ith inhomogeneous term is stored in position i*(d+1)+d.
    template<class XX> Polyhedron(uint d, uint nc, const XX* data);
  
    //! \brief Construct the polyhedron defined by the matrix equations \f$Ax\leq b\f$.
    explicit Polyhedron(const Matrix<Float>& A, const Vector<Float>& b);
  
    //! \brief Convert from a box. 
    explicit Polyhedron(const Box& bx);
  
    //! \brief Convert from a polytope. 
    explicit Polyhedron(const Polytope& p);
  
    //@}
  
  
    //@{
    //! \name Data access

    //! \brief The number of constraints. 
    size_t number_of_constraints() const { return this->b().size(); }
    //! \brief The matrix \f$A\f$ in the inequalities \f$Ax\leq b\f$. 
    Matrix<Float> A() const { return -dynamic_cast<const AffineFunction&>(this->function()).A(); }
    //! \brief The vector \f$b\f$ in the inequalities \f$Ax\leq b\f$. 
    Vector<Float> b() const { return dynamic_cast<const AffineFunction&>(this->function()).b(); }
    //@}
  
  
    //@{
    //! \name Geometric operations

    //! \brief Tests if the polyhedron is empty. (Not currently implemented.)
    virtual tribool empty() const;

    //! \brief Tests if the polyhedron is bounded. (Not currently implemented.)
    virtual tribool bounded() const;

    //! \brief The \a i<sup>th</sup> defining halfspace. 
    Polyhedron halfspace(size_t i) const;

    //! \brief Convert to a polytope. 
    operator Polytope () const ;

    //! \brief The intersection of two polyhedra. 
    friend Polyhedron intersection(const Polyhedron& p1, const Polyhedron& p2);

    //! \brief Convert a box to a polyhedron. 
    friend Polyhedron polyhedron(const Box& p);

    //! \brief Convert a polytope to a polyhedron. (Not currently implemented) 
    friend Polyhedron polyhedron(const Polytope& p);

    //! \brief Convert a polyhedron to a polytope. (Not currently implemented) 
    friend Polytope polytope(const Polyhedron& p);

    //@}

    //@{
    //! \name Input/output.

    //! \brief Write to an output stream. 
    std::ostream& write(std::ostream& os) const;
    //@}

};

  
inline std::ostream& operator<<(std::ostream& os, const Polyhedron& p) {
    return p.write(os);
}

std::ostream& operator<<(std::ostream& os, const Polytope& p);


} // namespace Ariadne

#endif // ARIADNE_POLYHEDRON_H
