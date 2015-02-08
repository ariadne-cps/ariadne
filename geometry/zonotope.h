/***************************************************************************
 *            zonotope.h
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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

/*! \file zonotope.h
 *  \brief Zonotopes in Euclidean space.
 */

#ifndef ARIADNE_ZONOTOPE_H
#define ARIADNE_ZONOTOPE_H

#include <iosfwd>

#include "utility/declarations.h"
#include "utility/tribool.h"

#include "algebra/vector.h"
#include "algebra/matrix.h"

#include "geometry/set_interface.h"

#include "output/graphics_interface.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;

class Zonotope;
template<class BS> class ListSet;

template<class X> class Affine;

class Figure;


/*!\brief A zonotope of arbitrary dimension.
 *
 * A zonotope is a set of the form \f$c+Ge\f$, where \f$||e||_{\infty}\leq1\f$.
 * The columns of the matrix \f$G\f$ are the <em>generators</em> of the
 * zonotope.
 *
 * Zonotopes are always bounded.
 * A zonotope always contains its centre point, so can never be empty.
 * However, it may not be regular.
 *
 *
 * The intersection and membership tests may be performed using algorithms from: <br>
 * Guibas, Leonidas J.; Nguyen, An; Zhang, Li, "Zonotopes as bounding volumes."
 * <i>Proceedings of the Fourteenth Annual ACM-SIAM Symposium on Discrete Algorithms</i>
 * (Baltimore, MD, 2003),  803--812, ACM, New York, 2003.
 *
 * \b Storage: A %Zonotope in dimension d with n generators is described by
 * d(n+1) real numbers. The ith component of the centre is given by the
 * ith element, and the ith component of the kth generator is given by the
 * (kd+i)th element.
 */


class Zonotope
    : public CompactSetInterface
    , public DrawableInterface
{
  private:
    Vector<Float64> _centre;
    Matrix<Float64> _generators;
    Vector<Float64> _error;
  public:
    //@{
    //! \name Constructors and destructors
    /*! \brief Virtual destructor. */
    virtual ~Zonotope();
    /*! \brief Default constructor yields a zonotope with dimension zero and no generators. */
    explicit Zonotope();
    /*! \brief Construct a zonotope of dimension \a d with no generators. */
    explicit Zonotope(Nat d);
    /*! \brief Construct a zonotope of dimension \a n with centre at the origin and \a m generators. */
    explicit Zonotope(Nat d, Nat m);

    /*! \brief Construct from centre, generators, and a uniform error term. */
    explicit Zonotope(const Vector<Float64>& c, const Matrix<Float64>& G, const Vector<Float64>& e);
    /*! \brief Construct from centre and generators. */
    explicit Zonotope(const Vector<Float64>& c, const Matrix<Float64>& G);
    /*! \brief Construct from interval centre and a generator matrix. */
    explicit Zonotope(const Vector<ExactInterval>& c, const Matrix<Float64>& G);
    /*! \brief Construct from centre and an interval generator matrix. */
    explicit Zonotope(const Vector<Float64>& c, const Matrix<ExactInterval>& G);
    /*! \brief Construct from an interval centre and an interval generator matrix. */
    explicit Zonotope(const Vector<ExactInterval>& c, const Matrix<ExactInterval>& G);


    /*! \brief Construct a zonotope of dimension \a d with centre at the origin and \a m generators from the data beginning at \a ptr. */
    template<class XX> explicit Zonotope(Nat d, Nat m, const XX* ptr);

    /*! \brief Construct a zonotope of dimension \a d with \a m generators from raw data.
     *  The data format is (c0,G00,G01,...,G0m,e0,c1,G10,...,G1m,e1,...).
     */
    explicit Zonotope(InitializerList< std::tuple<Float64,InitializerList<Float64>,Float64> > lst);


    /*! \brief Convert from a box. */
    Zonotope(const ExactBox& r);
    /*! \brief Copy constructor. */
    Zonotope(const Zonotope& z);
    /*! \brief Copy assignment operator. */
    Zonotope& operator=(const Zonotope& z);
    /*! \brief Cloning operator. */
    Zonotope* clone() const;
    //@}


    //@{
    //! \name Logical predicates
    /*! \brief Test for equality of representation. */
    friend Bool operator==(const Zonotope& z1, const Zonotope& z2);
    //@}

    //@{
    //! \name Data access
    /*! \brief The dimension of the zonotope. */
    Nat dimension() const;

    /*! \brief The number of generators of the zonotope. */
    Nat number_of_generators() const;

    /*! \brief The domain. */
    Vector<ExactInterval> domain() const;

    /*! \brief The centre. */
    const Vector<Float64>& centre() const;

    /*! \brief The matrix of principle directions. */
    const Matrix<Float64>& generators() const;

    /*! \brief The uniform error bound. */
    const Vector<Float64>& error() const;

    /*! \brief A bounding box for the set. */
    UpperBox bounding_box() const;

    /*! \brief The radius of the set in the supremum norm. */
    Float64 radius() const;

    /*! \brief Test if the set contains a point. */
    Tribool contains(const ExactPoint& pt) const;

    /*! \brief Test if the set is disjoint from a box. */
    Tribool separated(const ExactBox& bx) const;
    /*! \brief Test if the set is a inside of a box. */
    Tribool inside(const ExactBox& bx) const;

    //@}


    //@{
    //! \name Geometric binary predicates
    /*! \brief Tests disjointness of \a z and \a r. */
    friend Tribool separated(const Zonotope& z, const ExactBox& r);
    /*! \brief Tests if \a z and \a r intersect. */
    friend Tribool overlaps(const Zonotope& z, const ExactBox& r);
    /*! \brief Tests inclusion of \a z in \a r. */
    friend Tribool inside(const Zonotope& z, const ExactBox& r);
    /*! \brief Tests disjointness of \a r and \a z. */
    friend Tribool separated(const ExactBox& r, const Zonotope& z);
    //@}

    //@{
    //! \name Approximation operations.
    /*! \brief Compute an simplified approximation of the zonotope \a z. */
    friend Zonotope approximation(const Zonotope& z);
    /*! \brief Compute an over-approximation of the zonotope \a z. */
    friend Zonotope over_approximation(const Zonotope& z);
    /*! \brief Compute an over-approximation of the zonotope \a z without a uniform error term. */
    friend Zonotope error_free_over_approximation(const Zonotope&);
    /*! \brief Compute an over-approximation of the zonotope \a z by a non-coordinate aligned orthotope. */
    friend Zonotope orthogonal_over_approximation(const Zonotope&);
    /*! \brief Compute an over-approximation of a zonotope \a z with nonsingular generator matrix. */
    friend Zonotope nonsingular_over_approximation(const Zonotope&);
    /*! \brief Compute a cascade-over-approximation of the zonotope \a z with \a b blocks of \a d generators. */
    friend Zonotope cascade_over_approximation(const Zonotope& z, Nat b);
    //@}

    //@{
    //! \name Function operations.
    /*! \brief Compute the image of \a z under a function given by the affine form \a af. */
    friend Zonotope apply(const Vector<Affine<ExactInterval>>& af, const Zonotope& z);
    friend Zonotope apply(const VectorFunction<ValidatedTag>& f, const Zonotope& z);
    //@}

    //@{
    //! \name Input/output.
    /*! \brief Write to an output stream. */
    OutputStream& write(OutputStream& os) const;
    /*! \brief Draw on a canvas. */
    Void draw(CanvasInterface& c, const Projection2d& p) const;
    //@}
};


Tribool empty(const Zonotope& z);
Tribool bounded(const Zonotope& z);
Float64 radius(const Zonotope& z);
ExactBox bounding_box(const Zonotope& z);


Tribool contains(const Zonotope& z, const ExactPoint& pt);
Tribool separated(const Zonotope& z, const ExactBox& r);
Tribool overlaps(const Zonotope& z, const ExactBox& r);
Tribool inside(const Zonotope& z, const ExactBox& r);

Tribool separated(const Zonotope& z1, const Zonotope& z2);

ListSet<Zonotope> split(const Zonotope& z);

Zonotope approximation(const Zonotope& z);
Zonotope over_approximation(const Zonotope& z);
Zonotope error_free_over_approximation(const Zonotope&);
Zonotope orthogonal_over_approximation(const Zonotope&);
Zonotope nonsingular_over_approximation(const Zonotope&);
Zonotope cascade_over_approximation(const Zonotope& z, Nat b);
Zonotope orthogonal_approximation(const Zonotope& z);

Zonotope apply(const Affine<ExactInterval>& af, const Zonotope& z);
Zonotope apply(const VectorFunction<ValidatedTag>& f, const Zonotope& z);

OutputStream& operator<<(OutputStream& os, const Zonotope& z);
InputStream& operator>>(InputStream& is, Zonotope& z);


template<class X> inline
Zonotope::Zonotope(Nat d, Nat m, const X* ptr)
    : _centre(d,ptr), _generators(d,m,ptr+d), _error(d)
{
}


} // namespace Ariadne

#endif // ARIADNE_ZONOTOPE_H
