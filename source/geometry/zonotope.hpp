/***************************************************************************
 *            zonotope.hpp
 *
 *  Copyright 2008-17  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file zonotope.hpp
 *  \brief Zonotopes in Euclidean space.
 */

#ifndef ARIADNE_ZONOTOPE_HPP
#define ARIADNE_ZONOTOPE_HPP

#include <iosfwd>

#include "../utility/declarations.hpp"
#include "../utility/tribool.hpp"

#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"

#include "../geometry/set_interface.hpp"

#include "../output/graphics_interface.hpp"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;

template<class P, class X, class XE=X> class Zonotope;
template<class X, class XE=X> using ValidatedZonotope=Zonotope<ValidatedTag,X,XE>;

class ValidatedAffineConstrainedImageSet;
template<class BS> class ListSet;

template<class X> class Affine;

class Figure;


/*!\brief A zonotope of arbitrary dimension.
 *
 * A zonotope is a set of the form \f$c+Ge\f$, where \f$||e||_{\infty}\leq1\f$.
 * The columns of the matrix \f$G\f$ are the <em>generators</em> of the
 * zonotope.
 *
 * Zonotopes are always is_bounded.
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
template<class X, class XE> class Zonotope<ValidatedTag,X,XE>
    : public ValidatedCompactSetInterface
    , public DrawableInterface
{
    typedef ValidatedTag P;
    typedef typename X::PrecisionType PR;
    typedef typename XE::PrecisionType PRE;
  private:
    Vector<Value<X>> _centre;
    Matrix<Value<X>> _generators;
    Vector<Error<XE>> _error;
  public:
    //@{
    //! \name Constructors and destructors
    /*! \brief Virtual destructor. */
    virtual ~Zonotope();
    /*! \brief Default constructor yields a zonotope with dimension zero and no generators. */
    explicit Zonotope();
    /*! \brief Construct a zonotope of dimension \a d with no generators. */
    explicit Zonotope(DimensionType d);
    /*! \brief Construct a zonotope of dimension \a n with centre at the origin and \a m generators. */
    explicit Zonotope(DimensionType d, SizeType m);

    /*! \brief Construct from centre and generators. */
    explicit Zonotope(const Vector<X>& c, const Matrix<X>& G);
    /*! \brief Construct from centre, generators, and a uniform error term. */
    explicit Zonotope(const Vector<X>& c, const Matrix<X>& G, const Vector<XE>& e);

    /*! \brief Construct from centre, generators, and a uniform error term. */
    explicit Zonotope(const Vector<Value<X>>& c, const Matrix<Value<X>>& G, const Vector<Error<XE>>& e);

    /*! \brief Construct from centre and generators. */
    explicit Zonotope(const Vector<Value<X>>& c, const Matrix<Value<X>>& G);
    /*! \brief Construct from interval centre and a generator matrix. */
    explicit Zonotope(const Vector<Bounds<X>>& c, const Matrix<Value<X>>& G);
    /*! \brief Construct from centre and an interval generator matrix. */
    explicit Zonotope(const Vector<Value<X>>& c, const Matrix<Bounds<X>>& G);
    /*! \brief Construct from an interval centre and an interval generator matrix. */
    explicit Zonotope(const Vector<Bounds<X>>& c, const Matrix<Bounds<X>>& G);
    /*! \brief Construct from an interval centre, a generator matrix, and a uniform error term. */
    explicit Zonotope(const Vector<Bounds<X>>& c, const Matrix<Value<X>>& G, const Vector<Error<XE>>& e);


    /*! \brief Construct a zonotope of dimension \a d with \a m generators from raw data.
     *  The data format is (c0,G00,G01,...,G0m,e0,c1,G10,...,G1m,e1,...).
     */
    explicit Zonotope(InitializerList< std::tuple<X,InitializerList<X>,XE> > lst);


    /*! \brief Convert from a point. */
    Zonotope(const Point<Value<X>>& pt);
    Zonotope(const Point<Bounds<X>>& pt);
    /*! \brief Convert from a box. */
    Zonotope(const Box<Interval<UpperBound<X>>>& bx);
    /*! \brief Copy constructor. */
    Zonotope(const Zonotope& z);
    /*! \brief Copy assignment operator. */
    Zonotope& operator=(const Zonotope& z);
    /*! \brief Cloning operator. */
    Zonotope* clone() const;

    /*! \brief Construct a zonotope of dimension \a d with centre at the origin and \a m generators from the data beginning at \a ptr. */
    explicit operator ValidatedAffineConstrainedImageSet () const;

    //@}


    //@{
    //! \name Logical predicates
    /*! \brief Test for equality of representation. */
    friend Bool operator==(const Zonotope& z1, const Zonotope& z2);
    //@}

    //@{
    //! \name Data access
    /*! \brief The dimension of the zonotope. */
    DimensionType dimension() const;

    /*! \brief The number of generators of the zonotope. */
    SizeType number_of_generators() const;

    /*! \brief The domain. */
    Box<ExactIntervalType> domain() const;

    /*! \brief The centre. */
    const Vector<Value<X>>& centre() const;

    /*! \brief The matrix of principle directions. */
    const Matrix<Value<X>>& generators() const;

    /*! \brief The uniform error bound. */
    const Vector<Error<XE>>& error() const;

    /*! \brief A bounding box for the set. */
    UpperBoxType bounding_box() const;

    /*! \brief The radius of the set in the supremum norm. */
    PositiveFloatDPUpperBound radius() const;

    /*! \brief Test if the set contains a point. */
    ValidatedLowerKleenean contains(const ExactPoint& pt) const;

    /*! \brief Test if the set is disjoint from a box. */
    ValidatedLowerKleenean separated(const ExactBoxType& bx) const;
    /*! \brief Test if the set is a inside of a box. */
    ValidatedLowerKleenean inside(const ExactBoxType& bx) const;

    //@}


    //@{
    //! \name Geometric predicates

    /*! \brief Test if \a z is empty. Always returns \c false. */
    friend ValidatedKleenean empty(const Zonotope& z) { return false; }
    /*! \brief Test if \a z is bounded. Always returns \c true. */
    friend ValidatedKleenean is_bounded(const Zonotope& z) { return true; }
    /*! \brief An upper bound for the radius of \a z. */
    friend PositiveFloatDPUpperBound radius(const Zonotope& z) { return z.radius(); }
    /*! \brief A bounding box for the \a z. */
    friend UpperBoxType bounding_box(const Zonotope& z);

    /*! \brief Test if the z contains a point. */
    friend ValidatedLowerKleenean contains(const Zonotope& z, const ExactPoint& pt) { return z.contains(pt); }
    /*! \brief Tests inclusion of \a z in \a bx. */
    friend ValidatedLowerKleenean inside(const Zonotope& z, const ExactBoxType& bx) { return z.inside(bx); }
    /*! \brief Tests disjointness of \a z and \a bx. */
    friend ValidatedLowerKleenean separated(const Zonotope& z, const ExactBoxType& bx) { return z.separated(bx); }
    /*! \brief Tests disjointness of \a bx and \a z. */
    friend ValidatedLowerKleenean separated(const ExactBoxType& bx, const Zonotope& z) { return z.separated(bx); }
    /*! \brief Tests disjointness of \a bx and \a z. */
    friend ValidatedLowerKleenean separated(const Zonotope& z1, const Zonotope& z2) { return _separated(z1,z2); }
    //@}



    //@{
    //! \name Geometric operations

    /*! \brief The Cartesian product of two zonotopes. */
    friend Zonotope product(const Zonotope& z1, const Zonotope& z2) { return _product(z1,z2); }
    friend Zonotope product(const Zonotope& z1, const Interval<UpperBound<X>>& ivl2) { return _product(z1,ivl2); }
    //@}

    //@{
    //! \name Approximation operations.
    /*! \brief Compute an simplified approximation of the zonotope \a z. */
    friend Zonotope approximation(const Zonotope& z);
    /*! \brief Compute an over-approximation of the zonotope \a z. */
    friend Zonotope over_approximation(const Zonotope& z);
    /*! \brief Compute an over-approximation of the zonotope \a z without a uniform error term. */
    friend Zonotope error_free_over_approximation(const Zonotope& z) { return _error_free_over_approximation(z); }
    /*! \brief Compute an over-approximation of the zonotope \a z by a non-coordinate aligned orthotope. */
    friend Zonotope orthogonal_over_approximation(const Zonotope& z) { return _orthogonal_over_approximation(z); }
    /*! \brief Compute an over-approximation of a zonotope \a z with nonsingular generator matrix. */
    friend Zonotope nonsingular_over_approximation(const Zonotope& z) { return _nonsingular_over_approximation(z); }
    /*! \brief Compute a cascade-over-approximation of the zonotope \a z with \a blocks blocks of \a d generators. */
    friend Zonotope cascade_over_approximation(const Zonotope& z, SizeType blocks) { return _cascade_over_approximation(z,blocks); }

    /*! \brief Split into smaller subsets. */
    friend ListSet<Zonotope> split(const Zonotope& z);

    //@}

    //@{
    //! \name Function operations.
    /*! \brief Compute the image of \a z under a function given by the affine form \a af. */
    friend Zonotope apply(const Vector<Affine<ExactIntervalType>>& af, const Zonotope& z) { return _apply_affine(af,z); }
    friend Zonotope apply(const VectorMultivariateFunction<ValidatedTag>& f, const Zonotope& z) { return _apply(f,z); }
    friend Zonotope image(const Zonotope& z, const VectorMultivariateFunction<ValidatedTag>& f) { return _apply(f,z); }
    //@}

    //@{
    //! \name Input/output.
    /*! \brief Write to an output stream. */
    OutputStream& write(OutputStream& os) const;
    /*! \brief Write to an output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Zonotope& z) { return z.write(os); }
    /*! \brief Draw on a canvas. */
    Void draw(CanvasInterface& c, const Projection2d& p) const;
    //@}
  private:
    static Zonotope _apply_affine(const Vector<Affine<ExactIntervalType>>& f, const Zonotope& z);
    static Zonotope _apply(const VectorMultivariateFunction<ValidatedTag>& f, const Zonotope& z);

    static Zonotope _product(const Zonotope& z1, const Zonotope& z2);
    static Zonotope _product(const Zonotope& z1, const Interval<UpperBound<X>>& ivl2);

    static ListSet<Zonotope> _split(const Zonotope& z);

    static ValidatedLowerKleenean _separated(const Zonotope& z1, const Zonotope& z2);

    static Zonotope _approximation(const Zonotope& z);
    static Zonotope _over_approximation(const Zonotope& z);
    static Zonotope _error_free_over_approximation(const Zonotope&);
    static Zonotope _orthogonal_over_approximation(const Zonotope&);
    static Zonotope _nonsingular_over_approximation(const Zonotope&);
    static Zonotope _cascade_over_approximation(const Zonotope& z, SizeType blocks);
    static Zonotope _orthogonal_approximation(const Zonotope& z);
};






} // namespace Ariadne

#endif // ARIADNE_ZONOTOPE_HPP
