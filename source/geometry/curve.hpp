/***************************************************************************
 *            curve.hpp
 *
 *  Copyright  2007-8  Pieter Collins
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

/*! \file curve.hpp
 *  \brief A arbitraty curve in Euclidean space.
 */

#ifndef ARIADNE_CURVE_HPP
#define ARIADNE_CURVE_HPP

#include <memory>

#include "../utility/macros.hpp"
#include "../utility/stlio.hpp"
#include "../function/function.hpp"
#include "../output/graphics_interface.hpp"

namespace Ariadne {


template<class X> class Vector;

template<class X> class Point;
typedef Point<ExactNumericType> ExactPoint;
template<class IVL> class Box;
typedef Box<ExactIntervalType> ExactBoxType;
typedef Box<UpperIntervalType> UpperBoxType;


// Forward declarations for friends
class CurveInterface;

/*! \brief A curve in Euclidean space
 */
class CurveInterface
{
  public:
    typedef Dyadic GenericParameterType;
    typedef FloatDPValue ParameterType;
    typedef Point<FloatDPApproximation> PointType;
    typedef Vector<FloatDPApproximation> TangentVectorType;
  public:
    /*! \brief Destructor. */
    virtual ~CurveInterface() = default;
    /*! \brief Return a new dynamically-allocated copy of the curve. */
    virtual CurveInterface* clone() const = 0;
    /*! \brief The dimension of the space the curve lies in. */
    virtual DimensionType dimension() const = 0;
    /*! \brief The smoothness of the curve. */
    virtual DegreeType smoothness() const = 0;

    /*! \brief The point on the curve at a parameter value. */
    virtual PointType value(const ParameterType& s) const = 0;
    /*! \brief The tangent vector to the curve at a parameter value. */
    virtual TangentVectorType tangent(const ParameterType& s) const = 0;

    /*! \brief Write to an output stream. */
    virtual OutputStream& write(OutputStream& os) const = 0;
};

inline OutputStream& operator<<(OutputStream& os, const CurveInterface& c) {
    return c.write(os); }



/*! \brief A curve in Euclidean space
 */
class Curve
    : public CurveInterface
{
  public:
    typedef CurveInterface::GenericParameterType GenericParameterType;
    typedef CurveInterface::ParameterType ParameterType;
    typedef CurveInterface::PointType PointType;
    typedef CurveInterface::TangentVectorType TangentVectorType;
  public:
    /*! \brief Destructor. */
    virtual ~Curve() = default;
    /*! \brief Constructor. */
    Curve(const Function<EffectiveTag,IntervalDomainType,BoxDomainType>& f);
    /*! \brief Copy constructor. */
    Curve(const Curve& c);
    /*! \brief Return a new dynamically-allocated copy of the constraint. */
    virtual Curve* clone() const;
    /*! \brief The dimension of the set. */
    virtual DimensionType dimension() const;
    /*! \brief The smoothness of the curve. */
    virtual DegreeType smoothness() const;

    /*! \brief The value at a point. */
    virtual PointType value(const ParameterType& s) const;
    /*! \brief The tangent at a point. */
    virtual TangentVectorType tangent(const ParameterType& s) const;

    /*! \brief Write to an output stream. */
    virtual OutputStream& write(OutputStream& os) const;
  private:
    Function<EffectiveTag,IntervalDomainType,BoxDomainType> _function;
};



/*!\brief A line segment in Euclidean space. */
class InterpolatedCurve
    : public DrawableInterface
{
  public:
    typedef DoublePrecision PrecisionType;
    typedef CurveInterface::GenericParameterType GenericParameterType;
    typedef CurveInterface::ParameterType ParameterType;
    typedef CurveInterface::PointType PointType;
    typedef CurveInterface::TangentVectorType TangentVectorType;
  public:
    typedef std::map< ParameterType, PointType >::const_iterator ConstIterator;
  public:
    /*! \brief Create an empty curve. */
    InterpolatedCurve() : _points() { }
    /*! \brief Create a curve with a single point \a pt at parameter value 0. */
    InterpolatedCurve(const PointType& pt)
        : _points() { this->insert(0,pt); }
    /*! \brief Create a curve with a single point \a pt at parameter value \a s. */
    InterpolatedCurve(ParameterType s, const PointType& pt)
        : _points() { this->insert(s,pt); }
    InterpolatedCurve(GenericParameterType s, const PointType& pt)
        : _points() { PrecisionType pr; this->insert(ParameterType(s,pr),pt); }
    InterpolatedCurve(const RawFloatDP& s, const Vector<RawFloatDP>& pt)
        : _points() { this->insert(s,pt); }
    /*! \brief Create a segment from \a pt0 at parameter value 0 to \a pt1 at parameter value 1. */
    InterpolatedCurve(const PointType& pt0, const PointType& pt1)
        : _points() { this->insert(0,pt0); this->insert(1,pt1); }
    /*! \brief Insert a point with parameter value \a s and spacial value \a pt. */
    Void insert(const GenericParameterType& s, const PointType& pt);
    Void insert(const ParameterType& s, const PointType& pt);
    Void insert(const RawFloatDP& s, const Vector<RawFloatDP>& pt);

    /*! \brief The number of segments in the curve. */
    SizeType size() const { return this->_points.size(); }
    /*! \brief The dimension of the Euclidean space the line segment lies in. */
    DimensionType dimension() const { return this->_points.begin()->second.size(); }
    /*! \brief An Iterator to the first point in the curve. */
    ConstIterator begin() const { return this->_points.begin(); }
    /*! \brief An Iterator to the end point in the curve, NOT the one-past-the-end! */
    ConstIterator end() const { return --this->_points.end(); }

    /*! \brief A dynamically-allocated copy. */
    virtual InterpolatedCurve* clone() const;
    /*! \brief Draw on a two-dimensional canvas. */
    virtual Void draw(CanvasInterface& c, const Projection2d& p) const;
    /*! \brief A bounding box for the curve. */
    virtual UpperBoxType bounding_box() const;

    /*! \brief Write to an output stream. */
    virtual OutputStream& write(OutputStream& os) const;

  private:
    friend OutputStream& operator<<(OutputStream&, const InterpolatedCurve&);
  private:
    std::map< ParameterType, PointType > _points;
};

inline
OutputStream& operator<<(OutputStream& os, const InterpolatedCurve& curve) {
    return os << "InterpolatedCurve( size=" << curve.size() << ", points=" << curve._points << " )";
}



} // namespace Ariadne

#endif /* ARIADNE_CURVE_HPP */
