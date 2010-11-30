/***************************************************************************
 *            geometry2d.h
 *
 *  Copyright 2009-10  Davide Bresolin, Pieter Collins
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

/*! \file graphics_interface.h
 *  \brief Base graphics interface from which all plotting and drawing classes are inherited.
 */

#ifndef ARIADNE_GEOMETRY2D_H
#define ARIADNE_GEOMETRY2D_H

typedef unsigned int uint;

namespace Ariadne {

template<class R, class A> inline R numeric_cast(const A&);

struct Vector2d {
    double x,y; Vector2d(double xx, double yy) : x(xx), y(yy) { }
    template<class X, class Y> Vector2d(const X& xx, const Y& yy) : x(numeric_cast<double>(xx)), y(numeric_cast<double>(yy)) { }
};
inline Vector2d operator-(const Vector2d& v) { return Vector2d(-v.x,-v.y); }
inline Vector2d operator+(const Vector2d& v1, const Vector2d& v2) { return Vector2d(v1.x+v2.x,v1.y+v2.y); }
inline Vector2d operator-(const Vector2d& v1, const Vector2d& v2) { return Vector2d(v1.x-v2.x,v1.y-v2.y); }
inline Vector2d operator*(const double& s1, const Vector2d& v2) { return Vector2d(s1*v2.x,s1*v2.y); }
inline std::ostream& operator<<(std::ostream& os, const Vector2d& v) {
    return os << "["<<v.x<<","<<v.y<<"]"; }

struct Point2d {
    double x,y; Point2d(double xx, double yy) : x(xx), y(yy) { }
    template<class X, class Y> Point2d(const X& xx, const Y& yy) : x(numeric_cast<double>(xx)), y(numeric_cast<double>(yy)) { }
};
inline bool operator==(Point2d& pt1, const Point2d& pt2) { return pt1.x==pt2.x && pt1.y==pt2.y; }
inline Point2d& operator+=(Point2d& pt, const Vector2d& v) { pt.x+=v.x; pt.y+=v.y; return pt; }
inline Point2d& operator-=(Point2d& pt, const Vector2d& v) { pt.x-=v.x; pt.y-=v.y; return pt; }
inline std::ostream& operator<<(std::ostream& os, const Point2d& pt) {
    return os << "("<<pt.x<<","<<pt.y<<")"; }

struct Box2d {
    double xl,xu,yl,yu; Box2d() { }
    Box2d(double xxl, double xxu, double yyl, double yyu) : xl(xxl), xu(xxu), yl(yl), yu(yyu) { }
};
inline std::ostream& operator<<(std::ostream& os, const Box2d& bx) {
    return os << "["<<bx.xl<<","<<bx.xu<<"]x["<<bx.yl<<","<<bx.yu<<"]"; }


} // namespace Ariadne


#endif // ARIADNE_GEOMETRY2D_H
