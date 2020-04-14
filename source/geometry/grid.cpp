/***************************************************************************
 *            geometry/grid.cpp
 *
 *  Copyright  2008-20  Ivan S. Zapreev, Pieter Collins
 *
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

#include "../function/functional.hpp"
#include "../config.hpp"

#include <iostream>
#include <iomanip>

#include "../utility/macros.hpp"
#include "../utility/exceptions.hpp"
#include "../utility/stlio.hpp"
#include "../geometry/function_set.hpp"
#include "../geometry/list_set.hpp"
#include "../geometry/grid.hpp"


namespace Ariadne {

using ExactNumericType = Grid::ExactNumericType;
using Double = double;

struct Grid::Data
{
    Vector<FloatDP> _origin;
    Vector<FloatDP> _lengths;
};

Grid::~Grid()
{
}

Grid::Grid()
    : _data(new Data())
{
}

Grid::Grid(const Grid& gr)
    : _data(gr._data)
{
}

Grid& Grid::operator=(const Grid& gr)
{
    this->_data=gr._data; return *this;
}

Grid::Grid(Nat d)
    : _data(new Data())
{
    Vector<FloatDP> origin(d,FloatDP(0));
    Vector<FloatDP> lengths(d,FloatDP(1));
    this->_create(origin,lengths);
}

Grid::Grid(Nat d, FloatDP l)
    : _data(new Data())
{
    Vector<FloatDP> origin(d,FloatDP(0));
    Vector<FloatDP> lengths(d,l);
    this->_create(origin,lengths);
}

Grid::Grid(const Vector<FloatDP>& lengths)
    : _data(new Data())
{
    Vector<FloatDP> origin(lengths.size(),FloatDP(0));
    this->_create(origin,lengths);
}

Grid::Grid(const Vector<FloatDP>& origin, const Vector<FloatDP>& lengths)
    : _data(new Data())
{
    if(origin.size() != lengths.size()) {
        throw IncompatibleSizes(ARIADNE_PRETTY_FUNCTION);
    }
    this->_create(origin,lengths);
}

Void Grid::_create(const Vector<FloatDP>& origin, const Vector<FloatDP>& lengths)
{
    this->_data->_origin=origin;
    this->_data->_lengths=lengths;
}

DimensionType Grid::dimension() const
{
    return this->_data->_lengths.size();
}

const Vector<FloatDP>& Grid::origin() const
{
    return this->_data->_origin;
}

const Vector<FloatDP>& Grid::lengths() const
{
    return this->_data->_lengths;
}

ExactNumericType Grid::coordinate(Nat d, DyadicType x) const
{
    return ExactNumericType(add(approx,this->_data->_origin[d],mul(approx,this->_data->_lengths[d],x)));
}

ExactNumericType Grid::subdivision_coordinate(Nat d, DyadicType x) const
{
    return ExactNumericType(add(approx,this->_data->_origin[d],mul(approx,this->_data->_lengths[d],x)));
}

ExactNumericType Grid::subdivision_coordinate(Nat d, IntegerType n) const
{
    return ExactNumericType(add(approx,this->_data->_origin[d],mul(approx,this->_data->_lengths[d],n)));
}

Int Grid::subdivision_index(Nat d, const ExactNumericType& x) const
{
    FloatDP half=0.5;
    Int n=integer_cast<Int>(floor(add(approx,div(approx,sub(approx,x.raw(),this->_data->_origin[d]),this->_data->_lengths[d]),half)));
    FloatDP sc=add(approx,this->_data->_origin[d],mul(approx,this->_data->_lengths[d],n));
    if(sc == x.raw()) {
        return n;
    } else {
        ARIADNE_THROW(InvalidGridPosition,std::setprecision(20)<<"Grid::subdivision_index(Nat d,ExactNumericType x)","d="<<d<<", x="<<x<<", this->origin[d]="<<this->_data->_origin[d]<<", this->lengths[d]="<<this->_data->_lengths[d]<<" (closest value is "<<sc<<")");
    }
}

Int Grid::subdivision_lower_index(Nat d, const LowerNumericType& x) const
{
    Int n=integer_cast<Int>(floor(div(down,sub(down,x.raw(),this->_data->_origin[d]),this->_data->_lengths[d])));
    if(x.raw()>=add(approx,this->_data->_origin[d],mul(approx,this->_data->_lengths[d],(n+1)))) {
        return n+1;
    } else {
        return n;
    }
}

Int Grid::subdivision_upper_index(Nat d, const UpperNumericType& x) const
{
    Int n=integer_cast<Int>(ceil(div(up,sub(up,x.raw(),this->_data->_origin[d]),this->_data->_lengths[d])));
    if(x.raw()<=add(approx,this->_data->_origin[d],mul(approx,this->_data->_lengths[d],(n-1)))) {
        return n-1;
    } else {
        return n;
    }
}

Bool Grid::operator==(const Grid& g) const
{
    if(this->_data==g._data) {
        return true;
    } else {
        return this->_data->_origin==g._data->_origin && this->_data->_lengths==g._data->_lengths;
    }
}

Bool Grid::operator!=(const Grid& g) const
{
    return !(*this==g);
}

Array<Double> Grid::index(const ExactPointType& pt) const
{
    Array<double> res(pt.size());
    for(SizeType i=0; i!=res.size(); ++i) {
        res[i]=subdivision_index(i,pt[i]);
    }
    return res;
}

Array<Double> Grid::lower_index(const ExactBoxType& bx) const
{
    Array<double> res(bx.size());
    for(SizeType i=0; i!=res.size(); ++i) {
        res[i]=subdivision_lower_index(i,bx[i].lower());
    }
    return res;
}

Array<Double> Grid::upper_index(const ExactBoxType& bx) const
{
    Array<double> res(bx.size());
    for(SizeType i=0; i!=res.size(); ++i) {
        res[i]=subdivision_upper_index(i,bx[i].upper());
    }
    return res;
}

ExactPointType Grid::point(const Array<IntegerType>& a) const
{
    Vector<FloatDPValue> res(a.size());
    for(SizeType i=0; i!=res.size(); ++i) {
        res[i]=cast_exact(add(near,this->_data->_origin[i],mul(near,this->_data->_lengths[i],a[i])));
    }
    return res;
}

ExactPointType Grid::point(const Array<DyadicType>& a) const
{
    Vector<FloatDPValue> res(a.size());
    for(SizeType i=0; i!=res.size(); ++i) {
        res[i]=cast_exact(add(near,this->_data->_origin[i],mul(near,this->_data->_lengths[i],a[i])));
    }
    return res;
}

ExactBoxType Grid::box(const Array<DyadicType>& lower, const Array<DyadicType>& upper) const
{
    Vector<ExactIntervalType> res(lower.size());
    for(SizeType i=0; i!=res.size(); ++i) {
        res[i]=ExactIntervalType(this->subdivision_coordinate(i,lower[i]),
                        this->subdivision_coordinate(i,upper[i]));
    }
    return res;
}

Grid join(Grid const& g1, Grid const& g2) {
    return Grid(join(g1._data->_origin,g2._data->_origin),join(g1._data->_lengths,g2._data->_lengths));
}

OutputStream& operator<<(OutputStream& os, const Grid& gr)
{
    os << "Grid( ";
    os << "origin=" << gr.origin() << ", ";
    os << "lengths=" << gr.lengths() << " )";
    return os;
}

}
