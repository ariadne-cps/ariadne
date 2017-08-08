/***************************************************************************
 *            real_interface.hpp
 *
 *  Copyright 2013--17  Pieter Collins
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
 *  You should have received a copy of the GNU G3c767e04cec413f9afb4c30b521ca71ceb5b0409eneral Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iosfwd>

namespace Ariadne {

class DoublePrecision;
class MultiplePrecision;
class FloatMP;
class FloatDP;

template<class F> class Bounds;
using DyadicBounds = Bounds<Dyadic>;
using FloatDPBounds = Bounds<FloatDP>;
using FloatMPBounds = Bounds<FloatMP>;

using OutputStream = std::ostream;

class Real;

class RealInterface {
  public:
    virtual ~RealInterface() = default;
    virtual ValidatedReal _compute(Effort) const = 0;
    virtual FloatDPBounds _compute_get(DoublePrecision) const = 0;
    virtual FloatMPBounds _compute_get(MultiplePrecision) const = 0;
  public:
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

class ValidatedRealInterface {
  public:
    virtual ~ValidatedRealInterface() = default;
    virtual DyadicBounds _get() const = 0;
    virtual FloatDPBounds _get(DoublePrecision) const = 0;
    virtual FloatMPBounds _get(MultiplePrecision) const = 0;
  public:
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

} // namespace Ariadne
