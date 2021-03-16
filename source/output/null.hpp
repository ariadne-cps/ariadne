/***************************************************************************
 *            output/null.hpp
 *
 *  Copyright  2020-21  Mirko Albanese, Luca Geretti
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

/*! \file output/null.hpp
 *  \brief Null output in the case of missing graphics support.
 */

#include "output/graphics.hpp"

#include "config.hpp"

namespace Ariadne {

#if not(defined(HAVE_CAIRO_H)) || not(defined(HAVE_GNUPLOT_H))

class NullCanvas : public CanvasInterface
{
public:
    virtual Void initialise(StringType x, StringType y, StringType z, double lx, double ux, double ly, double uy, double lz, double uz) { };
    virtual Void initialise(StringType x, StringType y, double lx, double ux, double ly, double uy) { }
    virtual Void finalise() { }

    virtual Void write(const char* filename) const { }

    virtual Void move_to(double x, double y) { }
    virtual Void line_to(double x, double y) { }
    virtual Void circle(double x, double y, double r) { }
    virtual Void dot(double x, double y) { }
    virtual Void stroke() { }
    virtual Void fill() { }
    virtual Void set_line_width(double lw) { }
    virtual Void set_line_colour(double r, double g, double b) { }
    virtual Void set_fill_opacity(double fo) { }
    virtual Void set_fill_colour(double r, double g, double b) { }

    virtual Void plot_data(Array<double> data) { }
    virtual Void plot_bounds(Array<Array<double>> bounds) { }
    virtual Void plot_tensor_2d_image(Tensor<2, double> tensor) { }
    virtual Void plot_tensor_3d_image(Tensor<3, double> tensor) { }
    virtual Void plot_xz_projection(Tensor<3, double> tensor) { }
    virtual Void plot_xy_projection(Tensor<3, double> tensor) { }
    virtual Void plot_yz_projection(Tensor<3, double> tensor) { }

    virtual Vector2d scaling() const { return Vector2d(0,0); }
    virtual Box2d bounds() const { return Box2d(0,0,0,0); }
};

#endif

} // namespace Ariadne


