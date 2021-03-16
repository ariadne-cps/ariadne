/***************************************************************************
 *            output/gnuplot.hpp
 *
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



#include "../output/graphics.hpp"
#include "../config.hpp"
#include "../utility/string.hpp"
#include "algebra/tensor.hpp"
#include "numeric/float_bounds.hpp"
#include "geometry2d.hpp"
#include <string>
#include <ctype.h>

#ifdef HAVE_GNUPLOT_H

#include "gnuplot-iostream.hpp"

namespace Ariadne{

struct  _Range
{       
    double Xmin;
    double Xmax;
    double Ymin;
    double Ymax;
    double Zmin;
    double Zmax; 
};

struct  _Labels
{
    String xLabel = "";
    String yLabel = "";
    String zLabel = "";
};

class GnuplotCanvas : public CanvasInterface
{
friend class Figure;
private:
    Gnuplot *gnuplot;
    List<Point2d> geom;
    Colour lc, fc;
    double lw;
    _Range rng;
    Nat dim;
    Point2d Cpoint;
    double dr;
    bool isdot;
    Nat sizeX;
    Nat sizeY;
    bool isMultiplot;
    bool is2DPalette;
    bool is3DPalette;
    _Labels labels;
    GnuplotFileType fileType;

private:
    Void plot_2d(Array<double> data);

    Void plot_2d(Array<Array<double>> data);

    Void plot_3d(Array<Array<double>> data);

public:
    ~GnuplotCanvas();
    // Constructors - Create the canvas
    //Create canvas with dimensions
    GnuplotCanvas(String filename, GnuplotFileType fileType, Nat X = 800, Nat Y = 800);

    //CanvasInterface
    Void initialise(StringType x, StringType y, StringType z, double xl, double xu, double yl, double yu, double lz, double uz);
    Void initialise(StringType x, StringType y, double xl, double xu, double yl, double yu);
    Void write(const char* filename) const;
    Void finalise();
    Void move_to(double x, double y);
    Void line_to(double x, double y);
    Void circle(double x, double y, double r);
    Void dot(double x, double y);
    Void stroke();
    Void fill();
    Void set_dot_radius(double dr);
    Void set_line_width(double lw);
    Void set_line_colour(double r, double g, double b);
    Void set_fill_opacity(double o);
    Void set_fill_colour(double r, double g, double b);
    Vector2d scaling() const;
    Box2d bounds() const;

    Void set_3d_palette();

    Void plot_data(Array<double> data);
    Void plot_bounds(Array<Array<double>> bounds);
    Void plot_tensor_2d_image(Tensor<2, double> tensor);
    Void plot_tensor_3d_image(Tensor<3, double> tensor);
    Void plot_xy_projection(Tensor<3, double> tensor);
    Void plot_yz_projection(Tensor<3, double> tensor);
    Void plot_xz_projection(Tensor<3, double> tensor);
    //CanvasInterface

    //Set Multiplot - Multiple plot on same screen
    void set_multiplot(bool s);
    //Set Multiplot Layout
    void set_multiplot_layout(int nRow, int nCol, String Title);
    // Set X Label
    void set_x_label(String xLabel);
    // Set Y Label
    void set_y_label(String yLabel);
    // Set Z Label
    void set_z_label(String zLabel);
    // Set Title
    void set_title(String title);
    // Set Labels
    void set_xyz_label(String xLabel, String yLabel, String zLabel);
    // Set Labels and Title
    void set_labels(String xLabel, String yLabel, String zLabel, String title);
    // Set X, Y Range
    void set_range_2d(double minX, double maxX, double minY, double maxY);

    void set_range_3d(double minX, double maxX, double minY,  double maxY, double minZ, double maxZ);
    // Set X Log axis
    void set_x_log_axis();
    // Set Y Log axis
    void set_y_log_axis();
    // Set XY Log axis
    void set_xy_log_axis();
    // Set XZ Log axis
    void set_xz_log_axis();
    // Set YZ Log axis
    void set_yz_log_axis();
    // Set XYZ Log axis
    void set_xyz_log_axis();
    // Set Legend
    void set_legend();
    // Set View Projection of a 3D rapresentation
    void set_map();
    //Set 3D palette
    void set_3d_palette(double min, double max, double step);
    //Set 2D palette
    //void set2DPalette(Image2D& image, double min, double max, double step);
    void set_2d_palette(double min, double max, double step);
    //Unset colorbox
    void unset_color_box();

};

GnuplotCanvas::~GnuplotCanvas()
{
    delete gnuplot;
}

} // namespace Ariadne



#endif