/***************************************************************************
 *            io/tensor_drawing.hpp
 *
 *  Copyright  2026  Ariadne contributors
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

/*! \file io/tensor_drawing.hpp
 *  \brief Drawing adapters for algebraic Tensor data.
 */

#ifndef ARIADNE_IO_TENSOR_DRAWING_HPP
#define ARIADNE_IO_TENSOR_DRAWING_HPP

#include <limits>
#include <type_traits>

#include "algebra/tensor.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "io/figure.hpp"
#include "io/graphics_interface.hpp"

namespace Ariadne {

template<SizeType N, class PR>
class TensorDrawable
    : public LabelledDrawable2d3dInterface
    , public Drawable2d3dInterface
{
    using X = Float<PR>;
    Tensor<N,X> _tensor;

  public:
    explicit TensorDrawable(Tensor<N,X> const& tensor) : _tensor(tensor) { }

    TensorDrawable* clone() const override { return new TensorDrawable(*this); }
    TensorDrawable* clone2d3d() const override { return new TensorDrawable(*this); }

    DimensionType dimension() const override { return _tensor.rank(); }

    Void draw(CanvasInterface& canvas, Projection2d const& projection) const override {
        if constexpr (N == 2u) {
            this->_draw_rank2(canvas);
        } else if constexpr (N == 3u) {
            ARIADNE_ASSERT(projection.argument_size() == this->dimension());
            this->_draw_rank3_projection(canvas, projection.x_coordinate(), projection.y_coordinate());
        }
    }

    Void draw(CanvasInterface& canvas, Projection3d const&) const override {
        if constexpr (N == 3u) {
            this->_draw_rank3(canvas);
        }
    }

    Void draw(CanvasInterface& canvas, Variables2d const& variables) const override {
        if constexpr (N == 2u) {
            this->_draw_rank2(canvas);
        } else if constexpr (N == 3u) {
            // Preserve the previous MultiplePrecision labelled-2d behaviour.
            if constexpr (!std::is_same_v<PR,MultiplePrecision>) {
                this->_draw_rank3_projection(canvas, this->_coordinate(variables.x()), this->_coordinate(variables.y()));
            }
        }
    }

    Void draw(CanvasInterface& canvas, Variables3d const&) const override {
        if constexpr (N == 3u) {
            this->_draw_rank3(canvas);
        }
    }

  private:
    static DimensionType _coordinate(RealVariable const& variable) {
        if (variable == RealVariable("x")) { return 0u; }
        if (variable == RealVariable("y")) { return 1u; }
        if (variable == RealVariable("z")) { return 2u; }
        return 3u;
    }

    Void _draw_rank2(CanvasInterface& canvas) const {
        for (SizeType frame=0; frame!=_tensor.size(1); ++frame) {
            canvas.move_to(0.0, _tensor[{0,frame}].get_d());
            for (SizeType i=1; i!=_tensor.size(0); ++i) {
                canvas.line_to(numeric_cast<double>(i), _tensor[{i,frame}].get_d());
            }
            canvas.fill();
        }
    }

    Void _draw_rank3_projection(CanvasInterface& canvas, DimensionType ix, DimensionType iy) const {
        if (ix == 0u && iy == 1u) {
            canvas.set_heat_map(true);
            for (SizeType frame=0; frame!=_tensor.size(2); ++frame) {
                canvas.move_to(0.0, _tensor[{0,0,frame}].get_d());
                for (SizeType x2=0; x2!=_tensor.size(1); ++x2) {
                    for (SizeType x1=0; x1!=_tensor.size(0); ++x1) {
                        if (x2 == 0u && x1 == 0u) { ++x1; }
                        canvas.line_to(numeric_cast<double>(x1), _tensor[{x1,x2,frame}].get_d());
                    }
                    canvas.line_to(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max());
                }
                canvas.fill_3d();
            }
        } else if (ix == 1u && iy == 0u) {
            canvas.set_heat_map(true);
            for (SizeType frame=0; frame!=_tensor.size(2); ++frame) {
                canvas.move_to(0.0, _tensor[{0,0,frame}].get_d());
                for (SizeType x1=0; x1!=_tensor.size(1); ++x1) {
                    for (SizeType x2=0; x2!=_tensor.size(0); ++x2) {
                        if (x2 == 0u && x1 == 0u) { ++x2; }
                        canvas.line_to(numeric_cast<double>(x1), _tensor[{x1,x2,frame}].get_d());
                    }
                    canvas.line_to(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max());
                }
                canvas.fill_3d();
            }
        } else if (ix == 0u && iy == 2u) {
            for (SizeType frame=0; frame!=_tensor.size(2); ++frame) {
                for (SizeType x2=0; x2!=_tensor.size(1); ++x2) {
                    canvas.move_to(0.0, _tensor[{0,x2,frame}].get_d());
                    for (SizeType x1=1; x1!=_tensor.size(0); ++x1) {
                        canvas.line_to(numeric_cast<double>(x1), _tensor[{x1,x2,frame}].get_d());
                    }
                    canvas.fill();
                }
            }
        } else if (ix == 2u && iy == 0u) {
            for (SizeType frame=0; frame!=_tensor.size(2); ++frame) {
                for (SizeType x2=0; x2!=_tensor.size(1); ++x2) {
                    canvas.move_to(0.0, _tensor[{0,x2,frame}].get_d());
                    for (SizeType x1=1; x1!=_tensor.size(0); ++x1) {
                        canvas.line_to(_tensor[{x1,x2,frame}].get_d(), numeric_cast<double>(x1));
                    }
                    canvas.fill();
                }
            }
        } else if (ix == 1u && iy == 2u) {
            for (SizeType frame=0; frame!=_tensor.size(2); ++frame) {
                for (SizeType x1=0; x1!=_tensor.size(1); ++x1) {
                    canvas.move_to(0.0, _tensor[{x1,0,frame}].get_d());
                    for (SizeType x2=1; x2!=_tensor.size(0); ++x2) {
                        canvas.line_to(numeric_cast<double>(x1), _tensor[{x1,x2,frame}].get_d());
                    }
                    canvas.fill();
                }
            }
        } else if (ix == 2u && iy == 1u) {
            for (SizeType frame=0; frame!=_tensor.size(2); ++frame) {
                for (SizeType x1=0; x1!=_tensor.size(1); ++x1) {
                    canvas.move_to(0.0, _tensor[{x1,0,frame}].get_d());
                    for (SizeType x2=1; x2!=_tensor.size(0); ++x2) {
                        canvas.line_to(_tensor[{x1,x2,frame}].get_d(), numeric_cast<double>(x1));
                    }
                    canvas.fill();
                }
            }
        }
    }

    Void _draw_rank3(CanvasInterface& canvas) const {
        for (SizeType frame=0; frame!=_tensor.size(2); ++frame) {
            canvas.move_to(0.0, _tensor[{0,0,frame}].get_d());
            for (SizeType x2=0; x2!=_tensor.size(1); ++x2) {
                for (SizeType x1=0; x1!=_tensor.size(0); ++x1) {
                    if (x2 == 0u && x1 == 0u) { ++x1; }
                    canvas.line_to(numeric_cast<double>(x1), _tensor[{x1,x2,frame}].get_d());
                }
                canvas.line_to(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max());
            }
            canvas.fill_3d();
        }
    }
};

template<SizeType N, class PR>
TensorDrawable<N,PR> tensor_drawable(Tensor<N,Float<PR>> const& tensor) {
    return TensorDrawable<N,PR>(tensor);
}

} // namespace Ariadne

#endif // ARIADNE_IO_TENSOR_DRAWING_HPP
