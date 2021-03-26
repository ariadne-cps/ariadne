
#include "tensor.hpp"
#include "numeric/float.decl.hpp"
#include "algebra/algebra.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_value.hpp"

#include <limits.h>


namespace Ariadne {
//DRAWABLE INTERFACE
    template <SizeType N, class X>
    Void Tensor<N, X>::draw(CanvasInterface& canvas, const Projection2d& p) const {

    }

    template <SizeType N, class X>
    Void Tensor<N, X>::draw3d(CanvasInterface& canvas, const Projection3d& p, ProjType proj) const {

    }
template<>
Void Tensor<2ul, Ariadne::Vector<Ariadne::Bounds<Ariadne::FloatDP>>>::draw3d(Ariadne::CanvasInterface&, Ariadne::Variables3d const&, Ariadne::ProjType) const
{ }

    template <>
    Void Tensor<2, Value<RawFloatType<DoublePrecision>>>::draw(CanvasInterface& canvas, const Projection2d& p) const {
        //2D Drawing
        for(SizeType frame=0; frame!=_ns[1]; ++frame){
            SizeType index = _index({0, frame});
            canvas.move_to(0.0, _a[index].get_d());
            for(SizeType i=1; i!=_ns[0]; ++i){
                index = _index({i, frame});
                canvas.line_to(numeric_cast<double>(i), _a[index].get_d());
            }
            canvas.fill();       
        }
    }

    template <>
    Void Tensor<2, Value<RawFloatType<DoublePrecision>>>::draw3d(CanvasInterface& canvas, const Projection3d& p, ProjType proj) const {

    }

    template <>
    Void Tensor<2, Value<RawFloatType<MultiplePrecision>>>::draw(CanvasInterface& canvas, const Projection2d& p) const {
        //2D Drawing
        for(SizeType frame=0; frame!=_ns[1]; ++frame){
            SizeType index = _index({0, frame});
            canvas.move_to(0.0, _a[index].get_d());
            for(SizeType i=1; i!=_ns[0]; ++i){
                index = _index({i, frame});
                canvas.line_to(numeric_cast<double>(i), _a[index].get_d());
            }
            canvas.fill();       
        }
    }

    template <>
    Void Tensor<2, Value<RawFloatType<MultiplePrecision>>>::draw3d(CanvasInterface& canvas, const Projection3d& p, ProjType proj) const {

    }

    template <>
    Void Tensor<3, Value<RawFloatType<DoublePrecision>>>::draw(CanvasInterface& canvas, const Projection2d& p) const {

    }

    template <>
    Void Tensor<3, Value<RawFloatType<DoublePrecision>>>::draw3d(CanvasInterface& canvas, const Projection3d& p, ProjType proj) const {
        if(proj == ProjType::no_proj){
            for(SizeType frame=0; frame!=_ns[2]; ++frame){
                SizeType index = _index({0, 0, frame});
                canvas.move_to(0.0, _a[index].get_d());
                for(SizeType x2=0; x2!=_ns[1]; ++x2){
                    for(SizeType x1=0; x1!=_ns[0]; ++x1){
                        if (x2 == 0 && x1 == 0){
                            x1++;
                        }
                        index = _index({x1, x2, frame});
                        canvas.line_to(numeric_cast<double>(x1), _a[index].get_d());
                    }
                    canvas.line_to(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max());
                }
                canvas.fill3d();
            }
        }else if(proj ==ProjType::x1_proj){
            for(SizeType frame=0; frame!=_ns[2]; ++frame){
                for(SizeType x2=0; x2!=_ns[1]; ++x2){
                    SizeType index = _index({0, x2, frame});
                    canvas.move_to(0.0, _a[index].get_d());
                    for(SizeType x1=1; x1!=_ns[0]; ++x1){
                        index = _index({x1, x2, frame});
                        canvas.line_to(numeric_cast<double>(x1), _a[index].get_d());
                    }
                    canvas.fill();
                }
            } 
        }else if(proj == ProjType::x2_proj){
            for(SizeType frame=0; frame!=_ns[2]; ++frame){
                for(SizeType x1=0; x1!=_ns[1]; ++x1){
                    SizeType index = _index({x1, 0, frame});
                    canvas.move_to(0.0, _a[index].get_d());
                    for(SizeType x2=1; x2!=_ns[0]; ++x2){
                        index = _index({x1, x2, frame});
                        canvas.line_to(numeric_cast<double>(x1), _a[index].get_d());
                    }
                    canvas.fill();
                }
            } 
        }
    }

    template <>
    Void Tensor<3, Value<RawFloatType<MultiplePrecision>>>::draw(CanvasInterface& canvas, const Projection2d& p) const {

    }

    template <>
    Void Tensor<3, Value<RawFloatType<MultiplePrecision>>>::draw3d(CanvasInterface& canvas, const Projection3d& p, ProjType proj) const {
        if(proj == ProjType::no_proj){
            for(SizeType frame=0; frame!=_ns[2]; ++frame){
                SizeType index = _index({0, 0, frame});
                canvas.move_to(0.0, _a[index].get_d());
                for(SizeType x2=0; x2!=_ns[1]; ++x2){
                    for(SizeType x1=0; x1!=_ns[0]; ++x1){
                        if (x2 == 0 && x1 == 0){
                            x1++;
                        }
                        index = _index({x1, x2, frame});
                        canvas.line_to(numeric_cast<double>(x1), _a[index].get_d());
                    }
                    canvas.line_to(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max());
                }
                canvas.fill3d();
            }
        }else if(proj ==ProjType::x1_proj){
            for(SizeType frame=0; frame!=_ns[2]; ++frame){
                for(SizeType x2=0; x2!=_ns[1]; ++x2){
                    SizeType index = _index({0, x2, frame});
                    canvas.move_to(0.0, _a[index].get_d());
                    for(SizeType x1=1; x1!=_ns[0]; ++x1){
                        index = _index({x1, x2, frame});
                        canvas.line_to(numeric_cast<double>(x1), _a[index].get_d());
                    }
                    canvas.fill();
                }
            } 
        }else if(proj == ProjType::x2_proj){
            for(SizeType frame=0; frame!=_ns[2]; ++frame){
                for(SizeType x1=0; x1!=_ns[1]; ++x1){
                    SizeType index = _index({x1, 0, frame});
                    canvas.move_to(0.0, _a[index].get_d());
                    for(SizeType x2=1; x2!=_ns[0]; ++x2){
                        index = _index({x1, x2, frame});
                        canvas.line_to(numeric_cast<double>(x1), _a[index].get_d());
                    }
                    canvas.fill();
                }
            } 
        }
    }

//LABELLED DRAWABLE INTERFACE

    template <SizeType N, class X>
    Void Tensor<N, X>::draw(CanvasInterface& canvas, const Variables2d& p) const {

    }

    template <SizeType N, class X>
    Void Tensor<N, X>::draw3d(CanvasInterface& canvas, const Variables3d& p, ProjType proj) const {

    }


    template <>
    Void Tensor<2, Value<RawFloatType<DoublePrecision>>>::draw(CanvasInterface& canvas, const Variables2d& p) const {
        //2D Drawing
        for(SizeType frame=0; frame!=_ns[1]; ++frame){
            SizeType index = _index({0, frame});
            canvas.move_to(0.0, _a[index].get_d());
            for(SizeType i=1; i!=_ns[0]; ++i){
                index = _index({i, frame});
                canvas.line_to(numeric_cast<double>(i), _a[index].get_d());
            }
            canvas.fill();       
        }
    }

    template <>
    Void Tensor<2, Value<RawFloatType<DoublePrecision>>>::draw3d(CanvasInterface& canvas, const Variables3d& p, ProjType proj) const {

    }

    template <>
    Void Tensor<2, Value<RawFloatType<MultiplePrecision>>>::draw(CanvasInterface& canvas, const Variables2d& p) const {
        //2D Drawing
        for(SizeType frame=0; frame!=_ns[1]; ++frame){
            SizeType index = _index({0, frame});
            canvas.move_to(0.0, _a[index].get_d());
            for(SizeType i=1; i!=_ns[0]; ++i){
                index = _index({i, frame});
                canvas.line_to(numeric_cast<double>(i), _a[index].get_d());
            }
            canvas.fill();       
        }
    }

    template <>
    Void Tensor<2, Value<RawFloatType<MultiplePrecision>>>::draw3d(CanvasInterface& canvas, const Variables3d& p, ProjType proj) const {

    }

    template <>
    Void Tensor<3, Value<RawFloatType<DoublePrecision>>>::draw(CanvasInterface& canvas, const Variables2d& p) const {

    }

    template <>
    Void Tensor<3, Value<RawFloatType<DoublePrecision>>>::draw3d(CanvasInterface& canvas, const Variables3d& p, ProjType proj) const {
        if(proj == ProjType::no_proj){
            for(SizeType frame=0; frame!=_ns[2]; ++frame){
                SizeType index = _index({0, 0, frame});
                canvas.move_to(0.0, _a[index].get_d());
                for(SizeType x2=0; x2!=_ns[1]; ++x2){
                    for(SizeType x1=0; x1!=_ns[0]; ++x1){
                        if (x2 == 0 && x1 == 0){
                            x1++;
                        }
                        index = _index({x1, x2, frame});
                        canvas.line_to(numeric_cast<double>(x1), _a[index].get_d());
                    }
                    canvas.line_to(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max());
                }
                canvas.fill3d();
            }
        }else if(proj ==ProjType::x1_proj){
            for(SizeType frame=0; frame!=_ns[2]; ++frame){
                //SizeType index = _index({0, 0, frame});
                //canvas.move_to(0.0, _a[index].get_d());
                for(SizeType x2=0; x2!=_ns[1]; ++x2){
                    SizeType index = _index({0, x2, frame});
                    canvas.move_to(0.0, _a[index].get_d());
                    for(SizeType x1=1; x1!=_ns[0]; ++x1){
                        index = _index({x1, x2, frame});
                        canvas.line_to(numeric_cast<double>(x1), _a[index].get_d());
                    }
                    canvas.fill();
                }
            } 
        }else if(proj == ProjType::x2_proj){
            for(SizeType frame=0; frame!=_ns[2]; ++frame){
                //SizeType index = _index({0, 0, frame});
                //canvas.move_to(0.0, _a[index].get_d());
                for(SizeType x1=0; x1!=_ns[1]; ++x1){
                    SizeType index = _index({x1, 0, frame});
                    canvas.move_to(0.0, _a[index].get_d());
                    for(SizeType x2=1; x2!=_ns[0]; ++x2){
                        index = _index({x1, x2, frame});
                        canvas.line_to(numeric_cast<double>(x1), _a[index].get_d());
                    }
                    canvas.fill();
                }
            } 
        }
    }

    template <>
    Void Tensor<3, Value<RawFloatType<MultiplePrecision>>>::draw(CanvasInterface& canvas, const Variables2d& p) const {

    }

    template <>
    Void Tensor<3, Value<RawFloatType<MultiplePrecision>>>::draw3d(CanvasInterface& canvas, const Variables3d& p, ProjType proj) const {
        if(proj == ProjType::no_proj){
            for(SizeType frame=0; frame!=_ns[2]; ++frame){
                SizeType index = _index({0, 0, frame});
                canvas.move_to(0.0, _a[index].get_d());
                for(SizeType x2=0; x2!=_ns[1]; ++x2){
                    for(SizeType x1=0; x1!=_ns[0]; ++x1){
                        if (x2 == 0 && x1 == 0){
                            x1++;
                        }
                        index = _index({x1, x2, frame});
                        canvas.line_to(numeric_cast<double>(x1), _a[index].get_d());
                    }
                    canvas.line_to(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max());
                }
                canvas.fill3d();
            }
        }else if(proj ==ProjType::x1_proj){
            for(SizeType frame=0; frame!=_ns[2]; ++frame){
                //SizeType index = _index({0, 0, frame});
                //canvas.move_to(0.0, _a[index].get_d());
                for(SizeType x2=0; x2!=_ns[1]; ++x2){
                    SizeType index = _index({0, x2, frame});
                    canvas.move_to(0.0, _a[index].get_d());
                    for(SizeType x1=1; x1!=_ns[0]; ++x1){
                        index = _index({x1, x2, frame});
                        canvas.line_to(numeric_cast<double>(x1), _a[index].get_d());
                    }
                    canvas.fill();
                }
            } 
        }else if(proj == ProjType::x2_proj){
            for(SizeType frame=0; frame!=_ns[2]; ++frame){
                //SizeType index = _index({0, 0, frame});
                //canvas.move_to(0.0, _a[index].get_d());
                for(SizeType x1=0; x1!=_ns[1]; ++x1){
                    SizeType index = _index({x1, 0, frame});
                    canvas.move_to(0.0, _a[index].get_d());
                    for(SizeType x2=1; x2!=_ns[0]; ++x2){
                        index = _index({x1, x2, frame});
                        canvas.line_to(numeric_cast<double>(x1), _a[index].get_d());
                    }
                    canvas.fill();
                }
            } 
        }
    }
}