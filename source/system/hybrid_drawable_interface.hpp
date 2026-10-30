// Abstract drawing hook used by hybrid objects. No rendering dependency.
#ifndef ARIADNE_SYSTEM_HYBRID_DRAWABLE_INTERFACE_HPP
#define ARIADNE_SYSTEM_HYBRID_DRAWABLE_INTERFACE_HPP

#include "system/discrete_location.hpp"
#include "utility/container.hpp"

namespace Ariadne {

class CanvasInterface;
struct Variables2d;

//! \ingroup GraphicsModule
//! \brief Base interface for drawable objects
class HybridDrawableInterface {
  public:
    //! brief The type of data needed to project to a 2d image.
    typedef Variables2d ProjectionType;
    //! brief Virtual destructor.
    virtual ~HybridDrawableInterface() = default;
    //! brief Draw the object on the canvas \a c using line segments and fill/stroke commands.
    virtual Void draw(CanvasInterface& c, const Set<DiscreteLocation>& q, const Variables2d& v) const = 0;
};


} // namespace Ariadne

#endif // ARIADNE_SYSTEM_HYBRID_DRAWABLE_INTERFACE_HPP
