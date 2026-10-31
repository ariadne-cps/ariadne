#ifndef ARIADNE_ANALYSIS_HYBRID_STORAGE_ORBIT_HPP
#define ARIADNE_ANALYSIS_HYBRID_STORAGE_ORBIT_HPP

#include "analysis/hybrid_storage.hpp"
#include "evolution/hybrid_orbit.hpp"

namespace Ariadne {

template<>
class Orbit<HybridStorage>
    : public HybridDrawableInterface
{
    struct Data;
  public:
    typedef HybridStorage EnclosureListType;

    Orbit(const HybridStorage& initial_set);
    Orbit(const HybridStorage& initial_set, const HybridStorage& reach_set,
          const HybridStorage& intermediate_set, const HybridStorage& final_set);

    HybridStorage const& initial() const;
    HybridStorage const& reach() const;
    HybridStorage const& intermediate() const;
    HybridStorage const& final() const;

    Void draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d& v) const;
  private:
    SharedPointer<Data> _data;
};

} // namespace Ariadne

#endif // ARIADNE_ANALYSIS_HYBRID_STORAGE_ORBIT_HPP
