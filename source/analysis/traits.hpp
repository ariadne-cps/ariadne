#ifndef ARIADNE_ANALYSIS_TRAITS_HPP
#define ARIADNE_ANALYSIS_TRAITS_HPP

#include "evolution/traits.hpp"
#include "geometry/set_interface.hpp"
#include "hybrid/hybrid_set_interface.hpp"

// Storage and set-interface associations are analysis-level concerns.
// The generic space traits preserve the EuclideanSpace API during staging;
// hybrid specialisation moves those aliases out of the system-level HybridSpace.
namespace Ariadne {

class HybridSpace;
class LabelledStorage;
class HybridStorage;
class GridTreePaving;

template<class SYS> struct AnalysisTraits;

template<> struct AnalysisTraits<IteratedMap> {
    using StorageType = LabelledStorage;
};

template<> struct AnalysisTraits<VectorField> {
    using StorageType = LabelledStorage;
};

template<> struct AnalysisTraits<DifferentialInclusion> {
    using StorageType = LabelledStorage;
};

template<> struct AnalysisTraits<HybridAutomatonInterface> {
    using StorageType = HybridStorage;
};

template<class SPACE> struct SpaceAnalysisTraits {
    using BoundingDomainType = typename SPACE::BoundingDomainType;
    using BoundedSetInterfaceType = typename SPACE::BoundedSetInterfaceType;
    using OvertSetInterfaceType = typename SPACE::OvertSetInterfaceType;
    using OpenSetInterfaceType = typename SPACE::OpenSetInterfaceType;
    using CompactSetInterfaceType = typename SPACE::CompactSetInterfaceType;
    using LocatedSetInterfaceType = typename SPACE::LocatedSetInterfaceType;
    using RegularSetInterfaceType = typename SPACE::RegularSetInterfaceType;
    using RegularLocatedSetInterfaceType = typename SPACE::RegularLocatedSetInterfaceType;
    using SetApproximationType = typename SPACE::SetApproximationType;
};

// EuclideanSpace intentionally describes the geometry only; the
// discretisation used by reachability analysis is an analysis concern.
// Do not add SetApproximationType back to EuclideanSpace.
template<> struct SpaceAnalysisTraits<EuclideanSpace> {
    using BoundingDomainType = EuclideanSpace::BoundingDomainType;
    using BoundedSetInterfaceType = EuclideanSpace::BoundedSetInterfaceType;
    using OvertSetInterfaceType = EuclideanSpace::OvertSetInterfaceType;
    using OpenSetInterfaceType = EuclideanSpace::OpenSetInterfaceType;
    using CompactSetInterfaceType = EuclideanSpace::CompactSetInterfaceType;
    using LocatedSetInterfaceType = EuclideanSpace::LocatedSetInterfaceType;
    using RegularSetInterfaceType = EuclideanSpace::RegularSetInterfaceType;
    using RegularLocatedSetInterfaceType = EuclideanSpace::RegularLocatedSetInterfaceType;
    using SetApproximationType = GridTreePaving;
};

template<> struct SpaceAnalysisTraits<HybridSpace> {
    using BoundingDomainType = HybridExactBoxes;
    using BoundedSetInterfaceType = EffectiveHybridBoundedSetInterface;
    using OvertSetInterfaceType = EffectiveHybridOvertSetInterface;
    using OpenSetInterfaceType = EffectiveHybridOpenSetInterface;
    using CompactSetInterfaceType = EffectiveHybridCompactSetInterface;
    using LocatedSetInterfaceType = EffectiveHybridLocatedSetInterface;
    using RegularSetInterfaceType = EffectiveHybridRegularSetInterface;
    using RegularLocatedSetInterfaceType = EffectiveHybridRegularLocatedSetInterface;
    using SetApproximationType = HybridGridTreePaving;
};

} // namespace Ariadne
#endif // ARIADNE_ANALYSIS_TRAITS_HPP
