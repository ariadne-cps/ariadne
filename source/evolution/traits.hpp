#ifndef ARIADNE_EVOLUTION_TRAITS_HPP
#define ARIADNE_EVOLUTION_TRAITS_HPP

// Associations between systems and their evolution algorithms belong to evolution,
// never to the system-description classes.
namespace Ariadne {

class IteratedMap;
class VectorField;
class DifferentialInclusion;
class HybridAutomatonInterface;

class IteratedMapEvolver;
class VectorFieldEvolver;
class DifferentialInclusionEvolver;
class HybridEvolverInterface;

class LabelledEnclosure;
class HybridEnclosure;

template<class SYS> struct EvolutionTraits;

template<> struct EvolutionTraits<IteratedMap> {
    using EvolverType = IteratedMapEvolver;
    using EnclosureType = LabelledEnclosure;
};

template<> struct EvolutionTraits<VectorField> {
    using EvolverType = VectorFieldEvolver;
    using EnclosureType = LabelledEnclosure;
};

template<> struct EvolutionTraits<DifferentialInclusion> {
    using EvolverType = DifferentialInclusionEvolver;
    using EnclosureType = LabelledEnclosure;
};

template<> struct EvolutionTraits<HybridAutomatonInterface> {
    using EvolverType = HybridEvolverInterface;
    using EnclosureType = HybridEnclosure;
};

} // namespace Ariadne
#endif // ARIADNE_EVOLUTION_TRAITS_HPP
