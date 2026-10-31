// Adapters between evolved enclosures and analysis-level grid storage.
#include "function/functional.hpp"
#include "analysis/storage.hpp"
#include "evolution/enclosure.hpp"
#include "paving/grid_paving.hpp"
#include "symbolic/space.hpp"
#include "logging/logging.hpp"

namespace Ariadne {

Void Enclosure::adjoin_outer_approximation_to(Storage& storage, Nat fineness) const
{
    if (this->auxiliary_mapping().result_size()!=0 &&
            this->auxiliary_mapping().managed_pointer()==storage.auxiliary_mapping().managed_pointer())
    {
        ARIADNE_WARN("enclosure="<<*this<<", "
                    "enclosure.auxiliary_mapping()="<<this->auxiliary_mapping()<<", "
                    "storage.auxiliary_mapping()="<<storage.auxiliary_mapping());
    }
    this->configuration().paver().adjoin_outer_approximation(storage.state_set(),this->state_set(),fineness);
}


Storage Enclosure::outer_approximation(const Grid& grid, Nat fineness) const
{
    Storage paving(grid,this->auxiliary_mapping());
    this->adjoin_outer_approximation_to(paving,fineness);
    return paving;
}

LabelledStorage::LabelledStorage(Grid const& grid,
                                 RealSpace const& state_space,
                                 EffectiveVectorMultivariateFunction const& auxiliary_mapping,
                                 RealSpace const& auxiliary_space)
    : LabelledStorage(GridTreePaving(grid), state_space, auxiliary_mapping, auxiliary_space) { }

LabelledStorage::LabelledStorage(GridTreePaving const& paving,
                                 RealSpace const& state_space,
                                 EffectiveVectorMultivariateFunction const& auxiliary_mapping,
                                 RealSpace const& auxiliary_space)
    : Storage(paving,auxiliary_mapping)
    , _state_variables(state_space.variable_names())
    , _auxiliary_variables(auxiliary_space.variable_names())
{ }


LabelledGridTreePaving inner_approximation(EffectiveEuclideanSetInterface const& set, LabelledGrid const& grid, Nat fineness) {
    LabelledGridTreePaving paving(grid); paving.euclidean_set().adjoin_inner_approximation(set,fineness); return paving;
}

const RealSpace LabelledStorage::state_space() const {
    return RealSpace(this->_state_variables);
}

const RealSpace LabelledStorage::auxiliary_space() const {
    return RealSpace(this->_auxiliary_variables);
}

const RealSpace LabelledStorage::state_auxiliary_space() const {
    return join(this->state_space(),this->auxiliary_space());
}

const LabelledGrid LabelledStorage::grid() const {
    return LabelledGrid(this->euclidean_set().grid(),this->state_space());
}

const Mapping LabelledStorage::auxiliary_function() const {
    return this->Storage::auxiliary_mapping();
}

const LabelledMapping LabelledStorage::auxiliary_mapping() const {
    return LabelledMapping(this->state_space(), this->auxiliary_function(), this->auxiliary_space());
}

const LabelledMapping LabelledStorage::auxiliary_data() const {
    return this->auxiliary_mapping();
}

Void LabelledStorage::draw(CanvasInterface& canvas, const Variables2d& axes) const
{
    Projection2d proj=projection(this->state_auxiliary_space(),axes);
    this->euclidean_set().draw(canvas,proj);
}

} // namespace Ariadne
