#include "reasoning/storage_orbit.hpp"

namespace Ariadne {

Orbit<Storage>::
Orbit(const Storage& initial_set)
    : _data(new Orbit<Storage>::Data(initial_set.grid(),initial_set.auxiliary_mapping()))
{
    this->_data->initial=initial_set;
}


Orbit<Storage>::
Orbit(const Storage& initial_set,
      const Storage& reach_set,
      const Storage& intermediate_set,
      const Storage& final_set)
    : _data(new Data(initial_set.grid(),initial_set.auxiliary_mapping()))
{
    this->_data->initial=initial_set;
    this->_data->reach=reach_set;
    this->_data->intermediate=intermediate_set;
    this->_data->final=final_set;
}


Storage const&
Orbit<Storage>::
initial() const
{
    return this->_data->initial;
}

Storage const&
Orbit<Storage>::
reach() const
{
    return this->_data->reach;
}

Storage const&
Orbit<Storage>::
intermediate() const
{
    return this->_data->intermediate;
}

Storage const&
Orbit<Storage>::
final() const
{
    return this->_data->final;
}

} // namespace Ariadne
