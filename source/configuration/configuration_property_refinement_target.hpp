/***************************************************************************
 *            configuration/configuration_property_refinement_target.hpp
 *
 *  Copyright  2007-20  Luca Geretti
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

/*! \file configuration/configuration_property_refinement_target.hpp
 *  \brief A set of objectives for refining a configuration based on a task result.
 */

#ifndef ARIADNE_CONFIGURATION_PROPERTY_REFINEMENT_TARGET_HPP
#define ARIADNE_CONFIGURATION_PROPERTY_REFINEMENT_TARGET_HPP

#include <functional>
#include "configurable.hpp"
#include "configuration_property_interface.hpp"

namespace Ariadne {

template<class R> struct TaskInput;
template<class R> struct TaskOutput;

enum class OptimisationCriterion;

//! \brief A set of objectives for refinement of a property.
//! \details
template<class R> class ConfigurationPropertyRefinementTarget {
  public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    typedef TaskObjective<R> ObjectiveType;

    ConfigurationPropertyRefinementTarget(ConfigurationPropertyPath const& path, List<ObjectiveType> const& objectives)
        : _path(path), _objectives(objectives) { }
    ConfigurationPropertyPath const& path() const { return _path; }
    List<ObjectiveType> const& objectives() const { return _objectives; }

    //! \brief Comparison for set appending
    //! \details Only one rule for a specific property is expected to be used.
    Bool operator<(ConfigurationPropertyRefinementTarget<R> const& r) const { return _path < r._path; }
  private:
    ConfigurationPropertyPath const _path;
    List<ObjectiveType> const _objectives;
};


} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_PROPERTY_REFINEMENT_TARGET_HPP
