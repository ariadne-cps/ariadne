/***************************************************************************
 *            verification/configuration_property_refinement_rule.hpp
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

/*! \file configuration/configuration_refinement_rule.hpp
 *  \brief A rule for refining a configuration based on a task result.
 */

#ifndef ARIADNE_CONFIGURATION_PROPERTY_REFINEMENT_RULE_HPP
#define ARIADNE_CONFIGURATION_PROPERTY_REFINEMENT_RULE_HPP

#include <functional>
#include "configuration/configurable.hpp"
#include "configuration/configuration_property_interface.hpp"
#include "verification/safety_objective_measure.decl.hpp"

namespace Ariadne {

template<class R> struct TaskInput;
template<class R> struct TaskOutput;

//! \brief A rule for refinement of a property.
template<class R> class ConfigurationPropertyRefinementRule {
  public:
    typedef double RatioType;
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    typedef Configuration<R> ConfigurationType;
    typedef typename R::ConfigurationRefinementObjectiveType ObjectiveType;
    typedef Pair<ConfigurationPropertyPath,RatioType> ResultType;
    typedef std::function<RatioType(InputType const&, OutputType const&, ObjectiveType const&)> FunctionType;

    ConfigurationPropertyRefinementRule(ConfigurationPropertyPath const& path, FunctionType const& function) : _path(path), _ratio_function(function) { }
    ConfigurationPropertyPath const& path() const { return _path; }
    //! \brief Find the ratio for refinement
    ResultType get_ratio(InputType const& i, OutputType const& o, ObjectiveType const& obj) { return make_pair(_path,_ratio_function(i,o,obj)); }
    //! \brief Comparison for set appending
    //! \details Only one rule for a specific property is expected to be used.
    Bool operator<(ConfigurationPropertyRefinementRule<R> const& r) const { return _path < r._path; }
  private:
    ConfigurationPropertyPath const _path;
    FunctionType const _ratio_function;

};


} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_PROPERTY_REFINEMENT_RULE_HPP
